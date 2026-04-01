from io import BytesIO
import zipfile
from celery import shared_task
import time
import os
from google.cloud import storage, run_v2
import uuid
from datetime import timedelta

from rest_framework import status
from rest_framework.views import APIView
from rest_framework_simplejwt.views import TokenObtainPairView
from rest_framework.response import Response
from .serializer import VarSerializer, FreqSerializer, CustomTokenObtainPairSerializer, SubSerializer, CNVarSerializer, VarCountSerializer
from .models import Variant, Frequency, Submission, CNV, VarCounts, UploadJob
from rest_framework.permissions import AllowAny, IsAuthenticated
from django.conf import settings
from django.db import transaction, connection 
from django.core.exceptions import *
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from django.db.models import Q, F
from django.db.models.expressions import  RawSQL
from django.http import HttpResponse, StreamingHttpResponse
from channels.layers import get_channel_layer
from asgiref.sync import async_to_sync
import re, csv, json
import xml.etree.ElementTree as ET
import xmltodict
import requests
from drf_spectacular.utils import extend_schema, OpenApiExample, OpenApiResponse,OpenApiRequest

# Create your views here.


# Retrieves data from the SLGVD given a user query
def FilterArgs(query, db):

    print("In FilterArgs")

    args = {} # For the key:value pairs that makeup the database queries
    op = None # For the logical operators seperating the individual key:value pairs
    field = None

    # Identify fields specified by user
    field_pat = re.compile(r"(?:\[)(chr|pos|consequence|gene|start_pos|end_pos)(?:\])")
    match = re.findall(field_pat, query)

    field = match
    print(field)

    # Identify logical operators specified by user
    lop_pat = re.compile(r"AND|OR")
    match_1 = re.findall(lop_pat, query)

    op = match_1
    print(op)

    # Assign field values and store in args dict
    def FilterFieldVals(query,i):

                # Search variants by a single value for a specified field
                if "," not in query:
                    if 'pos' in field[i] :
                        query = query.replace("[pos]","")
                        args['position'] = query
                    
                    elif 'chr' in field[i]:
                        query = query.replace("[chr]","")
                        args['chromosome'] = query

                    elif 'consequence' in field[i]:
                        query = query.replace("[consequence]","")
                        args['consequence__iexact'] = query

                    elif 'gene' in field[i]:
                        query = query.replace("[gene]","")
                        args['gene_name__iexact'] = query

                # Search variants by multiple values for a specified field
                if "," in query:
                    
                    if 'pos' in field[i]:
                        query = query.replace("[pos]","")
                        pos = [val.strip() for val in query.split(",")]
                        args['position__in'] = pos

                    elif 'chr' in field[i]:
                        query = query.replace("[chr]","")
                        chr = [val.strip() for val in query.split(",")]
                        args['chromosome__in'] = chr

                    elif 'consequence' in field[i]:
                        query = query.replace("[consequence]","")
                        con = [val.strip() for val in query.split(",")]
                        args['consequence__in'] = con

                    elif 'gene' in field[i]:
                        query = query.replace("[gene]","")
                        genes = [val.strip() for val in query.split(",")]
                        args['gene_name__in'] = genes
        
    # Search variants by a positional range
    if ":" in query:
        
        range = [val.strip() for val in query.split(":")]
        args['chromosome__iexact'] = range[0]

        if db == 'ssv':
            args['position__range'] = (range[1],range[2])
        elif db == 'gcnv':
            args['start_pos__gte'] = range[1]
            args['end_pos__lte'] = range[2]

        op = ['AND']

    # General search
    else:

        # General search with fields specified
        if len(field) != 0:

            if len(field) == 1:

                FilterFieldVals(query,0)
                
            elif len(field) > 1:
                values_and_fields = [term for term in query.split() if term not in op]
                print(values_and_fields)

                for i,val in enumerate(field):
                    
                    FilterFieldVals(values_and_fields[i],i)
                    
        # General search without fields specified
        else:
            
            # Idenitfying chromosome in query
            chr_pattern = re.compile(r'^(\d{1,2}|[XY])$')
            chr = re.findall(chr_pattern, query)
            if len(chr) != 0: 
                args['chromosome__iexact'] = chr[0]

            # Identifying position in query
            pos_pattern = re.compile(r'^(\d+)$')
            pos = re.findall(pos_pattern, query)
            if len(pos) != 0: 
                args['position__in'] = pos

            # Identifying gene_name in query
            gen_pattern = re.compile(r'^(?:[A-Za-z][A-Za-z0-9-.]*)$')
            genes = re.findall(gen_pattern, query)
            if len(genes) != 0: 
                    args['gene_name__iexact'] = genes[0]

            # Identifying consequence in query
            con_pattern = re.compile(r'^(?:[A-Za-z _]*)$')
            con = re.findall(con_pattern, query)
            if len(con) != 0: 
                args['consequence__iexact'] = con[0]

            # Identifying variant id in query
            variant_pattern = re.compile(r'^(?:\d{1,2}|[XY])p\d+r\w+a\w+$')
            variants = re.findall(variant_pattern, query)
            if len(variants) != 0: 
                args['variation_id__iexact'] = variants[0]

            op = ['OR']
  
    
    return args,op

class RetrieveDBData(APIView):

    permission_classes = [AllowAny]
    
    # Details for the Swagger UI API documentation
    @extend_schema(

        summary = "Retrieve Variant Data",
        description = "Retrieve information on variants including the homozygous and heterozygous counts of each variant allele and total number of individuals whose variants are stored for calculation of allele frequencies.",

        request = {

            "properties":{
                "type":"object",
                "body":{
                    "type":"object",
                    "items": {"type":"string"}
                },
                "example":"{\"body\": \"{\\\"query\\\":\\\"hes4\\\"}\"}",
            }
            
        },

        responses = {
            200: OpenApiResponse(
                response = {
                    "type": "object",
                    "properties" : {
                        "var_data": {"type":"array", "items": {"type": "object"}},
                        "freq_data": {"type": "array", "items": {"type": "object"}},
                        "sub_data": {"type": "array", "items": {"type": "object"}},
                        "tot_individuals": {"type": "integer"}
                    }
                },
                examples = [
                    OpenApiExample(
                        "Example Response",
                        value = {
                            'var_data': [{'variation_id': '1p935222rCaA', 'chromosome': 1, 'position': 935222, 'ref_allele': 'C', 'alt_allele': 'A', 'gene_name': 'HES4', 'consequence': 'Missense', 'submission_id': 2}, {'variation_id': '1p935459rAaG', 'chromosome': 1, 'position': 935459, 'ref_allele': 'A', 'alt_allele': 'G', 'gene_name': 'HES4', 'consequence': '5 Prime UTR', 'submission_id': 2}],

                            'freq_data': [{'variation_id': '1p935222rCaA', 'homo_count': 1, 'het_count': 0, 'last_updated': {'hom': [2]}}, {'variation_id': '1p935459rAaG', 'homo_count': 1, 'het_count': 0, 'last_updated': {'hom': [2]}}],

                            'sub_data': [{'submission_id': 2, 'no_individuals': 1, 'submission_date': '2024-12-09', 'username_id': '1A'}],

                            'tot_individuals': 21,
                        }
                    )
                ]
            )
        }
    )

    def post(self, request, format = None):

        # Parsing the request body for the query and db values
        raw_body = request.body.decode('utf-8')
        data = json.loads(raw_body)
        print(data)

        query = json.loads(data['body'])['query']
        db = json.loads(data['body'])['db_type']
        load_more_times = json.loads(data['body'])['load_more_times']
        
        print(query)
        print(db)
        print(load_more_times)

        try:

            # Constructing the queries for the database abstraction API

            # Extracting the arguments and logical operators for the queries from the users submitted search keywords
            args,op = FilterArgs(query, db)

            # Limiting the number of variants retrieved
            limit = 2500
            offset = (load_more_times - 1) * limit

            print(args,op)

            if len(args) > 0:
                
                # If the number of search fields is > 1. For example for the search query, SAMD11[gene] AND 1[chr]
                if len(args) > 1:
                    query_0 = Q()

                    # If only one logical operator is used to seperate search keywords. For example, SAMD11[gene] AND 1[chr]
                    if len(op) == 1:
                        if op[0] == 'OR':
                            for field,val in args.items():
                                query_0 |= Q(**{field:val})

                        elif op[0] == 'AND':

                            for field,val in args.items():
                                query_0 &= Q(**{field:val})

                    # If more than one logical operator is used to seperate search keywords. For example, SAMD11[gene] AND 1[chr] OR 2[chr]
                    elif len(op) > 1:
                        z = 0 
                        for field,val in args.items():

                            if op[z] == 'OR':
                                query_0 |= Q(**{field:val})
                            
                            elif op[z] == 'AND':
                                query_0 &= Q(**{field:val})

                            z += 1

                    print(query_0)

                    
                    # If the user is searching for short sequence variants (variants of length < 50 bp), the Variant model object filters data from ls50_variants db_table based on the query contructed.
                    if db == 'ssv':
                        variants = Variant.objects.filter(query_0).order_by('variation_id')[offset: offset + 2500]
                    
                    # If the user is searching for germline CNVs, , the CNV model object filters data from mt50_variants db_table based on the query contructed.
                    elif db == 'gcnv':
                        variants = CNV.objects.filter(query_0).order_by('variation_id')[offset: offset + 2500]

                # If the number of search fields is 1.
                elif len(args) == 1:

                    if db == "ssv":
                        variants = Variant.objects.filter(**args).order_by('variation_id')[offset: offset + 2500]
                    elif db == "gcnv":
                        variants = CNV.objects.filter(**args).order_by('variation_id')[offset: offset + 2500]
                    

                print(variants)

                # If there are variants retrieved from the db.
                if len(variants) > 0:

                    if db == 'ssv':

                        # Serialize the variant data retrieved to translate it to json format using the VarSerializer.
                        var_data = VarSerializer(variants, many = True).data               

                        # Obtaining frequency data of each these variants
                        variant_ids = [var['variation_id'] for var in var_data]
                        #print(variant_ids)

                        query_1 = Q()

                        for variant in variant_ids:
                            query_1 |= Q(variation_id__exact = variant)

                        #print(query_1)
                        print(len(var_data))

                        frequencies = Frequency.objects.filter(query_1)
                        freq_data = FreqSerializer(frequencies, many = True).data

                        # Obtaining all submission to calculate total number of individuals
                        sub_ids = set(Variant.objects.values_list("submission_id", flat=True))

                        sub_data = SubSerializer(Submission.objects.filter(submission_id__in = sub_ids), many=True).data
                        tot_individuals = sum([sub['no_individuals'] for sub in sub_data])

                        print(len(var_data), len(freq_data), len(sub_data), tot_individuals)

                    elif db == 'gcnv':

                        # Serialize the variant data retrieved that is translate it to a workable data format using the VarSerializer.
                        var_data = CNVarSerializer(variants, many = True).data 

                        # Obtaining number of individuals whose CNV data are archived in SLGVD
                        sub_ids = list(set([var['submission_id'] for var in var_data]))
                        print(sub_ids)

                        sub_data = SubSerializer(Submission.objects.filter(submission_id__in = sub_ids), many = True).data
                        print(sub_data)

                        tot_individuals = sum([sub['no_individuals'] for sub in sub_data])

                    # Final response generated for SSVs
                    if db == 'ssv':
                        return Response({"var_data": var_data,
                                        "freq_data": freq_data,
                                        "sub_data":sub_data,
                                        "tot_individuals":tot_individuals},
                                        status = status.HTTP_200_OK)
                    
                    # Final response generated for gcnv
                    elif db == 'gcnv':
                        return Response({"var_data": var_data,
                                         "sub_data":sub_data,
                                         "tot_individuals":tot_individuals},
                                        status = status.HTTP_200_OK)

            return Response({'message':'Not Found'}, status=status.HTTP_404_NOT_FOUND)

        except Exception as e:
            print(f"An error occurred: {e}")
 
# Retrieves additional information on variants stored in the SLGVD 
# from the dbSNP database
class RetrieveExtData(APIView):

    permission_classes = [AllowAny]

    @extend_schema(exclude = True)
    def post(self, request, format = None):

        raw_body = request.body.decode('utf-8')
        data = json.loads(raw_body)
        print(data)
    
        try:

            # Posting request to Entrez system
            query_eutil = json.loads(data['body'])['query_eutil']
            print(query_eutil)

            db = 'snp'
            base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
            search_url = f"{base_url}esearch.fcgi?db={db}&term={query_eutil}&usehistory=y"

            output = requests.get(search_url)
            print(output.text)


            if output:
                
                root = ET.fromstring(output.text)
                webenv = root.findtext("WebEnv")
                query_key = root.findtext("QueryKey")

                fetch_url = f"{base_url}efetch.fcgi?db={db}&query_key={query_key}&WebEnv={webenv}"

                data_eutil_xml = f"<variant>{requests.get(fetch_url).text}</variant>"

                data_eutil_json = json.dumps(xmltodict.parse(data_eutil_xml))

                print(data_eutil_json)
                
                return Response({'data_eutil':data_eutil_json},
                                status=status.HTTP_200_OK)
            
            else:
                return Response('Data not found', status=status.HTTP_404_NOT_FOUND)


        except Exception as e:
            print(e)
            return Response(f"An error occurred :{e}", status=status.HTTP_400_BAD_REQUEST)

# Retrieves a summary of the data submitted by an authenticated 
# user
class RetrieveDataSumByUser(APIView):

    permission_classes = [IsAuthenticated]

    @extend_schema(exclude = True)
    def post(self, request, format = None):

        username = self.request.user.username

        try:

            # Retrieve total variant count by user
            count = Variant.objects.filter(submission__username = username).count()
            count += CNV.objects.filter(submission__username = username).count()

            # Retrieve last submission date by user
            last_submission = SubSerializer(Submission.objects.filter(username = username).order_by("-submission_date")[0]).data
            last_sub_date = last_submission['submission_date']

            if count and last_submission:

                return Response({"count":count,
                             "last_sub_date":last_sub_date},
                             status=status.HTTP_200_OK)

            
        except Exception as e:
            print(e)

            return Response({"count":0,
                             "last_sub_date": "No previous submissions"},
                             status=status.HTTP_200_OK)

# Retrieves total count of short sequence variants (SSVs) and gCNVs 
# and the total number of variants by type of consequence
class RetrieveDBSum(APIView):

    permission_classes = [AllowAny]

    def post(self, request, format = None):

        try:

            with transaction.atomic():
                
                var_counts = VarCountSerializer(VarCounts.objects.get(entry_id = 1)).data

                print(var_counts)
                return Response({
                    "var_counts": var_counts
                    },
                    status=status.HTTP_200_OK
                )

        except Exception as e:

            return Response(
                {"Error":f"An error occurred {e}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR)

class CustomTokenObtainPairView(TokenObtainPairView):
    serializer_class = CustomTokenObtainPairSerializer

def UpdateVarCounts(variant_objects, db_type, action_type):

    # Dictionary of Sequnce Ontologies (SO) for variant consequences
    # that have been grouped under their parent SOs
    parent_SO = {
        'Exonic Variants': ['Missense','Synonymous',
        '5 Prime UTR', '3 Prime UTR', 'Non-coding exon', 'Frameshift','Stop loss','Initiator codon','Stop gain','Inframe insertion','Disruptive inframe deletion','Inframe deletion','Disruptive inframe insertion','Start loss','Stop retained','5 Prime UTR premature start codon gain'],
        'Intronic Variants': ['Intronic',],
        'Intergenic Variants':['Upstream gene',
        'Downstream gene','Intergenic',],
        'Regulatory Site Variant': ['TF Binding site'],
        'Splice Variants':['Splice region','Splice donor','Exon loss','Splice acceptor'],
        'Feature Ablation':['TFBS Ablation'],
        'Intragenic Variants':['Intragenic']
    }

    var_counts = VarCounts.objects.get(entry_id = 1)

    if db_type == "ssv":

        ssv_count = len(variant_objects)
        print(ssv_count)

        # Count the number of each consequence type
        # in var_data
        consequences = [var['consequence'] for var in variant_objects]

        #print(consequences)

        # Method to count the number of each type 
        # of variant consequence 
        def count_con(consequences, parent_SO):

            con_counts = {}
            none_n = 0

            for con in consequences:
                if con == None:
                    none_n += 1
                    continue

                cons = con.split("/")
                cons = [v.replace("variant","").strip() for v in cons]

                # Here, the parent SO term for a given consequence is identified 
                # from the parent_SO dict. Parent SO term refers to the Sequence Ontology term under which that particular consequence is categorized.
                parent_cons = [k for k,v in parent_SO.items() if any(item in v for item in cons)]

                for parent_con in parent_cons:
                    if parent_con in con_counts.keys():
                        con_counts[parent_con] += 1
                    else:
                        con_counts[parent_con] = 1

            print(f"None: {none_n}")
            return con_counts

        con_counts = count_con(consequences, parent_SO)

        print(con_counts)

        # Maps each parent variant to its respective count field
        con_relation_mapping = {
            'Exonic Variants': 'exo_count',
            'Intronic Variants': 'intro_count' ,
            'Intergenic Variants': 'inter_count',
            'Regulatory Site Variant': 'reg_count',
            'Splice Variants': 'spli_count',
            'Intragenic Variants': 'intra_count'
        }

        if action_type == "upload":
            var_counts.ssv_count += ssv_count

            for con,count in con_counts.items():

                field = con_relation_mapping[con]
                print(field)

                setattr(var_counts, field, getattr(var_counts,field) + count)

        elif action_type == "remove":
            var_counts.ssv_count -= ssv_count

            for con,count in con_counts.items():

                field = con_relation_mapping[con]
                print(field)

                setattr(var_counts, field, getattr(var_counts,field) - count)
        
    if db_type == "gcnv":

        if action_type == "upload":
            var_counts.cnv_count += len(variant_objects)
        elif action_type == "remove":
            var_counts.cnv_count -= len(variant_objects)

    var_counts.save()


# Upload the VCF file to GCS
def upload_to_gcs(vcf_files, bucket_name = "slgvd-uploads"):

    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)

    file_paths = []
    
    for file in vcf_files:
        # Create a unique filename to prevent overwriting
        unique_filename = f"{uuid.uuid4()}_{file.name}"
        blob = bucket.blob(unique_filename)
        blob.upload_from_file(file)

        file_paths.append(f"gs://{bucket_name}/{unique_filename}")

    return file_paths

# Method to trigger vcf parser in management/commands
def trigger_vcf_parse_job(job_id, gcs_paths):

    client = run_v2.JobsClient()
    job_path = f"projects/hgu-variationdb/locations/asia-south1/jobs/vcf-parser"

    print(job_path)
    
    for path in gcs_paths:
        overrides = {
            "container_overrides": [
                {
                    "args": ["python", "manage.py", "parse_vcf", "--id", str(job_id),
                             "--path", path]
                }
            ]
        }

        # 2. Build the request object
        request = run_v2.RunJobRequest(
            name=job_path,
            overrides=overrides
        )

        # Trigger the start and return immediately
        cloud_operation = client.run_job(request = request)
        print(f"Triggered job for {path}. Operation: {cloud_operation.operation.name}")


# Cleans the text file containing the filenames of the log files 
# that were created during the previous submission
def clear_downloadable_file(LOGS_DIR):

    with open(f"{LOGS_DIR}/downloadable_files.txt", "w"):
        pass

# Check upload job status
class JobStatus(APIView):

    permission_classes = [IsAuthenticated]

    def post(self, request):
        job = UploadJob.objects.latest('id')

        return Response({
            "status": job.status
        })


# Uploads data to the database
class DataUpload(APIView):

    permission_classes = [IsAuthenticated]

    # Details for the Swagger UI API documentation
    @extend_schema(

        operation_id="file_upload",
        summary = "Upload and Process Variant Data Files",
        description = "This API allows users to parse .vcf/ .tsv files and upload variant data to the SLGVD. This API can only be accessed by authenticated users. Users need to provide authentication details and obtain a JWT token using which the user then can upload data to the database.",

        request = {

            "properties":{
                "type":"object",
                "properties":{
                    "files":{
                        "type": "array",
                        "items":{"type":"string", "format":"binary"}
                    }
                }
            }
            
        },

        examples = [
            OpenApiExample(
                name = "Example Request",
                description = "Example of file upload request",
                value = {"files": ["file1.tsv", "file2.tsv"]}
            )
        ],

        responses = {

            200: OpenApiResponse(
                response = {
                    "type":"object",
                    "properties":{
                        "message": {"type":"string"},
                        "log":{"type":"string"}
                    }
                },

                examples = [

                    OpenApiExample(
                            name = "Successful upload",
                            value = {
                                "message":"File processing completed",
                                "log": "Submission Id: 123 \nSubmission Date: 2024-12-09 \n\nFile: file1.tsv \nTotal number of variants uploaded: 10 \nTotal number of allele frequencies updated: 5 \nVariants whose allele frequencies were updated: \n['1p935222rCaA', '1p935459rAaG'] \n\nVariant data upload/update completed successfully",
                            }
                        )

                ]
                
            ),

            500: OpenApiResponse(

                response = {
                    "type":"object",
                    "properties":{
                        "message": {"type":"string"},
                        "log":{"type":"string"}
                    }
                },
                
                examples = [
                    
                    OpenApiExample(
                            name = "Error Response",
                            value = {
                                "error": "An error occurred during submission",
                                "log": "An error occurred during submission: <error details> \nVariant data upload/update unsuccessful \n\n"
                            }
                        )

                ]

            )
        },

    )

    def post(self, request, format=None):
        
        user = self.request.user
        files = request.FILES.getlist('files')
        fileType = request.POST.get('file-type')

        print(f"File type: {fileType}")

        # channel_layer = get_channel_layer()

        # async_to_sync(channel_layer.group_send)(
        #     "progress_progress_room",
        #     {
        #         "type":"send_progress",
        #         "progress_data":{"progress":0, "status":"Processing started"}
        #     }
        # )

        # Save file in Cloud storage
        gcs_paths = upload_to_gcs(files)

        # Clear the names of log files in the downloadable_file from previous submission
        # clear_downloadable_file(LOGS_DIR)

        group_name = f"upload_{user.username}"

        # Create a Job ID
        job = UploadJob.objects.create(status = "queued", progress = 0, file_path=gcs_paths, file_type = fileType, username = user.username)

        print(job)

        # result = process_upload.delay(job.id, group_name)
        # print("Queued Celery task id:", result.id)

        try:
            trigger_vcf_parse_job(job.id, gcs_paths)
        except Exception as e:
            print(e)
            job.status = "failed"
            job.save()
            
        return Response({'job_id': job.id, 'group': group_name}, status = status.HTTP_202_ACCEPTED)



# Enable users to download data
class Download(APIView):

    permission_classes = [AllowAny]

    @extend_schema(exclude = True)
    def post(self, request, format = None):

        raw_body = request.body.decode('utf-8')
        data = json.loads(raw_body)
        print(data)

        job_id = json.loads(data['body']).get('job_id')
        print(job_id)
        
        content = json.loads(data['body']).get('data')
        print(data)

        type = json.loads(data['body'])['type']
        print(type)

        db = json.loads(data['body']).get('db')
        print(db)


        # Allow downloading data that are displayed in the results table
        if type == "csv":

            response = HttpResponse(content_type='text/csv')
            response['Content-Disposition'] = 'attachment; filename="results.csv"'

            fieldnames = content[0].keys()
            writer = csv.DictWriter(response,fieldnames=fieldnames)

            writer.writeheader()
            writer.writerows(content)
        
        elif type == 'vcf':
            response = HttpResponse(content_type='text/vcf')
            response['Content-Disposition'] = 'attachment; filename="results.vcf"'

            response.write("##fileformat=VCFv4.2\n")

            if db == 'ssv':
                response.write("##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description=\"Variant Consequence\">\n")
                response.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n")

                response.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

                for row in content:
                    info = row.get('info','.')
                    line = f"{row['chromosome']}\t{row['position']}\t{row['variation_id']}\t{row['ref_allele']}\t{row['alt_allele']}\t.\t.\t{info}\n"
                    response.write(line)

                return response
            
            if db == 'gcnv':
                response.write("##INFO=<ID=END,Number=1,Type=String,Description=\"End Position\">\n")
                response.write("##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description=\"Variant Consequence\">\n")
                response.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Site Frequency\">\n")

                response.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

                for row in content:
                    info = row.get('info','.')
                    line = f"{row['chromosome']}\t{row['start_pos']}\t{row['variation_id']}\tN\t.\t.\t.\t{info}\n"
                    response.write(line)

                return response

        # Enable downloading of the submission log after a successful submission
        # of data to the database
        elif type == 'txt/plain':

            job = UploadJob.objects.get(id = job_id)

            stat = job.status

            if stat != 'done':

                return Response({'message': 'File parsing in progress'}, status = status.HTTP_202_ACCEPTED)
                
            elif stat == 'done': 

                job = UploadJob.objects.get(id = job_id)
                gcs_path = job.log_path[0]
    
                storage_client = storage.Client()
                blob_name = "/".join(gcs_path.replace("gs://","").split("/")[1:])
    
                bucket = storage_client.bucket("slgvd-uploads")
                blob = bucket.blob(blob_name)
    
                response = StreamingHttpResponse(
                    blob.open("rb"),
                    content_type = "text/plain"
                )
            
                # response = HttpResponse(content, content_type='text/plain')
                response['Content-Disposition'] = f'attachment; filename="Submission_Log"'

        return response

# Edits/ Updates fields of variant entries 
class Update(APIView):

    permission_classes = [IsAuthenticated]

    @extend_schema(exclude = True)
    def post(self, request, format = None):

        raw_body = request.body.decode('utf-8')
        data = json.loads(json.loads(raw_body)['body'])['data']
        print(data)

        fields = data[0]
        var_id = data[1]
        db = data[2]

        print(fields)
        print(var_id)
        print(db)

        try:
            with transaction.atomic():
                
                if db == "ssv":
                    Variant.objects.filter(**var_id).update(**fields)
                elif db == "gcnv":
                    CNV.objects.filter(**var_id).update(**fields)

        except Exception as e: 
            print(e)

        return Response({"Message": "Successfully Updated"}, status=status.HTTP_200_OK)

# Removes a specified variant/s / submission (and the entire collection of variants identified by that submission) from the database permanently 
class Remove(APIView):

    permission_classes = [IsAuthenticated]

    @extend_schema(exclude = True)
    def post(self, request, format =None):

        raw_body = request.body.decode('utf-8')
        data = json.loads(json.loads(raw_body)['body'])['data']
        print(data)

        print(data[0])

        if 'submission_id' in data[0].keys():
            sub_id = data[0]['submission_id']
            print(sub_id)
        
        var_id = data[1]
        db = data[2]
        
        print(var_id, db)

        try:

            with transaction.atomic():

                # When only the submission id and no variation id is specified, the backend will remove that specific submission entry and the entire collection of variants identified by that id. 
                if len(var_id) == 0 or var_id['variation_id__in'][0] == '':
                    
                    if Submission.objects.filter(submission_id = sub_id).exists():

                        if db == "ssv":

                            # Update the variant counts in the var_count relation 
                            var_objects = VarSerializer(Variant.objects.filter(submission_id = sub_id), many = True).data
                            UpdateVarCounts(var_objects, "ssv", "remove")

                            # In the event that submission didn't submit that particular variant but only updated it's genotype counts, then the program will decrement the count and remove the submission id from its last_updated field. 
                            def update_ssv_freqs(sub_id):  

                                freq_hom = Frequency.objects.annotate(

                                    json_check = RawSQL(
                                        "JSON_CONTAINS(last_updated->'$.hom', %s)",
                                        (f"[{sub_id}]",)
                                    )

                                ).filter(json_check = True)

                                freq_het = Frequency.objects.annotate(

                                    json_check = RawSQL(
                                        "JSON_CONTAINS(last_updated->'$.het', %s)",
                                        (f"[{sub_id}]",)
                                    )

                                ).filter(json_check = True)

                                print(freq_hom)
                                print(freq_het)

                                if freq_hom:

                                    for freq in freq_hom:
                                        freq.last_updated["hom"].remove(int(sub_id))
                                        freq.homo_count -= 1

                                        freq.save()
                                        
                                if freq_het:

                                    for freq in freq_het:
                                        freq.last_updated["het"].remove(int(sub_id))

                                        freq.het_count -= 1
                                        freq.save()


                            if len(var_objects) == 0:

                                # If no entries were submitted but some existing entries's allele counts were updated then
                                update_ssv_freqs(sub_id)

                                return Response({
                                    'message':f'Allele counts of certain variants were updated. However, no short sequence variants (SSVs) were submitted by the submission id:{sub_id}.'},
                                    status=status.HTTP_404_NOT_FOUND)


                            update_ssv_freqs(sub_id)

                        if db == "gcnv":

                            # Update the variant counts in the var_count relation 
                            var_objects = CNVarSerializer(CNV.objects.filter(submission_id = sub_id), many = True).data
                            UpdateVarCounts(var_objects, "gcnv", "remove")

                            if len(var_objects) == 0:
                                return Response({
                                    'message':f'No germline copy number variants (gCNVs) were submitted in the submission id:{sub_id}'},
                                    status=status.HTTP_404_NOT_FOUND)

                            # In the event that submission didn't submit that particular variant but only updated it's genotype counts, then the program will decrement the count and remove the submission id from its last_updated field.
                            freqs = CNV.objects.annotate(

                                json_check = RawSQL(
                                    "JSON_CONTAINS(last_updated, %s)",
                                    (f"[{sub_id}]",)
                                )

                            ).filter(json_check = True)

                            print(freqs)

                            for freq in freqs:
                                freq.site_count =- 1
                                print(freq.last_updated)
                                freq.last_updated.remove(int(sub_id))

                                freq.save()

                        Submission.objects.filter(submission_id = sub_id).delete()

                    else:

                        return Response({"message":f"No submission by the submission_id: {sub_id} exist."}, status=status.HTTP_200_OK)

                       
                    
                else:
                    # Deleting the variant entry completely from the database which also deletes its entry in the allele_count relation in the case of SSVs

                    if db == "ssv":
                        var_object = VarSerializer(Variant.objects.filter(**var_id), many = True).data
                        UpdateVarCounts(var_object, "ssv", "remove")

                        if len(var_object) == 0:
                            return Response({
                                'message':f'No short sequence variant (SSV) with the {var_id}was found.'},
                                status=status.HTTP_404_NOT_FOUND)

                        Variant.objects.filter(**var_id).delete()

                    elif db == "gcnv":
                        var_object = CNV.objects.filter(**var_id)
                        UpdateVarCounts(var_object, "gcnv", "remove")

                        if len(var_object) == 0:
                            return Response({
                                'message':f'No germline Copy Number Variant (gCNV with the {var_id} was found.'},
                                status=status.HTTP_404_NOT_FOUND)

                        CNV.objects.filter(**var_id).delete()

                return Response({"message":"Deleted Successfully"}, status=status.HTTP_200_OK)

        except Exception as e:
            print(e)







        










        














