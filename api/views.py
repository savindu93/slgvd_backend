from rest_framework import status
from rest_framework.views import APIView
from rest_framework_simplejwt.views import TokenObtainPairView
from rest_framework.response import Response
from .serializer import VarSerializer, FreqSerializer, CustomTokenObtainPairSerializer, SubSerializer, CNVarSerializer, VarCountSerializer
from .models import Variant, Frequency, Submission, CNV, VarCounts
from rest_framework.permissions import AllowAny, IsAuthenticated
from django.db import transaction, connection 
from django.core.exceptions import *
from django.db.models import Q, F
from django.db.models.expressions import  RawSQL
from django.http import HttpResponse
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

        channel_layer = get_channel_layer()

        async_to_sync(channel_layer.group_send)(
            "progress_progress_room",
            {
                "type":"send_progress",
                "progress_data":{"progress":0, "status":"Processing started"}
            }
        )

        try:

            with transaction.atomic():

                # Creating new submission entry
                no_individuals = len(files)
                username = user.username

                Submission.objects.create(no_individuals = no_individuals, username_id = username)
                last_entry = Submission.objects.last()
                sub_id = last_entry.submission_id
                sub_date = last_entry.submission_date
                print(sub_id)

                # Content for the log file summarizing the data uploaded/ changed during sumission
                log = f'Submission Id: {sub_id} \nSubmission Date: {sub_date} \n\n'

                if fileType == "tsv":
                    
                    # Iterate through each file
                    for z, file in enumerate(files):
                        decoded_file = file.read().decode("utf-8").splitlines()
                        reader = csv.reader(decoded_file, delimiter='\t')

                        tot_rows = sum(1 for row in reader)


                        # Resetting file pointer to the beginning
                        file.seek(0)
                        reader = csv.reader(decoded_file, delimiter='\t')
                    
                        # All objects to be added to db
                        variant_objects = []
                        freq_objects = []

                        # Variants whose allele frequencies were updated
                        no_freqs_updated = 0
                        freqs_updated = []

                        # Variants with consequence types that are not registered in the database
                        con_unregistered = []

                        # For throttling messages sent via redis server
                        prog_number = []

                        # Iterate through each line in the file
                        for i, row in enumerate(reader):
                            progress = int((i + 1) / tot_rows * 100)
                            print(row)
        
                            if row[0] == 'PASS':
                                # A single object to be added to db
                                variant_object = {}
                                freq_object = {}

                                # Identify chr
                                pat = re.compile(r"chr(\d{1,2}|X|Y)")
                                pat_match = pat.search(row[1]).group(1)

                                if pat_match not in ["X","Y"]:
                                    chr = int(pat.search(row[1]).group(1))
                                else:
                                    chr = pat_match

                                # Get position
                                pos = int(row[2])

                                # method to shorten alleles of length > 10 bp with run length encoding
                                def encoding_allele(allele):

                                    def shorten_allele(allele):

                                        seq_bases = []
                                        char_count = {}

                                        for i, char in enumerate(allele):

                                            base = list(char_count.keys())
                                            base_count = 0

                                            if len(base) != 0:
                                                base_count = char_count[base[0]]

                                            if char not in base:

                                                if base_count >= 3:
                                                    seq_bases.append(f"{base_count}{base[0]}")

                                                elif base_count <3 and base_count > 0:

                                                    if base_count == 2:
                                                        seq_bases.append(f"{base[0]*2}")
                                                    
                                                    else:
                                                        seq_bases.append(base[0])
                                                
                                                char_count.clear()

                                                char_count[char] = 1

                                                if i == len(allele) - 1:
                                                    seq_bases.append(list(char_count.keys())[0])
                                            
                                            elif char in base:
                                                    char_count[char] += 1

                                                    if i == len(allele) - 1:
                                                        base_count = char_count[base[0]]
                                                        # base = list(char_count.keys())[0]

                                                        if base_count == 2:
                                                            seq_bases.append(f"{char*2}")
                                                        
                                                        if base_count >= 3:
                                                            seq_bases.append(f"{base_count}{char}")
                                                    
                                        seq = "".join(seq_bases)

                                        return seq
                                    
                                    def shorten_dinucleotide(allele):

                                        # print(len(allele))
                                        seq_bases = []
                                        seq = allele

                                        for z in range(2):
                                            # print(z)

                                            if len(seq_bases) > 0:
                                                seq = "".join(seq_bases)
                                                seq_bases.clear() 

                                                if seq != allele:
                                                    # print(f"In break {z}")
                                                    break      

                                            dinucleotide_count = {}

                                            dinucleotides = [seq[i:i+2] for i in range(z, len(seq), 2)]
                                            # print(dinucleotides)

                                            if z == 1:
                                                seq_bases.append(seq[0])

                                            skip = False

                                            for i, char in enumerate(dinucleotides):

                                                dinucleotide = list(dinucleotide_count.keys())
                                                char_count = 0

                                                if len(dinucleotide) != 0:
                                                    char_count = dinucleotide_count[dinucleotide[0]]

                                                if char not in dinucleotide:

                                                    if char_count >= 3:
                                                        seq_bases.append(f"{char_count}({dinucleotide[0]})")

                                                    elif char_count <3 and char_count > 0:

                                                        if char_count == 2:
                                                            seq_bases.append(f"{dinucleotide[0]*2}")
                                                        
                                                        else:
                                                            seq_bases.append(dinucleotide[0])
                                                    
                                                    dinucleotide_count.clear()

                                                    if len(char) == 1:
                                                        seq_bases.append(char)
                                                        continue
                                                    
                                                    if skip == True:
                                                        seq_bases.append(char)
                                                        skip = False
                                                        continue

                                                    if char[0].isdigit() | char[1].isdigit():
                                                        seq_bases.append(char)

                                                        if char[1].isdigit():
                                                            skip = True

                                                        continue

                                                    dinucleotide_count[char] = 1

                                                    if i == len(dinucleotides) - 1:
                                                        seq_bases.append(char)
                                                
                                                elif char in dinucleotide:
                                                    dinucleotide_count[char] += 1

                                                    if i == len(dinucleotides) - 1:
                                                        char_count = dinucleotide_count[dinucleotide[0]]
                                                        # dinucleotide = list(dinucleotide_count.keys())[0]

                                                        if char_count == 2:
                                                            seq_bases.append(f"{char*2}")
                                                        
                                                        if char_count >= 3:
                                                            seq_bases.append(f"{char_count}({char})")

                                            # print(allele)  
                                            # print(seq_bases)
                                            seq = "".join(seq_bases)
                                            # print(seq)
                                            # print(len(seq))

                                        return seq
                                    
                                    s_allele = allele

                                    if len(s_allele) >= 10:
                                        s_allele = shorten_allele(s_allele)
                                        
                                    # print(s_ref)
                                    # print(len(ref))
                                    # print(len(s_ref))

                                    # print(s_alt)
                                    # print(len(alt))
                                    # print(len(s_alt))

                                    s_allele = shorten_dinucleotide(s_allele)

                                    return s_allele
                                
                                # Ref & Alt alleles
                                ref = row[3]
                                s_ref = encoding_allele(ref)

                                alt = row[4]

                                # In some cases where heterozygous genotypes are found, both copies can have different variant alleles. alt_1 is for 2nd alt allele.
                                alt_1 = ""

                                def encoding_allele(allele):

                                    def shorten_allele(allele):

                                        seq_bases = []
                                        char_count = {}

                                        for i, char in enumerate(allele):

                                            base = list(char_count.keys())
                                            base_count = 0

                                            if len(base) != 0:
                                                base_count = char_count[base[0]]

                                            if char not in base:

                                                if base_count >= 3:
                                                    seq_bases.append(f"{base_count}{base[0]}")

                                                elif base_count <3 and base_count > 0:

                                                    if base_count == 2:
                                                        seq_bases.append(f"{base[0]*2}")
                                                    
                                                    else:
                                                        seq_bases.append(base[0])
                                                
                                                char_count.clear()

                                                char_count[char] = 1

                                                if i == len(allele) - 1:
                                                    seq_bases.append(list(char_count.keys())[0])
                                            
                                            elif char in base:
                                                    char_count[char] += 1

                                                    if i == len(allele) - 1:
                                                        base_count = char_count[base[0]]
                                                        # base = list(char_count.keys())[0]

                                                        if base_count == 2:
                                                            seq_bases.append(f"{char*2}")
                                                        
                                                        if base_count >= 3:
                                                            seq_bases.append(f"{base_count}{char}")
                                                    
                                        seq = "".join(seq_bases)

                                        return seq
                                    
                                    def shorten_dinucleotide(allele):

                                        # print(len(allele))
                                        s_alleles = []
                                        seq_bases = []
                                        seq = allele

                                        for z in range(2):

                                            seq = allele    

                                            seq_bases.clear()
                                            dinucleotide_count = {}

                                            dinucleotides = [seq[i:i+2] for i in range(z, len(seq), 2)]
                                            print(dinucleotides)

                                            if z == 1:
                                                seq_bases.append(seq[0])

                                            skip = False

                                            for i, char in enumerate(dinucleotides):

                                                dinucleotide = list(dinucleotide_count.keys())
                                                char_count = 0

                                                if len(dinucleotide) != 0:
                                                    char_count = dinucleotide_count[dinucleotide[0]]

                                                if char not in dinucleotide:

                                                    if char_count >= 3:
                                                        seq_bases.append(f"{char_count}({dinucleotide[0]})")

                                                    elif char_count <3 and char_count > 0:

                                                        if char_count == 2:
                                                            seq_bases.append(f"{dinucleotide[0]*2}")
                                                        
                                                        else:
                                                            seq_bases.append(dinucleotide[0])
                                                    
                                                    dinucleotide_count.clear()

                                                    if len(char) == 1:
                                                        seq_bases.append(char)
                                                        continue
                                                    
                                                    if skip == True:
                                                        seq_bases.append(char)
                                                        skip = False
                                                        continue

                                                    if char[0].isdigit() | char[1].isdigit():
                                                        seq_bases.append(char)

                                                        if char[1].isdigit():
                                                            skip = True

                                                        continue

                                                    dinucleotide_count[char] = 1

                                                    if i == len(dinucleotides) - 1:
                                                        seq_bases.append(char)
                                                
                                                elif char in dinucleotide:
                                                    dinucleotide_count[char] += 1

                                                    if i == len(dinucleotides) - 1:
                                                        char_count = dinucleotide_count[dinucleotide[0]]
                                                        # dinucleotide = list(dinucleotide_count.keys())[0]

                                                        if char_count == 2:
                                                            seq_bases.append(f"{char*2}")
                                                        
                                                        if char_count >= 3:
                                                            seq_bases.append(f"{char_count}({char})")

                                            # print(allele)  
                                            # print(seq_bases)
                                            seq = "".join(seq_bases)
                                            s_alleles.append(seq)
                                            # print(seq)
                                            # print(len(seq))

                                        # print(len(s_alleles[0]))
                                        # print(len(s_alleles[1]))

                                        if len(s_alleles[0]) < len(s_alleles[1]):
                                            seq = s_alleles[0]
                                        else:
                                            seq = s_alleles[1]
                                        
                                        return seq
                                    
                                    s_allele = allele

                                    if len(s_allele) >= 10:
                                        s_allele = shorten_allele(s_allele)
                                        
                                    # print(s_ref)
                                    # print(len(ref))
                                    # print(len(s_ref))

                                    # print(s_alt)
                                    # print(len(alt))
                                    # print(len(s_alt))

                                    s_allele = shorten_dinucleotide(s_allele)

                                    return s_allele
                                
                                if "," in alt:
                                    alt_alleles = alt.split(",")
                                    alt = alt_alleles[0]
                                    s_alt = encoding_allele(alt_alleles[0])

                                    alt_1 = alt_alleles[1]
                                    s_alt_1 = encoding_allele(alt_alleles[1])

                                else:
                                    s_alt = encoding_allele(alt)

                                # Consequence
                                con_values = row[10].split("&")

                                # Dictionary of SO of variant consequences and their custom named consequences
                                consequences = {

                                    'missense_variant': 'Missense',
                                    'intron_variant': 'Intronic',
                                    'upstream_gene_variant': 'Upstream gene',
                                    'synonymous_variant': 'Synonymous',
                                    'downstream_gene_variant': 'Downstream gene',
                                    'intergenic_region': 'Intergenic',
                                    'TF_binding_site_variant': 'TF Binding site',
                                    '5_prime_UTR_variant' : '5 Prime UTR',
                                    '3_prime_UTR_variant' : '3 Prime UTR',
                                    'splice_region_variant':'Splice region',
                                    'non_coding_exon_variant': 'Non-coding exon',
                                    'frameshift_variant' : 'Frameshift', 
                                    'splice_donor_variant':'Splice donor region',
                                    'splice_acceptor_variant':'Splice acceptor', 
                                    'TFBS_ablation' :'TFBS Ablation', 
                                    'exon_loss_variant':'Exon loss',
                                    'stop_lost':'Stop loss', 
                                    'initiator_codon_variant':'Initiator codon', 
                                    'stop_gained':'Stop gain', 
                                    'inframe_insertion':'Inframe insertion',
                                    'start_lost':'Start loss',  
                                    'protein_protein_contact':'Protein binding site', 
                                    'intragenic_variant':'Intragenic',  
                                    'sequence_feature':'Not specified', 
                                    'non_coding_exon_variant':'Non-coding exon', 
                                    'disruptive_inframe_deletion':'Disruptive inframe deletion', 
                                    'inframe_deletion':'Inframe deletion',
                                    'disruptive_inframe_insertion':'Disruptive inframe insertion', 'stop_retained_variant': 'Stop retained', '5_prime_UTR_premature_start_codon_gain_variant' : '5 Prime UTR premature start codon gain'
                                }

                                con = "/ ".join([consequences.get(con_value,"Unknown") for con_value in con_values]) 

                                # Gene name
                                gen = row[12]

                                # Hom/ Het count
                                het = row[67]
                                hom = row[68]

                                if hom == 'true':
                                    homo_count = 1
                                    het_count = 0
                                    last_updated = {"hom":[sub_id]}

                                elif het == 'true':
                                    het_count = 1
                                    homo_count = 0
                                    last_updated = {"het":[sub_id]}
                                
                                # Create new entries or update existing entries
                                def create_update_entries(var_id, variant_object, freq_object, sub_id, progress, no_freqs_updated):

                                    # Check if a variant with an identical id exist
                                    entry = Variant.objects.filter(variation_id__exact = var_id)
                                    
                                    print(entry)

                                    # Check for any unregistered consequences in the con variable
                                    if "Unknown" in con:
                                        con_unregistered.append(var_id)

                                    # If no variant already exist new variant entries are added to the database
                                    if not entry:
                                        variant_objects.append(variant_object)
                                        freq_objects.append(freq_object)
                                    
                                    # If a variant with same ID exist the allele counts for that variant are updated
                                    else:

                                        async_to_sync(channel_layer.group_send)(
                                            "progress_progress_room",
                                            {
                                                "type":"send_progress",
                                                "progress_data":{"progress":progress, "status": f"Updating {var_id} count"}
                                            }
                                        )

                                        no_freqs_updated += 1
                                        freqs_updated.append(var_id)

                                        if hom == 'true':
                                            freq = Frequency.objects.get(variation_id__exact = var_id)

                                            freq.homo_count += 1

                                            if "hom" not in freq.last_updated.keys():
                                                freq.last_updated["hom"] = [sub_id]
                                            else:
                                                freq.last_updated["hom"].append(sub_id)

                                            freq.save()

                                        else:
                                            freq = Frequency.objects.get(variation_id__exact = var_id)

                                            freq.het_count += 1

                                            if "het" not in freq.last_updated.keys():
                                                freq.last_updated["het"] = [sub_id]
                                            else:
                                                freq.last_updated["het"].append(sub_id)
                                            
                                            freq.save()

                                    return no_freqs_updated
                                    
                                # variation_id
                                var_id = f"{chr}p{pos}r{s_ref}a{s_alt}"

                                print(var_id)
                                # Add object field:value pairs to object dict
                                # Variants
                                variant_object['variation_id'] = var_id
                                variant_object['chromosome'] = chr
                                variant_object['position'] = pos
                                variant_object['ref_allele'] = ref
                                variant_object['alt_allele'] = alt
                                variant_object['gene_name'] = gen
                                variant_object['consequence'] = con
                                variant_object['submission_id'] = sub_id

                                # Allele counts
                                freq_object['variation_id'] = var_id
                                freq_object['homo_count'] = homo_count
                                freq_object['het_count'] = het_count
                                freq_object['last_updated'] = last_updated

                                no_freqs_updated = create_update_entries(var_id, variant_object, freq_object, sub_id, progress, no_freqs_updated)

                                # To create/ update variant entries for the 2nd variant allele of the heterozygous genotypes, if present.
                                if alt_1 != "":
                                    variant_object_1 = {}
                                    freq_object_1 = {}

                                    # variation_id
                                    var_id = f"{chr}p{pos}r{s_ref}a{s_alt_1}"
                                    
                                    print(var_id)
                                    # Add object field:value pairs to object dict
                                    # Variants
                                    variant_object_1['variation_id'] = var_id
                                    variant_object_1['chromosome'] = chr
                                    variant_object_1['position'] = pos
                                    variant_object_1['ref_allele'] = ref
                                    variant_object_1['alt_allele'] = alt_1
                                    variant_object_1['gene_name'] = gen
                                    variant_object_1['consequence'] = con
                                    variant_object_1['submission_id'] = sub_id

                                    # Allele counts
                                    freq_object_1['variation_id'] = var_id
                                    freq_object_1['homo_count'] = homo_count
                                    freq_object_1['het_count'] = het_count
                                    freq_object_1['last_updated'] = last_updated

                                    no_freqs_updated = create_update_entries(var_id, variant_object_1, freq_object_1, sub_id, progress, no_freqs_updated)

                            if progress % 5 == 0 and progress not in prog_number:

                                prog_number.append(progress)

                                async_to_sync(channel_layer.group_send)(
                                    "progress_progress_room",
                                    {
                                        "type":"send_progress",
                                        "progress_data":{"progress":progress, "status": f"Processing File {z + 1}"}
                                    }
                                )

                        log += f'File: {file.name} \nTotal number of variants uploaded: {len(variant_objects)} \
                        \nTotal number of allele frequencies updated: {no_freqs_updated} \
                        \nVariants whose allele frequencies were updated: \n{freqs_updated} \n\n \
                        \nVariants with unregistered consequence types: \n{con_unregistered} \n \
                        **Note If variants with unregistered consequence types are present, please kindly \
                            contact administrators and inform them to include them in system. \n\n'

                        # Saving new variant entries to database
                        if len(variant_objects) != 0:

                            async_to_sync(channel_layer.group_send)(
                                "progress_progress_room",
                                {
                                    "type":"send_progress",
                                    "progress_data":{"progress":100, "status": f"Adding data to database"}
                                }
                            )

                            variant_serializer = VarSerializer(data = variant_objects, many=True)

                            
                            if variant_serializer.is_valid():
                                # print("Valid data")
                                print(f"Variant data: \n {variant_serializer.data}")
                                variant_serializer.save()

                                                                
                                freq_serializer = FreqSerializer(data = freq_objects, many = True)
                                print(f"Freq_initial : \n {freq_serializer.initial_data}")

                                if freq_serializer.is_valid():
                                    
                                    print(f"Freq data: \n {freq_serializer.data}")
                                    freq_serializer.save()

                                else:
                                    raise freq_serializer.errors 
                                
                                # Update total SSV variation counts in var_counts relation
                                UpdateVarCounts(variant_objects, "ssv", "upload")

                            
                            else :
                                raise variant_serializer.errors
                    
                if fileType == "gcnv":

                    for z, file in enumerate(files):
                        decoded_file = file.read().decode("utf-8").splitlines()
                        reader = csv.reader(decoded_file, delimiter='\t')

                        tot_rows = sum(1 for row in reader)

                        # Resetting file pointer to the beginning
                        file.seek(0)
                        reader = csv.reader(decoded_file, delimiter='\t')
                    
                        # All objects to be added to db
                        variant_objects = []

                        # CNVs whose site counts & frequencies were updated
                        no_freqs_updated = 0
                        freqs_updated = []

                        # For throttling messages sent via redis server
                        prog_number = []
                        
                        for i, row in enumerate(reader):
                            progress = int((i + 1) / tot_rows * 100)
        
                            if not row[0].startswith("#"):

                                # A single object to be added to db
                                variant_object = {}
                                freq_object = {}
                                # Consequence
                                genotype = int(row[9].split(":")[0])

                                if genotype != 0:

                                    # Identify chr
                                    pat = re.compile(r"chr(\d{1,2}|X|Y)")
                                    pat_match = pat.search(row[0]).group(1)

                                    if pat_match not in ["X","Y"]:
                                        chr = int(pat.search(row[0]).group(1))
                                    else:
                                        chr = pat_match

                                    # Get position
                                    start_pos = int(row[1])

                                    pat = re.compile(r"END=(\d+)")
                                    end_pos = int(pat.search(row[7]).group(1))

                                    if genotype == 1:
                                        con = "Copy loss"
                                        class_ = "DEL"
                                    elif genotype == 2:
                                        con = "Copy gain"
                                        class_ = "DUP"

                                    copy_state = int(row[9].split(":")[1])

                                    last_updated = [sub_id]

                                    # variation_id
                                    var_id = f"{class_}_{chr}p{start_pos}:{end_pos}"

                                    # Add object field:value pairs to object dict
                                    # Variants
                                    variant_object['variation_id'] = var_id
                                    variant_object['chromosome'] = chr
                                    variant_object['start_pos'] = start_pos
                                    variant_object['end_pos'] = end_pos
                                    variant_object['site_count'] = 1
                                    variant_object['consequence'] = con
                                    variant_object['submission_id'] = sub_id
                                    variant_object['last_updated'] = last_updated

                                    # Checking for existing variants with identical variation ids
                                    entry = CNV.objects.filter(variation_id__exact = var_id)
                                    
                                    print(entry)

                                    # If no such entries are found save CNV as new CNV
                                    if not entry:
                                        variant_objects.append(variant_object)
                                    
                                    # Or Update the site count
                                    else:
                                        
                                        print("In update")

                                        async_to_sync(channel_layer.group_send)(
                                            "progress_progress_room",
                                            {
                                                "type":"send_progress",
                                                "progress_data":{"progress":progress, "status": f"Updating {var_id} count"}
                                            }
                                        )

                                        no_freqs_updated += 1
                                        freqs_updated.append(var_id)

                                        freq = CNV.objects.get(variation_id__exact = var_id)

                                        freq.site_count += 1
                                        freq.last_updated.append(sub_id)

                                        freq.save()

                            if progress % 5 == 0 and progress not in prog_number:

                                prog_number.append(progress)

                                async_to_sync(channel_layer.group_send)(
                                    "progress_progress_room",
                                    {
                                        "type":"send_progress",
                                        "progress_data":{"progress":progress, "status": f"Processing File {z + 1}"}
                                    }
                                )


                        log += f'File: {file.name} \nTotal number of variants uploaded: {len(variant_objects)} \nTotal number of allele frequencies updated: {no_freqs_updated} \nVariants whose allele frequencies were updated: \n{freqs_updated} \n\n'

                        # Saving new variant entries to database
                        if len(variant_objects) != 0:

                            async_to_sync(channel_layer.group_send)(
                                "progress_progress_room",
                                {
                                    "type":"send_progress",
                                    "progress_data":{"progress":100, "status": f"Adding data to database"}
                                }
                            )

                            variant_serializer = CNVarSerializer(data = variant_objects, many=True)

                            
                            if variant_serializer.is_valid():
                                print("Valid data")
                                print(f"Variant data: \n {variant_serializer.data}")
                                variant_serializer.save()
                            
                            else :
                                raise variant_serializer.errors
                            
                            # Update total gCNV counts in var_counts relation
                            UpdateVarCounts(variant_objects, "gcnv", "upload")


            async_to_sync(channel_layer.group_send)(
                                "progress_progress_room",
                                {
                                    "type":"send_progress",
                                    "progress_data":{"progress":100, "status": f"File Processing Complete"}
                                }
                            )

            log += f"Variant data upload/update completed successfully"

            return Response({"message":"File processing completed", "log" : log},
                    status=status.HTTP_200_OK)

        except Exception as e:
            log = f'An error occurred during submission: {e} \
            Variant data upload/update unsuccessful '

            return Response({"error":str(e), "log" : log}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

# Enable users to download data
class Download(APIView):

    permission_classes = [AllowAny]

    @extend_schema(exclude = True)
    def post(self, request, format = None):

        raw_body = request.body.decode('utf-8')
        data = json.loads(raw_body)
        print(data)

        content = json.loads(data['body'])['data']
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

            response = HttpResponse(content, content_type='text/plain')
            response['Content-Disposition'] = 'attachment; filename="Submission_Log.txt"'

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







        










        




