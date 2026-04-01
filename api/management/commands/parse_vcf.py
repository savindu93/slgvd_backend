from django.core.management.base import BaseCommand
from api.models import UploadJob, Variant, Frequency, Submission, CNV, VarCounts
from api.serializer import VarSerializer, FreqSerializer, CustomTokenObtainPairSerializer, SubSerializer, CNVarSerializer, VarCountSerializer
from api.views import UpdateVarCounts
from google.cloud import storage
from channels.layers import get_channel_layer
from asgiref.sync import async_to_sync
from django.http import HttpResponse
from rest_framework import status
from rest_framework.response import Response
from django.db import transaction, connection 
from django.core.exceptions import *
import requests
import time
import re, csv, json
import os
import io
import uuid

channel_layer = get_channel_layer()

def send_progress(progress, status):
  async_to_sync(channel_layer.group_send)(
          "progress_progress_room",
          {
              "type":"send_progress",
              "progress_data":{"progress":progress, "status":status}
          }
      )

def upload_txt_to_gcs(content, filename, bucket_name = "slgvd-uploads"):

  storage_client = storage.Client()
  bucket = storage_client.bucket(bucket_name)

  unique_filename = f"{uuid.uuid4()}_{filename}"
  blob = bucket.blob(unique_filename)

  # Upload text to file
  blob.upload_from_string(content, content_type = "text/plain")

  return f"gs://{bucket_name}/{unique_filename}"
    
class Command(BaseCommand):

  def add_arguments(self, parser):
    parser.add_argument('--id', type = int)
    parser.add_argument('--path', type = str)


  def handle(self, *args, **options):

    job_id = options['id']
    gcs_path = options['path']

    job = UploadJob.objects.get(id = job_id)

    # Get UploadJobs entry relevant to the job_id and update status
    job.status = "processing"
    job.save(update_fields=['status'])
    send_progress(0, "Processing started")

    fileType = job.file_type

    print(job)
    print(gcs_path)
    
    try:

      # Download from GCS
      storage_client = storage.Client()
      bucket_name = "slgvd-uploads"
      blob_name = '/'.join(gcs_path.split('/')[3:])

      print(blob_name)
      
      bucket = storage_client.bucket(bucket_name)
      blob = bucket.blob(blob_name)
      

      # Parse VCF file and save variants in Cloud SQL
      print("Start parse")
  
      print(f"Got the job: {job}")
  
      with transaction.atomic():

          # Creating new submission entry
          # no_individuals = len(files)
          no_individuals = 1
          username = job.username

          Submission.objects.create(no_individuals = no_individuals, username_id = username)
          last_entry = Submission.objects.last()
          sub_id = last_entry.submission_id
          sub_date = last_entry.submission_date
          print(sub_id)

          # Content for the log file summarizing the data uploaded/ changed during sumission
          log = f'Submission Id: {sub_id} \nSubmission Date: {sub_date} \n\n'

          blob_reader = blob.open("r", encoding = "utf-8")
          tot_rows = sum(1 for line in blob_reader if not line.startswith('#'))
          blob_reader.close()
        
          if fileType == "tsv":

                            
              # Iterate through each file
              # for z, file in enumerate(files):
              with blob.open("r", encoding="utf-8") as file:
                
                  # decoded_file = file.read().decode("utf-8").splitlines()
                  decoded_lines = (line for line in file if not line.startswith('#'))
                  reader = csv.reader(decoded_lines, delimiter='\t')

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

                  # Measure time for sending progress at intervals
                  last_time = time.time()

                  # Iterate through each line in the file
                  for i, row in enumerate(reader):

                      now = time.time()
                      
                      progress = int((i + 1) / tot_rows * 100)
                      # print(row)
  
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
                          def create_update_entries(var_id, variant_object, freq_object, sub_id, progress, no_freqs_updated, last_time):
                            
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
                                  print(now - last_time)
                                  
                                  # if (now - last_time) >= 20:
                                  #     send_progress(progress, f"Updating {var_id} count")
                                  #     last_time = time.time()

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

                          no_freqs_updated = create_update_entries(var_id, variant_object, freq_object, sub_id, progress, no_freqs_updated, last_time)

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

                              no_freqs_updated = create_update_entries(var_id, variant_object_1, freq_object_1, sub_id, progress, no_freqs_updated, last_time)

                      if progress % 5 == 0 and progress not in prog_number and now - last_time >= 20:

                          prog_number.append(progress)
                          send_progress(progress, f"Processing File")
                          last_time = time.time()


                  log += f'File: {blob_name} \nTotal number of variants uploaded: {len(variant_objects)} \
                  \nTotal number of allele frequencies updated: {no_freqs_updated} \
                  \nVariants whose allele frequencies were updated: \n{freqs_updated} \n\n \
                  \nVariants with unregistered consequence types: \n{con_unregistered} \n \
                  **Note If variants with unregistered consequence types are present, please kindly \
                      contact administrators and inform them to include them in system. \n\n'

                  # Saving new variant entries to database
                  if len(variant_objects) != 0:
                      
                      send_progress(100, "Adding data to database")

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

              #for z, file in enumerate(files):
              with open(filePath, "r", encoding="utf-8") as file:
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

                                  send_progress(progress, f"Updating {var_id} count")

                                  no_freqs_updated += 1
                                  freqs_updated.append(var_id)

                                  freq = CNV.objects.get(variation_id__exact = var_id)

                                  freq.site_count += 1
                                  freq.last_updated.append(sub_id)

                                  freq.save()

                      if progress % 5 == 0 and progress not in prog_number:

                          prog_number.append(progress)

                          send_progress(progress, f"Processing File {z + 1}")

                  

                  log += f'File: {blob_name} \nTotal number of variants uploaded: {len(variant_objects)} \nTotal number of allele frequencies updated: {no_freqs_updated} \nVariants whose allele frequencies were updated: \n{freqs_updated} \n\n'

                  # Saving new variant entries to database
                  if len(variant_objects) != 0:

                      send_progress(100, "Adding data to database")

                      variant_serializer = CNVarSerializer(data = variant_objects, many=True)

                      
                      if variant_serializer.is_valid():
                          print("Valid data")
                          print(f"Variant data: \n {variant_serializer.data}")
                          variant_serializer.save()
                      
                      else :
                          raise variant_serializer.errors
                      
                      # Update total gCNV counts in var_counts relation
                      UpdateVarCounts(variant_objects, "gcnv", "upload")

      log += f"Variant data upload/update completed successfully"
    
      # Saving log file as a text file in the GCS bucket
      log_path = upload_txt_to_gcs(log, f"log_{blob_name}")

      print(log_path)
      
      job.status = "done"
      job.progress = 100
      job.log.append(log_path)
      job.save(update_fields = ['progress', 'log_path', 'status'])

      send_progress(job.progress, f'File Processing Complete')

    except Exception as e:
    
      job.status = "failed"
      job.error = str(e)
      job.save(update_fields = ['status', 'error'])
    
      send_progress(job.progress, f'Failed (e)')
    
      log = f'An error occurred during submission: {e} \
      Variant data upload/update unsuccessful '
    
      # with open(f"{LOGS_DIR}/downloadable_files.txt", "a") as f:
      #         f.write(f"{job.id}-{sub_id}.log\n")




