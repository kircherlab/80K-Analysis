import pysam
import pandas
import sys
import yaml

config_path = "/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/config/filter_bc_from_bam_config.yml"
with open(config_path) as conf:
    config = yaml.load(conf, Loader=yaml.FullLoader)
    conf.close()

# TODO: add options for calling from the command line
# TODO: add verbose option to only give warnings with verbose?

# helpful functions
def determine_MD_NM(differences,sequence,alnLength):
  NM = 0
  MD = ''
  lpos = -1
  for refPos,seqPos,edit in differences:
    if type(edit[1]) != type('T'): 
      NM+=edit[1]
      if edit[0] == 'INS': pass
      else:
        if lpos+1 <= refPos:
          MD += "%d^%s"%(refPos-(lpos+1),edit[2])
          lpos = refPos-1+edit[1]
        else:
          MD += "^%s"%(edit[2])
          lpos += edit[1]
    else:
      if lpos+1 <= refPos:
        MD += "%d%s"%(refPos-(lpos+1),edit[1])
      else:
        MD += "%s"%(edit[1])
      lpos = refPos+len(edit[1])-1
      NM += len(edit[1])
  if lpos+1 <= alnLength:
    MD += "%d"%(alnLength-(lpos+1))
  return MD,NM


def aln_length(cigarlist):
  tlength = 0
  for operation,length in cigarlist:
    if operation == 0 or operation == 2 or operation == 3 or operation >= 6: tlength += length
  return tlength


def expected_length_filter(read, expected_length):
  """Filter reads that are of an expected length (AS tag is used because match gives 1 => if AS >= expected_length, then read is of expected length and quality)"""
  if read.get_tag('AS') >= expected_length:
    return True
  return False


def calculate_sequence_identity(read):
  """Compute sequence identity form a read using cigar string, MD tag, and NM tag"""
  try: 
    cigar = read.cigarstring
    MD = read.get_tag('MD')
    NM = read.get_tag('NM')
    alnLength = aln_length(read.cigartuples)
  except KeyError:
    sys.stderr.write("Error: No MD or NM tag found")
    return -1
  return (alnLength - NM)/alnLength, NM


def get_number_of_matches_per_strand(read):
  """Get number of matches per forward and reverse strand assuming read has XA tag"""
  # split ";" on right side of XA tag
  xa_tag = read.get_tag('XA').rstrip(';')
  xa_list = xa_tag.split(';')
  forward_matches = 0
  reverse_matches = 0
  for xa in xa_list:
    if '+' in xa.split(",")[-3]: # because oligo name might include ";"
      forward_matches += 1
    elif '-' in xa.split(",")[-3]:
      reverse_matches += 1
    else:
      sys.stderr.write("Error: XA tag does not contain strand information")
      sys.exit()
  return forward_matches, reverse_matches


def get_XA_information(read):
  """Select the XA information from the alternative alignment on the positive strand"""
  try:
    prepare_xa = read.get_tag("XA").split(",") # oligo name might include ";"
    reference_name = prepare_xa[0]
    next_elem = prepare_xa[1:] # remove first name
  except KeyError:
    sys.stderr.write("Error: No XA tag found")
    return -1
  while len(next_elem) > 0:
      next_elem = ",".join(next_elem) # join back
      next_elem = next_elem.split(";") # join back, split on ; in order to find one element
      elem = next_elem[0]
      if "+" in elem.split(",")[-3]:
          return f"{reference_name},{elem}"
      # join again with one less
      next_elem = ";".join(next_elem[1:])
      prepare_xa = next_elem.split(",")
      reference_name = prepare_xa[0]
      next_elem = prepare_xa[1:] # remove first name
  sys.stderr.write("Error: No forward alignment found in XA tag")
  return -1


def print_read_information(read):
  """Debugging and showing read information to the user"""
  print("Reference: ", read.reference_name)
  print("Cigar info: ", read.cigarstring, read.cigartuples)
  sequence_identity, num_mismatches = calculate_sequence_identity(read)
  print("Identity and mismatch", sequence_identity, num_mismatches)
  print("flag: ", read.flag)
  print("reference start: ", read.reference_start)
  print("mapping quality: ", read.mapping_quality)
  print("Tag information", read.tags) 


def get_barcode(read):
  """The barcode is stored in XI tag:XI:Z:<barcode>,YI:I:<unknown>"""
  # check if it does not include "N" 
  barcode = read.get_tag("XI").split(",")[0]
  if "N" in barcode.upper():
      sys.stderr.write("Error: Barcode of {read.query_name}")
      return "failed"
  return barcode


def prepare_table_information(read, case="normal"):
  """Prepares the information of the association:
    @case: 
      - normal: the information from the alignment are taken
      - fix_mapping_quality: the mapping quality is set to 1 but all information is taken from the alignment
      - rescue: we rescue a better alignment from the XA tag
  """
  barcode = get_barcode(read)
  reference_name = read.reference_name
  position=read.reference_start # 0-based
  cigarstring=read.cigarstring
  nm=read.get_tag("NM")
  md=read.get_tag("MD")
  if case == "normal":
    mapping_quality = read.mapping_quality
  
  if case == "fix_mapping_quality":
    mapping_quality = 1
    
  if case == "rescue":
    # if NM is 0, then we can rescue the alignment from XA tag
    xa_info = get_XA_information(read)
    reference_name = xa_info.split(",")[0]
    position = int(xa_info.split(",")[1])
    cigarstring = xa_info.split(",")[2]
    nm = int(xa_info.split(",")[3])
    md = read.get_tag("MD")
    mapping_quality = 1
    if nm != 0:
      md = "unknown" # if NM != 0: we don't know MD tag
  return f"{barcode}\t{reference_name}\t{position};{cigarstring};NM:i:{nm};MD:Z:{md};{mapping_quality}"





# use split as test data 
identity_threshold = config["general"]["identity_threshold"]
mismatches_threshold = config["general"]["mismatches_threshold"]
use_expected_alignment_length = config["general"]["use_expected_alignment_length"]
expected_alignment_length = config["general"]["expected_alignment_length"]
# code from /data/gpfs-1/users/kisa11_c/work/coding/MPRA/bin/removeSequenceErrors.py
bamfile = config["files"]["merged_bam"]

input_file = pysam.Samfile( bamfile, "rb" )
# open output_bamfile for writing
# output_file = pysam.Samfile( output_bamfile, "wb", template=input_file )
count = 0
high_mismatches_count = 0
high_identity_count = 0
unmapped_count = 0
reversed_read_count = 0
rescued_read_count = 0
not_rescuable_count = 0
count_no_second_alignment = 0
count_low_identity_but_high_score = 0
no_xi_count = 0
change_quality_count = 0
high_quali_reversed_read_count = 0
high_quality_same_as_xs = 0
high_quality_alignment = 0
low_quality_alignment = 0
show_example = True
print("start reading bam file ...")
with open(config["files"]["assignment_tbl"], "w") as output_file:
  for read in input_file:
    count += 1
    # if count > 10: sys.exit() # debug
    # # check if NOS3 in reference name (debug)
    # if "NOS3" not in read.reference_name: continue
    #### Counting of specific read properties
    # check if all alignments have XI tag
    try:
      # print(read.get_tag("XI"))
      xi = read.get_tag("XI")
    except:
      no_xi_count += 1 # 0 xi count in sample data

    # skip reads with certain properties
    if read.is_unmapped:
      unmapped_count += 1
      continue # because of weird line (NB501960:812:HH53WAFX5:1:11101:13917:2374	4	*	0	0	None	*	0	0	GGTG)
    # check if read has low identity but high alignment score
    sequence_identity, num_mismatches = calculate_sequence_identity(read)
    if sequence_identity < identity_threshold:
      if read.get_tag('AS') > 265:
        count_low_identity_but_high_score += 1
    
    
    if not sequence_identity >= identity_threshold: continue
    # high alignment identity => potentially interesting alignments
    high_identity_count += 1

    if num_mismatches > mismatches_threshold: # throw warning
      high_mismatches_count += 1
      # print("WARNING: number of missmatches from %s %d"%(read.query_name, num_mismatches))
    
    # filter for expected sequence length
    if use_expected_alignment_length:
      if not expected_length_filter(read,expected_alignment_length): continue
    
    
    # modify mapping quality for reads on forward with AS > XS
    if read.mapping_quality < 1:
      low_quality_alignment += 1
      if read.flag == 0: # check if it is a best alignment
        if read.get_tag("AS") > read.get_tag("XS"):
          change_quality_count += 1
          output_file.write(prepare_table_information(read, case="fix_mapping_quality") + "\n")
          

      # rescue reads with reversed read (flag == 16) with AS = XS and check if only on other alignment on forward strand is given
      if read.flag == 16:
        reversed_read_count += 1
        if read.get_tag("AS") == read.get_tag("XS"):
          try: 
            read.get_tag("XA")
          except:
            # sys.stderr.write("WARNING: read (%s) can not be rescued because XA tag is not given\n"%(read.query_name))
            count_no_second_alignment += 1
            continue
          # check if only one alignment on forward strand is given
          forward_matches, reverse_matches = get_number_of_matches_per_strand(read)
          if forward_matches == 1: # best alignment found on forward strand
            # change alignment to forward strand with the information from XA tag
            output_file.write(prepare_table_information(read, case="rescue") + "\n")
            rescued_read_count += 1
          elif forward_matches > 1:
            not_rescuable_count += 1
            #Throw warning that the read can not be rescued
            # sys.stderr.write("WARNING: read (%s) can not be rescued because more than one alignment on forward strand is given\n"%(read.query_name))
    
    else: # mapping quality >= 1
      # normal case
      output_file.write(prepare_table_information(read, case="normal") + "\n")
      high_quality_alignment += 1
      # check if reversed read has AS = XS
      if read.flag == 16:
        high_quali_reversed_read_count += 1
        if read.get_tag("AS") == read.get_tag("XS"):
          # if get_tag("XA") is not None
          high_quality_same_as_xs += 1
          try:
            read.get_tag("XA")
          except:
            # sys.stderr.write("WARNING: read (%s) can not be rescued because XA tag is not given"%(read.query_name))
            count_no_second_alignment += 1
            continue  
        
input_file.close()
output_file.close()
print("unmapped count: ", unmapped_count)
print("reversed_read_count: ", reversed_read_count)
print("reversed rescued_read_count: ", rescued_read_count)
print("It is a best alignment: ", change_quality_count)
print("All rescued reads: ", change_quality_count + rescued_read_count)
print("not_rescuable_count: ", not_rescuable_count)
print("high_identity_count: ", high_identity_count)
print("high_mismatches_count: ", high_mismatches_count)
print("count: ", count)
print("proportion of high quality: ", high_identity_count/count)
print("high mapping quality bwa: ", high_quality_alignment)
print("low mapping quality bwa: ", low_quality_alignment)
print("count_no_second_alignment: ", count_no_second_alignment)
print("count_low_identity_but_high_score: ", count_low_identity_but_high_score)
print("No xi count: ", no_xi_count)
print("high_quali_reversed_read_count: ", high_quali_reversed_read_count)
print("high_quality_same_as_xs: ", high_quality_same_as_xs)