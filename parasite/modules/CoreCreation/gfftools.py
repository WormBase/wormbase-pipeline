import re
import sys
import pprint
import copy
import itertools

##

def parse_gff(input_gff, source_to_WB = False, seq_region_synonyms = False, has_parentage = True):

  if seq_region_synonyms:
    synonyms = _parse_synonyms_file(seq_region_synonyms)    

  gff_fh = open(input_gff)

  genes       = {}  
  transcripts = {}
  exons       = {}
  CDS         = {}

  for line in gff_fh:
    if re.match("#", line):
      continue
    line = line.rstrip()
    fields = line.split("\t")
    if len(fields) != 9:
      continue
    this_feature = {
      'scaffold' : fields[0],
      'source'   : fields[1],
      'type'     : fields[2],
      'start'    : fields[3],
      'end'      : fields[4],
      'score'    : fields[5],
      'strand'   : fields[6],
      'phase'    : fields[7],
    }

    if seq_region_synonyms:
      try:
        this_feature["scaffold"] = synonyms[fields[0]]    
      except KeyError:
        sys.exit("No synonym found for "+fields[0])

    attributes = fields[8].split(";")
    for attribute_pair in attributes:
      if re.match("^$", attribute_pair): 
        continue
      try:
        attribute, value = attribute_pair.split("=")
        this_feature[attribute] = value
      except:
        sys.exit("Couldn't parse col9 attributes:\n"+line)

    if this_feature["type"] == "three_prime_UTR" or this_feature["type"] == "five_prime_UTR" or this_feature["type"] == "region":
      continue

    if source_to_WB is True:
      this_feature["source"] = "WormBase_imported"

    if this_feature["type"] == 'gene' or this_feature["type"] == 'pseudogene':
      genes = _add_to_genes(genes, this_feature)

    elif this_feature["type"] == 'mRNA' or this_feature["type"] == 'tRNA' or this_feature["type"] == 'rRNA':
      transcripts = _add_to_transcripts(transcripts, this_feature)
      if has_parentage is True:
        genes = _add_children_to_genes(genes, this_feature)

    elif this_feature["type"] == 'exon':
      exons = _add_to_exons(exons, this_feature)

    elif this_feature["type"] == "CDS":
      CDS = _add_to_cds(CDS, this_feature, has_parentage = has_parentage)

    else: 
      sys.exit("Unexpected feature type\n:"+line)

  transcripts = _fix_pseudogenes(genes, transcripts)

  return genes, transcripts, exons, CDS

##

def extract_gene_names(gff):
  fh = open(gff)
  names = []
  for line in fh:
    if re.match("#", line):
      continue
    fields = line.split("\t")
    if len(fields) != 9:
      continue
    if fields[2] == 'gene':
      a = re.search("Name=([^;]+)",fields[8])
      if a:
        names.append(a.group(1))
      else:
        sys.exit("Couldn't extract gene name:\n"+line)
  return(names)


##

def generate_ids_from_parent_ids(features, suffix):
  for mrna, list_of_cds in features.items():
    for cds in list_of_cds:
      new_id = mrna+suffix
      cds["ID"] = new_id
  return features

##

def id_to_name(features):
  for feature in features:
    features[feature]["Name"] = features[feature]["ID"]
  return features

##

def locus_tag_to_name(features):
  for feature in features:
    features[feature]["Name"] = features[feature]["locus_tag"]
  return features

##

def transcript_names_from_gene_names(genes, transcripts, format = "e"):
  for gene in genes:
    n = 1
    for child in genes[gene]["children"]:
      if format == "e":
        transcripts[child]["Name"] = genes[gene]["Name"]+".t"+str(n)
      elif format == "apollo":
        padding = 5 - len(str(n))
        transcripts[child]["Name"] = genes[gene]["Name"]+ "-" + "".join( itertools.repeat(str("0"), padding) ) + str(n)
      else:
        sys.exit("Format "+ format + "not recognised")
      n = n + 1
  return transcripts

##

def transcripts_from_cds(transcripts, CDS):
  for ID in CDS:
    if len(CDS[ID]) == 1:
      transcripts[ID] = copy.deepcopy(CDS[ID][0])
      transcripts[ID]["type"] = "mRNA"
    else:
      sys.exit("Trying to make an mRNA from multiple CDS segments- you need to write some more code to do that")
  return transcripts

##

def genes_from_transcripts(transcripts):
  genes = copy.deepcopy(transcripts)
  for gene in genes:
    genes[gene]["type"] = "gene"
    genes[gene]["children"] = []
    genes[gene]["children"].append(genes[gene]["ID"])
    genes[gene]["ID"] = genes[gene]["ID"]+".gene"
    transcripts[gene]["Parent"] = genes[gene]["ID"]
  return(genes, transcripts)

##

def exons_from_transcripts(transcripts):
  exons = {}
  for transcript in transcripts:
    this_exon = copy.deepcopy(transcripts[transcript])
    this_exon["Parent"] = transcript
    this_exon["type"] = "exon"
    if transcript not in exons:
      exons[transcript] = []
    exons[transcript].append(this_exon)
  return(exons)

##

def _parse_synonyms_file(seq_region_synonyms):
  seq_fh = open(seq_region_synonyms, "r")
  scaffold_synonyms = {}
  for line in seq_fh:
    line = line.rstrip()
    fields = line.split("\t")
    if len(fields) != 4:
      sys.exit("Seq region synonyms file should have 4 columns\n:"+line)
    scaffold_synonyms[fields[2]] = fields[1]
  return(scaffold_synonyms)

##

def print_gff(out_file, genes, transcripts, exons, CDS):
  fh = open(out_file, "w")
  print("##gff-version 3", file = fh)
  for gene in genes:
    attributes = "ID="+genes[gene]["ID"]+";Name="+genes[gene]["Name"]
    _print_line(genes[gene], "gene", attributes, fh)
    for child in genes[gene]["children"]:
      attributes = "ID="+transcripts[child]["ID"]+";Name="+transcripts[child]["Name"]+";Parent="+transcripts[child]["Parent"]
      _print_line(transcripts[child], transcripts[child]["type"], attributes, fh)
      for exon in exons[child]:
        attributes = "Parent="+child
        _print_line(exon, "exon", attributes, fh)   
      if transcripts[child]["type"] == "mRNA":
        for this_CDS in CDS[child]:
          attributes = "ID="+this_CDS["ID"]+";Parent="+child
          _print_line(this_CDS, "CDS", attributes, fh)

##

def _print_line(data, feature_type, attributes, out_fh):
  print(data["scaffold"], data["source"], feature_type, data["start"], data["end"], data["score"], data["strand"], data["phase"], attributes, sep='\t' ,file=out_fh)

##

def _add_to_genes(genes, feature):
  if "ID" in feature:
    gene_id = feature["ID"]
  elif "Name" in feature:
    gene_id = feature["Name"]
  else:
    sys.exit("Gene has neither an ID nor a Name attribute:\n"+line)

  try: 
    genes[gene_id].update(feature)
  except KeyError:
    genes[gene_id] = feature
  
  return genes

##

def _add_to_transcripts(transcripts, feature):
  try:
    transcripts[feature["ID"]] = feature
  except KeyError:
    sys.exit("transcript doesn't have an ID attribute:\n"+line)
  return transcripts

##

##

def _add_children_to_genes(genes, feature):
  if "ID" not in feature:
    sys.exit("mRNA doesn't have an ID attribute:\n"+line)
  if "Parent" not in feature:
    sys.exit("mRNA doesn't have a parent attribute:\n"+line)
  if feature["Parent"] not in genes:
    genes[feature["Parent"]] = {}
  if "children" not in genes[feature["Parent"]]:
    genes[feature["Parent"]]["children"] = [] 
  genes[feature["Parent"]]["children"].append(feature["ID"])

  return genes

##

def _add_to_exons(exons, feature):
  if "Parent" not in feature:
    sys.exit("exon doesn't have a Parent attribute:\n"+line)
  if feature["Parent"] not in exons:
    exons[feature["Parent"]] = []
  exons[feature["Parent"]].append(feature)
  return exons

##

def _add_to_cds(CDS, feature, has_parentage=True):
  if has_parentage is True:
    if "Parent" not in feature:
      sys.exit("CDS doesn't have a Parent attribute:\n"+line)
    if feature["Parent"] not in CDS:
      CDS[feature["Parent"]] = []  
    CDS[feature["Parent"]].append(feature)
  else:
    if feature["ID"] not in CDS:
      CDS[feature["ID"]] = []
    CDS[feature["ID"]].append(feature)
  return CDS

##

def _fix_pseudogenes(genes, transcripts):
  for gene in genes:
    if genes[gene]["type"] == "pseudogene":
      for child in genes[gene]["children"]:
        transcripts[child]["type"] == "pseudogenic_transcript" 
  return transcripts  

##



