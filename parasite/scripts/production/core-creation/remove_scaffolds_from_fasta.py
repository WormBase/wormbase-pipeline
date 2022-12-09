from Bio import SeqIO
from ProductionUtils import *
import sys

scaffolds_tsv = sys.argv[1]
input_fasta = sys.argv[2]

blacklisted_scaffolds = csvlines2list(scaffolds_tsv)

ffile = SeqIO.parse(sys.argv[2], "fasta")
header_set = set(line.strip() for line in flatten(blacklisted_scaffolds))

for seq_record in ffile:
    try:
        header_set.remove(seq_record.name)
    except KeyError:
        print(seq_record.format("fasta"))
        continue
if len(header_set) != 0:
    print(len(header_set),'of the headers from list were not identified in the input fasta file.', file=sys.stderr)