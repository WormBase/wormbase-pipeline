#!/usr/bin/env python
# based on
#       * VCF v4.2 (http://samtools.github.io/hts-specs/VCFv4.2.pdf)
#       * AGR JSON 1.0.1.3 (https://docs.google.com/document/d/1yAECtOs1VCEs3mhplJMXqg1akBQJyJbjPnvG2RrE5Aw/edit)

import json
import datetime
import argparse
import operator
import re
import sys

def genotype_string(variation, allStrains, zygosity):
    variationStrains = set(variation["strains"])
    genotypes = []

    for s in allStrains:
        if s in variationStrains:
            if zygosity == 'homozygous':
                genotypes.append('1/1')
            else:
                genotypes.append('1/2')
        else:
            genotypes.append('./.')

    return "\t".join(genotypes)


def get_header_info(gff):
    chr_lengths = {}
    assembly = ''
    with open(gff, 'r') as f:
        line = f.readline()
        while line:
            if line.startswith('#'):
                columns = line.split()
                if line.startswith('##sequence-region'):
                    chr_lengths[columns[1]] = columns[3]
                elif line.startswith('#!assembly'):
                    assembly = columns[-1]
            else:
                return assembly, chr_lengths
            line = f.readline()
    return assembly, chr_length


def get_strains(variations):
    strains = set()
    for variation in variations:
        for strain in variation["strains"]:
            strains.add(strain)
    return tuple(sorted(strains))


chromosomes = ('I', 'II', 'III', 'IV', 'V', 'X', 'MtDNA')

chrom2ncbi = {
    'I': 'RefSeq:NC_003279.8',
    'II': 'RefSeq:NC_003280.10',
    'III': 'RefSeq:NC_003281.10',
    'IV': 'RefSeq:NC_003282.8',
    'V': 'RefSeq:NC_003283.11',
    'X': 'RefSeq:NC_003284.9',
    'MtDNA': 'RefSeq:NC_001328.1',
}

expand_iupac = {
    'R': 'A,G',
    'Y': 'C,T',
    'S': 'C,G',
    'W': 'A,T',
    'K': 'G,T',
    'M': 'A,C',
}

parser = argparse.ArgumentParser()
parser.add_argument("-j", "--json", help="JSON input file")
parser.add_argument("-g", "--gff", help="Corresponding GFF file")
parser.add_argument("-o", "--out", help="Output VCF file")

args = parser.parse_args()
assembly, chr_lengths = get_header_info(args.gff)

vcf_file = open(args.out, 'w')

vcf_file.write("##fileformat=VCFv4.2\n" +
               datetime.datetime.today().strftime("##fileDate=%Y%m%d") + "\n" +
               "##reference=" + assembly + "\n" +
               "##source=AllianceJSON\n")

for chr in chromosomes:
    vcf_file.write("##contig=<ID=" + chr + ",length=" + chr_lengths[chr] + ">\n")

vcf_file.write("##FORMAT<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

headers = ['#CHROM', 'POS', 'ID', 'REF',
           'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

nt_regex = re.compile('^[ACGT]$')

with open(args.json, 'r') as read_file:
    parsed = json.load(read_file)

    # get all strains for column headers
    strains = get_strains(parsed["data"])
    for s in strains:
        headers.append('WB:' + s)  # need curie form of strain
    vcf_file.write("\t".join(headers) + "\n")

    vcf_lines = []
    for v in (parsed["data"]):
        vcf_data = {}

        refSeq = "" if v["genomicReferenceSequence"] == "N/A" else v["genomicReferenceSequence"]
        varSeq = "" if v["genomicVariantSequence"] == "N/A" else v["genomicVariantSequence"]
        pos = int(v["start"])
        zygosity = 'homozygous'
        
        if 'paddedBase' in v.keys():
            if pos == 1:
                refSeq = refSeq+v["paddedBase"]
                varSeq = varSeq+v["paddedBase"]
            else:
                refSeq = v["paddedBase"]+refSeq
                varSeq = v["paddedBase"]+varSeq
                pos = pos-1  # include the padding base in POS
        vcf_data["chromosome"] = v["chromosome"]
        vcf_data["pos"] = pos

        if len(varSeq) == 1:
            if nt_regex.match(varSeq) is None:
                if varSeq in expand_iupac.keys():
                    varSeq = expand_iupac[varSeq]
                    zygosity = 'heterozygous'
                else:
                    print("Unrecognised alternative allele " + varSeq + " for " + v["alleleId"], file=sys.stderr)
                    next
        
        gtString = genotype_string(v, strains, zygosity)
        
        vcf_data["line"] = "\t".join([v["chromosome"], str(
            pos), v["alleleId"], refSeq, varSeq, '.', 'PASS', '.', 'GT', gtString])

        vcf_lines.append(vcf_data)

    for v in sorted(vcf_lines, key=operator.itemgetter('chromosome', 'pos')):
        vcf_file.write(v["line"] + "\n")

    vcf_file.close()
