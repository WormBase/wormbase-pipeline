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
import subprocess

def genotype_string(variation, allStrains):
    variationStrains = set(variation["strains"])
    genotypes = []

    for s in allStrains:
        if s in variationStrains:
            genotypes.append('1/1')
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


def get_refseq_from_fasta(v, fasta):
    faidx_call = subprocess.run(["samtools", "faidx", fasta, str(v["chromosome"]) + ':' + str(v["start"]) + '-' + str(v["end"])], stdout=subprocess.PIPE, text=True)
    faidx_lines = faidx_call.stdout.split("\n");
    faidx_lines.pop(0)
    return ''.join(faidx_lines).upper()

def get_padbase_from_fasta(v, refSeq, fasta):
    
    if refSeq == '':
        pos = int(v["start"])
    else:
        if v["start"] == 1:
            pos = int(v["end"]) + 1
        else:
            pos = int(v["start"]) - 1
    faidx_call = subprocess.run(["samtools", "faidx", fasta, str(v["chromosome"]) + ':' + str(pos) + '-' + str(pos)], stdout=subprocess.PIPE, text=True)
    faidx_lines = faidx_call.stdout.split("\n");
    faidx_lines.pop(0)
    return faidx_lines[0].upper()

def get_strains(variations):
    strains = set()
    for variation in variations:
        for strain in variation["strains"]:
            strains.add(strain)
    return tuple(sorted(strains))


chromosomes = ('I', 'II', 'III', 'IV', 'V', 'X', 'MtDNA')

chrom2ncbi = {
    'FB': {
	'2L': 'NT_033779.5',
	'2R': 'NT_033778.4',
	'3L': 'NT_037436.4',
	'3R': 'NT_033777.3',
	'4': 'NC_004353.4',
	'X': 'NC_004354.4',
	'Y': 'NC_024512.1',
	'mitochondrion_genome': 'NC_024511.2',
	'Unmapped_Scaffold_8_D1580_D1567': 'NW_007931083.1',
	'211000022278279': 'NW_007931104.1',
	'211000022278436': 'NW_001845431.1',
	'211000022278449': 'NW_001845819.1',
	'211000022278760': 'NW_001846712.1',
	'211000022279165': 'NW_001846812.1',
	'211000022279188': 'NW_001845284.1',
	'211000022279264': 'NW_001847227.1',
	'211000022279392': 'NW_001846198.1',
	'211000022279681': 'NW_001845031.1',
	'211000022280328': 'NW_001844935.1',
	'211000022280341': 'NW_001846187.1',
	'211000022280347': 'NW_001845870.1',
	'211000022280481': 'NW_001845220.1',
	'211000022280494': 'NW_001845164.1',
	'211000022280703': 'NW_001845199.1',
	'rDNA': 'NW_007931121.1',
    },
    'MGI': {
	'1': 'NC_000067.6',
	'2': 'NC_000068.7',
	'3': 'NC_000069.6',
	'4': 'NC_000070.6',
	'5': 'NC_000071.6',
	'6': 'NC_000072.6',
	'7': 'NC_000073.6',
	'8': 'NC_000074.6',
	'9': 'NC_000075.6',
	'10': 'NC_000076.6',
	'11': 'NC_000077.6',
	'12': 'NC_000078.6',
	'13': 'NC_000079.6',
	'14': 'NC_000080.6',
	'15': 'NC_000081.6',
	'16': 'NC_000082.6',
	'17': 'NC_000083.6',
	'18': 'NC_000084.6',
	'19': 'NC_000085.6',
	'X': 'NC_000086.7',
	'Y': 'NC_000087.7',
	'MT': 'NC_005089.1',
    },
    'RGD': {
	'1': 'NC_005100.4',
	'2': 'NC_005101.4',
	'3': 'NC_005102.4',
	'4': 'NC_005103.4',
	'5': 'NC_005104.4',
	'6': 'NC_005105.4',
	'7': 'NC_005106.4',
	'8': 'NC_005107.4',
	'9': 'NC_005108.4',
	'10': 'NC_005109.4',
	'11': 'NC_005110.4',
	'12': 'NC_005111.4',
	'13': 'NC_005112.4',
	'14': 'NC_005113.4',
	'15': 'NC_005114.4',
	'16': 'NC_005115.4',
	'17': 'NC_005116.4',
	'18': 'NC_005117.4',
	'19': 'NC_005118.4',
	'20': 'NC_005119.4',
	'X': 'NC_005120.4',
	'Y': 'NC_024475.1',
	'MT': 'NC_001665.2',
    },
    'SGD': {
	'chrI': 'NC_001133.9',
	'chrII': 'NC_001134.8',
	'chrIII': 'NC_001135.5',
	'chrIV': 'NC_001136.10',
	'chrV': 'NC_001137.3',
	'chrVI': 'NC_001138.5',
	'chrVII': 'NC_001139.9',
	'chrVIII': 'NC_001140.6',
	'chrIX': 'NC_001141.2',
	'chrX': 'NC_001142.9',
	'chrXI': 'NC_001143.9',
	'chrXII': 'NC_001144.5',
	'chrXIII': 'NC_001145.3',
	'chrXIV': 'NC_001146.8',
	'chrXV': 'NC_001147.6',
	'chrXVI': 'NC_001148.4',
	'chrmt': 'NC_001224.1',
    },
    'WB': {
        'I': 'RefSeq:NC_003279.8',
        'II': 'RefSeq:NC_003280.10',
        'III': 'RefSeq:NC_003281.10',
        'IV': 'RefSeq:NC_003282.8',
        'V': 'RefSeq:NC_003283.11',
        'X': 'RefSeq:NC_003284.9',
        'MtDNA': 'RefSeq:NC_001328.1',
    },
    'ZFIN': {
	'1': 'NC_007112.7',
	'2': 'NC_007113.7',
	'3': 'NC_007114.7',
	'4': 'NC_007115.7',
	'5': 'NC_007116.7',
	'6': 'NC_007117.7',
	'7': 'NC_007118.7',
	'8': 'NC_007119.7',
	'9': 'NC_007120.7',
	'10': 'NC_007121.7',
	'11': 'NC_007122.7',
	'12': 'NC_007123.7',
	'13': 'NC_007124.7',
	'14': 'NC_007125.7',
	'15': 'NC_007126.7',
	'16': 'NC_007127.7',
	'17': 'NC_007128.7',
	'18': 'NC_007129.7',
	'19': 'NC_007130.7',
	'20': 'NC_007131.7',
	'21': 'NC_007132.7',
	'22': 'NC_007133.7',
	'23': 'NC_007134.7',
	'24': 'NC_007135.7',
	'25': 'NC_007136.7',
	'MT': 'NC_002333.2',
    },
}

expand_iupac = {
    'R': 'A,G',
    'Y': 'C,T',
    'S': 'C,G',
    'W': 'A,T',
    'K': 'G,T',
    'M': 'A,C',
    'B': 'C,G,T',
    'D': 'A,G,T',
    'H': 'A,C,T',
    'V': 'A,C,G',
    'N': 'A,C,G,T'
}

parser = argparse.ArgumentParser()
parser.add_argument("-j", "--json", help="JSON input file")
parser.add_argument("-g", "--gff", help="Corresponding GFF file")
parser.add_argument("-o", "--out", help="Output VCF file")
parser.add_argument("-m", "--mod", help="Acronym for MOD")
parser.add_argument("-f", "--fasta", help="FASTA file")
parser.add_argument("-s", "--strains", default=False, help="Input includes strain data")

args = parser.parse_args()
assembly, chr_lengths = get_header_info(args.gff)

vcf_file = open(args.out, 'w')

vcf_file.write("##fileformat=VCFv4.2\n" +
               datetime.datetime.today().strftime("##fileDate=%Y%m%d") + "\n" +
               "##reference=" + assembly + "\n" +
               "##source=AllianceJSON\n")

for chr in chrom2ncbi[args.mod].keys():
    if chr in chr_lengths.keys():
        vcf_file.write("##contig=<ID=" + chr + ",accession=\"" + chrom2ncbi[args.mod][chr] + "\",length=" + chr_lengths[chr] + ">\n")
    else:
        vcf_file.write("##contig=<ID=" + chr + ",accession=\"" + chrom2ncbi[args.mod][chr] + "\">\n")

vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

headers = ['#CHROM', 'POS', 'ID', 'REF',
           'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

nt_regex = re.compile('^[ACGT]$')

with open(args.json, 'r') as read_file:
    parsed = json.load(read_file)

    # get all strains for column headers
    if args.strains:
        strains = get_strains(parsed["data"])
        for s in strains:
            headers.append('WB:' + s)  # need curie form of strain
    
    vcf_file.write("\t".join(headers) + "\n")

    vcf_lines = []
    for v in (parsed["data"]):
        vcf_data = {}

        if 'genomicReferenceSequence' not in v.keys():
            if 'chromosome' not in v.keys() or 'start' not in v.keys() or 'end' not in v.keys():
                print("Neither genomicVariantSequence nor chr/start/end coordinates specified for " + v["alleleId"], file=sys.stderr)
                continue
            refSeq = get_refseq_from_fasta(v, args.fasta)
        else:
            if v["genomicReferenceSequence"] == "N/A":
                refSeq = ""  
            else:
                if args.mod == 'WB':
                    # don't need to check against reference for WB - makes it very slow for HTP variants
                    refSeq = v["genomicReferenceSequence"]
                else:
                    refSeq = get_refseq_from_fasta(v, args.fasta)
                    if v["genomicReferenceSequence"].upper() != refSeq:
                        print("Specified genomic reference allele (" + v["genomicReferenceSequence"] + ") doesn't match reference sequence ("
                              + refSeq + ") at specified coordinates for " + v["alleleId"], file=sys.stderr)
                        continue

    
        if 'genomicVariantSequence' not in v.keys():
            print("Missing reference allele for " + v["alleleId"], file=sys.stderr)
            continue

        varSeq = "" if v["genomicVariantSequence"] == "N/A" else v["genomicVariantSequence"]
        pos = int(v["start"])
        
        if refSeq == '' or varSeq == '':
            if refSeq == '' and varSeq == '':
                # FB has transgenes with insertions of unknown length
                print("No reference or alternative allele specified for " + v["alleleId"] + ', skipping', file=sys.stderr)
                continue

            if 'paddedBase' in v.keys():
                padBase = v["paddedBase"]
            else:
                padBase = get_padbase_from_fasta(v, refSeq, args.fasta)
            if pos == 1:
                refSeq = refSeq+padBase
                varSeq = varSeq+padBase
            else:
                refSeq = padBase+refSeq
                varSeq = padBase+varSeq
                if refSeq == '':
                    pos = pos-1  # include the padding base in POS for deletions

        vcf_data["chromosome"] = v["chromosome"]
        vcf_data["pos"] = pos

        if len(varSeq) == 1:
            if nt_regex.match(varSeq) is None:
                if varSeq in expand_iupac.keys():
                    varSeq = expand_iupac[varSeq]
                    if args.strains:
                        # Don't know genotypes of strains where there are multiple alternative alleles, so skip
                        continue
                else:
                    print("Unrecognised alternative allele " + varSeq + " for " + v["alleleId"], file=sys.stderr)
                    continue
        
        if args.strains:
            gtString = genotype_string(v, strains)
        
            vcf_data["line"] = "\t".join([v["chromosome"], str(
                pos), v["alleleId"], refSeq, varSeq, '.', 'PASS', '.', 'GT', gtString])
        else:
            vcf_data["line"] = "\t".join([v["chromosome"], str(pos), v["alleleId"], refSeq, varSeq, '.', '.' ,'.'])

        vcf_lines.append(vcf_data)

    for v in sorted(vcf_lines, key=operator.itemgetter('chromosome', 'pos')):
        vcf_file.write(v["line"] + "\n")

    vcf_file.close()
