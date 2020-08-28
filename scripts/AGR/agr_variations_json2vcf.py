#!/usr/bin/env python
# based on 
#       * VCF v4.3 (http://samtools.github.io/hts-specs/VCFv4.3.pdf)
#       * AGR JSON 1.0.1.3 (https://docs.google.com/document/d/1yAECtOs1VCEs3mhplJMXqg1akBQJyJbjPnvG2RrE5Aw/edit)

import json
import datetime
import sys
import operator


def genotype_string(variation, allStrains):
    variationStrains = set(variation["strains"])
    genotypes = []

    for s in allStrains:
        if s in variationStrains:
            genotypes.append('1/1')
        else:
            genotypes.append('./.')

    return "\t".join(genotypes)
            
def get_strains(variations):
    strains = set()
    for variation in variations:
        for strain in variation["strains"]:
            strains.add(strain)
    return tuple(sorted(strains))


chrom2ncbi = {
	'I': 'RefSeq:NC_003279.8',
	'II': 'RefSeq:NC_003280.10',
	'III': 'RefSeq:NC_003281.10',
	'IV': 'RefSeq:NC_003282.8',
	'V': 'RefSeq:NC_003283.11',
	'X': 'RefSeq:NC_003284.9',
	'MtDNA': 'RefSeq:NC_001328.1',
}


print "##fileformat=VCFv4.3"
print datetime.datetime.today().strftime("##fileDate=%Y%m%d")
print "##source=AllianceJSON"
print "##INFO<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
headers = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']

with open(sys.argv[1], 'r') as read_file:
    parsed = json.load(read_file)

    # get all strains for column headers
    strains = get_strains(parsed["data"])
    for s in strains:
        headers.append('WB:' + s) # need curie form of strain
    print "\t".join(headers)
    
    for v in sorted(parsed["data"],key=operator.itemgetter('chromosome','start')):

        refSeq   = "" if v["genomicReferenceSequence"]=="N/A" else v["genomicReferenceSequence"]
        varSeq   = "" if v["genomicVariantSequence"]  =="N/A" else v["genomicVariantSequence"]
        pos      = int(v["start"])
        gtString = genotype_string(v, strains) 

        
        if 'paddedBase' in v.keys():
            if pos == 1:
                refSeq = refSeq+v["paddedBase"]
                varSeq = varSeq+v["paddedBase"]
            else:
                refSeq = v["paddedBase"]+refSeq
                varSeq = v["paddedBase"]+varSeq
                pos = pos-1 # include the padding base in POS
        print "\t".join([v["chromosome"], str(pos), v["alleleId"], refSeq, varSeq, '.', 'PASS', 'GT', gtString])


