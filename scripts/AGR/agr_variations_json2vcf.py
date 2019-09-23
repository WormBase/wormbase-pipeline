#!/usr/bin/env python
# based on 
#       * VCF v4.3 (http://samtools.github.io/hts-specs/VCFv4.3.pdf)
#       * AGR JSON 1.0.0.8 (https://github.com/alliance-genome/agr_schemas/blob/release-1.0.0.8/)

import json
import datetime
import sys
import operator

print "##fileformat=VCFv4.3"
print datetime.datetime.today().strftime("##fileDate=%Y%m%d")
print "##source=AllianceJSON"
print "##INFO<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
print "\t".join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'])

with open(sys.argv[1], 'r') as read_file:
    parsed = json.load(read_file)
    for v in sorted(parsed["data"],key=operator.itemgetter('chromosome','start')):

        refSeq = "" if v["genomicReferenceSequence"]=="N/A" else v["genomicReferenceSequence"]
        varSeq = "" if v["genomicVariantSequence"]  =="N/A" else v["genomicVariantSequence"]
        pos    = int(v["start"])

        if 'paddedBase' in v.keys():
            if pos == 1:
                refSeq = refSeq+v["paddedBase"]
                varSeq = varSeq+v["paddedBase"]
            else:
                refSeq = v["paddedBase"]+refSeq
                varSeq = v["paddedBase"]+varSeq
                pos = pos-1 # include the padding base in POS
        print "\t".join([v["chromosome"],str(pos),v["alleleId"],refSeq,varSeq,'.','PASS','DP=100'])
