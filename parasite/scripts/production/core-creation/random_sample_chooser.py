import re
import sys
import subprocess


#this script only picks out the first gene/transcript of a gff file to use has a placeholder example
#once Compara has been run, choose a more interesting example using the set_core_samples script (read "Creating Core Database" doc on Confluence for more details)

gff=open(sys.argv[1],"r")

#choose the first line that corresponds to a transcript
for line in gff:
	ifsample=re.search("mRNA",line)
	if ifsample:
		break

#in this line, look for the transcript name and parent (a gene ID)
sample=line

ifmrna=re.search("Name=(.*?)(;|$)",sample)
mrna=ifmrna.group(1)

findgene=re.search("Parent=(.*?);",sample)
geneid=findgene.group(1)

#look for the line corresponding to this gene ID in the gff
gff.seek(0)
for line in gff:
        ifgeneid=re.search("ID="+geneid+";",line)
        if ifgeneid:
                break

geneline=line

#in this line, look for the gene name
ifgene=re.search("Name=(.*?)(;|$)",geneline)
gene=ifgene.group(1)

print sample
#look for the seq and position of the transcript
seq=re.search("^(.*?)\t",sample)
pos=re.search("([0-9]+)\t([0-9]+)\t",sample)
region=seq.group(1) + ":" + pos.group(1) + "-" + pos.group(2)

print "sample.gene_param: " + gene
print "sample.gene_text: " + gene
print "sample.location_param: " + region
print "sample.location_text: " + region
print "sample.search_text: ribosomal"
print "sample.transcript_param: " + mrna
print "sample.transcript_text: " + mrna

