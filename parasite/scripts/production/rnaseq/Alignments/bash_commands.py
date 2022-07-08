import os

def star_alignment(star, reference, threads, fastqs, outFileNamePrefix,
                           limitBAMsortRAM, sjdbOverhang, gtf, is_paired=1, extra_params=""):
    readFilesIn_string = ",".join([x[0] for x in fastqs])
    if is_paired==1:
        readFilesIn_string+= " "+",".join([x[1] for x in fastqs])
    bash_command = star + " " + \
                   "--genomeDir " + reference + " " + \
                   "--runThreadN " + threads + " " + \
                   "--outBAMsortingThreadN " + threads + " " + \
                   "--readFilesIn " + readFilesIn_string + " " + \
                   "--outFileNamePrefix " + outFileNamePrefix + " " + \
                   "--limitBAMsortRAM " + limitBAMsortRAM + " " + \
                   "--sjdbOverhang " + sjdbOverhang + " " + \
                   "--sjdbGTFfile " + gtf + " " + \
                   extra_params + ";\n\n"
    return (bash_command)

def cram2bam(cram, bam, fasta):
    bash_command="samtools view -b " + \
                 "-T " + fasta + " " + \
                 "-o " + bam + " " + \
                 cram + ";\n\n"
    return (bash_command)

def namesort_bam(inbam, outbam):
    bash_command="samtools sort -n " + \
                 "-o " + outbam + " " + \
                 inbam + ";\n\n"
    return bash_command
def positionsort_bam(inbam, outbam):
    bash_command="samtools sort " + \
                 "-o " + outbam + " " + \
                 inbam + ";\n\n"
    return bash_command
def fixmate_bam(inbam, outbam):
    bash_command="samtools fixmate -m " + \
                 inbam + " " + \
                 outbam + ";\n\n"
    return bash_command
def markdup_bam(inbam, outbam):
    bash_command="samtools markdup -r " + \
                 inbam + " " + \
                 outbam + ";\n\n"
    return bash_command
def index_bam(inbam):
    bash_command="samtools index " + \
                 inbam + ";\n\n"
    return bash_command
def merge_bams(inbams, outbam, threads):
    if (type(inbams) == str):
        inbams = [inbams]
    bash_command="samtools merge -f " + \
        "--threads " + str(threads) + " " + \
        outbam + " " + \
        " ".join(inbams) + ";\n\n"
    return(bash_command)
def sortrefname_bam(inbam, outbam, tmpdir, sortsamrefname):
    bash_command=sortsamrefname + " " + \
                 "--samoutputformat BAM " + \
                 "-o " + outbam + " " + \
                 "--tmpDir " + tmpdir + " " + \
                 inbam + ";\n\n"
    return(bash_command)
def cap_bam(inbam, outbam, cap_reads, biostar154220):
    bash_command=biostar154220 + " " + \
                 "-n " + str(cap_reads) + "  " + \
                 "--samoutputformat BAM " + \
                 "-o " + outbam + " " + \
                 inbam + ";\n\n"
    return bash_command
def bam2bigwig(inbam, outbw, binSize, bamCoverage):
    bash_command=bamCoverage + " " + \
                 "--bam " + inbam + " " + \
                 "--outFileName " + outbw + " " + \
                 "--binSize " + binSize + ";\n\n"
    return(bash_command)