#!/software/bin/perl -w
#
# align_RNASeq.pl
#
# script to process the RNASeq data using Tophat and Cufflinks , etc.
#
# currently we have for elegans:
#
# Waterston elegans RNASeq data:
#
# SRX001872 SRX001875 SRX004865 SRX004868 SRX008138 SRX008144 SRX014007
# SRX014010 SRX036882 SRX036970 SRX037198 SRX037288 SRX001873 SRX004863
# SRX004866 SRX004869 SRX008139 SRX011569 SRX014008 SRX035162 SRX036967
# SRX037186 SRX037199 SRX001874 SRX004864 SRX004867 SRX008136 SRX008140
# SRX014006 SRX014009 SRX036881 SRX036969 SRX037197 SRX037200
#
#
## Andy Fraser's Lab 
#     
# SRX007069	     SRX007170	     SRX007171	     SRX007172	     SRX007173	     
#
## Andy Fire's Lab
#	     
# SRX028190	     SRX028191	     SRX028192	     SRX028193	     SRX028194	     
# SRX028195	     SRX028196	     SRX028197	     SRX028198	     SRX028199	     
# SRX028200	     SRX028202	     SRX028203	     SRX028204
#
#
# and for brugia:
#
# some data taken directly from Matt Berriman group's FTP server
#
#
# and for remanei:
#
# SRX052082 SRX052083
#
# and for briggsae:
#
# SRX052079 SRX052081 SRX053351
#
# and for japonica:
#
# SRX100095 SRX100094 SRX100093 SRX100092 SRX100091 SRX100090


# Perl should be set to:
# alias perl /software/bin/perl
# otherwise it can't find LSF.pm

# Run this after TSL Features have been mapped in the Build.
# Run on farm-login2 as:
# bsub -I -q basement align_RNASeq.pl -noalign
#
# omit the '-noalign' if this is a frozen release and a fresh alignment will be done
# add '-check' if the run died and is being resumed
# add '-species' to run on a non-elegans Build
#
# bsub -I -e /nfs/wormpub/BUILD/autoace/logs/align.err -o /nfs/wormpub/BUILD/autoace/logs/align.out -q basement  align_RNASeq.pl [-noalign] [-check] [-species $SPECIES]

#------------------------------------------------------------------------------------
# Notes on using GSNAP



# to build a GMAP genome database
# requires a machine with 4Gb of memory, so specify this in LSF
# /software/worm/gmap/bin/gmap_build 
# Usage: gmap_build [options...] -d <genomename> <fasta_files>
#
# Options:
# --dir=STRING        Destination directory for installation (defaults to gmapdb directory specified at configure time)
# --db=STRING         Genome name
# -T STRING               Temporary build directory
# --sort=STRING       Sort chromosomes using given method:
#			      none - use chromosomes as found in FASTA file(s)
#			      alpha - sort chromosomes alphabetically (chr10 before chr 1)
#			      numeric-alpha - chr1, chr1U, chr2, chrM, chrU, chrX, chrY
#			      chrom - chr1, chr2, chrM, chrX, chrY, chr1U, chrU
# --gunzip            Files are gzipped, so need to gunzip each file first
# if you made a database called <genome>, your directory structure should look like this
#
#    /path/to/gmapdb/<genome>/
#    /path/to/gmapdb/<genome>/<genome>.chromosome
#    /path/to/gmapdb/<genome>/<genome>.chromosome.iit
#    ...
#    /path/to/gmapdb/<genome>/<genome>.version
#
# e.g.
# bsub6 -Is /software/worm/gmap/bin/gmap_build --dir /nfs/wormpub/RNASeq/elegans/GSNAP/ --db elegans -T /nfs/wormpub/RNASeq/elegans/GSNAP/ --sort none -k 15 elegans.genome.fa
# (the above runs out of memory if you specify only 4 Gb of memory)


# GSNAP options we probably require:
# --dir=<directory containing the genome> - sets the /path/to/gmapdb/ directory containing the genome database
# --db=STRING  - name of genome database
# --novelsplicing=1  - Look for novel splicing
# --use-splicing=STRING  - Look for splicing involving known sites or known introns
#                                         (in <STRING>.iit), at short or long distances
#                                         See README instructions for the distinction between known sites
#                                         and known introns
# --localsplicedist=1000  - Definition of local novel splicing event (default 200000)
#                           Have a small penalty if the splice exceeds this length
# --pairmax-rna=1000  - Max total genomic length for RNA-Seq paired reads, or other reads
#                                   that could have a splice (default 200000).  Used if -N or -s is specified.
#                                   Should probably match the value for --localsplicedist.
# --quality-protocol=STRING  - Protocol for input quality scores.  Allowed values:
#                                   illumina (ASCII 64-126) (equivalent to -J 64 -j -31)
#                                   sanger   (ASCII 33-126) (equivalent to -J 33 -j 0)
#                                 Default is sanger (no quality print shift)
#
#  --fails-as-input               Print completely failed alignments as input FASTA or FASTQ format
#  --format=STRING            Another output format type, other than default.
#                                   Currently implemented: sam
#
# e.g.
# bsub4 -Is /software/worm/gmap/bin/gsnap --dir /nfs/wormpub/RNASeq/elegans/GSNAP/ --db elegans --novelsplicing=1 --localsplicedist=1000 --pairmax-rna=1000 --quality-protocol=sanger --fails-as-input --format=sam /lustre/scratch101/ensembl/wormpipe/RNASeq/elegans/SRA/SRX049745/SRX049745/SRR137926/SRR137926.fastq
# requires at least 4 Gb of memory to run, otherwise it fails quietly with a truncated output file.
# size of output when adding in various parameters
# --novelsplicing=1       "User defined signal 2" after 3.1 Gb output when using bsub4
# --localsplicedist=1000 
# --pairmax-rna=1000      

# known splice site file
# I would recommend using known splice sites, rather
# than known introns, unless you are certain that all alternative
# splicing events are known are represented in your file.
#
# GSNAP can tell the difference between known site-level and known
# intron-level splicing based on the format of the input file.  To
# perform known site-level splicing, you will need to create a file with
# the following format:

# >NM_004448.ERBB2.exon1 17:35110090..35110091 donor 6678
# >NM_004448.ERBB2.exon2 17:35116768..35116769 acceptor 6678
# >NM_004448.ERBB2.exon2 17:35116920..35116921 donor 1179
# >NM_004448.ERBB2.exon3 17:35118099..35118100 acceptor 1179
# >NM_004449.ERG.exon1 21:38955452..38955451 donor 783
# >NM_004449.ERG.exon2 21:38878740..38878739 acceptor 783
# >NM_004449.ERG.exon2 21:38878638..38878637 donor 360
# >NM_004449.ERG.exon3 21:38869542..38869541 acceptor 360

# Each line must start with a ">" character, then be followed by an
# identifier, which may have duplicates and can have any format, with
# the gene name or exon number shown here only as a suggestion.  Then
# there should be the chromosomal coordinates which straddle the
# exon-intron boundary, so one coordinate is on the exon and one is on
# the intron.  (Coordinates are all 1-based, so the first character of a
# chromosome is number 1.)  Finally, there should be the splice type:
# "donor" or "acceptor".  You may optionally store the intron distance
# at the end.  GSNAP can use this intron distance, if it is longer than
# its value for --localsplicedist, to look for long introns at that
# splice site.  The same splice site may have different intron distances
# in the database; GSNAP will use the longest intron distance reported
# in searching for long introns.
                                                                                                
# Note that the chromosomal coordinates are in the sense direction.
# Therefore, genes on the plus strand of the genome (like NM_004448) have
# the coordinates in ascending order (e.g., 35110090..35110091).
# Genes on the minus strand of the genome (like NM_004449) have the
# coordinates in descending order (e.g., 38955452..38955451).

# if you have a GTF file, you can use the included program
# gtf_splicesites like this:

#    cat <gtf file> | gtf_splicesites > foo.splicesites

# Once you have built this splicesites or introns file, you process it
# as a map file (see "Building map files" above), by doing

#    cat foo.splicesites | iit_store -o <splicesitesfile>

# you put it in the maps subdirectory by doing
                                                                                                
#    cp <splicesitesfile>.iit /path/to/gmapdb/<genome>/<genome>.maps
                                                                                                
# Then, you may use the file by doing this:
                                                                                                
#    gsnap -d <genome> -s <splicesitesfile> <shortreads>


#------------------------------------------------------------------------------------



# by Gary Williams
#
# Last updated by: $Author: gw3 $
# Last updated on: $Date: 2012-10-12 11:13:54 $

#################################################################################
# Initialise variables                                                          #
#################################################################################

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Storable;
use Log_files;
use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;
#use Ace;
use Coords_converter;
use Carp;
use Feature_mapper;



##############################
# command-line options       #
##############################


my ($test, $store, $debug, $species, $verbose, $expt, $noalign, $check, $solexa, $illumina, $database, $nogtf, $norawjuncs, $tsl, $paired, $keepfastq, $onestage, $getsrx, $ribosome_occupancy);
GetOptions (
	    "test"               => \$test,
	    "store:s"            => \$store,
	    "debug:s"            => \$debug,
	    "verbose"            => \$verbose,
	    "database:s"         => \$database, # acedb database to use for GTF, splice junctions and TSL Features
	    "species:s"		 => \$species,  # set species explicitly
	    "noalign"            => \$noalign,  # do not run tophat to align unless the old BAM file is missing
	    "check"              => \$check,    # in each experiment, check to see if tophat/cufflinks worked - if not repeat them
	    "nogtf"              => \$nogtf,    # don't use the GTF file of existing gene models to create the predicted gene models
	    "norawjuncs"         => \$norawjuncs,  # don't use the raw junctions file of splice hints (for example in a new brugia assembly)
	    "tsl"                => \$tsl,      # look for reads with a trans-splice leader sequence
	    "solexa"             => \$solexa,   # reads have solexa quality scores, not phred + 33 scores
	    "expt:s"             => \$expt,     # do the alignment etc. for this experiment
	    "illumina"           => \$illumina, # read have illumina GA-Pipeline 1.3 quality scores, not phred + 33 scores
	    "paired"             => \$paired,   # the reads are paired-end - we may have to split the file
	    "keepfastq"          => \$keepfastq,# don't delete the fastq files after use
	    "onestage"           => \$onestage, # for comparing the old one-stage alignment against the new (default) two-stage (transcriptome then genome) alignment
	    "getsrx:s"           => \$getsrx,   # if this is set to the name of a 'SRX' entry from the NCBI SRA, then this file will be retrieved and unpacked - no alignment is done - '-paired' will split the file into pairs
	    "ribosome_occupancy" => \$ribosome_occupancy # this is a project from the Andy Fire lab looking at which codons have how many ribosomes - used by me to find START codons
	   );

my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			     -organism => $species
			   );
}

my $log = Log_files->make_build_log( $wormbase );

$species = $wormbase->species;
$log->write_to("Processing RNASeq data for $species\n");

my $Software = "/nfs/panda/ensemblgenomes/wormbase/software/packages";

#################################################################
# Quick utility for getting a file from the Short Read Archive.
# If the -paired switch is also used, then the files will be split
# into the paired reads
#################################################################

if ($getsrx) {

  get_SRX_file($getsrx);

  # unpack any .lite.sra file to make the .fastq files
  # now unpack the sra.lite files
  $log->write_to("Unpack the $getsrx .lite.sra file to fastq files\n");
  foreach my $dir (glob("$getsrx/SRR*")) {
    chdir $dir;
    foreach my $srr (glob("*.lite.sra")) {
      $log->write_to("Unpack $srr\n");
      my $options = '';
      if ($paired && $srr !~ /_1\./ && $srr !~ /_2\./) {$options='--split-files'} # specify that the file contains paired-end reads and should be split
#      my $status = $wormbase->run_command("/software/worm/sratoolkit/fastq-dump $options $srr", $log);
      my $status = $wormbase->run_command("$Software/sratoolkit/fastq-dump $options $srr", $log);
      if ($status) {$log->log_and_die("Didn't unpack the fastq file successfully\n")}
      $status = $wormbase->run_command("rm -f $srr", $log);
    }
  }

  exit(0);
}

##########################################
# Data to be processed
##########################################

my $currentdb = $wormbase->database('current');
my $dbdir  = $wormbase->autoace;                                     # Database path
$database = $dbdir unless $database;
#$database = $currentdb unless $database;
$species = $wormbase->species;                          # set the species if it was given

# global variables
my $coords;
my $mapper;
my $status;

# SRP000401  Deep sequencing of the Caenorhabditis elegans transcriptome using RNA isolated from various developmental stages under various experimental conditions RW0001
#


# quality score description in: http://en.wikipedia.org/wiki/FASTQ_format
my %expts;

if ($species eq 'elegans') {

  %expts = ( # key= SRA 'SRX' experiment ID, values = [Analysis ID, quality score metric]

# Add new analysis/condition objects to ~wormpub/DATABASES/geneace

#	     SRX001872  => ["RNASeq_Hillier.L2_larva", 'phred', 'single'], # needs a lot of memory to run tophat
	     SRX001873  => ["RNASeq_Hillier.Young_Adult", 'phred', 'single'],
	     SRX001874  => ["RNASeq_Hillier.L4_larva", 'phred', 'single'],
	     SRX001875  => ["RNASeq_Hillier.L3_larva", 'phred', 'single'],

	     SRX004863  => ["RNASeq_Hillier.early_embryo_Replicate3", 'phred', 'single'], 
	     SRX004864  => ["RNASeq_Hillier.early_embryo", 'solexa', 'single'], # solexa
	     SRX004865  => ["RNASeq_Hillier.late_embryo", 'phred', 'single'],
	     SRX004866  => ["RNASeq_Hillier.late_embryo-Replicate1", 'solexa', 'single'], # SRR016678 is phred - the other three are solexa
	     SRX004867  => ["RNASeq_Hillier.L1_larva", 'solexa', 'paired-end'], # solexa - paired end reads
	     SRX004868  => ["RNASeq_Hillier.L4_larva_Male", 'solexa', 'single'], # solexa
	     SRX004869  => ["RNASeq_Hillier.L1_larva_lin-35", 'solexa', 'single'], # SRR016691 and SRR016690 are phred, the other two are solexa

 	     SRX008136  => ["RNASeq_Hillier.Adult_spe-9", 'phred', 'single'], # phred, despite what the docs say
 	     SRX008138  => ["RNASeq_Hillier.dauer_daf-2", 'phred', 'single'], # phred, despite what the docs say
 	     SRX008139  => ["RNASeq_Hillier.dauer_entry_daf-2", 'phred', 'single'], # phred, despite what the docs say
 	     SRX008140  => ["RNASeq_Hillier.dauer_exit_daf-2", 'phred', 'single'], # phred, despite what the docs say
 	     SRX008144  => ["RNASeq_Hillier.L4_Larva_Replicate1", 'phred', 'single'], # phred, despite what the docs say

 	     SRX011569  => ["RNASeq_Hillier.embryo_him-8", 'phred', 'single'], # phred, despite what the docs say

 	     SRX014006  => ["RNASeq_Hillier.young_adult_Harposporium", 'phred', 'single'], # phred, despite what the docs say
 	     SRX014007  => ["RNASeq_Hillier.young_adult_Harposporium_control", 'phred', 'single'], # phred, despite what the docs say
 	     SRX014008  => ["RNASeq_Hillier.young_adult_S_macescens", 'phred', 'single'], # phred, despite what the docs say
 	     SRX014009  => ["RNASeq_Hillier.young_adult_S_macescens_control", 'phred', 'single'], # phred, despite what the docs say
 	     SRX014010  => ["RNASeq_Hillier.L4_larva_JK1107", 'phred', 'single'], # phred, despite what the docs say

 	     SRX035162  => ["RNASeq_Hillier.adult_D_coniospora_12hrs", 'phred', 'single'], # 76 bp reads

 	     SRX036967  => ["RNASeq_Hillier.adult_D_coniospora_control", 'phred', 'single'], # 76 bp reads

 	     SRX036881  => ["RNASeq_Hillier.L3_larva_Replicate1", 'phred', 'single'], # 76 bp reads
 	     SRX036882  => ["RNASeq_Hillier.adult_D_coniospora_5hrs", 'phred', 'single'], # 76 bp reads

 	     SRX036969  => ["RNASeq_Hillier.adult_E_faecalis", 'phred', 'single'], # 76 bp reads
 	     SRX036970  => ["RNASeq_Hillier.adult_P_luminescens", 'phred', 'single'], # 76 bp reads

 	     SRX037186  => ["RNASeq_Hillier.early_embryo_Replicate2", 'phred', 'single'], # 76 bp reads
 	     SRX037197  => ["RNASeq_Hillier.early_embryo_Replicate1", 'phred', 'single'], # modENCODE_3526
 	     SRX037198  => ["RNASeq_Hillier.embryo_Male_him-8-Replicate1", 'phred', 'single'], # 76 bp reads
 	     SRX037199  => ["RNASeq_Hillier.dauer_exit_daf-2-Replicate1", 'phred', 'single'], # 76 bp reads
 	     SRX037200  => ["RNASeq_Hillier.L4_larva_JK1107_Replicate1", 'phred', 'single'], # 76 bp reads

 	     SRX037288  => ["RNASeq_Hillier.L1_larva-Replicate1", 'phred', 'single'], # 76 bp reads

 	     SRX047446  => ["RNASeq_Hillier.late_embryo_Replicate2", 'phred', 'single'], # 76 bp reads
 	     SRX047469  => ["RNASeq_Hillier.L4_larva_Male_Replicate1", 'phred', 'single'], # 76 bp reads
 	     SRX047470  => ["RNASeq_Hillier.dauer_entry_daf-2_Replicate1", 'phred', 'single'], # 76 bp reads

 	     SRX047653  => ["RNASeq_Hillier.L2_larva_Replicate1", 'phred', 'single'], # 76 bp reads
 	     SRX047787  => ["RNASeq_Hillier.Young_Adult_Replicate1", 'phred', 'single'], # 76 bp reads
	     
## Andy Fraser's Lab 
	     
	     SRX007069   => ["RNASeq_Fraser.L4_larva", 'phred', 'single'], # N2, L4 larva
	     SRX007170   => ["RNASeq_Fraser.adult", 'phred', 'single'], # N2, Young Adult
	     SRX007171   => ["RNASeq_Fraser.smg-1_L4_larva", 'phred', 'single'], # smg-1, L4 larva
	     SRX007172   => ["RNASeq_Fraser.smg-1_adult", 'phred', 'single'], # smg-1, Young Adult
	     SRX007173   => ["RNASeq_Fraser.all_stages", 'phred', 'paired-end'], # N2, all_stages, paired end reads
	     
## Andy Fire's Lab
	     
	     SRX028190   => ["RNASeq_Fire.all_stages_dsDNALigSeq", 'phred', 'single'], # GSM577107: N2_mixed-stage_dsDNALigSeq
	     SRX028191   => ["RNASeq_Fire.all_stages_ssRNALigSeq", 'phred', 'single'], # GSM577108: N2_mixed-stage_ssRNALigSeq
	     SRX028192   => ["RNASeq_Fire.fem-3_dsDNALigSeq", 'phred', 'single'],      # GSM577110: fem-3_dsDNALigSeq
	     SRX028193   => ["RNASeq_Fire.fem-1_dsDNALigSeq", 'phred', 'single'],      # GSM577111: fem-1_dsDNALigSeq
	     SRX028194   => ["RNASeq_Fire.him-8_dsDNALigSeq", 'phred', 'single'],      # GSM577112: him-8_dsDNALigSeq
	     SRX028195   => ["RNASeq_Fire.him-8_ssRNALigSeq", 'phred', 'single'],      # GSM577113: him-8_ssRNALigSeq
	     SRX028196   => ["RNASeq_Fire.rrf-3_him-8_ssRNALigSeq", 'phred', 'single'], # GSM577114: rrf-3_him-8_ssRNALigSeq
	     SRX028197   => ["RNASeq_Fire.L1_larva_ssRNALigSeq", 'phred', 'single'],   # GSM577115: N2_L1_ssRNALigSeq
	     SRX028198   => ["RNASeq_Fire.L2_larva_ssRNALigSeq", 'phred', 'single'],   # GSM577116: N2_L2_ssRNALigSeq
	     SRX028199   => ["RNASeq_Fire.L3_larva_ssRNALigSeq", 'phred', 'single'],   # GSM577117: N2_L3_ssRNALigSeq
	     SRX028200   => ["RNASeq_Fire.L4_larva_ssRNALigSeq", 'phred', 'single'],   # GSM577118: N2_L4_ssRNALigSeq
	     SRX028201   => ["RNASeq_Fire.L1_larva_CircLigSeq", 'phred', 'single'],    # GSM577119: N2_L1_CircLigSeq
	     SRX028202   => ["RNASeq_Fire.L2_larva_CircLigSeq", 'phred', 'single'],    # GSM577120: N2_L2_CircLigSeq
	     SRX028203   => ["RNASeq_Fire.L3_larva_CircLigSeq", 'phred', 'single'],    # GSM577121: N2_L3_CircLigSeq
	     SRX028204   => ["RNASeq_Fire.all_stages_polysomes", 'phred', 'single'],   # GSM577122: N2_mixed-stage_polysomes - these are spceial in that they show that this transcript is translated, not just transcribed

## new entries 7 Oct 2011
	    SRX017684 => ["RNASeq.Shin.L1_larva", 'phred', 'single'], # 454 sequencing http://www.biomedcentral.com/1741-7007/6/30

	    SRX006984 => ["RNASeq.Gerstein.embryo_Replicate1", 'phred', 'single'], # 
	    SRX006985 => ["RNASeq.Gerstein.embryo_Replicate2", 'phred', 'single'], # 
	    SRX006986 => ["RNASeq.Gerstein.embryo_Replicate3", 'phred', 'single'], # 
	    SRX006987 => ["RNASeq.Gerstein.L1_larva.starved_Replicate1", 'phred', 'single'], # 
	    SRX006988 => ["RNASeq.Gerstein.L1_larva.starved_Replicate2", 'phred', 'single'], # 
	    SRX006989 => ["RNASeq.Gerstein.L1_larva.starved_Replicate3", 'phred', 'single'], # 

	    SRX026728 => ["RNASeq.Martazavi.L3_larva", 'phred', 'paired-end'], # paired-end http://genome.cshlp.org/content/20/12/1740.long
	    SRX026729 => ["RNASeq.Martazavi.L1_larva", 'phred', 'paired-end'], # paired-end http://genome.cshlp.org/content/20/12/1740.long

	    SRX047635 => ["RNASeq_Hillier.adult_D_coniospora_12hrs_Replicate2", 'phred', 'single'], # 76 bp reads

	    SRX049268 => ["RNASeq_Hillier.L4_larva_cap4", 'phred', 'single'], # 76 bp reads. Array capture experiment to remove messages already detected in other RNAseq experiments on C. elegans N2 larval L4 library.
	    SRX049269 => ["RNASeq_Hillier.L4_larva_RRcap3", 'phred', 'single'], # 76 bp reads. Ribominus kit was used on the C. elegans N2 larval L4 library and subjected to array capture to remove messages already detected in other RNAseq experiments. 
	    SRX049270 => ["RNASeq_Hillier.late_embryo_cap6", 'phred', 'single'], # 76 bp reads. Array capture experiment on late embryo LE-2 library to remove messages already detected in other RNAseq experiments.

	    SRX049744 => ["RNASeq_Hillier.L1_larva_lin-35_cap1", 'phred', 'single'], # 	modENCODE_3521 76 bp reads. Array capture experiment to remove messages already detected in other RNAseq experiments on C. elegans N2 larva L1 (lin-35) library.
	    SRX049745 => ["RNASeq_Hillier.L4_larva_Male_cap2", 'phred', 'single'], #  modENCODE_3524 76 bp reads. Array capture experiment to remove messages already detected in other RNAseq experiments on C. elegans N2 larva L4 male library.

	    SRX050610 => ["RNASeq_Hillier.L1_larva_lin-35-cap1", 'phred', 'single'], # 76 bp reads. Array capture experiment to remove messages already detected in other RNAseq experiments on C. elegans N2 larva L1 (lin-35) library. 

	    SRX050630 => ["RNASeq_Hillier.L4_larva_Male_cap2_Replicate2", 'phred', 'single'], # NOT SURE IF THIS IS A REPLICATE - ASKED LADEANA. Array capture experiment to remove messages already detected in other RNAseq experiments on C. elegans N2 larva L4 male library.

	    SRX085111 => ["RNASeq_Hillier.blastula_embryo", 'phred', 'paired-end'], # 100 bp reads. paired-end. EE_50-60
	    SRX085112 => ["RNASeq_Hillier.blastula_embryo_Replicate2", 'phred', 'paired-end'], # 100 bp reads. paired-end. EE_50-60 NOT SURE IF THIS IS A REPLICATE

	    SRX085217 => ["RNASeq_Hillier.gastrulating_embryo", 'phred', 'paired-end'], # 76 bp reads. paired-end. EE_50-120
	    SRX085218 => ["RNASeq_Hillier.gastrulating_embryo_Replicate2", 'phred', 'paired-end'], # 100 bp reads. paired-end. EE_50-120
	    SRX085219 => ["RNASeq_Hillier.4-cell_embryo", 'phred', 'single'], # Illumina sequencing of C. elegans N2 hand-picked 4 cell embryos DSN-Negative, total RNA.
# get "terminate called after throwing an instance of 'std::bad_alloc'" in cufflinks with this one	    SRX085220 => ["RNASeq_Hillier.gastrulating_embryo_Replicate3", 'phred', 'single'], # 100 bp reads. paired-end. EE_50-120

	    SRX085286 => ["RNASeq_Hillier.Pharyngeal_muscle_late_embryo", 'phred', 'paired-end'], # 100 bp reads. paired-end. Pharyngeal muscle total RNA, late embryo. 
	    SRX085287 => ["RNASeq_Hillier.Pharyngeal_muscle_late_embryo_Replicate2", 'phred', 'paired-end'], # 76 bp reads. paired-end. Pharyngeal muscle total RNA, late embryo.

	    SRX092371 => ["RNASeq_Hillier.2-cell_embryo", 'phred', 'paired-end'], # 76 bp reads. paired-end. EE_50-30
	    SRX092372 => ["RNASeq_Hillier.2-cell_embryo_Replicate2", 'phred', 'paired-end'], # 100 bp reads. paired-end. EE_50-30

	    SRX092477 => ["RNASeq_Hillier.1-cell_embryo", 'phred', 'paired-end'], # 76 bp reads. paired-end. EE_50-0
	    SRX092478 => ["RNASeq_Hillier.1-cell_embryo_Replicate2", 'phred', 'paired-end'], # 100 bp reads. paired-end. EE_50-0
	    SRX092479 => ["RNASeq_Hillier.28-cell_embryo", 'phred', 'paired-end'], #  100 bp reads. paired-end. EE_50-90
	    SRX092480 => ["RNASeq_Hillier.28-cell_embryo_Replicate2", 'phred', 'paired-end'], # 76 bp reads. paired-end. EE_50-90

	    SRX099901 => ["RNASeq_Hillier.1-cell_embryo_Replicate3", 'phred', 'single'], # 100 bp reads. EE_50-0
	    SRX099902 => ["RNASeq_Hillier.1-cell_embryo_Replicate4", 'phred', 'single'], # 101 bp reads. EE_50-0

	    SRX099907 => ["RNASeq_Hillier.2-cell_embryo_Replicate3", 'phred', 'single'], # 100 bp reads. EE_50-30
	    SRX099908 => ["RNASeq_Hillier.2-cell_embryo_Replicate4", 'phred', 'single'], # 101 bp reads. EE_50-30

	    SRX099915 => ["RNASeq_Hillier.28-cell_embryo_Replicate3", 'phred', 'single'], # 101 bp reads. EE_50-90

	    );

  if ($ribosome_occupancy) {

    # these are the C. elegans data-sets for the ribosome occupancy
    # project "SRA049309.1" by Michael Stadler from Andy Fire's lab
    # From the Ribosome occupancy paper:
    # http://www.ncbi.nlm.nih.gov/pubmed/22045228
    
    %expts = ( # key= SRA 'SRX' experiment ID, values = [Analysis ID, quality score metric]
	      
	      SRX118140  => ["L4 ribosome footprints biological replicate 2", 'phred', 'single'], 
	      SRX118139  => ["L3 ribosome footprints biological replicate 2", 'phred', 'single'], 
	      SRX118118  => ["C. elegans L1 ribosome footprints", 'phred', 'single'], 
	      SRX116363  => ["C. elegans L1 mRNAseq", 'phred', 'single'], 
	      SRX116361  => ["C. elegans L1 stage mRNAseq", 'phred', 'single'], 
	      SRX118144  => ["L3 ribosome footprints prepared without cycloheximide, technical replicate 2", 'phred', 'single'], 
	      SRX118143  => ["L3 ribosome footprints prepared without cycloheximide, technical replicate 1", 'phred', 'single'], 
	      SRX118128  => ["L3 ribosome footprints biological replicate 1", 'phred', 'single'], 
	      SRX118127  => ["C. elegans L4 ribosome footprints biological replicate 1", 'phred', 'single'], 
	      SRX118126  => ["C. elegans L4 mRNAseq biological replicate 1", 'phred', 'single'], 
	      SRX118125  => ["C. elegans L3 mRNAseq biological replicate 1", 'phred', 'single'], 
	      SRX118124  => ["C. elegans L3 mRNAseq biological replicate 2", 'phred', 'single'], 
	      SRX118120  => ["C. elegans L2 mRNAseq biological replicate 2", 'phred', 'single'], 
	      SRX118119  => ["C. elegans L2 mRNAseq biological replicate 1", 'phred', 'single'], 
	      SRX118117  => ["C. elegans L1 ribosome footprints", 'phred', 'single'], 
	      SRX118116  => ["C. elegans L1 mRNAseq", 'phred', 'single'], 
	     );
  }

} elsif ($species eq 'brugia') {

# data from Matt Berriman's group.
# /nfs/disk69/ftp/pub4/pathogens/Brugia/malayi
#
# Seven libraries have thus far been sequenced :-
#
# Adult male
# Adult female
# mature microfillariae
# immature microfillariae
# L3 stage
# L4 stage
# eggs embryos
#
# The directory is structured as follows
#
# DATA : contains the raw sequence data in fastq files and is
# organised as follows
#
# - DATA/SLX/library_name/lane/lane_1.fastq - first of the pair
# - DATA/SLX/library_name/lane/lane_2.fastq - second of the pair

# The dataset names are what the Berriman group called them - these data are not downloaded from the SRA

  %expts = ( # key= SRA 'SRX' experiment ID, values = [Analysis ID, quality score metric]
# Add new analysis/condition objects to ~wormpub/DATABASES/brugia

	    Adult_female => ["RNASeq_brugia_Berriman.Adult_female", 'phred', 'paired-end'],
	    Adult_male => ["RNASeq_brugia_Berriman.Adult_male", 'phred', 'paired-end'],
	    BmL3_1361258 => ["RNASeq_brugia_Berriman.BmL3_1361258", 'phred', 'paired-end'],
	    eggs_embryos => ["RNASeq_brugia_Berriman.eggs_embryos", 'phred', 'paired-end'],
	    immature_female => ["RNASeq_brugia_Berriman.immature_female", 'phred', 'paired-end'],
	    L3_stage => ["RNASeq_brugia_Berriman.L3_stage", 'phred', 'paired-end'],
	    L4 => ["RNASeq_brugia_Berriman.L4", 'phred', 'paired-end'],
	    microfillariae => ["RNASeq_brugia_Berriman.microfillariae", 'phred', 'paired-end'],
	    mayhew => ["RNASeq_brugia_mayhew.mf.L2.L3", 'phred', 'single'],
	   );
	    
} elsif ($species eq 'remanei') {
# Add new analysis/condition  objects to ~wormpub/DATABASES/remanei

  %expts = ( 
	    # SRX076136 # 454 library - will have to deal with the 454 quality file?
	    SRX101895 => ["RNASeq.remanei.Adult_male.Replicate1", 'phred', 'paired-end'],
	    SRX101894 => ["RNASeq.remanei.Adult_male.Replicate2", 'phred', 'single'],
	    SRX101893 => ["RNASeq.remanei.Adult_male.Replicate3", 'phred', 'paired-end'],
	    SRX101892 => ["RNASeq.remanei.Adult_male.Replicate4", 'phred', 'single'],
	    SRX101891 => ["RNASeq.remanei.Adult_male.Replicate5", 'phred', 'paired-end'],
	    SRX101890 => ["RNASeq.remanei.Adult_female.Replicate1", 'phred', 'paired-end'],
	    SRX101889 => ["RNASeq.remanei.Adult_female.Replicate2", 'phred', 'single'],
	    SRX101888 => ["RNASeq.remanei.Adult_female.Replicate3", 'phred', 'paired-end'],
	    SRX101887 => ["RNASeq.remanei.Adult_female.Replicate4", 'phred', 'single'],
	    SRX101886 => ["RNASeq.remanei.Adult_female.Replicate5", 'phred', 'paired-end'],
	    SRX101885 => ["RNASeq.remanei.L2_larva.Replicate1", 'phred', 'single'],
	    SRX101884 => ["RNASeq.remanei.L2_larva.Replicate2", 'phred', 'paired-end'],
	    SRX101883 => ["RNASeq.remanei.L2_larva.Replicate3", 'phred', 'paired-end'],
	    SRX101882 => ["RNASeq.remanei.early_embryo.Replicate1", 'phred', 'paired-end'],
	    SRX101881 => ["RNASeq.remanei.early_embryo.Replicate2", 'phred', 'paired-end'],
	    SRX101880 => ["RNASeq.remanei.L4_larva.Replicate1", 'phred', 'paired-end'],
	    SRX101879 => ["RNASeq.remanei.L4_larva.Replicate2", 'phred', 'paired-end'],
	    SRX052083 => ["RNASeq.remanei.L4_larva", 'phred', 'single'],
	    SRX052082 => ["RNASeq.remanei.L2_larva", 'phred', 'paired-end'],
	     
	   );

} elsif ($species eq 'briggsae') {
# Add new analysis/condition  objects to ~wormpub/DATABASES/briggsae
  %expts = (
	    SRX103668 => ["RNASeq.briggsae.Adult.Replicate1", 'phred', 'single'],
	    SRX103667 => ["RNASeq.briggsae.Adult.Replicate2", 'phred', 'paired-end'],
	    SRX103666 => ["RNASeq.briggsae.L2_larva.Replicate1", 'phred', 'single'],
	    SRX103656 => ["RNASeq.briggsae.early_embryo.Replicate1", 'phred', 'single'],
	    SRX103655 => ["RNASeq.briggsae.early_embryo.Replicate2", 'phred', 'paired-end'],
	    SRX100089 => ["RNASeq.briggsae.L4_larva.Replicate1", 'phred', 'paired-end'],
	    SRX100088 => ["RNASeq.briggsae.L2_larva.Replicate2", 'phred', 'paired-end'],
	    SRX100087 => ["RNASeq.briggsae.early_embryo.Replicate3", 'phred', 'paired-end'],
	    SRX100086 => ["RNASeq.briggsae.Adult.Replicate3", 'phred', 'paired-end'],
	    SRX089120 => ["RNASeq.briggsae.Adult.Replicate4", 'phred', 'paired-end'],
	    SRX089119 => ["RNASeq.briggsae.Adult.Replicate5", 'phred', 'paired-end'],
	    SRX089070 => ["RNASeq.briggsae.early_embryo.Replicate4", 'phred', 'paired-end'],
	    SRX089069 => ["RNASeq.briggsae.early_embryo.Replicate5", 'phred', 'paired-end'],
	    SRX053351 => ["RNASeq.briggsae.all_stages", 'phred', 'single'],
	    SRX052081 => ["RNASeq.briggsae.L4_larva", 'phred', 'single'],
	    SRX052079 => ["RNASeq.briggsae.L2_larva", 'phred', 'single'],
	   );

} elsif ($species eq 'brenneri') {
# Add new analysis/condition  objects to ~wormpub/DATABASES/brenneri
  %expts = (
	    SRX103657 => ["RNASeq.brenneri.Adult_male.Replicate1", 'phred', 'paired-end'],
	    SRX103654 => ["RNASeq.brenneri.Adult_female.Replicate1", 'phred', 'paired-end'],
	    SRX100780 => ["RNASeq.brenneri.Adult_male.Replicate2", 'phred', 'paired-end'],
	    SRX100779 => ["RNASeq.brenneri.Adult_male.Replicate3", 'phred', 'single'],
	    SRX100778 => ["RNASeq.brenneri.Adult_male.Replicate4", 'phred', 'single'],
	    SRX100777 => ["RNASeq.brenneri.Adult_male.Replicate5", 'phred', 'paired-end'],
	    SRX100776 => ["RNASeq.brenneri.Adult_female.Replicate2", 'phred', 'paired-end'],
	    SRX100775 => ["RNASeq.brenneri.Adult_female.Replicate3", 'phred', 'single'],
	    SRX100774 => ["RNASeq.brenneri.Adult_female.Replicate4", 'phred', 'single'],
	    SRX100773 => ["RNASeq.brenneri.Adult_female.Replicate5", 'phred', 'paired-end'],
	    SRX100772 => ["RNASeq.brenneri.L4_larva.Replicate1", 'phred', 'paired-end'],
	    SRX100771 => ["RNASeq.brenneri.L4_larva.Replicate2", 'phred', 'paired-end'],
	    SRX100770 => ["RNASeq.brenneri.early_embryo.Replicate1", 'phred', 'paired-end'],
	    SRX100769 => ["RNASeq.brenneri.early_embryo.Replicate2", 'phred', 'paired-end'],
	   );


} elsif ($species eq 'japonica') {

  # NEED TO ADD THESE ANALYSIS OBJECTS TO ~wormpub/DATABASES/japonica WHEN IT HAS ITS models.wrm UPDATED

  %expts = ( SRX100095 => ["RNASeq_Hillier.japonica.early_embryo", 'phred', 'paired-end'],
	     SRX100094 => ["RNASeq_Hillier.japonica.adult_female", 'phred', 'paired-end'],
	     SRX100093 => ["RNASeq_Hillier.japonica.adult_male", 'phred', 'paired-end'],
	     SRX100092 => ["RNASeq_Hillier.japonica.L4_larva", 'phred', 'paired-end'],
	     SRX100091 => ["RNASeq_Hillier.japonica.L4_larva", 'phred', 'paired-end'],
	     SRX100090 => ["RNASeq_Hillier.japonica.L2_larva", 'phred', 'paired-end'],
	   );

} else {
  $log->log_and_die("Unkown species: $species\n");
}



##########################################
# Set up database paths                  #
##########################################

my $script = "align_RNASeq.pl";
my $lsf;
my $store_file;
my $scratch_dir;
my $job_name;
$wormbase->checkLSF;
$lsf = LSF::JobManager->new();
$store_file = $wormbase->build_store; # make the store file to use in all wormpub perl commands
$scratch_dir = $wormbase->logs;
$job_name = "worm_".$wormbase->species."_RNASeq";


#my $RNASeqDir   = $wormbase->rnaseq;
my $RNASeqBase   = "/nfs/nobackup/ensembl_genomes/wormbase/BUILD/RNASeq/$species";
my $RNASeqSRADir    = "$RNASeqBase/SRA";
my $RNASeqGenomeDir = "$RNASeqBase/Genome";
chdir $RNASeqSRADir;


my @SRX = keys %expts;
  
if (!$expt) {

  if (!$check) {

    # delete results from the previous run
    foreach my $SRX (@SRX) {
      chdir "$RNASeqSRADir/$SRX";
      $wormbase->run_command("rm -rf tophat_out/", $log) unless ($noalign);
      $wormbase->run_command("rm -rf cufflinks/genes.fpkm_tracking", $log);
      $wormbase->run_command("rm -rf cufflinks/isoforms.fpkm_tracking", $log);
      $wormbase->run_command("rm -rf TSL/TSL_evidence.ace", $log) if ($tsl);
      $wormbase->run_command("rm -rf Introns/Intron.ace", $log);
    }

    my @chrom_files;
    if ($wormbase->assembly_type eq 'contig') {
      @chrom_files = ("${species}_masked.dna");
    } else {
      @chrom_files = $wormbase->get_chromosome_names('-prefix' => 1, '-mito' => 1);
    }

    if (! $noalign) { # only create the genome sequence index if the assembly changed
      # Build the bowtie index for the reference sequences by:
      mkdir "$RNASeqGenomeDir", 0777;
      mkdir "$RNASeqGenomeDir/reference-indexes/", 0777;
      chdir "$RNASeqGenomeDir/reference-indexes/";
      my $G_species = $wormbase->full_name('-g_species' => 1);
      unlink glob("${G_species}*");
      
      foreach my $chrom_file (@chrom_files) {
	$wormbase->run_command("rm -f ./$chrom_file", $log);
	if ($wormbase->assembly_type eq 'chromosome') {$chrom_file .= '_masked.dna'} # changes the contents of @chrom_files
	my $source_file = "${database}/CHROMOSOMES/${chrom_file}";
	if (-e $source_file) {
	  my $copy_cmd = "cp $source_file .";
	  $status = $wormbase->run_command($copy_cmd, $log);

	} elsif (-e "${source_file}.gz") {
	  my $copy_cmd = "gunzip -c ${source_file}.gz > ./$chrom_file";
	  $status = $wormbase->run_command($copy_cmd, $log);

	} else {
	  $log->log_and_die("Can't locate the repeat-masked genome file: $source_file\n");
	}
      }
#      my $bowtie_cmd = "bsub -I -M 4000000 -R \"select[mem>4000] rusage[mem=4000]\" /software/worm/bowtie/bowtie-build " . (join ',', @chrom_files) . " $G_species";
      my $bowtie_cmd = "bsub -I -M 4000000 -R \"select[mem>4000] rusage[mem=4000]\" $Software/bowtie/bowtie-build " . (join ',', @chrom_files) . " $G_species";
      $status = $wormbase->run_command($bowtie_cmd, $log);
      if ($status != 0) {  $log->log_and_die("Didn't create the bowtie indexes $RNASeqGenomeDir/reference-indexes/${G_species}.*\n"); }
    }

    # make the file of splice junctions:    
    # CDS, Pseudogene and non-coding-transcript etc introns
    # allows:
    # Coding_transcript Transposon_CDS Pseudogene tRNAscan-SE-1.23 Non_coding_transcript ncRNA Confirmed_cDNA Confirmed_EST Confirmed_UTR
    # Note: the curated CDS model introns are also added to this file
    
    my $splice_juncs_file = "$RNASeqGenomeDir/reference-indexes/splice_juncs_file";
    unless ($norawjuncs) {
      unlink $splice_juncs_file;
      unlink "${splice_juncs_file}.tmp";
      my $splice_junk_cmd = "grep -h intron ${database}/GFF_SPLITS/*.gff | egrep 'curated|Coding_transcript|Transposon_CDS|Pseudogene|tRNAscan-SE-1.23|Non_coding_transcript|ncRNA Confirmed_cDNA|Confirmed_EST|Confirmed_UTR' | awk '{OFS=\"\t\"}{print \$1,\$4-2,\$5,\$7}' > ${splice_juncs_file}.tmp";
      $status = $wormbase->run_command($splice_junk_cmd, $log);
      
      $status = $wormbase->run_command("sort -u ${splice_juncs_file}.tmp > $splice_juncs_file", $log);
      if ($status != 0) {  $log->log_and_die("Didn't create the splice_juncs_file\n"); }
    }
    

    # Make a GTF file of current transcripts. Used by tophat and cufflinks.
    my $gtf_file = "$RNASeqGenomeDir/transcripts.gtf";
    unless ($nogtf) {
      unlink $gtf_file;
      my $scripts_dir = $ENV{'CVS_DIR'};
      my $err = "$scratch_dir/make_GTF_transcript.pl.lsf.err";
      my $out = "$scratch_dir/make_GTF_transcript.pl.lsf.out";
      my @bsub_options = (-e => "$err", -o => "$out");
      push @bsub_options, (-q =>  "long",
			   #-F =>  "100000000", # there is no file size limit in Sanger LSF - don't impose one - keep this commented out
			   -M =>  "14000000", 
			   -R => "\"select[mem>14000] rusage[mem=14000]\"", 
			   -J => $job_name);
      my $cmd = "make_GTF_transcript.pl -database $database -out $gtf_file -species $species";
      $log->write_to("$cmd\n");
      $cmd = $wormbase->build_cmd($cmd);
      $lsf->submit(@bsub_options, $cmd);
      $lsf->wait_all_children( history => 1 );
      $log->write_to("The GTF file is written.\n");
      for my $job ( $lsf->jobs ) {
	if ($job->history->exit_status ne '0') {
	  $log->write_to("Job $job (" . $job->history->command . ") exited non zero: " . $job->history->exit_status . "\n");
	}
      }
      $lsf->clear;

      if ($species eq 'elegans') {
	$wormbase->check_file($gtf_file, $log,
			      minsize => 29900000,
			     );
      }
    }



    # now index this transcriptome GTF file with an initial run of tophat and a dummy small fastq read file
    unless ($nogtf) {
      # we don't want this initial run to do anything more than just index the GTF file
      # so we we create a dummy fastq file with 60 reads in
      # 60 reads is about the minimum that tophat will run successfully with
      chdir "$RNASeqGenomeDir"; # make the dummy tophat output in this directory
      unlink glob("transcriptome-gtf*");
      $status = $wormbase->run_command("rm -rf tophat_out/", $log);
      my $dummy_file = "$RNASeqGenomeDir/dummy.fastq";
      open (DUMMY, "> $dummy_file") || die "";
      print DUMMY '@IL21_2586:2:1:5:1516/1
CTCATCAATNATACAATCAACATATAGTNNTNCAAANNTTNAATAANNNNNNNT
+
EEEEFCDE7%6EDCDCD6@B7;DDE776%%6%7;D7%%73%77B77%%%%%%%6
@IL21_2586:2:1:5:1477/1
GATCGGAAGNGCGGTTCAGCAGGAATGCNNANATCGNNAGNGCGGTNNNNNNNG
+
@>@<>>==1%5=>:=>;>:.>=06>=11%%7%5=/1%%5,%5(775%%%%%%%/
@IL21_2586:2:1:5:979/1
GATCGGAAGNGCGGTTCAGCAGGAATGCNNANATAGNNAGGGCGGTNNANNNNG
+
@<>4>>521%50<=<>2/5,56..2<1*%%%%,-%.%%*((%*(/5%%%%%%%,
@IL21_2586:2:1:5:1787/1
CGGTTCAGCNGGAATGCCGAGATCGGAANNGNGGTTNNGCAGGAATNNCNNNNA
+
>>>??=?=5%5==?>=<:9?=>><9=>5%%1%5)55%%11575/<5%%7%%%%,
@IL21_2586:2:1:5:415/1
ATCGCTGAANAAATACGTTTTGCTGTTGANTNATACANTGCCAGTANNCNNNCG
+
DD>D9DCC6%6CCDC8;DDDA>9C5@D=%%7%6D?61%6)603)14%%-%%%.%
@IL21_2586:2:1:5:1362/1
GTTGTATATNATTGCAAATGCCTCAAGAANTNAAAAANAATACTGCNNANNNTA
+
DDC@@DDD7%7>@@=BDDA=85?1><D6.%3%758:6%7:<78:3(%%7%%%6-
@IL21_2586:2:1:5:1313/1
GATCGGAAGNGCGGTTCAGCAGGAATGCCNANATCGGNAGAGCTGTNNANNNGT
+
@>>=>><=5%5=>>=>===2=<+8;8,5*%5%52-(5%5**,,&*5%%5%%%,%
@IL21_2586:2:1:5:712/1
CCTCTGGCTNCCGTAACTTTCGTTTTCCANANTTTCANGCTTTGAANNANNNAG
+
DDDDDCDC6%3>DCA>7ADA;<CDDB>A6%3%6DC66%)6ED@0<6%%7%%%17
@IL21_2586:2:1:5:1564/1
GATCGGAAGNGCGGTTCAGCAGGAATGCCNANATCGGAAGAGCGGTNNANNNGG
+
>>@=>>=>/%/9>>=?=>><>=>=>==0*%1%5</50(5.:5*&55%%1%%%.,
@IL21_2586:2:1:5:1136/1
CCCATCTTCNACTTCCAGAGGTTCAGCAGNTNTGATCCGAGCAAACNNANNNGT
+
>A@@A>AA6%6>?A><>>>=:>>>?>>>-%6%6==>;37>69>>>)%%6%%%67
@IL21_2586:2:1:5:1023/1
CCCGTGTTGNGTCAAATTAAGCCGCAGGCNCNACTCCTGGTGGTGCNNTNNNGT
+
>>>??8>?5%,?>>>:@?>?6<=>=</51%5%5<?54575>(4575%%5%%%5>
@IL21_2586:2:1:5:789/1
CAGAAAAGCNATAGCAGCAAATGCAGCAANANTATTAAGTTAAACTNTCNNTTT
+
DFEEFAD96%-D?B=>:7@>CA6666-66%4%6>C=;3=:@@99>6%67%%6@=
@IL21_2586:2:1:5:496/1
GATCGGAAGNGCGGTTCAGCAGGAATGCCNANATCGGAAGAGGGGTNCANNGGG
+
A=>2>>515%5-,><</::1+6=,33;,%%%%%0,,6778731.0.%%%%%%*(
@IL21_2586:2:1:5:1831/1
CCCGTCCGTNCCTCTTAACCATTATCTCANANCAGAAAACCAACAAAATNNAAC
+
DDDDDDA=6%6DDDDDDDBABDDDDDDB6%6%6@4@A=6-6?@=A?466%%67;
@IL21_2586:2:1:5:243/1
GACGTGGGCNGCATCAATCTGACGGATGANGNATGGCTCTTCAATCCCANNTGG
+
?><><==;,%16=><<=>0>=82<=5=;&%.%5=75-?5=5%,;5**,+%%5*,
@SRR350954.945 BB0BG1A
CCTTCATTCTCAATGATACGAGTCTTGATCTCTTCGTTAACTCGATCTTCCTCC
+SRR350954.945 BB0BG1A
IIIIIIIIIIIIIIIIIIIIIIIIIIIFIIIIIIIIIIIIIIIIIIIIIIIIII
@SRR350954.1072 BB0BG1
AAATATCTTGGCGGTTTTGAAGACACCGGCGTCGCATCACTGCTACGTGGCTTT
+SRR350954.1072 BB0BG1
IIIIIIIIIIIIIIIIIIIGIIHIFIIIIIIIIIHFHIIIGHIIIIIHHIIIII
@SRR350954.1265 BB0BG1
CTCTCCGCAAAGTTCCACAATACGACAAATTCATCACCGAGCGATTCGAGCGAT
+SRR350954.1265 BB0BG1
IIIIIIIIHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIGIGIDHIIIIIGIGI
';

      close(DUMMY);
      my $err = "$scratch_dir/align_RNASeq.pl.lsf.initial.err";
      my $out = "$scratch_dir/align_RNASeq.pl.lsf.initial.out";
      my @bsub_options = (-e => "$err", -o => "$out");
      push @bsub_options, (-q =>  "normal",
			   #-F =>  "100000000", # there is no file size limit in Sanger LSF - don't impose one - keep this commented out
			   -M =>  "4000000", 
			   -R => "\"select[mem>4000 && tmp>10000] rusage[mem=4000]\"", 
			   -J => "tophat_dummy_$species");
      my $gtf = "--GTF $RNASeqGenomeDir/transcripts.gtf";
      my $gtf_index = "--transcriptome-index $RNASeqGenomeDir/transcriptome-gtf";
      my $G_species = $wormbase->full_name('-g_species' => 1);
      my $cmd = "$Software/tophat/tophat $gtf $gtf_index $RNASeqGenomeDir/reference-indexes/$G_species $dummy_file";
      $log->write_to("$cmd\n");
      $lsf->submit(@bsub_options, $cmd);
      
      
      $lsf->wait_all_children( history => 1 );
      $log->write_to("The initial tophat job to index the GTF file has completed!\n");
      for my $job ( $lsf->jobs ) {
	if ($job->history->exit_status ne '0') {
	  $log->write_to("Job $job (" . $job->history->command . ") exited non zero: " . $job->history->exit_status . " Try making a larger dummy.fastq file?\n");
	}
      }
      $lsf->clear;
    }
  }

  # run tophat against a few experiments at a time - we don't want to
  # have too many fastq files ungzipped at a time otherwise we may run
  # out of quota filespace.

  my @args;
  my $count=0;
  foreach my $arg (@SRX) {

    if ($check) {

      # we have had a problem with incomplete intron.ace files
      if (-e "$arg/Introns/Intron.ace" && ! -z "$arg/Introns/Intron.ace") {
	my $last_line = `tail -1 "$arg/Introns/Intron.ace"`;
	my $last_char = chop $last_line;
	if ($last_char ne "\n") {
	  $wormbase->run_command("rm -rf $arg/Introns/Intron.ace", $log);
	  $log->write_to("Incomplete $arg intron file removed\n");
	}
      }

      if (-e "$RNASeqSRADir/$arg/tophat_out/accepted_hits.bam" &&
	  -e "$RNASeqSRADir/$arg/cufflinks/genes.fpkm_tracking" &&
	  (-e "$RNASeqSRADir/$arg/TSL/TSL_evidence.ace" || !$tsl) &&
	  -e "$RNASeqSRADir/$arg/Introns/Intron.ace") {next}
    }

    push @args, $arg;
    if (($count % 15) == 14) {
      $log->write_to("Running alignments on: @args\n");
      print "Running alignments on: @args\n";
      &run_align($check, $noalign, $tsl, @args);
      @args=();
      $count = 0;
    }
    $count++;
  }
  if (@args) {
    $log->write_to("Running alignments on: @args\n");
    print "Running alignments on: @args\n";
    &run_align($check, $noalign, $tsl, @args);
  }
  
  ####################################################################################
  # all analyses have now run - munge results and put them where they are expected
  ####################################################################################
  
  # now sleep for 5 minutes to give the lustre file-system a chance to
  # sort out the files we have just written - was having intermittant
  # problems with one or two of the intron.ace files as if it hadn't
  # finished flushing the file buffers before the file was read in the
  # next sections.
  sleep 300;

  # now write out a table of what has worked
  # and make the expresssion tarball for Wen to put into SPELL
  # and write out .ace files for the FPKM expression levels of genes, transcripts, pseudogenes and CDSs


  $log->write_to("\nResults\n");
  $log->write_to("-------\n\n");


  chdir $RNASeqSRADir;

  # get the Life_stages of the Condition objects (same name as the Analysis objects)
  my %life_stages = get_life_stage_from_conditions();

  # get the valid CDS and Pseudogenes IDs
  my %CDSs = get_CDSs();
  my %Pseudogenes = get_Pseudogenes();

  # open a .ace file to hold the FPKM expression levels of genes, transcripts and CDSs
  my $misc_dynamic = $wormbase->misc_dynamic;
  open (EXPRACE, "> $misc_dynamic/RNASeq_expression_levels_${species}.ace") || $log->log_and_die("Can't open $misc_dynamic/RNASeq_expression_levels_${species}.ace\n");

  foreach my $SRX (@SRX) {
    $log->write_to("$SRX");
    if (-e "$SRX/tophat_out/accepted_hits.bam") {$log->write_to("\ttophat OK");} else {{$log->write_to("\ttophat ERROR");}}
    if (-e "$SRX/cufflinks/genes.fpkm_tracking") {$log->write_to("\tcufflinks OK");} else {$log->write_to("\tcufflinks ERROR");}
    if ($tsl) {if (-e "$SRX/TSL/TSL_evidence.ace") {$log->write_to("\tTSL OK");} else {$log->write_to("\tTSL ERROR");}}

    # we have had a problem with incomplete intron.ace files
    if (-e "$SRX/Introns/Intron.ace" && ! -z "$SRX/Introns/Intron.ace") {
      my $last_line = `tail -1 "$SRX/Introns/Intron.ace"`;
      my $last_char = chop $last_line;
      if ($last_char ne "\n") {
	$wormbase->run_command("rm -rf $SRX/Introns/Intron.ace", $log);
	$log->write_to("\tincomplete intron file removed");
      }
    }

    if (-e "$SRX/Introns/Intron.ace") {$log->write_to("\tIntrons OK");} else {$log->write_to("\tintron ERROR");}
    $log->write_to("\n");

    my $analysis = $expts{$SRX}->[0];
    my $life_stage = $life_stages{$analysis};
    if (!defined $life_stage) {$log->write_to("WARNING: no Condition object found for Analysis object $analysis\n");}
    my %CDS_SRX;
    my %Pseudogene_SRX;

# Wen says: "Gary, the file is good, I just wonder what should we do
# with the "0" values. SPELL data are all log2 transformed. That is
# how they are stored in mysql and on the website there are also
# options for users to disply linear or log2 transformed data.  I will
# not be able to handle the "0" in SPELL. Can we replace "0" with a
# very small value, such as 0.0000000001?  If you have a lot tables,
# you tarzip them together and place them on a ftp site for me to
# download."
    
    if (-e "$SRX/cufflinks/genes.fpkm_tracking") {
      open (EXPR, "<$SRX/cufflinks/genes.fpkm_tracking") || $log->log_and_die("Can't open $SRX/cufflinks/genes.fpkm_tracking\n");
      open (EXPROUT, ">$SRX.out") || $log->log_and_die("Can't open $SRX.out\n");
      while (my $line = <EXPR>) {
	my @f = split /\s+/, $line;
	if ($f[0] eq 'tracking_id') {next;}
	if ($f[0] =~ /CUFF/) {$log->log_and_die("Cufflinks for $SRX has failed to put the Gene IDs in the output file - problem with the GTF file?\n");}
	my $gene_expr = $f[10];
	if ($f[9] == 0) {$gene_expr = 0.0000000001}
	if ($f[12] eq "OK" || $f[12] eq "LOWDATA") {
	  # print the file for Wen's SPELL data
	  print EXPROUT "$f[0]\t$gene_expr\n";
	  # and print to the Gene model ace file
	  if (defined $life_stage) {
	    print EXPRACE "\nGene : \"$f[0]\"\n";
	    print EXPRACE "RNASeq_FPKM  \"$life_stage\"  \"$f[9]\"  From_analysis \"$analysis\"\n";
	  }
	}
      }
      close(EXPROUT);
      close (EXPR);
    }

    my $number_of_complaints = 5; #complain only 5 times about Sequences names not in the normal format, otherwise we are swamped by warnings in remanei
    # now get the isoform (coding transcript and CDS) expression values and write out to the ace file
    open (EXPR, "<$SRX/cufflinks/isoforms.fpkm_tracking") || $log->log_and_die("Can't open $SRX/cufflinks/isoforms.fpkm_tracking\n");
    while (my $line = <EXPR>) {
      my @f = split /\s+/, $line;
      if ($f[0] eq 'tracking_id') {next;}
      if ($f[0] =~ /CUFF/) {$log->log_and_die("Cufflinks for $SRX has failed to put the Transcript IDs in the output file - problem with the GTF file?\n");}
      if ($f[12] eq "OK" || $f[12] eq "LOWDATA") {
	if (defined $life_stage) {
	  my ($sequence_name) = ($f[0] =~ /(^\S+?\.\d+[a-z]?)/);
	  if (!defined $sequence_name) {
	    # complain up to 5 times if it is not even like 'T27B1.t1' - a tRNA gene
	    if ($f[0] !~ /(^\S+?\.t\d+)/ && $number_of_complaints-- > 0) {print "Sequence name $f[0] is not in the normal format\n";} 
	    $sequence_name = $f[0]; # use the name as given
	  }

	  # Pseudogenes don't have Transcript objects in acedb, only do CDSs
	  if (exists $CDSs{$sequence_name}) {
	    print EXPRACE "\nTranscript : \"$f[0]\"\n";
	    print EXPRACE "RNASeq_FPKM  \"$life_stage\"  \"$f[9]\"  From_analysis \"$analysis\"\n";
	  }

	  # sum the Pseudogene and CDS scores from their isoforms
	  if (exists $CDSs{$sequence_name}) {
	    $CDS_SRX{$sequence_name} += $f[9]; # sum the FPKMs for this CDS
	  } elsif (exists $Pseudogenes{$sequence_name}) {
	    $Pseudogene_SRX{$sequence_name} += $f[9]; # sum the FPKMs for this Pseudogene
	  }
	}
      }
    }
    close (EXPR);

    # now get the CDSs and Pseudogenes that were seen and output their FPKMs to the ace file
    foreach my $CDS (keys %CDS_SRX) {
      print EXPRACE "\nCDS : \"$CDS\"\n";
      print EXPRACE "RNASeq_FPKM  \"$life_stage\"  \"$CDS_SRX{$CDS}\"  From_analysis \"$analysis\"\n";
    }
    foreach my $Pseudogene (keys %Pseudogene_SRX) {
      print EXPRACE "\nPseudogene : \"$Pseudogene\"\n";
      print EXPRACE "RNASeq_FPKM  \"$life_stage\"  \"$Pseudogene_SRX{$Pseudogene}\"  From_analysis \"$analysis\"\n";
    }


  } # foreach SRX
  $log->write_to("\n");

  close(EXPRACE);

  # make a tarball
  $status = $wormbase->run_command("tar cf expr.tar *.out", $log);
  unlink "expr.tar.gz";
  $status = $wormbase->run_command("gzip -f expr.tar", $log);
  my $autoace = $wormbase->autoace;
  $status = $wormbase->run_command("cp expr.tar.gz $autoace", $log); # this will probably be changed to the autoace/OUTPUT directory soon

  
  # make the ace file of RNASeq spanned introns to load into acedb
  my $splice_file = $wormbase->misc_dynamic."/RNASeq_splice_${species}.ace";
  chdir $RNASeqSRADir;
  $status = $wormbase->run_command("rm -f $splice_file", $log);
  $status = $wormbase->run_command("cat */Introns/virtual_objects.${species}.RNASeq.ace > $splice_file", $log);
  $status = $wormbase->run_script("acezip.pl -file $splice_file", $log);
  $status = $wormbase->run_command("cat */Introns/Intron.ace > ${splice_file}.tmp", $log);
  $status = $wormbase->run_script("acezip.pl -file ${splice_file}.tmp", $log);
  # flatten the results of all libraries at a position into one entry
  open (FEAT, "< ${splice_file}.tmp") || $log->log_and_die("Can't open file ${splice_file}.tmp\n");
  open (FLAT, ">> $splice_file") || $log->log_and_die("Can't open file $splice_file\n");
  my %splice;
  while (my $line = <FEAT>) {
    if ($line =~ /^Feature_data/ || $line =~ /^\s*$/) { # new clone
      foreach my $start (keys %splice) {
	foreach my $end (keys %{$splice{$start}}) {
	  my $total= 0;
	  my $string = "";
	  foreach my $library (keys %{$splice{$start}{$end}}) {
	    my $value = $splice{$start}{$end}{$library};
	    $total += $value;
	    $string .= "$library $value "
	  }
	  # filter out any spurious introns with only 1 or 2 reads
	  if ($total > 2) {print FLAT "Feature RNASeq_splice $start $end $total \"$string\"\n";}
	}
      }
      # reset things for the new clone
      %splice = ();
      print FLAT $line;
    } else {
      my @feat = split /\s+/, $line;
      $splice{$feat[2]}{$feat[3]}{$feat[5]} = $feat[4];
    }
  }
  # and do the last clone
  foreach my $start (keys %splice) {
    foreach my $end (keys %{$splice{$start}}) {
      my $total= 0;
      my $string = "";
      foreach my $library (keys %{$splice{$start}{$end}}) {
	my $value = $splice{$start}{$end}{$library};
	$total += $value;
	$string .= "$library $value "
      }
	  # filter out any spurious introns with only 1 or 2 reads
      if ($total > 2) {print FLAT "Feature RNASeq_splice $start $end $total \"$string\"\n";}
    }
  }
  close(FLAT);
  close(FEAT);
  $status = $wormbase->run_command("rm -f ${splice_file}.tmp", $log);


} else { # we have a -expt parameter
  
  &run_tophat($check, $noalign, $tsl, $expt, $solexa, $illumina, $paired);
}


$log->mail();
print "Finished.\n" if ($verbose);

exit(0);


###################################
##########  SUBROUTINES  ##########
###################################

# submit this script to LSF on the farm2 head node for each experiment
# and wait for the jobs to complete.  Check that the job worked.

sub run_align {
  my ($check, $noalign, $tsl, @args) = @_;

  foreach my $arg (@args) {

    # pull over the SRA files and unpack them to make a fastq file
    # the aspera commmand to pull across the files only works on the farm2-login head node
    if ((!$check && !$noalign) || !-e "tophat_out/accepted_hits.bam") {
      get_SRA_files($arg);
    }

    # "the normal queue can deal with memory requests of up to 15 Gb, but 14 Gb is better" - Peter Clapham, ISG
    my $err = "$scratch_dir/align_RNASeq.pl.lsf.${arg}.err";
    my $out = "$scratch_dir/align_RNASeq.pl.lsf.${arg}.out";
    my @bsub_options = (-e => "$err", -o => "$out");
    push @bsub_options, (-q =>  "long",
                         #-F =>  "100000000", # there is no file size limit in Sanger LSF - don't impose one - keep this commented out
			 -M =>  "14000000", 
			 -R => "\"select[mem>14000] rusage[mem=14000]\"",
			 -J => $job_name);
    my $cmd = "$script -expt $arg";      # -expt is the parameter to make the script run an alignment and analysis on a dataset $arg
    if ($check) {$cmd .= " -check";}
    if ($noalign) {$cmd .= " -noalign";}
    if ($nogtf) {$cmd .= " -nogtf";}
    if ($tsl) {$cmd .= " -tsl";}
    if ($onestage) {$cmd .= " -onestage";}
    if ($database) {$cmd .= " -database $database";}
    if ($species) {$cmd .= " -species $species";}
    if ($norawjuncs) {$cmd .= " -norawjuncs";}
    if ($expts{$arg}[1] eq 'solexa') {$cmd .= " -solexa";}
    if ($expts{$arg}[1] eq 'illumina1.3') {$cmd .= " -illumina";}
    if ($expts{$arg}[2] eq 'paired-end') {$cmd .= " -paired";}
    $log->write_to("$cmd\n");
    $cmd = $wormbase->build_cmd($cmd);
    $lsf->submit(@bsub_options, $cmd);
  }
  
  $lsf->wait_all_children( history => 1 );
  $log->write_to("This set of Tophat jobs have completed!\n");
  for my $job ( $lsf->jobs ) {
    if ($job->history->exit_status ne '0') {
      $log->write_to("Job $job (" . $job->history->command . ") exited non zero: " . $job->history->exit_status . "\n");
    }
  }
  $lsf->clear;

  # clear up the fastq files
  if (!$keepfastq) {
    foreach my $arg (@args) {
      # if it is a SRA file, then we can safely  delete it and pull it over again next time
      if ($arg =~ /SRX/) {my $status = $wormbase->run_command("rm -rf $RNASeqSRADir/$arg/SRR", $log);}
    }
  }


}

#####################################################################
# run tophat on one experiment to create the accepted_hits.bam file
# this can fail through lack of memory if that happens, try running it
# again on the hugemem queue 

# ARG: list of names of SRX* directories in the $wormbase->rnaseq
# directory

# probably best to only give about five/ten input directories at a
# time otherwise your filesize quota may be exceeded

sub run_tophat {
  
  my ($check, $noalign, $tsl, $arg, $solexa, $illumina, $paired) = @_;
  
  chdir "$RNASeqSRADir/$arg";
  my $G_species = $wormbase->full_name('-g_species' => 1);

  my $cmd_extra = "";
  if ($solexa)   {$cmd_extra = "--solexa-quals"} 
  if ($illumina) {$cmd_extra = "--solexa1.3-quals"} 

  if ((!$check && !$noalign) || !-e "$RNASeqSRADir/$arg/tophat_out/accepted_hits.bam") {
    if (glob("$RNASeqSRADir/$arg/SRR/SRR*/*.lite.sra") ) {
      # unpack any .lite.sra file to make the .fastq files
      # now unpack the sra.lite files
      $log->write_to("Unpack the $arg .lite.sra file to fastq files\n");
      foreach my $dir (glob("$RNASeqSRADir/$arg/SRR/SRR*")) {
	chdir $dir;
	foreach my $srr (glob("*.lite.sra")) {
	  $log->write_to("Unpack $srr\n");
	  my $options = '';
	  if ($paired && $srr !~ /_1\./ && $srr !~ /_2\./) {$options='--split-files'} # specify that the file contains paired-end reads and should be split
#	  $status = $wormbase->run_command("/software/worm/sratoolkit/fastq-dump $options $srr", $log);
	  $status = $wormbase->run_command("$Software/sratoolkit/fastq-dump $options $srr", $log);
	  if ($status) {$log->log_and_die("Didn't unpack the fastq file successfully\n")}
	  $status = $wormbase->run_command("rm -f $srr", $log);
	}
      }
    }



    $wormbase->run_command("rm -rf tophat_out/", $log);

    chdir "$RNASeqSRADir/$arg";
    my @files = glob("SRR/*/*.fastq");
    my $joined_file = join ",", @files;
    $log->write_to("Have fastq files: @files\n");
    $log->write_to("Check read length in file: $files[0]\n");
    my $seq_length = get_read_length($files[0]);
    
    # do we have paired reads?
    my @files1 = sort glob("SRR/*/*_1.fastq"); # sort to ensure the two sets of files are in the same order
    my @files2 = sort glob("SRR/*/*_2.fastq");
    if ((@files1 == @files2) && @files2 > 0) {
      $log->write_to("Have paired-end files.\n");
      my $joined1 = join ",", @files1;
      my $joined2 = join ",", @files2;
      $joined_file = "$joined1 $joined2";
      print "Made paired-read joined files: $joined_file\n";
      # no longer need to set the inner-distance -r parameter
      # It is now set to 50 by default.
      #my $inner = 50;
      #$cmd_extra .= " -r $inner ";
    }

    # Is the read length less than 50?  
    # If so then we should set --segment-length to be less than half the
    # read length, otherwise the alignments are not done properly and we
    # don't get a junctions.bed file output
    my $segment_length = 25; # the default value
    if ($seq_length < 50) {
      $segment_length = int($seq_length / 2);
      $segment_length--; # we want it to be less than half the length
      $cmd_extra .= " --segment-length $segment_length ";
      # and the documentation now says to set the --segment-mismatches to 1 or 0
      $cmd_extra .= " --segment-mismatches 0 ";
    }


    $log->write_to("run tophat $joined_file\n");
    my $gtf_index = ''; # use the GTF index unless we are not using the GTF file
    $gtf_index = "--transcriptome-index $RNASeqGenomeDir/transcriptome-gtf" unless ($nogtf || $onestage);
    my $raw_juncs = ''; # use the EST raw junctions hint file unless we specify otherwise
    $raw_juncs = "--raw-juncs $RNASeqGenomeDir/reference-indexes/splice_juncs_file" unless $norawjuncs;
#    $status = $wormbase->run_command("/software/worm/tophat/tophat $cmd_extra --min-intron-length 30 --max-intron-length 5000 $gtf_index $raw_juncs $RNASeqGenomeDir/reference-indexes/$G_species $joined_file", $log);
    $status = $wormbase->run_command("$Software/tophat/tophat $cmd_extra --min-intron-length 30 --max-intron-length 5000 $gtf_index $raw_juncs $RNASeqGenomeDir/reference-indexes/$G_species $joined_file", $log);
    if ($status != 0) {  $log->log_and_die("Didn't run tophat to do the alignment successfully\n"); } # only exit on error after gzipping the files

  } else {
    $log->write_to("Check tophat_out/accepted_hits.bam: already done\n");
  }

##############################################################################
# now get the RPKM values
##############################################################################

# cufflinks
# use a GTF of our known transcripts - CDS, Pseudogene and non-coding-transcript
# to restrict the created gene structures to our own


# now run cufflinks
  if (!$check || !-e "cufflinks/genes.fpkm_tracking") {
    $log->write_to("run cufflinks\n");
    chdir "$RNASeqSRADir/$arg";
    mkdir "cufflinks", 0777;
    chdir "cufflinks";
    my $gtf = ''; # use the existing gene models unless we specify otherwise
    $gtf = "--GTF $RNASeqGenomeDir/transcripts.gtf" unless $nogtf;
    $status = $wormbase->run_command("$Software/cufflinks/cufflinks $gtf ../tophat_out/accepted_hits.bam", $log);
    if ($status != 0) {  $log->log_and_die("Didn't run cufflinks to get the isoform/gene expression successfully\n"); }
  } else {
    $log->write_to("Check cufflinks/genes.fpkm_tracking: already done\n");
  }


#############################################################
# now get the intron spans
#############################################################

  if (!$check || !-e "Introns/Intron.ace") {
    $log->write_to("run Introns\n");
    chdir "$RNASeqSRADir/$arg";
    mkdir "Introns", 0777;
    my $analysis = $expts{$arg}[0];
    my $status = &Intron_stuff($cmd_extra, $analysis);

    if ($status != 0) {  $log->log_and_die("Didn't run the Intron stuff successfully\n"); }
  } else {
    $log->write_to("Check Introns/Intron.ace: already done\n");
  }

#############################################################
# now get the TSL sites
#############################################################

# now run TSL stuff
  if ($tsl) {
    if (!$check || !-e "TSL/TSL_evidence.ace") {
      $log->write_to("run TSL\n");
      chdir "$RNASeqSRADir/$arg";
      mkdir "TSL", 0777;
      my $analysis = $expts{$arg}[0];
      $status = &TSL_stuff($analysis);
      
      if ($status != 0) {  $log->log_and_die("Didn't run the TSL stuff successfully\n"); }
    } else {
      $log->write_to("Check TSL/TSL_evidence.ace: already done\n");
    }
  }
}

#################################################################################################################
# take a quick look in the read files to determine the length of the reads
# read a hundred entries because the TSL stuff deals with reads of varying lengths, so we want a good sample
sub get_read_length {
  my ($filename) = @_;

  my $count = 100; # number to look at
  my $min_length = 1000000;

  open (PEEK, "<$filename") || $log->log_and_die("get_read_length() : cant open $filename\n");
  while (my $line = <PEEK>) {    # skip the title line
    $line = <PEEK>; # this is the sequence line
    chomp $line;
    if (length $line < $min_length) {$min_length = length $line}
    $line = <PEEK>; # second title line
    $line = <PEEK>; # quality line
    if (!$count--) {last}
  }
  close(PEEK);
  return $min_length;
}

#################################################################################################################
# pull over the SRR files from the Short Read Archive

sub get_SRA_files {
  my ($arg) = @_;

  # is this a SRA ID?
  if ($arg !~ /^SRX/) {return}

  my $status = $wormbase->run_command("rm -rf $RNASeqSRADir/$arg/SRR", $log);

  chdir "$RNASeqSRADir";

  get_SRX_file($arg);

  chdir "$RNASeqSRADir/$arg";
  mkdir "SRR", 0777;

  foreach my $dir (glob("SRR*")) {
    if ($dir eq "SRR") {next}
    $status = $wormbase->run_command("mv -f $dir $RNASeqSRADir/$arg/SRR", $log);
  }

}

#################################################################################################################
# this does the work of pulling the NCBI SRA file across

sub get_SRX_file {
  my ($arg) = @_;

  # get the name of the subdirectory in the SRA
  my ($dirbit) = ($arg =~ /SRX(\d\d\d)/);

  # the Aspera file transfer system is error prone and (in Sanger) only works when run on the farm2-login head node
  $log->write_to("Get the SRA files for $arg\n");
  my $failed = 1;
  while ($failed) {
    my $cmd = "~gw3/.aspera/connect/bin/ascp -k 1 -l 300M -QTr -i /net/nas17b/vol_homes/homes/gw3/.aspera/connect/etc/asperaweb_id_dsa.putty anonftp\@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByExp/litesra/SRX/SRX${dirbit}/$arg ./";
    print "cmd = $cmd\n";
    my $return_status = system($cmd);
    if ( ( $return_status >> 8 ) == 0 )  {$failed=0}
  }
}


#################################################################################################################
# gets the reads that don't align and which have a bit of TSL sequence
# on them - runs tophat using this library and parses out the TSL
# sites from the bam file then matches to the known TSL sites and
# writes evidence for the TSL Features.
#
# Assumes we have a TSL directory set up to get the results and that
# the current working directory is "$RNASeqSRADir/$arg"


# TSL sequences from 
# PLOS Genetics
# Nov 2006
# Operon Conservation and the Evolution of trans-Splicing in the Phylum Nematoda
# David B. Guiliano, Mark L. Blaxter
# http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.0020198
# Figure 4



sub TSL_stuff {

  my ($analysis) = @_;

  $coords = Coords_converter->invoke($database, 0, $wormbase);
  $mapper = Feature_mapper->new($database, undef, $wormbase);

  mkdir "TSL", 0777;
  mkdir "TSL-tmp", 0777;


  my $outfile="TSL1";
  foreach my $file (glob "SRR/*/*.fastq") {
    $log->write_to("Running get_TSL_RNASeq_reads.pl on $file\n");
    $wormbase->run_script("get_TSL_RNASeq_reads.pl -infile $file -outfile TSL-tmp/$outfile",$log);
    $outfile++;
  }

  # now concatenate all the TSL files to make one
  my $tsl_fastq = "TSL/sl.fq";
  $wormbase->run_command("cat TSL-tmp/* > $tsl_fastq", $log);
  $wormbase->run_command("rm -rf TSL-tmp", $log);


  # Is the read length of the TSL-stripped reads less than 50?  
  # If so then we should set --segment-length to be less than half the
  # read length, otherwise the alignments are not done properly and we
  # don't get a junctions.bed file output
  my $segment_length = 25; # the default value
  my $seq_length = get_read_length($tsl_fastq);
  my $cmd_extra='';
  if ($seq_length < 50) {
    $segment_length = int($seq_length / 2);
    $segment_length--; # we want it to be less than half the length
    $cmd_extra .= " --segment-length $segment_length ";
    # and the documentation now says to set the --segment-mismatches to 1 or 0
    $cmd_extra .= " --segment-mismatches 0 ";
  }

  # now run tophat on this set of TSL reads
  my $G_species = $wormbase->full_name('-g_species' => 1);

  $log->write_to("run tophat with $tsl_fastq\n");

  $status = $wormbase->run_command("$Software/tophat/tophat $cmd_extra --output-dir TSL $RNASeqGenomeDir/reference-indexes/$G_species $tsl_fastq", $log);

#  $log->write_to("remove TSL fastq files\n");
#  $wormbase->run_command("rm -f $tsl_fastq", $log);
  if ($status != 0) { return $status;  }

  # now parse the results looking for TSL sites
  my $TSLfile = "TSL/accepted_hits.bam";
  my $samtools = "$Software/samtools/samtools view $TSLfile";

  my %results; # hash of TSL sites found by RNAseq alignment
  $log->write_to("get the TSL sites found by alignment\n");
  open (HITS, "$samtools |") || $log->log_and_die("can't run samtools\n");
  while (my $line = <HITS>) {
    my $sense = '+';
    my ($id, $flags, $chrom, $pos, $seq) = ($line =~ /^(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/);
    my ($tsl) = ($id =~ /(SL\d)/);
    if ($flags & 0x10) { # find and deal with reverse alignments
      $sense = '-';
      $pos = $pos + (length $seq) - 2; # go past the 'A' of the 'AG' we added on to ensure we hit a splice site
    } else {
      $pos = $pos + 1; # go past the 'A' of the 'AG' we added on to ensure we hit a splice site
    }
    if (!defined $chrom) {print "undefined chrom in ID $id\n";}
    if (!defined $pos) {print "undefined pos in ID $id\n";}
    if (!defined $sense) {print "undefined sense in ID $id\n";}
    if (!defined $tsl) {print "undefined tsl in ID $id\n";}
    $results{$chrom}{$pos}{$sense}{$tsl}++;
    #print "store: ${chrom} ${pos} ${sense} ${tsl} value: $results{$chrom}{$pos}{$sense}{$tsl} $seq\n";
  }
  close(HITS);

  # now check that the TSL sites match known TSL Feature objects and add evidence to the Features
  $log->write_to("check that the TSL sites match known TSL Feature objects\n");
  my %found; # hash of TSL Feature object IDs with evidence from RNASeq - holds no. of supporting reads
  my %not_found; # hash of TSL Feature object IDs with no evidence from this RNASeq experiment

  print "Using database $database to look up Features\n";
  my $table_def = &write_feature_def;
  my $table_query = $wormbase->table_maker_query($database, $table_def);
  while(<$table_query>) {
    chomp;
    s/\"//g;  #remove "
    next if (/acedb/ or /\/\//);
    my @data = split("\t",$_);
    my ($feature, $clone, $start, $end, $method) = @data;
    if (!defined $feature) {next}
    if (!defined $clone) {$log->write_to("Feature $feature does not have a position mapped\n"); next}
    if (!defined $start) {$log->write_to("Feature $feature does not have a position mapped\n"); next}
    my ($TSL_chrom, $TSL_coord);
    ($TSL_chrom, $TSL_coord) = $coords->Coords_2chrom_coords( $clone, $start );
    if ($start < $end) { # forward sense
      print "looking for: ${TSL_chrom} ${TSL_coord} + ${method}\n";

      if (exists $results{$TSL_chrom}{$TSL_coord}{'+'}{$method}) {
	$found{$feature} = $results{$TSL_chrom}{$TSL_coord}{'+'}{$method};
	delete $results{$TSL_chrom}{$TSL_coord}{'+'}{$method};
	print "found $feature +\n";
      } else {
	$not_found{$feature} = 1;
	print "NOT found $feature +\n";
      }
    } else { # reverse sense
      print "looking for: ${TSL_chrom} ${TSL_coord} - ${method}\n";
      if (exists $results{$TSL_chrom}{$TSL_coord}{'-'}{$method}) {
	$found{$feature} = $results{$TSL_chrom}{$TSL_coord}{'-'}{$method};
	delete $results{$TSL_chrom}{$TSL_coord}{'-'}{$method};
	print "found $feature -\n";
      } else {
	$not_found{$feature} = 1;
	print "NOT found $feature -\n";
      }
    }    
  }


  # the only data left in %results now are the matches to the genome
  # where there is no defined TSL Feature object, so these are all new
  # TSL sites that we could define a new Feature for. This has been
  # started, but needs more work to get the flanking sequences and the
  # clone name.

  my $total_matches = keys %results; # so we can work out the no. of alignments per million (or billion?) matches
  my $full_species = $wormbase->full_name;
  # find all RNASeq evidence for TSL sites that do not have Features
  my $newaceout = "TSL/TSL_new_features.ace";

  open (NEWACE, ">$newaceout") || $log->log_and_die("Can't open ace file $newaceout\n");
  $log->write_to("find all RNASeq evidence for TSL sites that do not have Features\n");
  foreach my $TSL_chrom (keys %results) {
    foreach my $TSL_coord (keys %{$results{$TSL_chrom}}) {
      foreach my $sense (keys %{$results{$TSL_chrom}{$TSL_coord}}) {

	# get the flanking sequences and write the Feature_object coords for this new Feature
	my ($feat_clone, $clone_start, $clone_end, $left_flank, $right_flank) = get_feature_flanking_sequences($TSL_chrom, $TSL_coord, $TSL_coord+1, $sense);
	if (defined $feat_clone) {
	    
	  foreach my $method (keys %{$results{$TSL_chrom}{$TSL_coord}{$sense}}) {
	    if ($left_flank !~ /AGTTTGAG$/ || $method ne 'SL1') { # check the splice sequence isn't in the genome next to it
	      if (!defined $results{$TSL_chrom}{$TSL_coord}{$sense}{$method}) {next} # ignore the deleted entries
	      my $ft_id = "WBsf#${TSL_chrom}#${TSL_coord}#${sense}#${method}"; # totally bogus but unique Feature ID, needs to be changed before adding to Build
	      print NEWACE "\n\nFeature : \"$ft_id\"\n";
	      print NEWACE "Sequence \"$feat_clone\"\n";
	      print NEWACE "Flanking_sequences \"$feat_clone\" \"$left_flank\" \"$right_flank\"\n";
	      print NEWACE "Species \"$full_species\"\n";
	      print NEWACE "Description \"$method trans-splice leader acceptor site\"\n";
	      print NEWACE "SO_term SO:0000706\n";
	      print NEWACE "Method $method\n";
	      print NEWACE "Defined_by_analysis $analysis\n";
	      my $reads =  $results{$TSL_chrom}{$TSL_coord}{$sense}{$method}; # reads
	      print NEWACE "Remark \"Defined by RNASeq data from $analysis with $reads reads\"\n";
	      print NEWACE "\n\n";
	      print NEWACE "Sequence : \"$feat_clone\"\n";
	      print NEWACE "Feature_object \"$ft_id\" $clone_start $clone_end\n\n";
	      
	      #$log->write_to("No existing Feature object for: chrom: $TSL_chrom pos: $TSL_coord sense: $sense type: $method, reads= $results{$TSL_chrom}{$TSL_coord}{$sense}{$method}\n");
	    } else {
	      $log->write_to("Ignoring site chromosome: $TSL_chrom, position: $TSL_coord sense: $sense because it is just downstream of the last 8 bases of the standard SL1 and so it probably a spurious hit and not trans-spliced\n");
	    }
	  }
	} else {
	  $log->write_to("Couldn't find flanking sequence for chromosome: $TSL_chrom, position: $TSL_coord sense: $sense\n");
	}
      }
    }
  }
  close(NEWACE);

  # show all Feature objects that do not have evidence from this RNASeq data-set
#  $log->write_to("show all Feature objects that do not have evidence from this RNASeq data-set\n");
#  foreach my $feature (keys %not_found) {
#    $log->write_to("Not found: $feature\n");
#  }

  # write out evidence for matched existing Feature objects
  my $aceout = "TSL/TSL_evidence.ace";
  open (ACE, ">$aceout") || $log->log_and_die("Can't open ace file $aceout\n");
  foreach my $feature (keys %found) {
    print ACE "\n\nFeature : \"$feature\"\n";
    print ACE "Defined_by_analysis $analysis\n";
    my $reads = $found{$feature}; # reads
    print ACE "Remark \"Defined by RNASeq data from $analysis with $reads reads\"\n";
  }
  close(ACE);

  $log->write_to("TSL stuff done\n");



  return $status;
}

######################################
# finds matches to the TSL sequences #
######################################

# Returns name of type of TSL and length of TSL at start of sequence
# this routine is not used - using the script get_TSL_RNASeq_reads.pl instead

sub match_tsl {

  my ($seq, $sense, $SL_hashref)=@_;

  # searches for n-mers of the TSL (max 18, not 23 so that we have a decent 16 bases left to search with)
  # min 8 so that we are fairly confident that we have a TSL
  if ($sense eq '+') {
    for (my $i=18; $i>8; $i--) {
      
      # loops through the TSL sequences 
      foreach my $slname (keys %{$SL_hashref}){
	
	#next if $i > length $SL_hashref->{$slname};
	my $nmer=substr($SL_hashref->{$slname}, -$i);
	
	if ($seq=~ /^${nmer}/) { # both sequences are already in upercase
	  return ($slname, $i);
	}
      }
    }
  }
  return undef;
}


############################################################################
# this will write out an acedb tablemaker defn to a temp file
############################################################################

sub write_feature_def {
  my $def = "/tmp/Features_$$.def";
  open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
  my $species = $wormbase->full_name;
  my $txt = <<END;

Sortcolumn 1

Colonne 1
//Subtitle Feature
Width 20
Optional
Visible
Class
Class Feature
From 1
Condition ((Method = "SL1") OR (Method = "SL2")) AND (Species = "$species")

Colonne 2
//Subtitle Clone
Width 12
Mandatory
Visible
Class
Class Sequence
From 1
Tag Sequence

Colonne 3
//Subtitle Feature2
Width 12
Optional
Hidden
Class
Class Feature
From 2
Tag Feature_object
Condition IS \\%1

Colonne 4
//Subtitle Start_pos
Width 12
Optional
Visible
Integer
Right_of 3
Tag HERE

Colonne 5
//Subtitle End_pos
Width 12
Optional
Visible
Integer
Right_of 4
Tag HERE

Colonne 6
//Subtitle Type
Width 12
Optional
Visible
Class
Class Method
From 1
Tag Method

END

  print TMP $txt;
  close TMP;
  return $def;
}

############################################################################
# Reads the Junction file to create confirmed an intron ACE file
############################################################################


sub Intron_stuff {

  my ($cmd_extra, $analysis) = @_;

  $log->write_to("Reading splice junctions BED file\n");

  my $status = 0;

# BED format provides a flexible way to define the data lines that are
# displayed in an annotation track. BED lines have three required fields
# and nine additional optional fields. The number of fields per line
# must be consistent throughout any single set of data in an annotation
# track. The order of the optional fields is binding: lower-numbered
# fields must always be populated if higher-numbered fields are used.

# The first three required BED fields are:

#    1. chrom - The name of the chromosome (e.g. chr3, chrY,
#       chr2_random) or scaffold (e.g. scaffold10671).
#    2. chromStart - The starting position of the feature in the
#       chromosome or scaffold. The first base in a chromosome is
#       numbered 0.
#    3. chromEnd - The ending position of the feature in the chromosome
#       or scaffold. The chromEnd base is not included in the display of
#       the feature. For example, the first 100 bases of a chromosome
#       are defined as chromStart=0, chromEnd=100, and span the bases
#       numbered 0-99.

# The 9 additional optional BED fields are:

#    4. name - Defines the name of the BED line. This label is displayed
#       to the left of the BED line in the Genome Browser window when
#       the track is open to full display mode or directly to the left
#       of the item in pack mode.
#    5. score - A score between 0 and 1000. If the track line useScore
#       attribute is set to 1 for this annotation data set, the score
#       value will determine the level of gray in which this feature is
#       displayed (higher numbers = darker gray).
#    6. strand - Defines the strand - either '+' or '-'.
#    7. thickStart - The starting position at which the feature is drawn
#       thickly (for example, the start codon in gene displays).
#    8. thickEnd - The ending position at which the feature is drawn
#       thickly (for example, the stop codon in gene displays).
#    9. itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the
#       track line itemRgb attribute is set to "On", this RBG value will
#       determine the display color of the data contained in this BED
#       line. NOTE: It is recommended that a simple color scheme (eight
#       colors or less) be used with this attribute to avoid
#       overwhelming the color resources of the Genome Browser and your
#       Internet browser.
#   10. blockCount - The number of blocks (exons) in the BED line.
#   11. blockSizes - A comma-separated list of the block sizes. The
#       number of items in this list should correspond to blockCount.
#   12. blockStarts - A comma-separated list of block starts. All of the
#       blockStart positions should be calculated relative to
#       chromStart. The number of items in this list should correspond
#       to blockCount.

  my %seqlength;
  my %virtuals;

  $coords = Coords_converter->invoke($database, 0, $wormbase);

  my $output = "Introns/Intron.ace";
  my $junctions = "tophat_out/junctions.bed";
  my $old_virtual = "";
  open (ACE, ">$output") || $log->log_and_die("Can't open the file $output\n");
  open(BED, "<$junctions") || $log->log_and_die("Can't open the file $junctions\n");
  while (my $line = <BED>) {
    if ($line =~ /^track/) {next}
    my @cols = split /\s+/, $line;
    my $chrom = $cols[0];
    my $start = $cols[1] + 1;
    my $end = $cols[2];
    my $reads = $cols[4];
    my $sense = $cols[5];
    my @blocksizes = split /,/, $cols[10];
    my $splice_5 = $start + $blocksizes[0]; # first base of intron
    my $splice_3 = $end - $blocksizes[1]; # last base of intron

    # get the clone that this intron is on
    my ($clone, $clone_start, $clone_end) = $coords->LocateSpan($chrom, $splice_5, $splice_3);

    
    if (not exists $seqlength{$clone}) {
      $seqlength{$clone} = $coords->Superlink_length($clone);
    }

    my $virtual = "${clone}:Confirmed_intron_RNASeq";
      
    if ($old_virtual ne $virtual) {
      print ACE "\nFeature_data : \"$virtual\"\n";
      $old_virtual = $virtual
    }

    if ($sense eq '-') {
      ($clone_end, $clone_start) = ($clone_start, $clone_end)
    }

# when we have changes to acedb that can deal with a Feature_data Confirmed_intron from RNASeq, then do this:
#    print ACE "Confirmed_intron $clone_start $clone_end RNASeq $analysis $reads\n";
# until then we store it as a Feature_data Feature, which works quite well:
    print ACE "Feature RNASeq_splice $clone_start $clone_end $reads $analysis\n";

  }
  close(BED);
  close(ACE);

  
  # now add the Feature_data objects to the chromosome objects
  my $vfile = "Introns/virtual_objects." . $wormbase->species . ".RNASeq.ace";
  open(my $vfh, ">$vfile") or $log->log_and_die("Could not open $vfile for writing\n");
  foreach my $clone (keys %seqlength) {
    print $vfh "\nSequence : \"$clone\"\n";
    my $virtual = "${clone}:Confirmed_intron_RNASeq";
    printf $vfh "S_Child Feature_data $virtual 1 $seqlength{$clone}\n\n";
  }
  close($vfh);


  return $status;
}

############################################################################

sub get_life_stage_from_conditions {

  my %life_stages;
  my $table_def = &write_life_stage_def;
  my $table_query = $wormbase->table_maker_query($database, $table_def);
  while(<$table_query>) {
    chomp;
    s/\"//g;  #remove "
    next if (/acedb/ or /\/\//);
    my @data = split("\t",$_);
    my ($condition, $life_stage) = @data;
    if (!defined $condition || !defined $life_stage) {next}
    $life_stages{$condition} = $life_stage; # multiple life-stages get overwritten leaving just the last one
  }

  return %life_stages;
}
############################################################################
# this will write out an acedb tablemaker defn to a temp file
############################################################################

sub write_life_stage_def {
  my $def = "/tmp/Life_stages_$$.def";
  open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
  my $species = $wormbase->full_name;
  my $txt = <<END2;

Sortcolumn 1

Colonne 1 
Subtitle Analysis 
Width 50 
Optional 
Visible 
Class 
Class Analysis 
From 1 
Condition Sample
 
Colonne 2 
Subtitle Condition 
Width 50 
Mandatory 
Hidden
Class 
Class Condition 
From 1 
Tag Sample 
 
Colonne 3 
Width 40 
Mandatory 
Visible 
Class 
Class Life_stage 
From 2 
Tag Life_stage

END2

  print TMP $txt;
  close TMP;
  return $def;
}
############################################################################

sub get_CDSs {

  my %CDS;
  my $table_def = &write_CDS_def;
  my $table_query = $wormbase->table_maker_query($database, $table_def);
  while(<$table_query>) {
    chomp;
    s/\"//g;  #remove "
    next if (/acedb/ or /\/\//);
    my @data = split("\t",$_);
    my ($CDS) = @data;
    if (!defined $CDS) {next}
    $CDS{$CDS} = -1;
  }

  return %CDS;
}
############################################################################
# this will write out an acedb tablemaker defn to a temp file
############################################################################

sub write_CDS_def {
  my $def = "/tmp/CDS_$$.def";
  open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
  my $species = $wormbase->full_name;
  my $txt = <<END3;

Sortcolumn 1

Colonne 1
//Subtitle CDS
Width 80
Optional
Visible
Class
Class CDS
From 1

END3

  print TMP $txt;
  close TMP;
  return $def;
}
############################################################################

sub get_Pseudogenes {

  my %pseudogenes;
  my $table_def = &write_Pseudogenes_def;
  my $table_query = $wormbase->table_maker_query($database, $table_def);
  while(<$table_query>) {
    chomp;
    s/\"//g;  #remove "
    next if (/acedb/ or /\/\//);
    my @data = split("\t",$_);
    my ($Pseudogene) = @data;
    if (!defined $Pseudogene) {next}
    $pseudogenes{$Pseudogene} = -1;
  }

  return %pseudogenes;
}
############################################################################
# this will write out an acedb tablemaker defn to a temp file
############################################################################

sub write_Pseudogenes_def {
  my $def = "/tmp/Life_stages_$$.def";
  open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
  my $species = $wormbase->full_name;
  my $txt = <<END4;

Sortcolumn 1

Colonne 1
//Subtitle Pseudogene
Width 80
Optional
Visible
Class
Class Pseudogene
From 1

END4

  print TMP $txt;
  close TMP;
  return $def;
}



############################################################################
# get the clone coords and the flanking sequences


sub get_feature_flanking_sequences {
  my ($chromosome, $start, $end, $sense) = @_;

  # returned data
  my ($left_flank, $right_flank);

  my @clone_coords = $coords->LocateSpan($chromosome, $start, $end ); 
  my ($clone, $clone_start, $clone_end) = @clone_coords;
    
  # use reverse sense coords
  if ($sense eq '-') {
    ($clone_end, $clone_start) = ($clone_start, $clone_end);
    $clone_end--;
    $clone_start--;
  }

  # now get the left and right flanking sequences for this position
  ($left_flank, $right_flank) = $mapper->get_flanking_sequence_for_feature($clone, $clone_start, $clone_end, 1, 30); # fourth argument says this is a zero-length feature

   
  if (! defined $left_flank) {
    # try again on a larger sequence
    my @clone_coords = $coords->LocateSpan($chromosome, $start-5000, $end+5000 ); 
    $clone_start = $clone_coords[1]+5000;       # the start position of feature in clone coordinates
    $clone_end = $clone_coords[2]-5000;         # the end position of feature in clone coordinates
    $clone = $clone_coords[0];

    ($left_flank, $right_flank) = $mapper->get_flanking_sequence_for_feature($clone, $clone_start, $clone_end, 1, 30); # fourth argument says this is a zero-length feature

    if (! defined $left_flank) {
      # give up
      $log->write_to( "OUT of sequence when trying to find flanking regions in $clone, $clone_start, $clone_end\n");
      return (undef, undef, undef);
    }
  }
  
  return ($clone, $clone_start, $clone_end, $left_flank, $right_flank);
}

############################################################################


__END__

=pod

=head1 NAME - align_RNASeq.pl

=head2 DESCRIPTION

This aligns RNASeq data to the current genome. It requires the Transcript data to be available in the Build.

Run on the 'long' queue as:

bsub -I -e align.err -o align.out -q long  align_RNASeq.pl [-check]


=over 4

=item *

-verbose - output verbose chatter

=back

=over 4

=item *

-check - don't run everything, resume running and only run the things that appear to have failed or not run yet.

=back

=over 4

=item *

-noalign - don't do a short-read alignment against the genome - use the results remaining from the previous alignment in the rest of the analyses

=back

=over 4

=item *

-species $SPECIES - specify the species to use

=back

=over 4

=item *

-expt <name of SRX data> - only process the specified SRX set of data - this is usually only used by this script to LSF submit another instance of this script to run the mapping of one of a set of SRX data.

=back

=over 4

=item *

-species - set the species explicitly

=back

=over 4

=item *

-database - set the data to read, the default is autoace.

=back

=cut


