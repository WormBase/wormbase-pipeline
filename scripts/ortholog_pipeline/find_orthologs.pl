#!/usr/local/bin/perl -w
#files

use lib "/wormsrv2/scripts/ortholog_pipeline";

my $ver = shift;
my $ORTHO = "/nfs/disk100/wormpub/ORTHOLOGS";
my $BEST_MUTUALS_SEG_OFF = "$ORTHO/orthologs.best_mutuals.seg-off";
my $ORTHOLOG_SPANS="${ORTHO}/orthologs.ortholog.spans";
my $ORTHOLOGS_100K="${ORTHO}/orthologs.100k_spans";
my $ORTHOLOGS_200K="${ORTHO}/orthologs.200k_spans";
my $ORTHOLOGS_ALL="${ORTHO}/orthologs.all";

# Cutoffs, thresholds, etc
my $MAX_EVAL=1e-10;
my $MEX_EVAL_SEG_OFF=1e-20;
my $CUTOFF=100000;

print "copying files to $ORTHO";

#copy all files needed to ~wormpub/ORTHOLOGS
#system("scp acari:/acari/work2a/wormpipe/sort_dump/worm_brigpep.srt ${ORTHO}/briggsae_blastp");
#system("scp acari:/acari/work2a/wormpipe/sort_dump/worm_pep.srt ${ORTHO}/elegans_blastp");
system("scp wormsrv2:/wormsrv2/WORMPEP/wormpep$ver/wormpep$ver ${ORTHO}/wormpep");

print "slimming blastp results";
#slim down the blast results - we only want briggsae and elegans ones
system("perl clean_blast_results.pl -wormpep ${ORTHO}/wormpep -elegans ${ORTHO}/elegans_blastp -briggsae ${ORTHO}/briggsae_blastp");

print "converting to genes";
#convert to genes names rather than proteins.
system("cat ${ORTHO}/elegans_blastp_clean | convert_peps2genes.pl > ${ORTHO}/elegans_blastp_clean_gene");
system("cat ${ORTHO}/briggsae_blastp_clean | convert_peps2genes.pl > ${ORTHO}/briggsae_blastp_clean_genes");

print "writing elegans_positions.out";
#Create elegans_positions.out -req wormsrv2 and to query autoace
system("write_elegans_positions.pl > ${ORTHO}/elegans_postions.pl");

print "Finding best mutuals from the seg-off blastp data";
#run best_mutal_blastp
system("best_mutuals.pl -blast1 ${ORTHO}/elegans_blastp_clean_gene -blast2 ${ORTHO}/elegans_blastp_clean_gene -max_eval $MAX_EVAL_SEG_OFF -method seg-off -cutoff 100000 > $BEST_MUTUALS_SEG_OFF");

#run synteny assignments
print "Assigning orthologs by ortholog span";
system("find_orthos_by_synteny --max_eval $MAX_EVAL --cutoff $CUTOFF -method syteny_ortholog_spans --orthologs $BEST_MUTUALS_SEG_OFF > $ORTHOLOG_SPANS");
system("cat $BEST_MUTUALS_SEG_OFF $ORTHOLOG_SPANS > $ORTHOLOGS_ALL");

print "Assigning orthologs by range span: 100kb";
system("find_orthos_by_synteny --max_eval $MAX_EVAL --cutoff $CUTOFF -method syteny_100k -range 100000 --orthologs $ORTHOLOGS_ALL > $ORTHOLOGS_100K");
system("cat $ORTHOLOG_100K >> $ORTHOLOGS_ALL");

print "Assigning orthologs by range span: 200kb";
system("find_orthos_by_synteny --max_eval $MAX_EVAL --cutoff $CUTOFF -method syteny_200k -range 200000 --orthologs $ORTHOLOGS_ALL > $ORTHOLOGS_200K");
system("cat $ORTHOLOG_200K >> $ORTHOLOGS_ALL");

print "finished...";
