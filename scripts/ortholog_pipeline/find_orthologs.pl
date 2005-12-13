#!/usr/local/bin/perl -w
#files

use lib "/nfs/farm/Worms/Ensembl/bioperl-live";
use Getopt::Long;
my( $copy, $slim, $convert, $positions, $mutual, $span, $k_100, $k_200 );

GetOptions ( 'copy' => \$copy,
	     'slim' => \$slim,
	     'convert' => \$convert,
	     'positions' => \$positions,
	     'mutual' => \$mutual,
	     'span' => \$span,
	     '100k' => \$k_100,
	     '200k' => \$k_200
	   );

my $all = 1 unless ($copy or $slim or $convert or $positions or $mutual or $span or $k_100 or $k_200 );


my $scripts_dir = (-e "/wormsrv2/scripts") ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
$scripts_dir .= "/ortholog_pipeline";

die "no scripts\n" unless ( -e "$scripts_dir" );

my $ver = shift;
die "please give WS version no. eq $0 123\n" unless $ver;
my $ORTHO = "/nfs/disk100/wormpub/ORTHOLOGS";
my $BEST_MUTUALS_SEG_OFF = "$ORTHO/orthologs.best_mutuals.seg-off";
my $ORTHOLOG_SPANS="${ORTHO}/orthologs.ortholog.spans";
my $ORTHOLOGS_100K="${ORTHO}/orthologs.100k_spans";
my $ORTHOLOGS_200K="${ORTHO}/orthologs.200k_spans";
my $ORTHOLOGS_ALL="${ORTHO}/orthologs.all";

# Cutoffs, thresholds, etc
my $MAX_EVAL=1e-10;
my $MAX_EVAL_SEG_OFF=1e-20;
my $CUTOFF=100000;

#copy all files needed to ~wormpub/ORTHOLOGS
print "copying files to $ORTHO\n" if ( $all or $copy );
system("scp acari:/acari/work2a/wormpipe/sort_dump/worm_brigpep.srt ${ORTHO}/briggsae_blastp") if ( $all or $copy );
system("scp acari:/acari/work2a/wormpipe/sort_dump/worm_pep.srt ${ORTHO}/elegans_blastp") if ( $all or $copy );
system("scp wormsrv2:/wormsrv2/WORMPEP/wormpep$ver/wormpep$ver ${ORTHO}/wormpep") if ( $all or $copy );

#slim down the blast results - we only want briggsae and elegans ones
print "slimming blastp results\n" if ( $all or $slim );
system("perl $scripts_dir/clean_blast_results.pl -wormpep ${ORTHO}/wormpep -elegans ${ORTHO}/elegans_blastp -briggsae ${ORTHO}/briggsae_blastp") if ( $all or $slim );

print "converting to genes\n" if ( $all or $convert );
#convert to genes names rather than proteins.
system("cat ${ORTHO}/elegans_blastp_clean | $scripts_dir/convert_peps2genes.pl > ${ORTHO}/elegans_blastp_clean_gene") if ( $all or $convert );
system("cat ${ORTHO}/briggsae_blastp_clean | $scripts_dir/convert_peps2genes.pl > ${ORTHO}/briggsae_blastp_clean_genes") if ( $all or $convert );

print "writing elegans_positions.out\n" if ( $all or $positions );
#Create elegans_positions.out -req wormsrv2 and to query autoace
system("perl5.6.1 $scripts_dir/write_elegans_positions.pl -database ~wormpub/DATABASES/current_DB -version $ver > ${ORTHO}/elegans_postions.out") if ( $all or $positions );

print "Finding best mutuals from the seg-off blastp data\n" if ( $all or $mutual );
#run best_mutal_blastp
system("$scripts_dir/best_mutuals.pl -blast1 ${ORTHO}/elegans_blastp_clean_gene -blast2 ${ORTHO}/elegans_blastp_clean_gene -max_eval $MAX_EVAL_SEG_OFF -method seg-off -cutoff 100000 > $BEST_MUTUALS_SEG_OFF") if ( $all or $mutual );

#run synteny assignments
print "Assigning orthologs by ortholog span\n" if ( $all or $span );
system("$scripts_dir/find_orthos_by_synteny --max_eval $MAX_EVAL_SEG_OFF --cutoff $CUTOFF -method syteny_ortholog_spans --orthologs $BEST_MUTUALS_SEG_OFF > $ORTHOLOG_SPANS") if ( $all or $span );
system("cat $BEST_MUTUALS_SEG_OFF $ORTHOLOG_SPANS > $ORTHOLOGS_ALL") if $all;

print "Assigning orthologs by range span: 100kb\n" if ( $all or $k_100 );
system("$scripts_dir/find_orthos_by_synteny --max_eval $MAX_EVAL --cutoff $CUTOFF -method syteny_100k -range 100000 --orthologs $ORTHOLOGS_ALL > $ORTHOLOGS_100K") if ( $all or $k_100 );
system("cat $ORTHOLOGS_100K >> $ORTHOLOGS_ALL") if $all;

print "Assigning orthologs by range span: 200kb\n" if ( $all or $k_200 );
system("$scripts_dir/find_orthos_by_synteny --max_eval $MAX_EVAL --cutoff $CUTOFF -method syteny_200k -range 200000 --orthologs $ORTHOLOGS_ALL > $ORTHOLOGS_200K") if ( $all or $k_200 );
system("cat $ORTHOLOGS_200K >> $ORTHOLOGS_ALL") if $all;

print "finished...\n";
