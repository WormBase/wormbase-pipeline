#!/software/bin/perl -w


use strict;
use Getopt::Long;

$| = 1;

if (!exists($ENV{SORT_OPTS})) {
  die "You must set SORT_OPTS (to something like '-k2,2 -k8,8n -k10,10nr')\n";
}

my $host = 'ia64d';
my $user = 'wormro';
my $db   = 'worm_ensembl_elegans';
my $port = 3306;
my $start = undef;
my $count = undef;

my $outfile = undef;

&GetOptions(
  'host:s'       => \$host,
  'user:s'       => \$user,
  'db:s'         => \$db,
  'port:n'       => \$port,
  'start:n'      => \$start,
  'count:n'      => \$count,
  'outfile:s'    => \$outfile,
);

if (!defined($start) || !defined($count) || !defined($outfile)) {
  die "Must supply start, count and outfile\n";
}

my $comstr = "mysql -u $user -h $host -D$db -P $port -N -B -e 'select protein_feature_id, stable_id, seq_start, seq_end, hit_start, hit_end, hit_id, analysis_id, score, -log10(evalue), perc_ident from protein_feature p,translation_stable_id t where p.translation_id=t.translation_id and analysis_id in (select analysis_id from analysis where module=\"BlastPep\") limit $start,$count' | sort $ENV{SORT_OPTS} -S 1G -o $outfile";

system($comstr);
