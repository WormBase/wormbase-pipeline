#!/software/bin/perl -w


use strict;
use Getopt::Long;

if (!exists($ENV{SORT_OPTS})) {
  die "You must set SORT_OPTS (to something like '-k2,2 -k8,8n -k10,10nr')\n";
}

my $host = $ENV{'WORM_DBHOST'};
my $user = 'wormro';
my $db   = 'worm_ensembl_elegans';
my $port = $ENV{'WORM_DBPORT'};
my ($start, $count, $analysis, $outfile);


&GetOptions(
            'host:s'       => \$host,
            'user:s'       => \$user,
            'db:s'         => \$db,
            'port:n'       => \$port,
            'start:n'      => \$start,
            'count:n'      => \$count,
            'outfile:s'    => \$outfile,
            'analysis:i'   => \$analysis,

);

my $single = " ";
if($analysis){
    $single = "a.analysis_id = $analysis AND";
}


if (!defined($start) || !defined($count) || !defined($outfile)) {
  die "Must supply start, count and outfile\n";
}

my $comstr = "mysql -u $user -h $host -D$db -P $port -N -B -e 'SELECT protein_feature_id, stable_id, p.seq_start, p.seq_end, hit_start, hit_end, hit_name, logic_name, score, -log10(evalue), perc_ident FROM protein_feature p, translation t , analysis a WHERE $single p.translation_id=t.translation_id AND a.analysis_id=p.analysis_id AND program=\"blastp\" limit $start,$count' | sort $ENV{SORT_OPTS} -S 1G -o $outfile";

system($comstr);
