#!/nfs/disk100/wormpub/bin/perl -w
#
# $Id: process_sage_gff.pl,v 1.3 2006-11-30 14:35:58 gw3 Exp $;
#
# process the raw Sanger GFF dump to add data to SAGE tags
# Sheldon McKay <mckays@cshl.edu>
#

use lib $ENV{'CVS_DIR'};
use strict;
use Ace;
use Storable;
use Wormbase;
use Log_files;
use Getopt::Long;

my ($help, $debug, $test, $verbose, $store, $wormbase);
my $gff;

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "gff:s"      => \$gff
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}


# establish log file.
my $log = Log_files->make_build_log($wormbase);

my $db = Ace->connect( -path => $wormbase->database('current') );

my $c_reg = qr/CHROMOSOME_/;
my $s_reg = qr/SAGE_tag:(SAGE:[a-z]+)/;


my $go;

my @chroms = $wormbase->get_chromosome_names(-mito =>1, -prefix => 1);
my $gff_path = $wormbase->chromosomes;

foreach my $chrom (@chroms){
    my $file = $gff_path."/${chrom}.gff";
    my $tmp_file = "$file.tmp";
    open (GFF,"<$file") or $log->log_and_die("Cant open $file: !$\n");
    open (TMP,">$tmp_file") or  $log->log_and_die("Cant open $tmp_file: !$\n");
    while (<GFF>) {
	if( !$_ || /^\#/ || !/SAGE_tag/) {
	    print TMP;
	    next;
	}
	
	my ($ref,$src,$met,$start,$stop,$scr,$strand,$phase,$att) = split "\t";

	# infer strandedness from targets
	if ($att =~ /Target\s+\S+\s+(\d+)\s+(\d+)/) {
	    if ($2 < $1) {
		$strand = '-';
	    }
	    else {
		$strand = '+';
	    }
	}

	#$ref =~ s/$c_reg//;

	my ($tag) = $att =~ /$s_reg/;
	$strand =~ tr/\./+/;
	$tag or print TMP;

	# only mark up unambiguously mapped tags
	my @atts = &markup($tag);
	my $atts = join ';', "Sequence $tag", @atts;
	print TMP join ("\t",$ref,$src,$met,$start,$stop,$scr,$strand,$phase,$atts), "\n";
    }
    close TMP or  $log->error("Problem closing $tmp_file: !$\n");
    $wormbase->run_command("mv -f $tmp_file $file", $log);
}

$log->mail;
exit;

sub markup {
  my $tag = shift;
  my $src = shift;
  my $r = $db->fetch(SAGE_tag => $tag);
  my @genes    = ($r->Gene,$r->Transcript,$r->Pseudogene,$r->Predicted_CDS);
  my @unambig  = grep { scalar $_->get('Unambiguously_mapped')} @genes;
  my @three_p  = grep { scalar $_->get('Most_three_prime') } @genes;
  if (!@unambig && !@three_p) {
    my @unique = grep {$_->class eq 'Gene' || $_->class eq 'Pseudogene'} @genes;
    @unambig = @unique if @unique == 1;
  }
  
  my @count = $r->Results;
  my $count;
  for my $exp (@count) {
    my $cnt = $exp->get('Frequency')->right;
    $count += $cnt;
  }
  for (@unambig, @three_p) {
    my $name = eval{ $_->CGC_name} || eval{$_->Sequence_name} || $_->name;
    $_ = join(' ',$_->class,$name);
  }
  return "count $count", @unambig, @three_p;
}
