#!/usr/local/bin/perl5.6.1 -w

use DBI;
use strict;
use lib "/wormsrv2/scripts/";
use Getopt::Long;
use Wormbase;

#######################################
# command-line options                #
#######################################
my ($test, $debug, $help, $all);

GetOptions ("debug"   => \$debug,
	    "test"    => \$test,
	    "help"    => \$help,
	    "all"     => \$all
           );


my $dump_dir = "/wormsrv2/wormbase/ensembl_dumps";
my $acedb_database;
my $output = "$dump_dir/repeats.ace";
die &help if $help;

$output .= "_test" if $test;

open (OUT,">$output") or die "cant open $output\n";

# retrieve hash of acc2clone
my %acc2clone;
&FetchData("acc2clone",\%acc2clone);

# mysql database parameters
my $dbhost = "ecs1f";
my $dbuser = "wormro";
my $dbname = "wormrepeats";
my $dbpass = "";

my $wormrepeats_DB = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to db, $DBI::errstr";

# build hash of internal_id to clone via EMBL accession
my $sth_f = $wormrepeats_DB->prepare ( q{SELECT
					 internal_id , id
					 FROM clone
	  	  	     } );

my %int_id2clone;

$sth_f->execute;
my $ref_results = $sth_f->fetchall_arrayref;

foreach my $pair (@$ref_results) {
  my $int_id = $$pair[0];
  my $acc = $$pair[1];
  my $clone = $acc2clone{$acc};

  if ($clone) {
    $int_id2clone{$int_id} = $clone;
  }
  else {
    print "no clone found for acc $acc : internal_id $int_id\n";
  }
}

# now retrieve the data from wormrepeats database
$sth_f = $wormrepeats_DB->prepare ( q{SELECT
				      contig, seq_start, seq_end,
				      score, strand,
				      hstart,hend,hid
				      FROM repeat_feature
				      WHERE analysis = 1
				      AND score > 225
				      ORDER BY contig
				     }
				  );

$sth_f->execute;
$ref_results = $sth_f->fetchall_arrayref;
my ($contig, $seq_start, $seq_end, $score, $strand, $hstart, $hend, $hid);
my $current_clone = "";
my $clone;
foreach my $repeat (@$ref_results) {
  ($contig, $seq_start, $seq_end, $score, $strand, $hstart, $hend, $hid) = @$repeat;
  
  if($int_id2clone{$contig}) {
     $clone = $int_id2clone{$contig};
    if( "$clone" ne "$current_clone" ) {
      $current_clone = $clone;
      print OUT "\nSequence : \"$clone\"\n";
    }
    print OUT "Motif_homol $hid RepeatMasker $score $seq_start $seq_end";
    if( $strand == 1 ){
      print OUT " $hstart $hend\n";
    }
    else {
      print OUT " $hend $hstart\n";
    }
  }
  else {
    undef $clone;
  }
}

close OUT;

exit(0);  

sub help
  {
    print "===============================================\n$0\n
Extract RepeatMasker data from wormrepeats on ecs1f\n\n
Writes ace file $output\n
\t\t-test appends _test to output filename.

Must be able to access Wormbase.pm\
Takes about 2 mins to dump whole genome\n\n
================================================\n\n";
  }

