#!/usr/bin/perl5.8.0 -w

use strict;
use Getopt::Long;
use DBI;
use lib "/wormsrv2/scripts";
use Wormbase;

my ($debug,$fasta);

GetOptions (
	    "fasta=s" => \$fasta,
	    "debug=s" => \$debug
	   );

die "cant find $fasta\n" unless ( -e "$fasta" );

my $log = Log_files->make_build_log("$debug");

# mysql database parameters
my $dbhost = "ecs1f";
my $dbuser = "wormadmin";
my $dbname = "wormprot";
my $dbpass = "worms";

my %ws_id2length;
my %ws_id2seq;
my %sql_id2length;
my %sql_id2seq;

# get wormpep versions
&read_fasta($fasta,\%ws_id2seq);

foreach (keys %ws_id2seq ){
  $ws_id2length{$_} = length ($ws_id2seq{$_});
}

#get wormprot data and compare
my $wormprot = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
  || die "cannot connect to db, $DBI::errstr";

my $query = "select protein.proteinId, protein.length, peptide.peptide from protein, peptide where peptide.proteinId = protein.proteinId";

my $sth = $wormprot->prepare( "$query" );
$sth->execute();

my $ref_results = $sth->fetchall_arrayref;

foreach my $record (@$ref_results) {
  my ($protein_id, $length, $seq) = @$record;
  $sql_id2length{$protein_id} = $length;
  $sql_id2seq{$protein_id} = $seq;

  $log->write_to("ERROR:\tDifferent lengths for $protein_id\n") if( $length != $ws_id2length{$protein_id} );
  $log->write_to("ERROR:\tDifferent sequence for $protein_id\n") if( $seq ne $ws_id2seq{$protein_id} );
  
}


#check whats not in wormprot;
open (MISSING,">proteins missing from wormprot") or die "cant write new missing file\n";
foreach my $wp_protein (sort keys %ws_id2seq ) {
  if( !($sql_id2seq{$wp_protein}) ) {
    $log->write_to("MISSING : $wp_protein missing from $dbname\n");
    print MISSING "\>$wp_protein\n$ws_id2seq{$wp_protein}\n";
  }
}

$log->mail($debug);

exit (0);

sub read_fasta
  {
    my $fasta = shift;
    my $id2seq = shift;
    my $count;
    open (WS, "<$fasta") or die "$fasta\n";
    my ($id,$seq);
    my %id2seq;
    while (<WS>) {
      chomp;
      if( /\>(\w+)/ ) {
	if( defined $id ) {
	  warn "WARNING \$seq is NULL" unless $seq;
	  $$id2seq{$id} = $seq;
	  undef $id;
	  $seq = "";
	}
	$id = $1;
      }
      else {
	$seq .= $_;
      }
    }
    $$id2seq{$id} = $seq;
  }


=pod

  Checks that the contents of the wormprot database are in agreement with whats in the wp.fasta file for the current build.

=head1 OPTIONS

  -fasta   The wp.fasta file for the current version of wormpep.
  -debug   Who gets mailed.

=cut
