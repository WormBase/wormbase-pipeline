#!/usr/local/ensembl/bin/perl5.8.0 -w

use lib $ENV{'CVS_DIR'};
use lib "/nfs/farm/Worms/Ensembl/ensembl-pipeline/modules";

use strict;
use Getopt::Long;
use DBI;
use Wormbase;
use Bio::PrimarySeq;
use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::DBSQL::Protein::DBAdaptor;

my ($debug,$fasta,$fix);

GetOptions (
	    "fasta=s" => \$fasta,
	    "debug=s" => \$debug,
	    "fix"     => \$fix
	   );

die "cant find $fasta\n" unless ( -e "$fasta" );

my $log = Log_files->make_build_log("$debug");

# mysql database parameters
my $dbhost = "ecs1f";
my $dbuser = "wormadmin";
my $dbname = "worm_pep";
my $dbpass = "worms";


my $msg = "checking $dbname against $fasta\n";
$msg .= "will fix any errors found\n";
$log->write_to("$msg");

my %ws_id2length;
my %fasta_pep;
my %sql_id2length;
my %sql_id2seq;

#prepare file handles
my $fix_file;
my $miss_file;
open( $fix_file,">mysql_fix") or die "cant open fix\n";
open( $miss_file,">mysql_miss") or die "cant open miss\n";

# get wormpep versions
&read_fasta($fasta,\%fasta_pep);

my $db = Bio::EnsEMBL::Pipeline::DBSQL::Protein::DBAdaptor->new ( -host   => $dbhost,
                                                                  -dbname => $dbname,
                                                                  -user   => $dbuser,
                                                                  -pass   => $dbpass,

                                                        );

my $proteinAdaptor = $db->get_ProteinAdaptor;


foreach my $id ( keys %fasta_pep ){
  my $pep = $fasta_pep{$id};
  my $dbpep = $proteinAdaptor->fetch_Peptide_by_dbid($pep->primary_id);

  unless( $dbpep ){
    $log->write_to($pep->primary_id." not in database\n");
    &fix($pep, $miss_file);
    next;
  }
  if( $pep->seq ne $dbpep->seq ){
    $log->write_to($pep->primary_id." has discrepant sequences\n");
    &fix($pep, $fix_file);
    next;
  }
}

close $fix_file;
close $miss_file;


$log->mail();

exit (0);

sub read_fasta
  {
    my $fasta = shift;
    my $id2seq = shift;
    my $count;

    my $stream = Bio::SeqIO->new ( -file   => $fasta,
				   -format => 'fasta',
				 );
    $stream->alphabet('protein');
    while (my $pep = $stream->next_seq) {
      if( $pep->description =~ /^(CE\d{5})/ ) {
	my $CE = $1;
	$pep->primary_id($CE);
	$pep->primary_seq->primary_id($CE);
      }
      $$id2seq{$pep->primary_id} = $pep;
    }
  }

sub fix {
  my $pep = shift;
  my $file = shift;
  print $file $pep->primary_id,"\n",$pep->seq,"\n";
}

=pod

  Checks that the contents of the wormprot database are in agreement with whats in the wp.fasta file for the current build.

=head1 OPTIONS

  -fasta   The wp.fasta file for the current version of wormpep.
  -debug   Who gets mailed.

=cut
