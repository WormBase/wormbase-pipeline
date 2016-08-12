#!/usr/bin/env perl
#
# script to dump the EnsEMBL toplevel sequences as
# softmasked and hardmasked FASTA files
#

use strict;                                      
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use IO::File;
use Bio::SeqIO;

my $dbhost = "mysql-wormbase-pipelines";
my $dbuser = "wormro";
my $dbport = "4331";
my $dbpass = "";


my ($database, $mask, $softmask, $outfile, $ebi_header_prefix, $species_string);

GetOptions (
  'host=s'           => \$dbhost,
  'port=s'           => \$dbport,
  'user=s'           => \$dbuser,
  'pass=s'           => \$dbpass,
  'dbname=s'         => \$database,
  "softmask"         => \$softmask,
  "mask"             => \$mask, 
  "ebiblastheader=s" => \$ebi_header_prefix,
  "outfile:s"        => \$outfile,
	    );
my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -user   => $dbuser, 
  -host   => $dbhost,
  -port   => $dbport,
  -pass   => $dbpass,
  -dbname => $database,
    ) or die "Can't connect to Database $database";


if (defined $ebi_header_prefix) {
  $species_string = ucfirst( $db->get_MetaContainer->get_production_name );
}

my $seqio = (defined $outfile) 
    ? Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile")
    : Bio::SeqIO->new(-format => 'fasta', -fh => \*STDOUT);


foreach my $seq ( sort { $b->length <=> $a->length } @{$db->get_SliceAdaptor->fetch_all('toplevel')}) {
  &print_seq($seqio, $seq);
}

$seqio->close();

exit(0);


##############################################################
#
# Subroutines
#
##############################################################

sub print_seq {
  my ($seqio,$ens_slice) = @_;
 
  my $seq_name = $ens_slice->seq_region_name;

  if ($softmask) {
    $ens_slice = $ens_slice->get_repeatmasked_seq(undef, 1) ;
  } elsif ($mask) {
    $ens_slice = $ens_slice->get_repeatmasked_seq;
  }

  my ($id, $desc);
  if (defined $ebi_header_prefix) {
    $id = join(":", $ebi_header_prefix, $species_string, $seq_name);
    $desc = sprintf("dna:%s %s species:%s",
                    $ens_slice->coord_system->name,
                    $ens_slice->name,
                    $species_string);

  } else {
    $id = $seq_name;
    $desc = "1-" . $ens_slice->seq_region_length;
  }

  my $seqobj = Bio::PrimarySeq->new(-seq => $ens_slice->seq,
                                    -id => $id,
                                    -desc => $desc);

  $seqio->write_seq($seqobj);
}



##########################################


# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.

__END__

=pod

=head2 NAME - get_repeatmasked_chroms.pl

=head3 EXAMPLE USAGE 

get_repeatmasked_chrom -species brugia -soft

=head2 DESCRIPTION

This script "does exactly what it says on the tin". ie it dumps fasta sequences of the toplevel sequences in the EnsEMBL databases in a hardmasked/softmasked or raw FASTA format.

=head2 OPTIONS:

=over 4

B<-help>

print this perldoc

B<-debug NAME>

set a custom email target

B<-test>

use the test locations for the BUILD

B<-verbose>

be more verbose in the output

B<-store FILENAME>

pass a specific storable

B<-output FILENAME>

output filename

B<-database MYSQL_DATABASE_NAME>

specify a mysql database

B<-species WORMBASE_SPECIES>

specify the wormbase species

B<-softmask>

lowercase repetitive sequence in the FASTA instead of using Ns

=back

=cut
