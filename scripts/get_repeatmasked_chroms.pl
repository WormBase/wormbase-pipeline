#!/usr/local/bin/perl5.8.0 -w
#
# script to dump the EnsEMBL toplevel sequences as
# softmasked and hardmasked FASTA files


use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::PrimarySeq;
use Bio::SeqIO;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $species,$wormbase,$dbname, $dbhost, $dbport, $softmask, $nomask, $out_file);

GetOptions (
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "output:s"   => \$out_file,
	    "species:s"  => \$species,
            "softmask"   => \$softmask,
            "nomask"     => \$nomask,
	    "dbname:s"   => \$dbname,
            "dbhost:s"   => \$dbhost,
            "dbport:s"   => \$dbport,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species,
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

$species = $wormbase->species;
if (defined $dbname and $dbname !~ /$species/) {
  $log->write_to("WARNING: db $dbname does not look like it belongs to species $species. Proceeding anyway");
}

$dbhost = $ENV{WORM_DBHOST} if not defined $dbhost;
$dbport = $ENV{WORM_DBPORT} if not defined $dbport;
$dbname = "worm_ensembl_${species}" if not defined $dbname;

my $dbobj = Bio::EnsEMBL::DBSQL::DBAdaptor->
    new(
        '-user'   => "wormro", 
        '-host'   => $dbhost,
        '-port'   => $dbport,
        '-dbname' => $dbname,
        ) or die "Can't connect to Database $dbname";

$log->write_to("Building chromosomes\n");

if (not defined $out_file) {
  if ($softmask) {
    $out_file = $wormbase->softmasked_genome_seq;
  } elsif ($nomask) {
    $out_file = $wormbase->genome_seq;
  } else {
    $out_file = $wormbase->masked_genome_seq;
  }
}


my $seqio = Bio::SeqIO->new(-format => 'fasta',
                            -file => ">$out_file");

foreach my $seq ( sort { $b->length <=> $a->length } @{$dbobj->get_SliceAdaptor->fetch_all('toplevel')}) {
  &print_seq($seqio, $seq);
}

$seqio->close();

$log->write_to("Done\n");
$log->mail();
print "Finished.\n" if ($verbose);
exit(0);


##############################################################
#
# Subroutines
#
##############################################################

sub print_seq {
  my ($seqio,$ens_slice) = @_;
 
  my $seq_name = $ens_slice->seq_region_name;

  $log->write_to("\twriting chromosome $seq_name\n");

  if ($softmask) {
    $ens_slice = $ens_slice->get_repeatmasked_seq(undef, 1) ;
  } elsif (not $nomask) {
    $ens_slice = $ens_slice->get_repeatmasked_seq;
  }

  my $seqobj = Bio::PrimarySeq->new(-seq => $ens_slice->seq,
                                    -id => $seq_name,
                                    -desc => "1 " . $ens_slice->seq_region_length);

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
