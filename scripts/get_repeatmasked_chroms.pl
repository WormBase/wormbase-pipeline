#!/usr/local/bin/perl5.8.0 -w
#


use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

use lib '/software/worm/lib/bioperl-live';
use lib '/software/worm/ensembl/ensembl/modules';
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::PrimarySeq;
use Bio::SeqIO;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $species,$wormbase,$database, $softmask, $out_file);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "output:s"   => \$out_file,
	    "database:s" => \$database,
	    "species:s"  => \$species,
            "softmask"   => \$softmask,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species,
			     );
}
# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);
}
# establish log file.
my $log = Log_files->make_build_log($wormbase);

$species = $wormbase->species;
unless ($database =~ /$species/) {
  $log->log_and_die("are you sure you have the right species / database combo\n Species :$species\nDatabase : $database\n");
}

my $chr_assembly = ($wormbase->assembly_type eq 'contig') ? 0 : 1;

$log->write_to("Connecting to worm_dna\n");

my $dbobj = Bio::EnsEMBL::DBSQL::DBAdaptor->
    new(
        '-host'   => 'farmdb1',
        '-user'   => 'wormro',
        '-dbname' => $database
        )
    or die "Can't connect to Database $database";

$log->write_to("Building chromosomes\n");

if (not defined $out_file) {
  $out_file = ($softmask) ? $wormbase->softmasked_genome_seq : $wormbase->masked_genome_seq;
}


my $seqio = Bio::SeqIO->new(-format => 'fasta',
                         -file => ">$out_file");

foreach my $seq ( @{$dbobj->get_SliceAdaptor->fetch_all('toplevel')}) {
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

  $ens_slice = $softmask
      ? $ens_slice->get_repeatmasked_seq(undef, 1) 
      : $ens_slice->get_repeatmasked_seq;

  my $seqobj = Bio::PrimarySeq->new(-seq => $ens_slice->seq,
                                    -id => $seq_name,
                                    -desc => "1 " . $ens_slice->seq_region_length);

  $seqio->write_seq($seqobj);
}

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################




# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__


