#!/usr/local/bin/perl5.8.0 -w
#
# map_microarray.pl
#
# Add information to Microarray_results objects based on overlaps in GFF files 
#
# by Dan
#
# Last updated by: $Author: klh $                      
# Last updated on: $Date: 2013-10-14 10:16:24 $        

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Ace;

######################################
# variables and command-line options #
######################################
my $maintainers = "All";
my $help;       # Help perldoc
my $test;       # Test mode
my $debug;      # Debug mode, verbose output to user running script
my $noload;     
my $store;	# specify a frozen configuration file
my $species;
my $outfile;

GetOptions ("debug=s"   => \$debug,
	    "test"      => \$test,
	    "noload"    => \$noload,
            "species=s" => \$species,
	    "acefile=s"	=> \$outfile,
	    'store=s'	=> \$store
    );


############################
# recreate configuration   #
############################
my $wb;
if ($store){
  $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n");
} else {
  $wb = Wormbase->new(-debug    => $debug,
                      -test     => $test,
                      -organism => $species,
      );
}

###########################################
# Variables Part II (depending on $wb)    #
########################################### 

my $dbdir= $wb->autoace;                    # Database path
$outfile = $outfile ? $outfile : $wb->acefiles."/microarray_mappings.ace";
my $tace = $wb->tace;                                  # tace executable path

my $log = Log_files->make_build_log($wb);

my $db = Ace->connect(-path=>$dbdir,
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; $log->log_and_die();};

if ($debug) {
  my $count = $db->fetch(-query=> 'find PCR_product where Microarray_results');
  $log->write_to("checking $count PCR_products...\n");

  $count = $db->fetch(-query=> 'find Oligo_set where Microarray_results');
  $log->write_to("checking $count Oligo_sets\n");
}


open (OUTPUT, ">$outfile") or $log->log_and_die("Can't open the output file $outfile\n");

###########################################
# PCR_products and SMD_microarray results #
###########################################

$log->write_to("Making CDS/Pseudogene/Transcript connections to Microarray_results objects based on PCR_products\n");

my $i = $db->fetch_many(-query=> 'find PCR_product WHERE Microarray_results AND (Overlaps_CDS || Overlaps_pseudogene || Overlaps_transcript)');  
while (my $object = $i->next) {
  &process_object($object);
  $object->DESTROY();
}

###########################################
# Oligo_sets and Aff_microarray results   #
###########################################

$log->write_to("Making CDS/Pseudogene/Transcript connections to Microarray_results objects based on Oligo_set objects\n");

$i = $db->fetch_many(-query=> 'find Oligo_set WHERE Microarray_results AND (Overlaps_CDS || Overlaps_transcript || Overlaps_pseudogene)');  

while (my $object = $i->next) {
  
  &process_object($object);
  $object->DESTROY();
}


close OUTPUT;
$db->close;

unless ($noload) {
  $log->write_to("Loading file to autoace\n");
  $wb->load_to_database($wb->autoace, $outfile, 'microarray_mappings', $log);
}



$log->mail();

exit(0);

##################################
sub process_object {
  my ($obj) = @_;
    
  # Microarray_results
  
  my $microarray_results = $obj->Microarray_results;  
  my @CDSs               = $obj->Overlaps_CDS;
  my @pseudogenes        = $obj->Overlaps_pseudogene;
  my @transcripts        = $obj->Overlaps_transcript;

  my %genes;

  print OUTPUT "\nMicroarray_results : \"$microarray_results\"\n";
  
  if (@CDSs) {
    foreach my $cds (@CDSs) {
      my $gene = $obj->Overlaps_CDS->Gene;

      print OUTPUT "CDS $cds\n";
      if (defined $gene and not exists $genes{$gene}) {
        print OUTPUT "Gene $gene\n";
        $genes{$gene} = 1;
      }
    }    
  }

  if (@pseudogenes) {
    foreach my $pseudogene (@pseudogenes) {
      my $gene = $obj->Overlaps_pseudogene->Gene;
      print OUTPUT "Pseudogene $pseudogene\n";
      if (defined $gene and not exists $genes{$gene}) {
        print OUTPUT "Gene $gene\n";
        $genes{$gene} = 1;
      }
    }    
  }

  if (@transcripts) {
    foreach my $transcript (@transcripts) {
      my $gene = $obj->Overlaps_transcript->Gene;
      $gene =  $obj->Overlaps_transcript->Corresponding_CDS->Gene unless ($gene);
      print OUTPUT "Transcript $transcript\n";
      if (defined $gene and not exists $genes{$gene}) {
        print OUTPUT "Gene $gene\n" if (defined $gene);
        $genes{$gene} = 1;
      }
    }    
  }
  
  print OUTPUT "\n";

}


############################################

__END__

=pod

=head2 NAME - map_microarray.pl

=head1 USAGE

=over 4

=item map_microarray.pl [-options]

=back

map_microarray.pl connects PCR_products and micro_array_results  objects to CDS/transcript / Pseudogene /gene objects
via the defined_by tag in the database.

mandatory arguments:

=over 4

=item none

=back

optional arguments:

=over 4

=item -debug username, Debug mode

=item -test, Test mode, generate the acefile but do not upload them 

=item -help, Help pages

=item -outace specify location for the output Ace (if you want to use it outside of the build)

=item -store specifiy a configuration file

=back

=cut

