#!/usr/local/bin/perl5.8.0 -w
#
# map_microarray.pl
#
# Add information to Microarray_results objects based on overlaps in GFF files 
#
# by Anon
#
# Last updated by: $Author: mh6 $                      
# Last updated on: $Date: 2005-12-20 14:10:17 $        

use strict;
use Wormbase;
use lib $ENV{'CVS_DIR'};
use Getopt::Long;
use Ace;

######################################
# variables and command-line options #
######################################
my $maintainers = "All";
my $help;       # Help perdoc
my $test;       # Test mode
my $debug;      # Debug mode, verbose output to user running script
my $load;       # load file to autoace
my $store;	# specify a frozen configuration file

my $outfile;

GetOptions ("debug=s"   => \$debug,
	    "test"      => \$test,
	    "load"      => \$load,
            "help"      => \$help,
	    "outace=s"	=> \$outfile,
	    'store=s'	=> \$store
    );

# Display help if required
&usage("Help") if ($help);

############################
# recreate configuration   #
############################
my $wb;
if ($store){$wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n")}
else {$wb = Wormbase->new(-debug => $debug,-test => $test,)}

###########################################
# Variables Part II (depending on $wb)    #
########################################### 
$test  = $wb->test  if $wb->test;     # Test mode
$debug = $wb->debug if $wb->debug;    # Debug mode, output only goes to one user


#further variables
my $dbdir= $wb->autoace;                    # Database path
$outfile = $outfile ? $outfile : "$dbdir/acefiles/microarray_mappings.ace";
my $tace        = $wb->tace;                                  # tace executable path

# create log
my $log = Log_files->make_build_log($wb);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

# connect to database
print  "Opening database ..\n" if ($debug);
my $db = Ace->connect(-path=>$dbdir,
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};


if ($debug) {
    my $count = $db->fetch(-query=> 'find PCR_product where Microarray_results');
    print "checking $count PCR_products\n\n";

    $count = $db->fetch(-query=> 'find Oligo_set where Microarray_results');
    print "checking $count Oligo_sets\n\n";
}


open (OUTPUT, ">$outfile") or die "Can't open the output file $outfile\n";

###########################################
# PCR_products and SMD_microarray results #
###########################################

$log->write_to("Making CDS/Pseudogene/Transcript connections to Microarray_results objects based on PCR_products\n");

my $i = $db->fetch_many(-query=> 'find PCR_product WHERE Microarray_results AND (Overlaps_CDS || Overlaps_pseudogene || Overlaps_transcript)');  
while (my $obj = $i->next) {
    
  print "$obj\t" if ($debug);

  # Microarray_results
  
  my $microarray_results = $obj->Microarray_results;  
  my @CDSs               = $obj->Overlaps_CDS;
  my @pseudogenes        = $obj->Overlaps_pseudogene;
  my @transcripts        = $obj->Overlaps_transcript;

  
  if (@CDSs) {
    foreach my $cds (@CDSs) {
      my $gene = $obj->Overlaps_CDS->Gene;
      print OUTPUT "\nMicroarray_results : \"$microarray_results\"\n";
      print OUTPUT "CDS \"$cds\"\n";
      print OUTPUT "Gene $gene\n" if (defined $gene);
    }    
    print OUTPUT "\n";
  }

  if (@pseudogenes) {
    foreach my $pseudogene (@pseudogenes) {
      my $gene = $obj->Overlaps_pseudogene->Gene;
      print OUTPUT "\nMicroarray_results : \"$microarray_results\"\n";
      print OUTPUT "Pseudogene \"$pseudogene\"\n";
      print OUTPUT "Gene $gene\n" if (defined $gene);
    }    
    print OUTPUT "\n";
  }

  if (@transcripts) {
    foreach my $transcript (@transcripts) {
      my $gene = $obj->Overlaps_transcript->Gene;
      print OUTPUT "\nMicroarray_results : \"$microarray_results\"\n";
      print OUTPUT "Transcript \"$transcript\"\n";
      print OUTPUT "Gene $gene\n" if (defined $gene);
    }    
    print OUTPUT "\n";
  }
  
  $obj->DESTROY();
}

###########################################
# Oligo_sets and Aff_microarray results   #
###########################################

$log->write_to("Making CDS/Pseudogene/Transcript connections to Microarray_results objects based on Oligo_set objects\n");

$i = $db->fetch_many(-query=> 'find Oligo_set WHERE Microarray_results AND (Overlaps_CDS || Overlaps_transcript || Overlaps_pseudogene)');  

while (my $obj = $i->next) {
  
  print "$obj\t" if ($debug);

  # Microarray_results
    
  my $microarray_results = $obj->Microarray_results;  
  my @CDSs               = $obj->Overlaps_CDS;
  my @pseudogenes        = $obj->Overlaps_pseudogene;
  my @transcripts        = $obj->Overlaps_transcript;

  if (@CDSs) {
    foreach my $cds (@CDSs) {
      my $gene = $obj->Overlaps_CDS->Gene;
      print OUTPUT "\nMicroarray_results : \"$microarray_results\"\n";
      print OUTPUT "CDS \"$cds\"\n";
      print OUTPUT "Gene $gene\n" if (defined $gene);
    }    
    print OUTPUT "\n";
  }

  if (@pseudogenes) {
    foreach my $pseudogene (@pseudogenes) {
      my $gene = $obj->Overlaps_pseudogene->Gene;
      print OUTPUT "\nMicroarray_results : \"$microarray_results\"\n";
      print OUTPUT "Pseudogene \"$pseudogene\"\n";
      print OUTPUT "Gene $gene\n" if (defined $gene);
    }    
    print OUTPUT "\n";
  }

  if (@transcripts) {
    foreach my $transcript (@transcripts) {
      my $gene = $obj->Overlaps_transcript->Gene;
      print OUTPUT "\nMicroarray_results : \"$microarray_results\"\n";
      print OUTPUT "Transcript \"$transcript\"\n";
      print OUTPUT "Gene $gene\n" if (defined $gene);
    }    
    print OUTPUT "\n";
  }  

  $obj->DESTROY();
}


close OUTPUT;
$db->close;

if($load){
  $log->write_to("Loading file to autoace\n");
  my $command = "autoace_builder.pl -load $dbdir/acefiles/microarray_mappings.ace -tsuser microarray_mappings";
                                                                                   
  my $status = system($command);
  if(($status >>8) != 0){
    $log->write_to("ERROR: Loading microarray_mappings.ace file failed \$\? = $status\n");
  }
}



$log->mail("$maintainers","BUILD REPORT: $0");

exit(0);
sub usage {
    my $error = shift;

    if ( $error eq "Help" ) {

        # Normal help menu
        system( 'perldoc', $0 );
        exit(0);
    }
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

=item -load, loads file to autoace

=item -help, Help pages

=item -outace specify location for the output Ace (if you want to use it outside of the build)

=item -store specifiy a configuration file

=back

=cut

