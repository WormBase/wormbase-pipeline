#!/usr/local/bin/perl5.8.0 -w
#
# map_microarray.pl
#
# Add information to Microarray_results objects based on overlaps in GFF files 
#
# by Anon
#
# Last updated by: $Author: krb $                      
# Last updated on: $Date: 2004-10-20 12:03:04 $        


use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Cwd;
use Ace;


######################################
# variables and command-line options #
######################################

my $tace        = &tace;                                  # tace executable path
my $dbdir       = "/wormsrv2/autoace";                    # Database path

my $maintainers = "All";
my $rundate = &rundate;
my $runtime = &runtime;
my $help;       # Help perdoc
my $test;       # Test mode
my $debug;      # Debug mode, verbose output to user running script
my $log = Log_files->make_build_log();
my $verbose;    # verbose mode, extra output to screen
my $load;       # load file to autoace

my $outfile = "$dbdir/acefiles/microarray_mappings.ace";

GetOptions ("debug=s"   => \$debug,
	    "test"      => \$test,
	    "load"      => \$load,
            "help"      => \$help);


# Display help if required
&usage("Help") if ($help);

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
# Oligo_sets and Aff_microarray results #
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
  my $command = "autoace_minder.pl -load $dbdir/acefiles/microarray_mappings.ace -tsuser
RNAi_mappings";
                                                                                   
  my $status = system($command);
  if(($status >>8) != 0){
    $log->write_to("ERROR: Loading microarray_mappings.ace file failed \$\? = $status\n");
  }
}



$log->mail("$maintainers","BUILD REPORT: $0");

exit(0);

