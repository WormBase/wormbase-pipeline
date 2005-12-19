#!/usr/local/bin/perl5.8.0 -w
#
# molecular_names_for_genes.pl
# 
# by Keith Bradnam
#
# quick script to populate Molecular_name tag in ?Gene model during build
#
# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2005-12-19 17:42:26 $

#################################################################################
# Initialise variables                                                          #
#################################################################################

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;

##############################
# command-line options       #
##############################
my $test;                # test mode, uses ~wormpub/TEST_BUILD
my $debug;               # debug mode, email log file only goes to debugger
my $maintainers = "All"; # log file recipients
my $store;		# to specify any storable configuration files

GetOptions('debug=s'       => \$debug,
	   'test'          => \$test,
	   'store=s'	=> \$store
   	);
	
############################
# recreate configuration   #
############################
my $wb;
if ($store) { $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wb = Wormbase->new( -debug => $debug,-test=>$test) }
$debug = $wb->debug if $wb->debug;    # Debug mode, output only goes to one

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  $maintainers = $debug.'\@sanger.ac.uk';
}
my $tace = $wb->tace;
my $log = Log_files->make_build_log($wb);

# Set up top level base directory which is different if in test mode
# Make all other directories relative to this
my $basedir   = $wb->basedir;

##########################################################
#
# Main part of script
#
###########################################################
# simple counted for number of names assigned
my $counter = 0; 

# output file
open (OUT, ">$basedir/autoace/acefiles/molecular_names_for_genes.ace") or die "Can't write ace file";

# open tace pipe, and connect to AceDB using TableMaker
my $command="Table-maker -p $basedir/autoace/wquery/molecular_names_for_genes.def\nquit\n";
$log->write_to("Finding molecular names, using Table-maker...\n");

open (TACE, "echo '$command' | $tace $basedir/autoace |");
while (<TACE>) {
  chomp;
  next if ($_ eq "");
  last if (/\/\//);
  next if ($_ =~ m/^acedb>/);

  # get rid of quote marks
  s/\"//g;

  # split the line into various fields
  my ($gene,$cds,$transcript,$pseudogene,$protein, $coding_transcript) = split(/\t/, $_) ;

  # write output file


  if($cds){
    print OUT "Gene : \"$gene\"\n";
    print OUT "Molecular_name \"$cds\"\n\n";
    $counter++;
  }
  if($transcript){
    print OUT "Gene : \"$gene\"\n";
    print OUT "Molecular_name \"$transcript\"\n\n";
    $counter++;
  }
  if($pseudogene){
    print OUT "Gene : \"$gene\"\n";
    print OUT "Molecular_name \"$pseudogene\"\n\n";
    $counter++;
  }
  if($protein){
    print OUT "Gene : \"$gene\"\n";
    print OUT "Molecular_name \"$coding_transcript\n" if ($coding_transcript =~ /\w+/);
    print OUT "Molecular_name \"$protein\"\n";
    $counter++; # coding_transcripts will make this wrong as they get outputted on the same line as the protein.  Ho hum- Im sure no-one will notice!

    # also capture version without WP: prefix
    $protein =~ s/^WP\://;
    print OUT "Molecular_name \"$protein\"\n\n";
    $counter++;
  }
}
close TACE; 
close OUT;

$log->write_to("Found $counter molecular names for genes\n");

# load file to autoace using autoace_minder.pl -load
$log->write_to("Loading $basedir/autoace/molecular_names_for_genes.ace to autoace\n");
$command = "autoace_builder.pl -load $basedir/autoace/acefiles/molecular_names_for_genes.ace -tsuser molecular_names -database $basedir/autoace";

my $status = system($command);
if(($status >>8) != 0){
  $log->write_to("ERROR: Loading failed \$\? = $status\n");
}

# tidy up and exit
$log->mail("$maintainers");
exit(0);

__END__

=pod

=head1 NAME - molecular_names_for_genes.pl

=head2 DESCRIPTION


Simple file that uses a table-maker definition file to find all CDS, Transcript, and Pseudogene names
for a given ?Gene object and adds these to the Molecular_name tag in the ?Gene model.  It also adds
WormPep object names, and the raw wormpep accession (i.e. it strips off the WP: part).

This data is then loaded to autoace using autoace_minder.pl -load


Mandatory arguments:
                                                                                                              
=over 4
                                                                                                              
=item none                                                                                                              
=back
                                                                                                              
autoace_minder.pl OPTIONAL arguments:
                                                                                                              
=over 4
                                                                                                              
=item -debug <user>, email logfile goes to user rather than everyone
                                                                                                              
=item -test, uses test environment under ~wormpub/TEST_BUILD rather than autoace
                                                                                                              

=head1 AUTHOR Keith Bradnam (krb@sanger.ac.uk)
 
=back
 
=cut
