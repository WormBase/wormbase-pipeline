#!/usr/local/bin/perl5.8.0 -w
#
# molecular_names_for_genes.pl
# 
# by Keith Bradnam
#
# quick script to populate Molecular_name tag in ?Gene model during build
#
# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2012-05-11 15:29:59 $

#################################################################################
# Initialise variables                                                          #
#################################################################################

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

##############################
# command-line options       #
##############################
my ($help, $debug, $species, $test, $verbose, $store, $wormbase, $outfile, $no_load);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"      => \$store,
	    "species:s"	=> \$species,
            "outfile:s"  => \$outfile,
            "noload"    => \$no_load,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species
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




############################
# Get paths
############################

my $tace = $wormbase->tace;

# Set up top level base directory which is different if in test mode
# Make all other directories relative to this
my $basedir   = $wormbase->basedir;
my $db_dir    = $wormbase->autoace;

if (not defined $outfile) {
  $outfile = $wormbase->acefiles . "/molecular_names_for_genes.ace";
}

##########################################################
#
# Main part of script
#
###########################################################
# simple counted for number of names assigned
my $counter = 0; 

# output file
open (OUT, ">$outfile") or $log->log_and_die("Can't write ace file $outfile");

# open tace pipe, and connect to AceDB using TableMaker
my $command="Table-maker -p $db_dir/wquery/molecular_names_for_genes.def\nquit\n";
$log->write_to("Finding molecular names, using Table-maker...\n");

my %added_for_gene;

open (TACE, "echo '$command' | $tace $db_dir |");
while (<TACE>) {
  chomp;
  next if ($_ eq "");
  last if (/\/\//);
  next if ($_ =~ m/^acedb>/);

  # get rid of quote marks
  s/\"//g;

  # split the line into various fields
  my ($gene,$cds,$transcript,$pseudogene,$protein) = split(/\t/, $_) ;

  # write output file
  $added_for_gene{$gene} = {} if not exists $added_for_gene{$gene};

  if($cds and not exists $added_for_gene{$gene}->{$cds}){
    print OUT "Gene : \"$gene\"\n";
    print OUT "Molecular_name \"$cds\"\n\n";
    $added_for_gene{$gene}->{$cds} = 1;
    $counter++;
  }
  if($transcript and not exists $added_for_gene{$gene}->{$transcript}){
    print OUT "Gene : \"$gene\"\n";
    print OUT "Molecular_name \"$transcript\"\n\n";
    $added_for_gene{$gene}->{$transcript} = 1;
    $counter++;
  }
  if($pseudogene and not exists $added_for_gene{$gene}->{$pseudogene}){
    print OUT "Gene : \"$gene\"\n";
    print OUT "Molecular_name \"$pseudogene\"\n\n";
    $added_for_gene{$gene}->{$pseudogene} = 1;
    $counter++;
  }
  if($protein and not exists $added_for_gene{$gene}->{$protein}){
    print OUT "Gene : \"$gene\"\n";
    print OUT "Molecular_name \"$protein\"\n";
    $added_for_gene{$gene}->{$protein} = 1;
    $counter++;

    # also capture version without WP: prefix
    $protein =~ s/^WP\://;
    if (not exists $added_for_gene{$gene}->{$protein}) {
      print OUT "Molecular_name \"$protein\"\n";
      $added_for_gene{$gene}->{$protein} = 1;
      $counter++;
    }
    print OUT "\n";
  }
}
close TACE; 
close OUT;

$log->write_to("Found $counter molecular names for genes\n");

if (not $no_load) {
  # load file to autoace
  $wormbase->load_to_database($db_dir, $outfile, "molecular_names", $log);
}

# tidy up and exit
$log->mail();
print "Finished.\n" if ($verbose);
exit(0);




##############################################################
#
# Subroutines
#
##############################################################



##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################



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
