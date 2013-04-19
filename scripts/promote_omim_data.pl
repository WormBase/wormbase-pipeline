#!/software/bin/perl -w
#
# promote_omim_data.pl                           
# 
# by Paul Davis                         
#
# This script promoted the OMIM disease data to the level of the gene.
# Script extended to also populate the Potential_model_for Human disease data utilising the DO ontology.
#
# Last updated by: $Author: pad $     
# Last updated on: $Date: 2013-04-19 15:34:38 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $noload, $alternativedb, $file, $obo,);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "database:s" => \$alternativedb,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "noload"     => \$noload,
	    "file:s"     => \$file,
	    "obo:s"      => \$obo,
	   );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n";
  print "No FTP files will be created\n";
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

#################################
# Set up some useful paths      #
#################################
my $ace_dir;
if ($alternativedb){
  $ace_dir         = $alternativedb;     # DATABASE DIR
}
else {
  # Set up top level base directories (these are different if in test mode)
  $ace_dir         = $wormbase->autoace;     # AUTOACE DATABASE DIR
}

my $logs_dir        = $wormbase->logs;        # AUTOACE LOGS
my $tace            = $wormbase->tace;        # TACE PATH
my $giface          = $wormbase->giface;      # GIFACE PATH
my $WS_name = $wormbase->get_wormbase_version_name();
my %species = ($wormbase->species_accessors);
$species{$wormbase->species} = $wormbase;


##########################
# MAIN BODY OF SCRIPT
##########################
my ($ace_object,%objects,%acedata,%gene2species);
my ($qspecies,$gene,$database,$id,$type,$outfile);
if ($file) {
  $outfile = $file;
}
else {
  $outfile = $ace_dir."/acefiles/omim_db_data.ace";
}

unless ($test) {
  foreach my $organism(keys %species) {
    
    
  }
}

my $count = "0";
my $countgenes = "0";
my $omdiscount = "0";
my $dodiscount = "0";
open (OUT, ">$outfile") or $log->log_and_die("Can't write ace file $outfile");

# Additional information is now required to add the Disease Ontology terms to genes where an OMIM disease ID has 
# been identified by human protein orthology.
my %omim2do;
&gatherDOdata;

my $def = '/tmp/omim.def';
open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
    my $txt = <<END;


Sortcolumn 1

Colonne 1 
Width 16 
Optional 
Visible 
Class 
Class Gene
From 1
 
Colonne 2 
Width 40 
Optional 
Hidden 
Class 
Class Protein 
From 1 
Tag Ortholog_other  
Condition *ENSP*
 
Colonne 3 
Width 12 
Mandatory 
Visible 
Class 
Class Database 
From 2 
Tag Database  
Condition "OMIM"
 
Colonne 4 
Width 12 
Mandatory 
Visible 
Class 
Class Database_field 
Right_of 3 
Tag HERE   
 
Colonne 5 
Width 12 
Optional 
Visible 
Class 
Class Accession_number 
Right_of 4 
Tag HERE   

Colonne 6 
Width 30 
Optional 
Visible 
Class 
Class Species 
From 1 
Tag Species

// End of these definitions
END
print TMP $txt;
close TMP;


my $command = "Table-maker -p $def\nquit\n";
$log->write_to("\nRetrieving OMIM, using Table-maker and query ${def}...\n");
  
open (TACE, "echo '$command' | $tace $ace_dir | ") || die "Cannot query acedb. $command  $tace\n";
while (<TACE>) {
  $count++;
  chomp;
  s/\"//g;
  if ($debug) {print "$_\n";}
  if (/(WBGene\d+)\s+(OMIM)\s+(\S+)\s+(\d+)\s+(\S+\s\S+)/){
    $gene = $1;
    $database = $2;
    $type = $3;
    $id = $4;
    $ace_object = $1;
    $qspecies = $5;
    if ($type eq "gene") {
      $objects{$ace_object}++;
      push (@{$acedata{$ace_object}{gene}},$id);
    }
    elsif ($type eq "disease") {
      $objects{$ace_object}++;
      push (@{$acedata{$ace_object}{disease}},$id);
    }
    else {
      print "WARNING $_ contains unexpected Type data \"$type\".......\n";
      $log->write_to("WARNING $_ contains unexpected Type data \"$type\".......\n");
    }
    if ($qspecies eq "Caenorhabditis elegans") {
      $gene2species{$gene} = $qspecies; # save all the names and methods
    }
  }
  else {print "$_ line does not contain OMIM data.\n" if ($debug);}
}
foreach my $obj (keys %objects) {
  if (defined $gene2species{$obj}) {
    print OUT "\n// $obj $gene2species{$obj} \n\n";
  }
  else {
    print OUT "\n// $obj \n\n";
  }
  my $line;
  print OUT "Gene : \"$obj\"\n";
  foreach $line (@{$acedata{$obj}{gene}}) {
    print OUT "Database OMIM gene $line\n";
    # Now do a DOID lookup.
    my $aceomimid = "OMIM:$line";
    if (defined $omim2do{$aceomimid}) {
      my $val;
      if ($gene2species{$obj} eq "Caenorhabditis elegans") {
	foreach $val (@{$omim2do{$aceomimid}}){
	  print OUT "Potential_model  $val \"Homo sapiens\" Inferred_automatically \"Inferred from Human protein orthology and OMIM::DO_term conversion ($aceomimid).\"\n";
	}
      }
    }
  }
  # Now print out just the elegans disease data
  foreach $line (@{$acedata{$obj}{disease}}) {
    if ($gene2species{$obj} eq "Caenorhabditis elegans") {
      my $aceomimid = "OMIM:$line";
      # Now do a DOID lookup.
      if (defined $omim2do{$aceomimid}) {
	print OUT "Database OMIM disease $line\n";
	$omdiscount ++;
	# Now itterate throught the DOIDs connected to this OMIM.
	my $val;
	foreach $val (@{$omim2do{$aceomimid}}){
	  $dodiscount++;
	  print OUT "Potential_model  $val \"Homo sapiens\" Inferred_automatically \"Inferred from Human protein orthology and OMIM::DO_term conversion ($aceomimid).\"\n";
	}
      }
      else {
	# print OMIM disease data where the DOID mapping isn't present
	print OUT "Database OMIM disease $line\n";
      }
    }
  }
}
$countgenes = (keys %objects);
close OUT;

if ($noload){
  $log->write_to("Output NOT loaded into ".$wormbase->autoace."\n");
}

unless ($noload) {
  $log->write_to("loading $outfile to ".$wormbase->autoace."\n");
  $wormbase->load_to_database($wormbase->autoace,$outfile,'promote_omim_data.pl', $log);
}

# Close log files and exit
$log->write_to("\nOutput can be found here $outfile\n");
$log->write_to("----------\n");
$log->write_to("Processed $count lines from Table-maker touching $countgenes genes\n");
$log->write_to("$omdiscount OMIM diseases with $dodiscount DO_term connections\n");

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);

##############################################################
#
# Subroutines
#
##############################################################

sub gatherDOdata {
  my $obo_file;
  if (defined $obo){
    $obo_file = $obo;
  }
  else {
    $obo_file = $wormbase->primaries . "/citace/temp_unpack_dir/home/citace/Data_for_${WS_name}/Data_for_Ontology/disease_ontology.${WS_name}.obo";
  }
  open (OBO, "<$obo_file")  || die "Can't open OBO file: $obo_file\n\n";
  my $doid;
  my $omimid;
  $log->write_to("Retrieving DO_term data, using the citace obo file...\n");
  while (<OBO>) {
    # id: DOID:0050631
    if (/^id:\s+(DOID:\d+)/) {
      $doid = $1;
    }
    #xref: OMIM:203100
    if (/^xref:\s+(OMIM:\d+)/) {
      $omimid = $1;
      $omim2do{$omimid} ||= [];
      push @{$omim2do{$omimid}},$doid;
    }
  }
  my $countomim2do = (keys %omim2do);
  $log->write_to("Collected $countomim2do  OMIM::DO_terms\nFinished getting DO data...");
}


sub usage {
  my $error = shift;
  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

__END__

=pod

=head2 NAME - promote_omim_data.pl

=head1 USAGE

=over 4

=item promote_omim_data.pl  [-options]

=back

This script queries for C. elegans genes that have a Ortholog_other that is an ensembl protein with attached OMIM disease IDs and promotes them to the Database line of the elegans gene.

promote_omim_data.pl MANDATORY arguments:

=over 4

=item None at present.

=back

promote_omim_data.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script in the test build env.

=back

=over 4
    
=item -verbose, prints to screen additional information.

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Paul Davis (paul.davis@ebi.ac.uk)

=back

=cut
