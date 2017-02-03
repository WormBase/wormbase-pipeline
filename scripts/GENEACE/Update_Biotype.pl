#!/software/bin/perl -w
#
# Update_Biotype.pl                           
# 
# by pad                         
#
# This is a example of a good script templateScript to calculate and refresh the BioType of all genes.
#

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
#use Ace;
#use Sequence_extract;
#use Coords_converter;§

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $single, $build, $species, $outfile, $report);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "single:s"   => \$single,
	    "build"      => \$build,
	    "outfile:s"  => \$outfile, #full path to output file
	    "report:s"   => \$report,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -species => $species,
			     );
}

my $output;

# establish log file.
my $log = Log_files->make_build_log($wormbase);

if ($outfile) {
  $output = $outfile;
}
else {
  my $datestring = localtime();
  $output = "/tmp/${species}_$datestring";
}
open (ACE,">$output") or $log->log_and_die("can't write output: $!\n");

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  $log->write_to("In test mode\n") if ($verbose);

}



#################################
# Set up some useful paths      #
#################################

# Set up top level base directories (these are different if in test mode)
my $basedir         = $wormbase->basedir;     # BASE DIR
my $ace_dir         = $wormbase->autoace;     # AUTOACE DATABASE DIR
my $wormpep_dir     = $wormbase->wormpep;     # CURRENT WORMPEP
my $wormrna_dir     = $wormbase->wormrna;     # CURRENT WORMRNA
my $common_data_dir = $wormbase->common_data; # AUTOACE COMMON_DATA
my $chromosomes_dir = $wormbase->chromosomes; # AUTOACE CHROMSOMES
my $reports_dir     = $wormbase->reports;     # AUTOACE REPORTS
my $gff_dir         = $wormbase->gff;         # AUTOACE GFF
my $gff_splits_dir  = $wormbase->gff_splits;  # AUTOACE GFF SPLIT
my $logs_dir        = $wormbase->logs;        # AUTOACE LOGS

# some database paths
my ($seqdb, $geneace);
if (defined $build) {
  $seqdb = $wormbase->autoace;
  $geneace = $wormbase->autoace;
}
else {
  if ($species eq 'elegans') {
    $seqdb = $wormbase->database('camace');
  }
  else {
    $seqdb = $wormbase->database($species);
  }
  $geneace   = $wormbase->database('geneace');
}

# other paths
my $tace            = $wormbase->tace;        # TACE PATH
my $giface          = $wormbase->giface;      # GIFACE PATH

#ncRNA BioType mapping

my %rna2SOlookup = (
'asRNA'   => 'SO:0000655',  # ncRNA_gene
'lincRNA' => 'SO:0001641',  # lincRNA_gene
'miRNA'   => 'SO:0001265',  # miRNA_gene
'ncRNA'   => 'SO:0000655',  # ncRNA_gene
'piRNA'   => 'SO:0001638',  # piRNA_gene
'rRNA'    => 'SO:0001637',  # rRNA_gene
'scRNA'   => 'SO:0001266',  # scRNA_gene
'snoRNA'  => 'SO:0001267',  # snoRNA_gene
'snRNA'   => 'SO:0001268',  # snRNA_gene
'tRNA'    => 'SO:0001272',  # tRNA_gene
);


##########################
# MAIN BODY OF SCRIPT
##########################

# main stuff goes here
my $db = Ace->connect(-path=>$seqdb) or  $log->log_and_die("Couldn't connect to $seqdb\n". Ace->error);
my $gdb =  Ace->connect(-path=>$geneace) or  $log->log_and_die("Couldn't connect to $geneace\n". Ace->error);

$log->write_to("Using : $db for annotation data\n Using : $gdb for primary gene data.\n\n");

# example of running anther script
#$wormbase->run_script("other_script -options", $log);
my $command;
$command .= "query find Genes_$species\n";
$command .= "quit\n";


my $ccount;
my @SeqGenes;
if ($single) {
@SeqGenes = $db->fetch (-query => "FIND Gene $single");
}
else {
  @SeqGenes = $db->fetch (-query => "FIND Gene");
}
my @Generef;
my %biotype;
foreach my $SeqGene(@SeqGenes) {
  $ccount ++;
  print "Checking".$SeqGene->name."\n" if ($verbose);
  push (@Generef,"$SeqGene->name");
  if ($SeqGene->Corresponding_CDS){
    $biotype{$SeqGene} = 'SO:0001217';
    print "$SeqGene - Coding\n" if ($verbose);
    next;
  }
  if ($SeqGene->Corresponding_Pseudogene){
    $biotype{$SeqGene} = 'SO:0000336';
    print "$SeqGene - Pseudogene\n" if ($verbose);
    next;
  }
  if ($SeqGene->Corresponding_Transcript){
    if (defined $biotype{$SeqGene}){
      $log->write_to("WARNING: $SeqGene Already defined as coding\n");
      next;
    }
    else {
      my $Transcript = $SeqGene->Corresponding_Transcript->name; 
      my @Transcript_obj = $db->fetch (-query => "FIND Transcript $Transcript");
      my $type = $Transcript_obj[0]->Transcript->name;
      if (defined $rna2SOlookup{$type}) {
	$biotype{$SeqGene} = $rna2SOlookup{$type};
	print "$SeqGene - nonCoding with $rna2SOlookup{$type} ($type)\n" if ($verbose);
      }
      else {
	$biotype{$SeqGene} = 'SO:0001263';
	print "$SeqGene - nonCoding\n" if ($verbose);
      }
      next
    }
  }
  else {
    $log->write_to("Warning $SeqGene - No live annotation attached to this gene.\n") if ($verbose);
    if ($SeqGene->Corresponding_CDS_history) { 
      $log->write_to("Checked $SeqGene - Has History annotation associated with it.\n") if ($verbose);
    }
    else {$log->write_to("Error - $SeqGene appears to be dangling, why is it in the database $db?\n");}
    next;
  }
}


#Now query Geneace
my $gcount = 0;
my $mcount = 0;
my $bcount = 0;
my $qcount = 0;
my $ucount = 0;
my $missedcount = 0;
my @WBGenes;
if ($single) {
  @WBGenes = $gdb->fetch (-query => "FIND Gene WBGene00014836");
}
else {
  @WBGenes = $gdb->fetch (-query => "FIND Genes_$species WHERE Live");
}
foreach my $WBGene(@WBGenes) {
  $gcount ++;
  my $finalbio;
  my $WBGeneBio;
  if (defined $biotype{$WBGene}) {
    #I have a calculated biotype so I'm going to use it#
    $finalbio = $biotype{$WBGene};
    if ($WBGene->Biotype) {
      $WBGeneBio = $WBGene->Biotype;
      #test just for the hell of it#
      if ($WBGeneBio eq $finalbio) {
	$mcount++;
	$log->write_to("$WBGene - MATCH\n");
      }
      elsif ($WBGeneBio ne $finalbio) {
	$log->write_to("$WBGene - NO-MATCH Calculated $finalbio Geneace:$WBGeneBio\n");
	print ACE "Gene : $WBGene\nBiotype $finalbio\n\n";
	$bcount++;
      }
      else {
	$log->write_to("Warning missed round 1: $WBGene\n");
	$missedcount++;
      }
    }
    else { 
      $log->write_to("Warning missed round 2: $WBGene\n");
      # Need a biotype
      print ACE "Gene : \"$WBGene\"\nBiotype \"$biotype{$WBGene}\"\n\n";
      $missedcount++;






    }
  }
  else {
  #  no calculated biotype 
    if ($WBGene->Biotype) {
      $log->write_to("Warning $WBGene previously had a biotype of ".$WBGene->Biotype."\n");
      $qcount++;
    }
    else {
      $ucount++;
      $finalbio = "SO:0001867";
    }
  }
  #print "Gene : $WBGene\nBiotype $finalbio\n";
  next;
}


# Close log files and exit
$log->write_to("\n\nStatistics\n");
$log->write_to("----------\n\n");
$log->write_to("$ccount objects retrieved from $seqdb\n");
$log->write_to("$gcount $species genes retrieved from $geneace\n");
$log->write_to("$mcount genes have matching biotypes\n");
$log->write_to("$bcount genes have mis-matched biotypes\n");
$log->write_to("$qcount genes previously had a biotype, but not now\n");
$log->write_to("$ucount genes appear not to be cloned\n");
$log->write_to("$missedcount genes not dealt with\n\n");
$log->mail();
close(ACE);
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




# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - script_template.pl

=head1 USAGE

=over 4

=item script_template.pl  [-options]

=back

This script does...blah blah blah

script_template.pl MANDATORY arguments:

=over 4

=item None at present.

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item xxx (xxx@sanger.ac.uk)

=back

=cut
