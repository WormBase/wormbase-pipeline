#!/usr/local/bin/perl5.8.0 -w
#
# GeneID_updater.pl
# 
# by Paul Davis
#
# Script to refresh the CDS->WBGene connections in a chosen database from a chosen reference database.
# Script also refreshes Protein_IDs in the chosen database from the latest build.
#
# Last updated by: $Author: pad $
# Last updated on: $Date: 2005-09-13 12:47:03 $

use strict;                                      
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $geneID, $targetDB, $sourceDB, $fileout, $update, $public, $sourceDB2, $proteinID);

GetOptions (
	    'help'         => \$help,
            'debug=s'      => \$debug,
	    'geneID'       => \$geneID,
	    'sourceDB=s'   => \$sourceDB,
    	    'sourceDB2=s'  => \$sourceDB2,
	    'targetDB=s'   => \$targetDB,
	    'fileout=s'    => \$fileout,
	    'proteinID'    => \$proteinID,
	    'update'       => \$update,
	    'public'       => \$public,
	    );

my $maintainers = "All";
my $log;
my $output_file;
my $output_file2;
my %models2geneID;

# tace executable path
our $tace = &tace; 

# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if($debug){
    ($maintainers = $debug . '\@sanger.ac.uk');
}

&create_log_files;

##########################
# MAIN BODY OF SCRIPT
##########################

# Specify Paths and files
my $def_dir = "/nfs/disk100/wormpub/DATABASES/current_DB/wquery";
my $tablemaker_query =  "${def_dir}/SCRIPT:GeneID_updater.pl.def";
my $tablemaker_query2 = "${def_dir}/SCRIPT:HXcds2protID.def";

# Select and check source databases.
if ($sourceDB) {
  $sourceDB = $sourceDB;
}
elsif (!$sourceDB) {
  $sourceDB = "/wormsrv1/geneace";
}
if ($targetDB) {
  $targetDB = $targetDB;
}
elsif (!$targetDB) {
#$targetDB = "/nfs/disk100/wormpub/camace_pad";
$targetDB = "/wormsrv1/camace";
}
if ($proteinID) {
  $sourceDB2 ="/nfs/disk100/wormpub/DATABASES/current_DB";
}

#Check and set output file
if ($fileout) {
  $output_file = $fileout;
}
elsif (!$fileout) {
  $output_file = "/nfs/disk100/wormpub/camace_orig/updated_geneIDs.ace";
}

########################################################################
# make database connection for extracting WBGene<=>Sequence_name data  #
# Method   : tace/tablemaker                                           #
# SourceDB : /wormsrv1/geneace                                         # 
########################################################################

if ($geneID) {
print LOG "-geneID option selected, therefore Gene connections will be updated in $targetDB\n";
print "\n\nSOURCE Database for Gene IDs: $sourceDB.\n\n";
print "TARGET Database: $targetDB\n\n";
print "\nYou are creating the file $output_file.\n\n" if ($geneID);
print LOG "// Gathering data from $sourceDB; Object Class and GeneID\n";
    
# connect to AceDB using TableMaker, but use /wormsrv1/geneace for Table-maker definition
my $command = "Table-maker -p $tablemaker_query\nquit\n";
my $gene;
my $name;

open (TACE, "echo '$command' | $tace $sourceDB |");
while (<TACE>) {
    chomp;
    next if ($_ eq "");
    next if (/acedb\>/);
    #last if (/\/\//);
    
    # get rid of quote marks
    s/\"//g;

    next unless  (/^(\S+)\t(\S+)/);
    
    # split the line into various fields
    ($gene,$name) = (/^(\S+)\t(\S+)/);

    # add to hash. CDS, Pseudogene, or Transcript name is key, gene ID is value
    $models2geneID{$name} = $gene;

}
close TACE;

####################################################################
# make database connections for looping through elegans subclasses #
# Method   : AcePerl                                               # 
# TargetDB : /wormsrv1/camace                                      # 
####################################################################

print LOG "// Gathering data from $targetDB; Object Class and Name\n\n";
open (OUT, ">$output_file") or die "Can't write output .acefile: $output_file\n";

# Retrieve model names and assign to a hash#
my @classes = ('All_genes');

#my $history;
my $lookupname;
my $geneID;
my $query;

my $db = Ace->connect(-path => "$targetDB",
		      -program => $tace) || do { print LOG "Connection failed to $targetDB: ",Ace->error; die();};

foreach my $class (@classes) {

    $query = "find $class";
    my $i = $db->fetch_many(-query => $query);

    while (my $obj = $i->next) {
	
	$name = $obj->name;
	my $history;
        # histories AH6.1:wp999
	if ($name =~ /(\S+.+)\:(\S+.+)/) {
	    $lookupname = $1;
	    $history    = 1;
	}
	else {
	    $lookupname = $name;
	}

        # remove isoform letters from models with a cosmid.no name
	$lookupname =~ s/[a-z]$// unless ($lookupname =~ /\S+\.\D/);

        # lookup GeneID keyed on lookupname
	$geneID = $models2geneID{$lookupname};
		
        # Print out ace file full model name and modified class name.
	my $Tag = $obj->class;

#	print "// $name \t$lookupname\n";             # verbose debug line
	if ($Tag ne "Transposon")	{
	  print OUT "$Tag : \"$name\"\n";
	  if (!$history) {
	    print OUT "Gene\t\"$geneID\"\n\n";
	  }
	  elsif (defined $history) {
	  #print OUT "Gene_history\t\"$geneID\"\n\n" if (defined ($geneID));
	    print OUT "Gene_history\t\"$geneID\"\n\n" unless (!defined ($geneID));
	  }
	  print OUT "\n" if (!defined ($geneID));
	  print LOG "$name does not have a geneID please investigate.\n" if (!defined ($geneID));
	  $obj->DESTROY();
	}
      }
  }
$db->close;

close OUT;
}

########################
#  Refresh ProteinIDs  #
########################
if ($proteinID) {
  $output_file2 = "/nfs/disk100/wormpub/camace_orig/updated_proteinIDs.ace";
  print LOG "\n\n//-proteinID option selected, therefore Protein_ID connections will be updated in $targetDB\n";
  print "SOURCE Database for Protein IDs: $sourceDB2.\n\n";
  print LOG "\n//You are creating $output_file2\n"; 
  print LOG "// Gathering Protein_IDs from $sourceDB2\n";

# connect to AceDB using TableMaker
my $command = "query find curated_CDS where From_laboratory = HX\nshow -t Protein_id -a -f $output_file2\nquit\n";
system ("echo \"$command\" | $tace $sourceDB2");

##Hard Coded exceptions to add on the end of the ace file
open (OUT2, ">>$output_file2") or die "Can't write output .acefile: $output_file2\n";
print OUT2 "\nCDS	:	\"MTCE.3\"\n";
print OUT2 "protein_id	\"MTCE\"  \"CAA38153.1\"\n";
print OUT2 "\nCDS	:	\"MTCE.12\"\n";
print OUT2 "protein_id	\"MTCE\"  \"CAA38154.1\"\n";
print OUT2 "\nCDS	:	\"MTCE.16\"\n";
print OUT2 "protein_id	\"MTCE\"  \"CAA38155.1\"\n";
print OUT2 "\nCDS	:	\"MTCE.21\"\n";
print OUT2 "protein_id	\"MTCE\"  \"CAA38156.1\"\n";
print OUT2 "\nCDS	:	\"MTCE.23\"\n";
print OUT2 "protein_id	\"MTCE\"  \"CAA38157.1\"\n";
print OUT2 "\nCDS	:	\"MTCE.25\" \n";
print OUT2 "protein_id	\"MTCE\"  \"CAA38158.1\"\n";
print OUT2 "\nCDS	:	\"MTCE.26\"\n";
print OUT2 "protein_id	\"MTCE\"  \"CAA38159.1\"\n";
print OUT2 "\nCDS	:	\"MTCE.31\"\n";
print OUT2 "protein_id	\"MTCE\"  \"CAA38160.1\"\n";
print OUT2 "\nCDS	:	\"MTCE.34\"\n";
print OUT2 "protein_id	\"MTCE\"  \"CAA38161.1\"\n";
print OUT2 "\nCDS	:	\"MTCE.35\"\n";
print OUT2 "protein_id	\"MTCE\"  \"CAA38162.1\"\n";
#  $db->close;

  close OUT2;
}

##Logging##
print LOG "// upload file(s) complete\n";

if (!defined ($update)) {
  print LOG "// **You will have to manually load this data into /wormsrv1/camace.**\n";
  print LOG "// **Also load /nfs/disk100/wormpub/camace_orig/geneID_patch.ace for the exceptions**\n\n";
}

###############################################################
# Delete targetDB ID info & load new ones into $targetdb      #
###############################################################

if (($update) && ($geneID)) {
  &load_geneIDs;
}
elsif (($update) && ($proteinID)) {
  &load_protein_IDs;
}
elsif ($update) {
  &load_geneIDs;
  &load_protein_IDs;
}

print "Diaskeda same Poli\n"; #we had alot of fun#

&mail_maintainer("GeneID_updater.pl",$maintainers,$log);
close(LOG);

exit(0);


##############################################################
# Subroutines
##############################################################

sub create_log_files {

    # Create history logfile for script activity analysis
    $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

    # create main log file using script name for
    my $script_name = $1;
    $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
    my $rundate     = `date +%y%m%d`; chomp $rundate;
    $log        = "/wormsrv2/logs/$script_name.$rundate.$$";

    open (LOG, ">$log") or die "cant open $log";
    print LOG "$script_name\n";
    print LOG "started at ",`date`,"\n";
    print LOG "=============================================\n";
    print LOG "\n";
}

sub load_geneIDs {
  my $command = "query find curated_CDS\nedit -D Gene\nclear\n";
  $command   .= "query find elegans_pseudogenes\nedit -D Gene\nclear\n";
  $command   .= "query find elegans_RNA_genes\nedit -D Gene\nclear\n";
  $command   .= "pparse $output_file\n";
  $command   .= "pparse /nfs/disk100/wormpub/camace_orig/geneID_patch.ace\n";
  $command   .= "save\nquit\n";
  
  open (DB, "| $tace $targetDB -tsuser merge_split_update_geneID |") || die "Couldn't open $targetDB\n";
  print DB $command;
  close DB;
  
  print LOG "// updated Gene data in $targetDB\n";
}

sub load_protein_IDs {
    my $command = "query find curated_CDS\nedit -D Protein_ID\nclear\n";
    $command   .= "pparse /nfs/disk100/wormpub/camace_orig/updated_proteinIDs.ace\n";
    $command   .= "save\nquit\n";
    
    open (DB, "| $tace $targetDB -tsuser merge_split_update_proteinID |") || die "Couldn't open $targetDB\n";
    print DB $command;
    close DB;

    print LOG "// updated protein_ID data in $targetDB\n";
  }

sub usage 
  {
    my $error = shift;
    if ($error eq "Help") {
      # Help menu
      exec ('perldoc',$0);
    }
  }

__END__

=pod

=head2 NAME - GeneID_updater.pl

=head1 USAGE

=over 4

=item GeneID_updater.pl [-options]

=back

This script removes all prediction->WBGene id connections from a chosen target database and re-synchronises these connections with a chosen reference database.  Defaults are in place for routine updating of /wormsrv1/camace with data from wormsrv1/geneace.  The script can also populate the Public_name tags in the newly syncronised WBGene objects to aid curators when chosing which CDS should be renamed following a gene split event.  The script has also been modified to retrieve protein ids from the previous build and update these in a chosen target database ready for emnl dumping in the next build.


script_template.pl MANDATORY arguments:

=over 4

=item None,

=back

Update_checker.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=item -dubug, supply your user ID any you will be the only person to recieve the Log email.

=item -geneID, specifies that gene ids are to be updated.

=item -proteinID, specifies that protein ids are to be updated.

=item -update, Load the data automatically into targetDB.

=item -sourceDB, database from which you wish to retrieve gene_id data.

=item -sourceDB2, database from which you wish to retrieve protein data.

=item -targetDB, database you wish to synchronise.

=public -public, populates the public name data tag in gene objects in target database.

=back

=head1 REQUIREMENTS

=over 4

=item This script needs to run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Paul Davis (pad@sanger.ac.uk)

=back

=cut

