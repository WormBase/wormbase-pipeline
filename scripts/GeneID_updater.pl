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
# Last updated on: $Date: 2006-03-21 10:35:11 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Storable;

######################################
# variables and command-line options #
######################################

my ($help, $debug, $geneID, $targetDB, $sourceDB, $fileout, $update, $public, $sourceDB2, $proteinID, $version, $test, $store, $wormbase);

GetOptions (
	    'help'         => \$help,       #help documentation.
	    'test'         => \$test,       #test build
            'debug=s'      => \$debug,      #debug option for email
	    'geneID'       => \$geneID,     #update gene id's
	    'sourceDB=s'   => \$sourceDB,   #source for gene id info
    	    'sourceDB2=s'  => \$sourceDB2,  #source for protein id's
	    'targetDB=s'   => \$targetDB,   #target db
	    'fileout=s'    => \$fileout,    #specify an output file
	    'proteinID'    => \$proteinID,  #update protein id's option
	    'update'       => \$update,     #load ace files automatically into target db.
	    'public'       => \$public,     #retrieve public name info ....future
	    'version=s'    => \$version,    #version number for properly directing out files in future
	    'store:s'      => \$store,      #storable object
	    );


if ($store) {
  $wormbase = retrieve($store) or croak ("Can't restore wormbase from $store\n");
} 
else {
  $wormbase = Wormbase->new( -debug => $debug,
			     -test => $test,
			     );
}

my $output_file;
my $output_file2;
my %models2geneID;
my $version_no = ($version + 1);

# tace executable path
my $tace = $wormbase->tace;

# Display help if required
&usage("Help") if ($help);

# establish log file.
my $log = Log_files->make_build_log($wormbase);

##########################
# MAIN BODY OF SCRIPT
##########################

# Specify Paths and files
my $def_dir = "/nfs/disk100/wormpub/DATABASES/current_DB/wquery";
my $tablemaker_query =  "${def_dir}/SCRIPT:GeneID_updater.pl.def";
my $tablemaker_query2 = "${def_dir}/SCRIPT:HXcds2protID.def";
my $canonical = "/nfs/disk100/wormpub/DATABASES/camace";

# Select and check source databases.
if ($sourceDB) {
  $sourceDB = $sourceDB;
}
elsif (!$sourceDB) {
  $sourceDB = "/nfs/disk100/wormpub/DATABASES/geneace";
}
if ($targetDB) {
  $targetDB = $targetDB;
}
elsif (!$targetDB) {
#$targetDB = "/nfs/disk100/wormpub/camace_pad";
$targetDB = "$canonical";
}
if ($proteinID) {
  $sourceDB2 ="/nfs/disk100/wormpub/DATABASES/current_DB";
}

#Check and set output file
if ($fileout) {
  $output_file = $fileout;
}
elsif (!$fileout) {
  $output_file = "/nfs/disk100/wormpub/camace_orig/acefiles/updated_geneIDs_WS${version_no}.ace";
}

########################################################################
# make database connection for extracting WBGene<=>Sequence_name data  #
# Method   : tace/tablemaker                                           #
# SourceDB : /nfs/disk100/wormpub/DATABASES/geneace                    # 
########################################################################

if ($geneID) {
$log->write_to("\n==============================================================================================
geneID option selected, updating WBGeneID connections in $targetDB
----------------------------------------------------------------------------------------------\n");
$log->write_to("\nSOURCE Database for Gene IDs: $sourceDB.\n");
$log->write_to("TARGET Database: $targetDB\n");
$log->write_to("OUTPUT FILE:$output_file.\n\n") if ($geneID);
$log->write_to("Opening database connection.....\n/1/ Gathering data from $sourceDB; Object Class and GeneID\n");
    
# connect to AceDB using TableMaker, but use /nfs/disk100/wormpub/DATABASES/current_DB/wquery for Table-maker definition
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
# TargetDB : /nfs/disk100/wormpub/DATABASES/camace                 # 
####################################################################

$log->write_to("/2/ Gathering data from $targetDB; Object Class and Name\n\n");
open (OUT, ">$output_file") or die "Can't write output .acefile: $output_file\n";

# Retrieve model names and assign to a hash#
my @classes = ('All_genes');

#my $history;
my $lookupname;
my $geneID;
my $query;

my $db = Ace->connect(-path => "$targetDB",
		      -program => $tace) || do { 
$log->write_to("Connection failed to $targetDB: ",Ace->error);
die();
};

$log->write_to("==============================================================================================\nERROR TABLE\n==============================================================================================\n");
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
	  if (!defined ($geneID)) {$log->write_to("ERROR:$name does not have a geneID please investigate.\n");}
	  $obj->DESTROY();
	}
	#$log->write_to("Known errors\n<==========>C54G4.7:yk713c3.mRNA:wp126\nC54G4.7:yk728f5.mRNA:wp126\nF28A8.9:wp149\nY105E8A.30a:wp128\nY105E8A.30b:wp128\nZK228.9:wp144\n");
      }
  }
$log->write_to("==============================================================================================\nEXCEPTIONS to be ignored from above table.....\n==============================================================================================\nC54G4.7:yk713c3.mRNA:wp126\nC54G4.7:yk728f5.mRNA:wp126\nF28A8.9:wp149\nY105E8A.30a:wp128\nY105E8A.30b:wp128\nZK228.9:wp144\n==============================================================================================");
$db->close;

close OUT;
}

########################
#  Refresh ProteinIDs  #
########################
if ($proteinID) {
  $output_file2 = "/nfs/disk100/wormpub/camace_orig/acefiles/updated_proteinIDs_WS${version_no}.ace";
  $log->write_to("\n\n==============================================================================================\nproteinID option selected, updating WP:Protein_ID connections in $targetDB\n----------------------------------------------------------------------------------------------\n\n");
  $log->write_to("SOURCE Database for Protein IDs: $sourceDB2.\n");
  $log->write_to("OUTPUT FILE: $output_file2\n\n"); 
  $log->write_to("Opening database connection.....\n/3/ Gathering Protein_IDs from $sourceDB2\n\n");

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
$log->write_to("Upload file(s) completed.......\n");

if (!defined ($update)) {
  $log->write_to("\n\n==============================================================================================\n");
  $log->write_to("ACTION\n");
  $log->write_to("==============================================================================================\n");
  $log->write_to("\n\t**You will have to manually load:\n");
  $log->write_to("\t$output_file\n") if $geneID;
  $log->write_to("\t/nfs/disk100/wormpub/camace_orig/acefiles/geneID_patch.ace\n") if $geneID;
  $log->write_to("\t$output_file2") if $proteinID;
  $log->write_to("\n\n\tinto $canonical\n\n");
  $log->write_to("\tCOMMANDS:\n\ttace $canonical -tsuser merge_split\n");
  $log->write_to("\tpparse $output_file\n") if $geneID;;
  $log->write_to("\tpparse /nfs/disk100/wormpub/camace_orig/acefiles/geneID_patch.ace\n") if $geneID;;
  $log->write_to("\tpparse $output_file2\n") if $proteinID;
  $log->write_to("\tsave\n\tquit\n\n");
}

###############################################################
# Delete targetDB ID info & load new ones into $targetdb      #
###############################################################

if ($update) {
  $log->write_to("==============================================================================================");
  $log->write_to("ACTION");
  $log->write_to("==============================================================================================");
  $log->write_to("\nLoading files........\n");
  &load_geneIDs     if $geneID;
  &load_protein_IDs if $proteinID;
}

$log->mail();
print "Diaskeda same Poli\n"; #we had alot of fun#
exit(0);


##############################################################
# Subroutines
##############################################################

sub load_geneIDs {
  my $command = "query find curated_CDS\nedit -D Gene\nclear\n";
  $command   .= "query find elegans_pseudogenes\nedit -D Gene\nclear\n";
  $command   .= "query find elegans_RNA_genes\nedit -D Gene\nclear\n";
  $command   .= "pparse $output_file\n";
  $command   .= "pparse /nfs/disk100/wormpub/camace_orig/acefiles/geneID_patch.ace\n";
  $command   .= "save\nquit\n";
  
  open (DB, "| $tace $targetDB -tsuser merge_split_update_geneID |") || die "Couldn't open $targetDB\n";
  print DB $command;
  close DB;
  $log->write_to("\nLoaded Files:\n$output_file\n/nfs/disk100/wormpub/camace_orig/acefiles/geneID_patch.ace\n");
  $log->write_to("Updated Gene data in $targetDB\n");
}

sub load_protein_IDs {
    my $command = "query find curated_CDS\nedit -D Protein_ID\nclear\n";
    $command   .= "pparse $output_file2\n";
    $command   .= "save\nquit\n";
    open (DB, "| $tace $targetDB -tsuser merge_split_update_proteinID |") || die "Couldn't open $targetDB\n";
    print DB $command;
    close DB;
    $log->write_to("\nLoaded File:\n$output_file2\n");
    $log->write_to("Updated protein_ID data in $targetDB\n");
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

This script removes all prediction->WBGene id connections from a chosen target database and re-synchronises these connections with a chosen reference database.  Defaults are in place for routine updating of /nfs/disk100/wormpub/DATABASES/camace with data from /nfs/disk100/wormpub/DATABASES/geneace.  The script can also populate the Public_name tags in the newly syncronised WBGene objects to aid curators when chosing which CDS should be renamed following a gene split event.  The script has also been modified to retrieve protein ids from the previous build and update these in a chosen target database ready for emnl dumping in the next build.


=head2 GeneID_updater.pl MANDATORY arguments:

=over 4

=item None,

=back

GeneID_updater.pl  OPTIONAL arguments:

=over 4

=item -h, Help.

=item -dubug, supply user ID to limit logs email distribution.

=item -geneID, specifies that gene ids are to be updated.

=item -proteinID, specifies that protein ids are to be updated.

=item -update, Load the data automatically into targetDB.

=item -sourceDB, database from which you wish to retrieve gene_id data.

=item -sourceDB2, database from which you wish to retrieve protein data.

=item -targetDB, database you wish to synchronise.

=public -public, populates the public name data tag in gene objects in target database.

=back

=head1 REQUIREMENTS

=item Script no longer requires /wormsrv2.

=back

=head1 AUTHOR

=over 4

=item Paul Davis (pad@sanger.ac.uk)

=back

=cut

