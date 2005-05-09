#!/usr/local/bin/perl5.8.0 -w
#
# GeneID_updater.pl
# 
# by Paul Davis
#
# Script to do soemthing else quite different indeed.
#
# Last updated by: $Author: dl1 $
# Last updated on: $Date: 2005-05-09 14:15:13 $

use strict;                                      
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $targetDB, $sourceDB, $update, $all);
my $maintainers = "All";
my $log;
#our $logdir  = "/nfs/team71/worm/pad/Scripts";


GetOptions (
	    'help'         => \$help,
            'debug=s'      => \$debug,
	    'sourceDB=s'   => \$sourceDB,
	    'targetDB=s'   => \$targetDB,
	    'update'       => \$update,
	    );	     

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

# default databases
$sourceDB = "/wormsrv1/geneace" if (!$sourceDB);
$targetDB = "/wormsrv1/camace"  if (!$targetDB);

# tace executable path
our $tace = &tace; 


# paths and files
my $tablemaker_query =  "/nfs/disk100/wormpub/DATABASES/current_DB/wquery/SCRIPT:GeneID_updater.pl.def";
my $output_file      =  "/nfs/disk100/wormpub/camace_orig/updated_geneIDs.ace";

########################################################################
# make database connection for extracting WBGene<=>Sequence_name data  #
# Method   : tace/tablemaker                                           #
# SourceDB : /wormsrv1/geneace                                         # 
########################################################################

my %models2geneID; 

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
    last if (/\/\//);
    
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

print LOG "// Gathering data from $targetDB; Object Class and Name\n";

open (OUT, ">$output_file") or die "Can't write output .acefile: $output_file\n";

# Retrieve model names and assign to a hash#
my @classes = ('elegans_CDS', 'elegans_pseudogenes', 'elegans_RNA_genes'); 

my %tags = (
	    'elegans_CDS'         => 'CDS',
	    'elegans_pseudogenes' => 'Pseudogene',
	    'elegans_RNA_genes'   => 'Transcript'
	    );

my $history;
my $lookupname;
my $geneID;
my $object;
my $query;

my $db = Ace->connect(-path => "$targetDB",
		      -program => $tace) || do { print LOG "Connection failed to $targetDB: ",Ace->error; die();};

foreach my $class (@classes) {

    $query = "find $class";
    my $i = $db->fetch_many(-query => $query);

    while (my $obj = $i->next) {
	
	$name = $obj;

        # histories AH6.1:wp999
	if ($name =~ /(\S+.+)\:(\S+.+)/) {
	    $lookupname = $1;                 
	    $history    = 1;
	}
	else {
	    $lookupname = $name;
	}

        # remove isoform letters from model names
	$lookupname =~ s/[a-z]$//;

        # lookup GeneID keyed on lookupname
	$geneID = $models2geneID{$lookupname};
		
        # Print out ace file full model name and modified class name.

#	print "// $name \t$lookupname\n";             # verbose debug line
	
	print OUT "$tags{$class} : \"$name\"\n";
	unless ($history) {
	    print OUT "Gene\t\"$geneID\"\n\n";
	}
	else {
	    print OUT "Gene_history\t\"$geneID\"\n\n";
	}
	$obj->DESTROY();
    }
}
$db->close;

close OUT;

print LOG "// upload file complete\n";

###############################################################
# Delete targetDB Gene ID info & load new ones into $targetdb #
###############################################################

if ($update) {

    my $command = "query find elegans_CDS\nedit -D Gene\nclear\n";
    $command   .= "query find elegans_pseudogenes\nedit -D Gene\nclear\n";
    $command   .= "query find elegans_RNA_genes\nedit -D Gene\nclear\n";
    $command   .= "pparse $output_file\n";
    $command   .= "save\nquit\n";
    
    open (DB, "| $tace $targetDB |") || die "Couldn't open $targetDB\n";
    print DB $command;
    close DB;
    
    print LOG "// synchronised data in $targetDB\n";
}


###############
# hasta luego #
###############

&mail_maintainer("GeneID_updater.pl",$maintainers,$log);
close(LOG);

exit(0);


##############################################################
#
# Subroutines
#
##############################################################

sub create_log_files{

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







__END__

    =pod

    =head2 NAME - GeneID_updater.pl

    =head1 USAGE

    =over 4

    =item GeneID_updater.pl [-options]

    =back

    This script does...blah blah blah

    script_template.pl MANDATORY arguments:

    =over 4

    =item None,

    =back

    Update_checker.pl  OPTIONAL arguments:

    =over 4

    =item -h, Help

    =item -sourceDB, database from which you wish to synchronise

    =item -targetDB, database you wish to synchronise

    =item -update updates target database with new gene IDs

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
