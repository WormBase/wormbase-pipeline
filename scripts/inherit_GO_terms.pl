#!/usr/local/bin/perl -w
#
# inherit_GO_terms.pl
#
# map GO_terms to ?Sequence objects from ?Motif and ?Phenotype
#


$|=1;
use IO::Handle;
use Getopt::Std;
use Ace;
use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;

##############################
# Script variables (run)     #
##############################

my $maintainers = "All";
my $rundate = `date +%y%m%d`; chomp $rundate;
my $runtime = `date +%H:%M:%S`; chomp $runtime;

# grab version number from cvs
my $version = &get_cvs_version($0);


##############################
# command-line options       #
##############################
#
# -h : help
# -d : debug
# -m : ?Motif->?Protein->?Sequence
# -p : ?Phenotype->?RNAi->?Sequence
#

use vars qw/ $opt_d $opt_h $opt_m $opt_p/;
getopts ("hdmp");
&usage(0) if ($opt_h);
my $debug = $opt_d;

# only tell Dan if running debug mode
$maintainers = "dl1\@sanger.ac.uk" if ($debug);

##############################
# Paths etc                  #
##############################

my $tace      = "/nfs/disk100/wormpub/acedb/ace4/bin.ALPHA_4/tace");      # tace executable path
my $dbpath    = "/wormsrv2/autoace";                                      # Database path

########################################
# Open logfile                         #
########################################

my $log="/wormsrv2/logs/inherit_GO_terms.$rundate";

open (LOG,">$log");
LOG->autoflush();

print LOG "# inherit_GO_terms\n";     
print LOG "# version        : $version\n";
print LOG "# run details    : $rundate $runtime\n";
print LOG "\n";


my $out="/wormsrv2/wormbasr/misc_inherit_GO_term.ace";
open (OUT,">$out");
OUT->autoflush();


########################################
# Connect with acedb server            #
########################################

my $db = Ace->connect(-path=>$dbpath,
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

print LOG "inherit_GO_terms run STARTED at $runtime\n\n";

&motif($db) if ($opt_m);
&phenotype if ($opt_p);

$runtime = `date +%H:%M:%S`; chomp $runtime;
print LOG "\ninherit_GO_terms run ENDED at $runtime\n\n";
close LOG;

close OUT;

##############################
# mail $maintainer report    #
##############################

&mail_maintainer("inherit_GO_terms Report:",$maintainers,$log);

##############################
# hasta luego                #
##############################

exit(0);

########################################################################################
####################################   Subroutines   ###################################
########################################################################################

########################################################################################
# motif to sequence mappings                                                           #
########################################################################################

sub motif {
    my $db = shift;
    
    my ($motif,$obj,$term,$protein,$match,$pepobj) = "";
    my (@GO_terms,@pep_homols,@CDS) = "";


    my $i = $db->fetch_many(-query=> 'find Motif "INTERPRO*"');  
    while ($obj = $i->next) {
	$motif = $obj;
	@GO_terms = $obj->GO_term;
	@pep_homols = $obj->Pep_homol;
	
	print "\nMotif : $motif\n";
	foreach $term (@GO_terms) {
	    print "contains GO_term : $term with " . scalar (@pep_homols) . " attached Protein objects\n" if ($debug);
	
	    foreach $protein (@pep_homols) {
		print "maps to Protein: $protein " if ($debug);
		my $pepobj = $db->fetch(Protein=>$protein);
		@CDS = $pepobj->Corresponding_DNA;
		
		foreach $match (@CDS) {
		    print "== $match\n" if ($debug);
		    print OUT "\nSequence : \"$match\"\nGO_term $term\n";
		} # Sequence
	    }     # Protein
	}         # GO_term
    }             # Motif object

}


########################################################################################
# phenotype to sequence mappings                                                       #
########################################################################################

sub phenotype {
}
