#!/usr/local/bin/perl5.6.1 -w
#
# inherit_GO_terms.pl
#
# map GO_terms to ?Sequence objects from ?Motif and ?Phenotype
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2003-03-14 09:29:03 $      

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Ace;


$|=1;

##############################
# Script variables (run)     #
##############################

my $maintainers = "All";
my $rundate = `date +%y%m%d`; chomp $rundate;
my $runtime = `date +%H:%M:%S`; chomp $runtime;
our ($help, $debug, $motif, $phenotype, $log);


##############################
# command-line options       #
##############################

GetOptions ("help"      => \$help,
            "debug=s"   => \$debug,
	    "phenotype" => \$phenotype,
	    "motif"     => \$motif);


# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

&create_log_files;



##############################
# Paths etc                  #
##############################

my $tace      = &tace;      # tace executable path
my $dbpath    = "/wormsrv2/autoace";                                      # Database path


my $out="/wormsrv2/wormbase/misc/misc_inherit_GO_term.ace";
open (OUT,">$out");
OUT->autoflush();


########################################
# Connect with acedb server            #
########################################

my $db = Ace->connect(-path=>$dbpath,
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

print LOG "inherit_GO_terms run STARTED at $runtime\n\n";

&motif($db) if ($motif);
&phenotype if ($phenotype);

##############################
# read acefiles into autoace #
##############################

my $command =<<END;
pparse /wormsrv2/wormbase/misc/misc_inherit_GO_term.ace
save
quit
END

open (TACE,"| $tace -tsuser inherit_GO_terms $dbpath") || die "Couldn't open tace connection to $dbpath\n";
print TACE $command;
close (TACE);

print LOG "uploaded results into autoace\n\n";

$runtime = `date +%H:%M:%S`; chomp $runtime;
print LOG "\ninherit_GO_terms run ENDED at $runtime\n\n";
close LOG;

close OUT;

##############################
# mail $maintainer report    #
##############################

&mail_maintainer("inherit_GO_terms.pl Report:",$maintainers,$log);

##############################
# hasta luego                #
##############################

$db->close;
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
		    print OUT "\nSequence : \"$match\"\nGO_term $term GO_inference_type IEA\n";
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


######################################################################

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

##########################################

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

=head2   NAME - inherit_GO_terms.pl


=head1 USAGE

=over 4

=item inherit_GO_terms.pl [-options]

=back

inherit_GO_terms.pl assigns GO terms to sequences based on Interpro motifs
and RNAi phenotypes.

inherit_GO_terms.pl mandatory arguments:

=over 4

=item none, (but it won\'t do anything)

=back

inherit_GO_terms.pl OPTIONAL arguments:

=over 4

=item -motif, parse Interpro motif data

=item -phenotype, parse phenotype data

=item -debug, debug

=item -help, help

=back

=cut
