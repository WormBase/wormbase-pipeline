#!/usr/local/bin/perl5.8.0 -w
#
# inherit_GO_terms.pl
#
# map GO_terms to ?Sequence objects from ?Motif and ?Phenotype
#
# Last updated by: $Author: ar2 $     
# Last updated on: $Date: 2006-02-23 17:29:05 $      

use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Ace;

##############################
# Script variables (run)     #
##############################

my ($help, $debug, $motif, $phenotype,$store);
my $verbose;             # for toggling extra output
my $maintainers = "All"; # who receives emails from script
my $noload;              # generate results but do not load to autoace

##############################
# command-line options       #
##############################

GetOptions ("help"      => \$help,
            "debug=s"   => \$debug,
	    "phenotype" => \$phenotype,
	    "motif"     => \$motif,
	    "noload"    => \$noload,
    	    "store:s"   => \$store
    	);

# Display help if required
&usage("Help") if ($help);

# recreate configuration 
my $wb;
if ($store) { $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wb = Wormbase->new( -debug => $debug, -test => $debug, ) }

# Variables Part II (depending on $wb) 
$debug = $wb->debug if $wb->debug;    # Debug mode, output only goes to one user

my $log=Log_files->make_build_log($wb);

##############################
# Paths etc                  #
##############################

my $tace      = $wb->tace;      # tace executable path
my $dbpath    = $wb->autoace;                                      # Database path

my $out=$wb->acefiles."/inherited_GO_terms.ace";
open (OUT,">$out");
OUT->autoflush();

########################################
# Connect with acedb server            #
########################################

my $db = Ace->connect(-path=>$dbpath,
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};



&motif($db)     if ($motif);
&phenotype($db) if ($phenotype);


##############################
# read acefiles into autoace #
##############################

unless ($noload || $debug) {

  my $command = "pparse $out\nsave\nquit\n";
    
  open (TACE,"| $tace -tsuser inherit_GO_terms $dbpath") || die "Couldn't open tace connection to $dbpath\n";
  print TACE $command;
  close (TACE);  

  $log->write_to("uploaded results into autoace\n\n");
    
}

##############################
# mail $maintainer report    #
##############################
$log->mail();

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
		@CDS = $pepobj->Corresponding_CDS;
		
		foreach $match (@CDS) {
		    print "== $match\n" if ($debug);
		    print OUT "\nCDS : \"$match\"\nGO_term $term IEA Inferred_automatically\n";
		} # CDS
	    }     # Protein
	}         # GO_term
    }             # Motif object
}

########################################################################################
# phenotype to sequence mappings                                                       #
########################################################################################

sub phenotype {
  my $db = shift;
  
  my ($obj,$term,$rnai,$match,$rnaiobj) = "";
  my (@GO_terms,@rnai,@CDS) = "";
  
  my $i = $db->fetch_many(-query=> 'find Phenotype "*"');  
  while ($obj = $i->next) {
    $motif = $obj;
    @GO_terms = $obj->GO_term;
    @rnai = $obj->RNAi;
    
    print "RNA : $motif\n" if ($verbose);
    foreach $term (@GO_terms) {
      print "contains GO_term : $term with " . scalar (@rnai) . " attached RNAi objects\n" if ($debug);
      
      foreach $rnai (@rnai) {
	print "maps to RNAi $rnai " if ($debug);
	my $rnaiobj = $db->fetch(RNAi=>$rnai);
	@CDS = $rnaiobj->Predicted_gene;
	
	foreach $match (@CDS) {
	  print "== $match\n" if ($debug);
	  print OUT "\nCDS : \"$match\"\nGO_term $term IMP Inferred_automatically\n";
	} # CDS
      }   # RNAi
    }     # GO_term
  }       # Phenotype object
}


######################################################################


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
and RNAi phenotypes. The resulting acefile will be loaded into the database
as part of the script run (default database is autoace).

inherit_GO_terms.pl mandatory arguments:

=over 4

=item none, (but it won\'t do anything)

=back

inherit_GO_terms.pl OPTIONAL arguments:

=over 4

=item -motif, parse Interpro motif data

=item -phenotype, parse phenotype data

=item -noload, do not upload results to autoace

=item -debug, debug (results not loaded into autoace)

=item -help, help

=item -verbose, toggle extra output to screen

=item -store <storable_file>, specifiy stored commandline options

=back

=cut
