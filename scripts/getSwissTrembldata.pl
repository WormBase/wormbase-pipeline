#!/usr/local/bin/perl5.6.1 -w
#
# getSwissTrembldata.pl
#
# Originally crafted by Dan Lawson
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2003-01-02 16:16:59 $
#
# see pod documentation (i.e. 'perldoc getSwissTrembldata.pl') for more information.
#
#####################################################################################################

use lib "/wormsrv2/scripts/";
use Wormbase;
use Getopt::Long;
use strict;


 ###################
 # Paths and stuff #
 ###################

my $tace      = &tace;                         # tace executable path
my $dbdir     = "/wormsrv2/autoace";
my $outdir    = "/nfs/disk100/wormpub/analysis/SWALL";
$ENV{'ACEDB'} = $dbdir;

our %databases = (
	      'SW' => 'SWISSPROT',
	      'TR' => 'TREMBL',
	      'TN' => 'TREMBLNEW'
	      );

 ##############################
 # command-line options       #
 ##############################

my $help;       # Help perdoc
my $debug;      # Debug mode, verbose output to dl1 only

GetOptions (
            "debug"     => \$debug,
            "help"      => \$help
            );

# help page
&usage("Help") if ($help);

 ##############################
 # output files               #
 ##############################
$outdir = "/nfs/disk100/wormpub/analysis/SWALL_DEBUG" if (defined $debug);

open (ACE_GSC,  ">$outdir/output_stl.ace") or die "Failed to open output file: $outdir/output_stl.ace\n";
open (ACE_WTSI, ">$outdir/output_cam.ace") or die "Failed to open output file: $outdir/output_cam.ace\n";
open (OUTPUT,   ">$outdir/output_autoace") or die "Failed to open output file: $outdir/output_autoace\n";

# Grab accessions and sequence versions from autoace
my $command = "Table-maker -p '$dbdir/wquery/accession2clone.def'\nquit\n";

open (TACE, "echo '$command' | $tace | ") || die "Could not open pipe to tace\n";
while (<TACE>) {
    print;
    chomp;
    s/acedb\> //g;      # only need this is using 4_9i code, bug fixed in 4_9k onward (should be redundant)
    next if ($_ eq "");
    next if (/\/\//);
    s/\"//g;
    
    my ($acename,$acc) = (/^(\S+)\s+(\S+)/);
    #    print "\n// Parsing genome sequence $acename [$acc]\n";
    my $GSC = "";
    my $CDS_xref_count = 0;
    my $CDS_found_count = 0;
    
    
    my ($CDS_dbxref,$CDS_gene,$CDS_protein,$CDS_prod,$CDS_name,$pid,$ver,$db,$EMBL_acc);
    my ($CDS_on,$CDS_dbxref_ac,$CDS_dbxref_id,$CDS_dbxref_db);
    my $carryover = 0;
    open (PFETCH, "/usr/local/pubseq/bin/pfetch -F $acc |");
    while (<PFETCH>) {
	chomp;
	if (/^SV\s+(\S+)\.\d+/) {$EMBL_acc = $1; next;}
	if (/^DR/)              {$CDS_xref_count++; next;}
	#	if (/^DR/)              {print "DR lines $_\n"; $CDS_xref_count++; next;}
	if (/^FH\s+Key/)        {next;}
	#	if (/^FH\s+Key/)        {print "\nExpecting $CDS_xref_count CDS in this entry\n\n";next;}
	
	# begin CDS loop
	if (/^FT\s+CDS/) {
	    $CDS_on = 1;
	    ($CDS_dbxref,$CDS_gene,$CDS_protein,$CDS_prod,$CDS_name) = "";
	    ($pid,$ver,$db,$acc) = "";
	    next;
	}
	
	if (/^FT\s+\/gene=\"(\S+)\"/) {
	    $CDS_gene = $1;
	    next;
	}
	
	if (/^FT\s+\/standard_name=\"(\S+)\"/) {
	  $CDS_name = $1;
	}
	
	if (/^FT\s+\/product=\"(\S+.+)/) {
	    $CDS_prod = $1;
	    
	    # Set file on
	    $GSC = 1;
	    
	    # line ends in a speech mark (i.e. full entry)
	    if ($CDS_prod =~ /\"/) {
		chop $CDS_prod;
		next;
	    }
	    else {
		$carryover = 1;
		next;
	    }
	}
	if ($carryover == 1) {
	    (/^FT\s+(\S+.+)/);
	    $CDS_prod .= " $1";
	    if ($CDS_prod =~ /\"/) {
		chop $CDS_prod;
		$carryover = 0;
		next;
	    }
	    else {
		$carryover = 1;
	    }
	}
	
	if (/^FT\s+\/protein_id=\"(\S+)\"/) {
	    $CDS_protein = $1;
	    next;
	}
	
	if (/^FT\s+\/translation/) {
	    if (substr($CDS_prod,-1) eq ")") {chop $CDS_prod};
	    
	    ($CDS_name) = $CDS_prod =~ (/(\S+)$/) if !(defined $CDS_name); # should be filled thru "standard_name" field
	    
	    # check the entry from the protein_id directly
	    
	    #	    print "Checking with protein_id [$CDS_protein]\n";
	    ($CDS_dbxref_ac,$CDS_dbxref_id,$CDS_dbxref_db) = &get_from_protein_id($CDS_protein);
	    
	    ($pid,$ver) = split(/\./,$CDS_protein);
	    
	    if ($GSC) {
		print OUTPUT "EMBL [$acename|$EMBL_acc] CDS [gene:$CDS_name protein_id:$CDS_protein DB_xref:$CDS_dbxref_ac]\n";
		print ACE_GSC "\nSequence : \"$CDS_name\"\nProtein_id $acename $pid $ver\n";
		print ACE_GSC "Database $databases{$CDS_dbxref_db} $CDS_dbxref_id $CDS_dbxref_ac\n";
	    }
	    else {
		print OUTPUT "EMBL [$acename|$EMBL_acc] CDS [gene:$CDS_gene protein_id:$CDS_protein DB_xref:$CDS_dbxref_ac]\n";
		print ACE_WTSI "\nSequence : \"$CDS_gene\"\nProtein_id $acename $pid $ver\n";
		print ACE_WTSI "Database $databases{$CDS_dbxref_db} $CDS_dbxref_id $CDS_dbxref_ac\n";
	    }
	    $CDS_found_count++;
	    next;
	}
	
	if (/^SQ/) {
	    #	    print "\nFound $CDS_found_count CDS in this entry\n\n";
	    
	    next;
	}
	
    }
    close PFETCH;
}
close TACE;
close ACE_GSC;
close ACE_WTSI;
close OUTPUT;

 ##############################
 # a extremidade              #
 ##############################

exit(0);
 
sub get_from_protein_id {
    
    my $protein_id = shift;
    my $acc; 
    my $id; 
    my $db;
    
    open (LOOK, "/usr/local/pubseq/bin/pfetch -F $protein_id |");
    while (<LOOK>) {
	($acc = $1) if (/^AC\s+(\S+);/);
	($id  = $1) if (/^ID\s+(\S+)/);
    }
    close LOOK;
    
    
    # assign database based on length of accession and SWISSPROT species
    if ((length ($acc)) == 8) {
	$db = "TN";
    }
    elsif ($id =~ /_CAEEL/) {
	$db = "SW";
    }
    else {
	$db = "TR";
    }
    
    return ($acc,$id,$db);
}

###############################
# Prints help and disappears  #
###############################

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

=head2 NAME 

getSwissTrembldata.pl

=head2 USAGE

getSwissTrembldata.pl [-options]

getSwissTrembldata.pl runs a Tablemaker query on autoace to return the set of genome 
sequences and their accessions. Each accession is retrieved from EMBL and the CDS 
features (with protein_id) noted. Each protein_id is used to retrieve the SWALL 
flatfile and the ID and AC noted. Finally all of this data is reported in a flatfile
output which is used later in the build.

update_totals mandatory arguments:

=over 4

=item none

=back

update_totals OPTIONAL arguments:

=over 4

=item B<-help>, Help

=back

=head2 Dependencies

no dependencies.

=head1 Called by:
    
=over 4

=item autoace_minder

=back

=head2 AUTHOR

Dan Lawson (B<dl1@sanger.ac.uk>)

=cut
