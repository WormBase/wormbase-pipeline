#!/usr/local/bin/perl5.8.0 -w
# 
# correlate_Protein_IDs.pl
# v0.3
# dl
#
# 000619 : dl : PP version
#
# Usage : correlate_Protein_IDs.pl <file1> <file2>
#
# where <file1> is the Protein_ID file mailed by the EBI each monday
# and   <file2> is a tablemaker list from camace 
#

# 000711 : dan : PP version
# 000711 : dan : Use getz rather than efetch (to sequentially query EMBLNEW > EMBL
# 000712 : dan : revamp to do a lot more database calls - makes it's own flatfiles

use lib "/wormsrv2/scripts/";
use Wormbase;

#####################################################################################################
# get flatfile names from the command line                                                          #
#####################################################################################################

my $file1 = shift;
my $CDS = "";                          # CDS name
my $sequence = "";                     # Parent ?Sequence object
my $accession = "";                    # EMBL accession for parent sequence
my $PID_root = "";                     # Protein_ID root    (e.g. CAB02726)
my $PID_ver = "";                      # Protein_ID version (e.g. 1)
my $EMBL_id = "";
my $EMBL_ac = "";
my $tace = &tace; # tace executable path
my $db = glob("/wormsrv2/camace");                            # Database path
my $exec="$tace $db";     

my $debug = 1;

#####################################################################################################
# Main Loop                                                                                         # 
#####################################################################################################

   ##################################################################################################
   # open output file (.ace for uploading into camace                                               #
   ##################################################################################################

open (OUT, ">update_protein_ID.ace");

   ##################################################################################################
   # make a flatfile relating ACeDB_id, NDBL_ID & NDB_AC                                           #
   ##################################################################################################

print "make accession list from camace\n" if ($debug); 
&make_accession_list;

   ##################################################################################################
   # make a flatfile for Protein_id's in ACeDB                                                      #
   ##################################################################################################

print "collate Protein_ID's from camace\n" if ($debug); 
&protein_ID_list_from_camace;     

   ##################################################################################################
   # main loop : cycle through camace tablemaker output for predicted genes                         #
   ##################################################################################################

open (FILE2, "</tmp/protein_ID_list.camace");
while (<FILE2>) {
    select((select(STDOUT),$|=1)[0]);
    if (/^(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)/) {
	($acc_camace,$protein_ID_camace,$protein_version_camace,$gene) = ($1,$2,$3,$4);
	$PID_camace = $protein_ID_camace . "." . $protein_version_camace;
	print "camace  : $gene   \t'$PID_camace'\t[$acc_camace]\n";

	&get_EMBL;
	print "ID_list : $gene   \t'$PID_EMBL'\t[$acc_camace]\n";

	&get_protein_id;
	print "getz    : $gene   \t'$PID_EMBL2'\t[$acc_camace]\n";

	if ($PID_EMBL eq $PID_camace) {
	    print "protein_IDs match -> no change\n\n";
	}
	else {
	    if ($PID_EMBL ne "") {
		print "disparity in the protein_IDs -> update $gene\n\n";
		print OUT "\n// inherit protein_id data from EMBL database\n";
		print OUT "# Sequence : \"$gene\"\n";
		print OUT "# Protein_id $clone_camace $protein_ID_EMBL $protein_version_EMBL\n\n";
	    }
	    else {

		if ($PID_EMBL2 eq $PID_camace) {
		    print "protein_IDs match (old EMBL flatfile) -> no change\n\n";
		}
		else {
		    print "disparity in the protein_IDs -> update $gene\n\n";
		    print OUT "\n// retain protein_id data from ACeDB database\n";
		    print OUT "# Sequence : \"$gene\"\n";
		    print OUT "# Protein_id $clone_camace $protein_ID_camace $protein_version_camace\n\n";
		}
	    }
	}
	
	($protein_ID_EMBL,$protein_version_EMBL,$protein_ID_camace,$protein_version_camace) = "";
	($PID_camace,$PID_EMBL,$acc_EMBL,$clone_camace) = "";


    }
}
close (FILE2);
close (OUT);

exit(0);

#####################################################################################################
############################################ SUBROUTINES ############################################
#####################################################################################################


#####################################################################################################
# Make Protein_id list from camace                                                                  # 
#####################################################################################################
#
# not local
#

sub protein_ID_list_from_camace {

    my $command=<<EOF;
find elegans_CDS
show -a Protein_id
quit
EOF

    open (OUT, ">/tmp/protein_ID_list.camace");
    open(TACE, "echo '$command' | $exec -| ");
    while (<TACE>) {
	chomp;
	next if ($_ eq "");
	next if (/acedb/);
	next if (/\/\//);
	s/\"//g;
	
	if (/^Sequence : (\S+)/) {
	    $CDS = $1;
	    next;
	}
	if (/Protein_id\s+(\S+)\s+(\S+)\s+(\d+)/) {
	    ($sequence,$PID_root,$PID_ver) = ($1,$2,$3);
	    &get_accession;
	    select((select(STDOUT),$|=1)[0]);
	#    print "-> parse Protein_ID's for gene '$CDS' in entry '$sequence'\n" if ($debug);
	    print OUT "$EMBL_ac 00 $PID_root $PID_ver 0000000000 $CDS\n";
	}
    }    
    close (TACE);     
    close (OUT);
    
} 




sub get_EMBL {

    open (FILE1, "/usr/local/bin/agrep -w '$gene' $file1 |");
    while (<FILE1>) {
	chomp;
	if (/^$acc_camace\s+\d+\s+(\S+)\s+(\d+)/) {
	($protein_ID_EMBL,$protein_version_EMBL) = ($1,$2);
	    $PID_EMBL = $protein_ID_EMBL . "." . $protein_version_EMBL;
	    last;
	}
    }
    close (FILE1);
    return ($PID_EMBL,$protein_ID_EMBL,$protein_version_EMBL);
}    


sub get_protein_id {
    
    $next_protein_ID = "";
    open (GETZ, "getz -e  \"([{embl emblnew}-acc:$acc_camace]\!EMBL<EMBLNEW)\" |");

    while (<GETZ>) {
	if (/\/gene=\"$gene\"/) {
	    $next_protein_ID = 1;
	}
	if ((/\/protein_id=\"(\S+)\"/) && ($next_protein_ID == 1)) {
	    $PID_EMBL2 = $1;
	    ($protein_ID_EMBL2,$protein_version_EMBL2) = split (/\./, $PID_EMBL2);
	    last;
	}
    }
    close (GETZ);
   
    return ($PID_EMBL2,$protein_ID_EMBL2,$protein_version_EMBL2);
}


sub get_accession {

    $EMBL_ac = "";
    open (AC, "/usr/local/bin/agrep -w '$sequence' /tmp/camace_accessions |");
    while (<AC>) {
	($EMBL_ac) = (/\S+\s+\S+\s+(\S+)/);
	last;
    }    
    close (AC);

    if ($EMBL_ac eq "") {
	print "\n$sequence not found in EMBL\n";
    }
    return ($EMBL_ac);
}


sub get_accession2 {

    $EMBL_ac = "";

my $command2=<<EOF;
find Genome_Sequence $sequence
show -a Database
quit
EOF

    open(TACE2, "echo '$command2' | $exec -| ");
    while (<TACE2>) {
	chomp;
	next if ($_ eq "");
	next if (/acedb/);
	next if (/\/\//);
	s/\"//g;
	
	#print "$_";
	if (/Database\s+EMBL\s+(\S+)\s+(\S+)/) {
	    ($EMBL_id,$EMBL_ac) = ($1,$2);
	}
    }
    close (TACE2);     
    
    return ($EMBL_ac);
}


sub make_accession_list {

    my $command3=<<EOF;
find Genome_Sequence
show -a Database
quit
EOF

    open (OUT2, ">/tmp/camace_accessions") || die "failed to open file : '/tmp/camace_accessions'\n\n";
    open(TACE3, "echo '$command3' | $exec -| ");
    while (<TACE3>) {
	chomp;
	next if ($_ eq "");
	next if (/acedb/);
	next if (/\/\//);
	s/\"//g;
	
	if (/Sequence\s+\S+\s+(\S+)/) {
	    $sequence = $1;
	}
	if (/Database\s+EMBL\s+(\S+)\s+(\S+)/) {
	    select((select(STDOUT),$|=1)[0]);
	  #  print "-> Dump accession for $sequence\n" if ($debug);
	    print OUT2 "$sequence $1 $2\n";
	}
    }
    close (TACE3);     
    close (OUT2);
}


__END__

=pod

=head1 NAME - correlate_Protein_IDs.pl

=head2 AUTHOR

Dan Lawson dl1@sanger.ac.uk

=head2 USAGE

correlate_Protein_IDs.pl <file>

correlate_Protein_IDs.pl will query an ACeDB database
to produce two flatfiles:

/tmp/camace_accessions        : NDB_ID & NDBL_AC for all Genome_Sequences
/tmp/protein_ID_list.camace   : Protein_IDs for all C. elegans CDSs

The list of Protein_IDs is correlated with the input flatfile (weekly dump
from EBI), and the latest version of the flatfiles via getz. 

The script will report objects with mismatched Protein_IDs between ACeDB
and either of the EMBL views. It will attempt to create an .ace patch file
to update the database (manual uploading)

=head2 MANDATORY ARGUMENTS

=over 2

=item *
none

=back
=cut

=head2 REVISIONS

=over 2

=item *
v0.1
2000-07-11 : dan : PP version
2000-07-11 : dan : Use getz rather than efetch (to sequentially query EMBLNEW > EMBL

v0.3
2000-07-12 : dan : revamp to do a lot more database calls - makes it's own flatfiles

=back
=cut

