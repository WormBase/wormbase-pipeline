#!/usr/local/bin/perl5.6.1 -w

# krb 020829

# converts gadflyX.pep file with fly peptide entries into ace file
# Adds flybase gene name, and annotation IDs as well as separate Database
# lines for both Flybase genes and Gadfly

use strict;
#example gadfly3.pep header line
#>CG3645-PA translation from_gene[CG3645 CG3645 FBgn0031238] seq_release:3  mol_weight=58068  cds_boundaries:(2L:252,666..254,379[-]) aa_len
#gth:505 transcript_info:[CG3645-RA  seq_release:3] gene_info:[gene symbol:CG3645 FBgn0031238 gene_boundaries:(2L:252,589..271,732[-]) ]

my $dir = glob("~wormpipe/BlastDB");
my $acefile = "$dir/fly";
my $DBfile = "$dir/test_gadfly.pep";

my $pepfile = shift;

open (DB,">$DBfile") or die "cant open $DBfile\n";
open (ACE,">$acefile") or die "cant open $acefile\n";
open (PEP,"<$pepfile") or die "cant open $pepfile\n";

my $record_count = 0;
my $problem_count =0;
my ($gadID, $FBname, $FBgn);
my $count;
while (<PEP>)
  {
    chomp;
    if( /^>/ )
      {
	$record_count++;
	undef $FBgn;
	if( /^>(\S+).*from_gene\[\S+\s+(\S+)\s+(\S+)\]/ ) 
	  {
	    $count++;
	    #$1 = Gadfly id  $2 = Flybase gene name   $3 = Flybase id
	    ($gadID, $FBname, $FBgn)  = ($1, $2, $3);
	    unless( $FBgn =~ m/FBgn/ ) {
	      $FBgn = "non-assigned";
	    }
	    
	    #print ace file
	    print ACE "\n\nProtein : \"GADFLY:$gadID\"\n";
	    print ACE "Peptide \"GADFLY:$gadID\"\n";
	    print ACE "Species \"Drosophila melanogaster\"\n";
	    print ACE "DB_remark \"Flybase gene name is $FBname\"\n";
	    print ACE "Database \"Flybase\" \"$FBname\" \"$FBgn\"\n";
	    print ACE "Database \"Gadfly\" \"$gadID\" \"$gadID\"\n\n";
	
	    print ACE "Peptide : \"GADFLY:$gadID\"\n";

	    #write database file
	    print DB ">GADFLY:$gadID\n";
#	    if( $count > 20 ) { last;  }
	  }
#	else{
#	  print "PROBLEM : $_\n";
#	}
      }
    else
      {
	unless ( $FBgn ){$problem_count++;next;}
	if( "$FBgn" ne "non-assigned" )
	  {
	    print ACE "$_\n";
	    print DB "$_\n";
	  }
      }
  }
close PEP;
close DB;
close ACE;

print "\n\nabout to copy (scp) $acefile to /wormsrv2/wormbase/ensembl_dumps/\n";
system ("scp -r $acefile wormpub\@wormsrv2:/wormsrv2/wormbase/ensembl_dumps/") and warn "copy $acefile failed\n";

print "\n\nMaking gadfly BLASTable database \n";
system ("setdb $DBfile") and warn "cant setdb on $DBfile\n";

print "$record_count proteins - $problem_count problems\n";
exit (0);

__END__

