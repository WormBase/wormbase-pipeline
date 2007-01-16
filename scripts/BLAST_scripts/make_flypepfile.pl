#!/usr/local/ensembl/bin/perl -w

# krb 020829

# converts gadflyX.pep file with fly peptide entries into ace file
# Adds flybase gene name, and annotation IDs as well as separate Database
# lines for both Flybase genes and Gadfly

use strict;
use Getopt::Long;

my $version;
my $verbose;

GetOptions (
        	    "version:s"  => \$version,
	   		 "verbose"    => \$verbose
            );


#example gadfly3.pep header line
#>CG3645-PA translation from_gene[CG3645 CG3645 FBgn0031238] seq_release:3  mol_weight=58068  cds_boundaries:(2L:252,666..254,379[-]) aa_len
#gth:505 transcript_info:[CG3645-RA  seq_release:3] gene_info:[gene symbol:CG3645 FBgn0031238 gene_boundaries:(2L:252,589..271,732[-]) ]

#>CG2671-PA type=protein; loc=2L:complement(11218..11344,11410..11518,11779..12221,12286..12928,13520..13625,13683..14874,14933..15711,17053..17136); ID=l(2)gl-PA; name=l(2)gl-PA; db_xref='CG2671,FlyBase:FBgn0002121'; /gene=l(2)gl; len=1161

#>nop5-PA type=protein; loc=2L:complement(6916901..6918611); name=nop5-PA; 
#dbxref=FlyBase:FBpp0078997,GB_protein:AAF52455.2,FlyBase_Annotation_IDs:CG10206-PA,FlyBase:FBgn0026196;
# MD5=927705a12a71536536d1eea6c514361a; parent=FBtr0079369; release=r4.3; species=Dmel; length=511;

#>FBpp0072157 type=protein; loc=2R:join(19970258..19970592,19970994..19971107,19971175..19971808,19971862..19972342,19972455..19973683); ID=FBpp0072157; name=CG16786-PB; parent=FBgn0034974,FBtr0072248; dbxref=GB_protein:AAM68305.1,FlyBase:FBpp0072157,FlyBase_Annotation_IDs:CG16786-PB,GB_protein:AAM68305; MD5=6356cb499a107083f645b73322b5bd50; release=r5.1; species=Dmel; length=779;


my $blastdir    = "/nfs/acari/wormpipe/BlastDB";
my $acedir      = "/nfs/acari/wormpipe/ace_files";
my $source_file = "$blastdir/gadfly${version}.pep";
my $acefile     = "$acedir/flybase.ace";
# output initally goes to tmp file
my $pepfile  = "$blastdir/gadfly${version}.pep.tmp"; 


open (PEP,">$pepfile") or die "cant open $pepfile $!\n";
open (ACE,">$acefile") or die "cant open $acefile\n";
open (SOURCE,"<$source_file") or die "cant open $source_file\n";

my $record_count = 0;
my $problem_count =0;
my ($gadID, $FBname, $FBgn);
my $count;
while (<SOURCE>){
  chomp;
  if( /^>/ ){
    undef $FBgn;
    undef $gadID;
    undef $FBname;

    s/;//g; # get rid of ;'s
    $count++;
    $record_count++;

    my %fields = /(\w+)=(\S+)/g;

    $FBname = $fields{'name'} if $fields{'name'} ;
    ($FBgn) = $fields{'parent'} =~ /(FBgn\d+)/;
    foreach ( split(/,/,$fields{'dbxref'}) ) {
    	my($key, $value) = split(/:/);
 		if( $key eq 'FlyBase_Annotation_IDs') {
    		$gadID = $value;
    	}
    }
	        
    # some old style names still exist eg pp-CT*****.  In these cases
    # we need to use the 1st field of the "from_gene" fields.

    if( $gadID ){

      #print ace file
      print ACE "\n\nProtein : \"FLYBASE:$gadID\"\n";
      print ACE "Peptide \"FLYBASE:$gadID\"\n";
      print ACE "Species \"Drosophila melanogaster\"\n";

      #FlyBase_gn	      #Gadfly_ID
      print ACE "Gene_name \"$FBname\"\n" if $FBname;
      print ACE "Database \"FlyBase\" FlyBase_gn \"$FBgn\"\n" if ($FBgn);
      print ACE "Database \"FlyBase\" FlyBase_ID \"$gadID\"\n" if $gadID;
      print ACE "Description \"Flybase gene name is $FBname\"\n" if $FBname;
      print ACE "\nPeptide : \"FLYBASE:$gadID\"\n";

      #write database file
      print PEP ">$gadID\n";

    }
    else {
      # report problems?
      print "PROBLEM : $_\n" if ($verbose);
    }
  }
  else{
    if( $gadID ){
      print ACE "$_\n";
      print PEP "$_\n";
    }
  }
}
close SOURCE;
close PEP;
close ACE;

# Now overwrite source file with newly formatted file
system("mv $pepfile $source_file") && die "Couldn't overwrite original peptide file\n";

print "$record_count proteins\n" if ($verbose);

exit (0);


__END__

