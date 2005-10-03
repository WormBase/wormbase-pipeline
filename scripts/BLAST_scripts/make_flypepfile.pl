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

    $gadID = $fields{'ID'} if $fields{'ID'} ;
    $FBname = $fields{'name'};

    ($FBgn) = /FlyBase:(FBgn\d+)/;

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

# copy acefile across to /wormsrv2/wormbase/ensembl_dumps
system("scp $acefile wormsrv2:/wormsrv2/wormbase/ensembl_dumps")  && die "Couldn't copy acefile\n";


exit (0);


__END__

