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

my $blastdir    = "/nfs/acari/wormpipe/BlastDB";
my $acedir      = "/nfs/acari/wormpipe/ace_files";
my $source_file = "$blastdir/gadfly${version}.pep";
my $acefile     = "$acedir/gadfly.ace";
# output initally goes to tmp file
my $pepfile  = "$blastdir/gadfly${version}.pep.tmp"; 


open (PEP,">$pepfile") or die "cant open $pepfile $!\n";
open (ACE,">$acefile") or die "cant open $acefile\n";
open (SOURCE,"<$source_file") or die "cant open $source_file\n";

my $record_count = 0;
my $problem_count =0;
my ($gadID, $otherGadID, $FBname, $FBgn);
my $count;
while (<SOURCE>){
  chomp;
  if( /^>/ ){
    undef $FBgn;
    undef $gadID;
    undef $otherGadID;
    undef $FBname;
    
    $count++;
    $record_count++;
    
    if( /^>(\S+).*from_gene\[(\S+)\s+(\S+)\s+(\S+)\]/ ) {
      #$1 = Gadfly id  $2 = $3 = Flybase gene name   $4 = Flybase id
      ($gadID, $otherGadID, $FBname, $FBgn)  = ($1, $2, $3, $4);
      
    }
    else{
      # some are just >CG32015-PD translation from_gene[CG32015 CG32015 ]
      if( $_ =~ /from_gene\[(\S+) (\S+)/ ) {
	$gadID = $1;
	$FBname = $1;
      }
    }
    
    if( !($FBgn) ||          # no FBgn
	!( $FBgn =~ m/FBgn/ )# invalid FBgn
	) {
      $FBgn = $1 if ($_ =~ /(FBgn\d+)/);  #last ditch attemp to grap if from the header
    }

    # some old style names still exist eg pp-CT*****.  In these cases
    # we need to use the 1st field of the "from_gene" fields.

    if( $gadID ){
      if ("$gadID" =~ /^pp-/) {
	$gadID = $otherGadID;
      }


      #print ace file
      print ACE "\n\nProtein : \"GADFLY:$gadID\"\n";
      print ACE "Peptide \"GADFLY:$gadID\"\n";
      print ACE "Species \"Drosophila melanogaster\"\n";

      #FlyBase_gn	      #Gadfly_ID
      print ACE "Gene_name \"$FBname\"\n" if $FBname;
      print ACE "Database \"Flybase\" FlyBase_gn \"$FBgn\"\n" if ($FBgn);
      print ACE "Database \"Gadfly\" Gadfly_ID \"$gadID\"\n" if $gadID;
      print ACE "DB_remark \"Flybase gene name is $FBname\"\n" if $FBname;
      print ACE "Peptide : \"GADFLY:$gadID\"\n";

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
system("scp $acedir/gadfly.ace wormsrv2:/wormsrv2/wormbase/ensembl_dumps")  && die "Couldn't copy acefile\n";


exit (0);


__END__

