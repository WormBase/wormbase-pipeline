#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: ar2 $
# Last updated on: $Date: 2005-04-14 12:59:53 $ 

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use GENEACE::Geneace;
use Getopt::Long;

# ----- command line options
my ($input, $debug);
GetOptions ("input=s"  => \$input);

# ----- warn
if (!$input){
  print "\n\nUSAGE: -input path/filename\n\n";
  exit(0);
}

# ----- check user priviledges
my $user = `whoami`; chomp $user;
if ($user ne "wormpub"){
  print "\n\nYou need to be wormpub to do this!\n\n";
  exit(0);
}

# ----- global variables
my $ga = init Geneace();
my $database = $ga->curr_db;
my $tace = &tace;
my $allele_dir = "/wormsrv1/geneace/ALLELE_DATA/JAPANESE_KNOCKOUTS";


##########################################
# open logfile and various output streams
##########################################

my $rundate = &rundate;
my $log = "/wormsrv2/logs/parse_NBP_alleles.$rundate.$$";

open(LOG,    ">$log")                                              || die $!;
open(ACE,    ">$allele_dir/NBP_alleles.$rundate.ace")              || die $!;
open(DELETE, ">$allele_dir/NBP_alleles.$rundate.delete.ace")       || die $!;
open(CHECK,  ">$allele_dir/NBP_alleles_to_check_by_hand.$rundate") || die $!;


# ----- prepare a locus to gene_id conversion hash

my %Gene_info = $ga->gene_info();

# ----- prepare clone coords from GFF file
my %clone_info = $ga->get_clone_chrom_coords(); # key is clone, values: chrom, start_coord, end_coord

# ----- parse NBP flatfile
my (%NBP_info, %chrom_NBP_allele);

# only work with certain fields in their input file
my @NBP = `cut -f 1-6,9,13 $input`;


print LOG "\nProcessing \'$input\' input file\n";

foreach(@NBP){
  chomp;
  
  # skip header line
  next if m/^allele/;

  my($allele, $locus, $gene, $clone, $indel, $pheno, $primers, $mapper)=split(/\t/,$_);

  if ($allele =~ /^tm0.+/){
    $allele =~ /(tm)0(\d+)/;  # remove leading 0 in allele name
    $allele = $1.$2;	
  }

  # remove any whitespaces surrounding gene ID
  $gene =~ s/\s+//g;

  $pheno    = "NA" if !$pheno; $pheno =~ s/\"//g; # get rid of " from strange flatfile exported from filemaker
  $mapper   = "Mitani S" if $mapper;
  $mapper   = "NA" if !$mapper;

  # initialization
  my $L_clone=(); ; my $R_clone=(); my $site_L =(); my $site_R =(); my $remark_1 =(); my $insert=(); my @indel_info=();

  #----- process coordinates info in $indel

  $indel =~ s/\"//g; # get rid of hidden " from the strange flatfile

  if ($indel =~ /(\[(.+)\]|([A-Z0-9]+):)*\s*\d+\s*\/\s*(\d+)\s*-\s*(\[(.+)\]|([A-Z0-9]+):)*\s*(\d+)\s*\/\s*\d+\s*\((.+)\)/){
    if($9){
      $L_clone  = $2 if $2;
      $L_clone  = $3 if $3;
      $site_L   = $4;
      $R_clone  = $6 if $6;
      $R_clone  = $7 if $7;
      $site_R   = $8;
      $remark_1 = $9;
    }
    else{
      $site_L   = $1;
      $site_R   = $2;
      $remark_1 = $3;
    }
  }
  elsif ($indel =~ /(\[(.+)\]|([A-Z0-9]+):)*\s*\d+\s*\/\s*(\d+)\s*[-|+]\s*([\d\+atcgnACGTN\?]+)\s*[-|+]\s*(\[(.+)\]|([A-Z0-9]+):)*\s*(\d+)\s*\/\s*\d+\s*\((.+)\)/){
    if ($10){
      $L_clone  = $2 if $2;
      $L_clone  = $3 if $3;
      $site_L   = $4;
      $insert   = $5;
      $R_clone  = $7 if $7;
      $R_clone  = $8 if $8;
      $site_R   = $9;
      $remark_1 = $10;
    }
    else{
      $site_L   = $1;
      $insert   = $2;
      $site_R   = $3;
      $remark_1 = $4;
    }
  }
  else {    
    print CHECK "$allele -> $indel  \n" if $indel ne "deletion_site";  # this needs hand check as the format is quite versatile
  }

  $insert  = "NA" if !$insert;
  $L_clone = "NA" if !$L_clone;
  $R_clone = "NA" if !$R_clone;
  $site_L  = "NA" if !$site_L;
  $site_R  = "NA" if !$site_R;

  # convert CDS name to clone name by removing trailing digit
  $clone =~ s/\..+//;


  #----- process primers info in $primer (4 primers)
  my ($ext_f, $int_b, $ext_b, $int_f) = split(/,/, $primers);
  my @primers;

  if(!defined($ext_f) || !defined($int_b) || !defined($ext_b) || !defined($int_f)){
    print "ERROR: primers are missing for $allele\n";
    @primers = "NA"
  }
  else{
    # ----- get rid of hidden/strange space and double quotes
    $ext_f =~ s/\s|\"//g;
    $ext_b =~ s/\s|\"//g;
    $int_f =~ s/\s|\"//g;
    $int_b =~ s/\s|\"//g;
    

    if ($ext_f && $int_f && $int_b && $ext_b){
      push(@primers, $ext_f, $int_f, $int_b, $ext_b);
    }
  }

  push(@indel_info, $site_L, $site_R, $remark_1, $insert, $L_clone, $R_clone, $indel);
  push(@{$NBP_info{$allele}}, $locus, $gene, $clone, \@indel_info, $pheno, \@primers, $mapper);
  push(@{$chrom_NBP_allele{$clone_info{$clone}->[0]}}, $allele); # key: chromosome letter

}


print LOG "Finding flanking sequences for alleles\n";

&get_30_bp_flanks($database);


# ----- mail notice

my $recipients = "mt3\@sanger.ac.uk";
&mail_maintainer("NBP allele update", $recipients, $log);


# ----- The End

print LOG "\nScript finished\n";
close(LOG);
close(CHECK);
close(ACE);
close(DELETE);
exit(0);

#-------------------------------------
#       s u b r o u t i n e s
#-------------------------------------

sub get_30_bp_flanks {

  my $database = shift;

  my @LGs = qw (I II III IV V X);

  foreach my $chrom (@LGs){

    # -- prepare DNA files
    my $dna_file = "$database/CHROMOSOMES/CHROMOSOME_".$chrom.".dna";
    my @line = `egrep "[atcg]" $dna_file`;
    my $line;

    foreach (@line){chomp; $line .= $_}

    my $superlink;

    foreach my $allele (@{$chrom_NBP_allele{$chrom}}){

      print ACE "\nVariation : \"$allele\"\n";
      print ACE "Allele\n";
      print DELETE "\nVariation : \"$allele\"\n";
      print DELETE "-D Allele\n";
      
      my @indel_info = @{$NBP_info{$allele}->[3]};
      my $L_clone = $indel_info[4];
      my $R_clone = $indel_info[5];

      # depending on information present, determine left and right flanking sequences
      # initially set to NA, and if they are still 'NA' when you write the ace files
      # then you can't print any flanking sequences!
      my $DNA_L = "NA";
      my $DNA_R = "NA";

#      print "indel: $indel_info[0] $indel_info[1] $indel_info[2] $indel_info[3] $L_clone $R_clone\n";

      if($indel_info[0] eq "NA" || $indel_info[1] eq "NA"){
	# this is for when column 5 of input file cannot be parsed by regular expression
	# these are probably complex deletion + insertion alleles
	print "CHECK: $allele - flanking sequences need to be determined by hand\n";
      }
      else{

	if($L_clone eq "NA"){
	  $DNA_L =substr($line, $clone_info{$NBP_info{$allele}->[2]}->[1] + $indel_info[0]-1 -31, 30);
	}
	else{
	  $DNA_L =substr($line, $clone_info{$L_clone}->[1] + $indel_info[0]-1 -31, 30);
	}

	if($R_clone eq "NA"){
	  $DNA_R =substr($line, $clone_info{$NBP_info{$allele}->[2]}->[1] + $indel_info[1]-1 -31, 30);
	}
	else{
	  $DNA_R =substr($line, $clone_info{$L_clone}->[1] + $indel_info[1]-1 -31, 30);
	}
      }

      print ACE "Sequence \"$NBP_info{$allele}->[2]\"\n";
      print DELETE "-D Sequence \"$NBP_info{$allele}->[2]\"\n";

      # use Gene IDs from their input file
      if ( $NBP_info{$allele}->[1] =~ /^WBGene\d{8}$/ ){
	my $gene = $NBP_info{$allele}->[1]; 
	print ACE "Gene \"$gene\"  \/\/$NBP_info{$allele}->[0]\n";
	print DELETE "-D Gene \"$gene\"  \/\/$NBP_info{$allele}->[0]\n";
      }
      else{
	print "ERROR: $NBP_info{$allele}->[1] is a non standard Gene identifier\n";
      }

      # ----- remove tm and use the rest in as allele id to link back to NBP allele webpage
      my $allele_id = $allele;
      $allele_id =~ s/tm//;

      # ----- acefile
      print ACE "PCR_product     \"$allele"."_external\"\n";
      print ACE "PCR_product     \"$allele"."_internal\"\n";
      print ACE "Author \"$NBP_info{$allele}->[6]\"\n" if $NBP_info{$allele}->[6] ne "NA";
      if($DNA_L eq "NA" || $DNA_R eq "NA"){
#	print "ERROR: Can't determine flanking sequences for $allele. Please add by hand\n";
      }
      else{
	print ACE "Flanking_sequences \"$DNA_L\" \"$DNA_R\"\n";
      }
      print ACE "Deletion\n" if $indel_info[3] eq "NA";
      print DELETE "-D PCR_product     \"$allele"."_external\"\n";
      print DELETE "-D PCR_product     \"$allele"."_internal\"\n";
      print DELETE "-D Author \"$NBP_info{$allele}->[6]\"\n" if $NBP_info{$allele}->[6] ne "NA";
 
      if($DNA_L eq "NA" || $DNA_R eq "NA"){
	print DELETE "-D Flanking_sequences\n";
      }
      else{
	print DELETE "-D Flanking_sequences \"$DNA_L\" \"$DNA_R\"\n";
      }

      print DELETE "-D Deletion\n" if $indel_info[3] eq "NA";

      if ($indel_info[3] ne "NA"){
	print ACE "Insertion \"", lc($indel_info[3]),"\"\n";
	print ACE "Deletion\n";
	print ACE "Method \"Deletion_and_insertion_allele\"\n";
	print DELETE "-D Insertion \"", lc($indel_info[3]),"\"\n";
	print DELETE "-D Deletion\n";
	print DELETE "-D Method \"Deletion_and_insertion_allele\"\n";
      }
      else {
	print ACE "Method \"Deletion_allele\"\n";
	print DELETE "-D Method \"Deletion_allele\"\n";
      }
      print ACE "NBP_allele\n";
      print ACE "Species \"Caenorhabditis elegans\"\n";
      print ACE "Phenotype \"$NBP_info{$allele}->[4]\"\n" unless ($NBP_info{$allele}->[4] eq "NA");
      print ACE "Mutagen \"TMP\/UV\"\n";
      print ACE "Laboratory \"FX\"\n";
      print ACE "Map \"$chrom\"\n"; 
      print ACE "Remark \"Mutations at cosmid coordinates: $indel_info[6]\"\n";
      print ACE "Remark \"<A href='http:\\/\/www.shigen.nig.ac.jp\/c.elegans\/mutants/DetailsSearch?lang=english&seq=$allele_id' target=_new> more on $allele...<\/A>\"\n";
      print DELETE "-D NBP_allele\n";
      print DELETE "-D Species \"Caenorhabditis elegans\"\n";
      print DELETE "-D Phenotype \"$NBP_info{$allele}->[4]\"\n" unless ($NBP_info{$allele}->[4] eq "NA");
      print DELETE "-D Mutagen \"TMP\/UV\"\n";
      print DELETE "-D Laboratory \"FX\"\n";
      print DELETE "-D Map \"$chrom\"\n"; 
      print DELETE "-D Remark \"Mutations at cosmid coordinates: $indel_info[6]\"\n";
      print DELETE "-D Remark \"<A href='http:\\/\/www.shigen.nig.ac.jp\/c.elegans\/mutants/DetailsSearch?lang=english&seq=$allele_id' target=_new> more on $allele...<\/A>\"\n";


      # ----- only for those tm allele having primer information
      if (scalar @{$NBP_info{$allele}->[5]} != 1){

	print ACE "\nPCR_product : $allele"."_external\n";
	print ACE "Oligo $allele"."_external_f\n";
	print ACE "Oligo $allele"."_external_b\n";
	print ACE "Variation \"$allele\"\n";
	print ACE "Method Variation_PCR\n";
	
	print ACE "\nOligo : \"$allele"."_external_f\"\n";
	print ACE "Sequence \"$NBP_info{$allele}->[5]->[0]\"\n";
	print ACE "Length ", length($NBP_info{$allele}->[5]->[0]), "\n";
        print ACE "PCR_product \"$allele"."_external\"\n";

      	print ACE "\nOligo : \"$allele"."_external_b\"\n";
	print ACE "Sequence \"$NBP_info{$allele}->[5]->[3]\"\n";
	print ACE "Length ", length($NBP_info{$allele}->[5]->[3]), "\n";
        print ACE "PCR_product \"$allele"."_external\"\n";

        print ACE "\nPCR_product : $allele"."_internal\n";
        print ACE "Oligo $allele"."_internal_f\n";
        print ACE "Oligo $allele"."_internal_b\n";
        print ACE "Variation \"$allele\"\n";

        print ACE "\nOligo : \"$allele"."_internal_f\"\n";
	print ACE "Sequence \"$NBP_info{$allele}->[5]->[1]\"\n";
	print ACE "Length ", length($NBP_info{$allele}->[5]->[1]), "\n";
        print ACE "PCR_product \"$allele"."_internal\"\n";

      	print ACE "\nOligo : \"$allele"."_internal_b\"\n";
	print ACE "Sequence \"$NBP_info{$allele}->[5]->[2]\"\n";
	print ACE "Length ", length($NBP_info{$allele}->[5]->[2]), "\n";
        print ACE "PCR_product \"$allele"."_internal\"\n";

	print DELETE "\nPCR_product : $allele"."_external\n";
	print DELETE "-D Oligo $allele"."_external_f\n";
	print DELETE "-D Oligo $allele"."_external_b\n";
	print DELETE "-D Variation \"$allele\"\n";
	
	print DELETE "\nOligo : \"$allele"."_external_f\"\n";
	print DELETE "-D Sequence \"$NBP_info{$allele}->[5]->[0]\"\n";
	print DELETE "-D Length ", length($NBP_info{$allele}->[5]->[0]), "\n";
        print DELETE "-D PCR_product \"$allele"."_external\"\n";

      	print DELETE "\nOligo : \"$allele"."_external_b\"\n";
	print DELETE "-D Sequence \"$NBP_info{$allele}->[5]->[3]\"\n";
	print DELETE "-D Length ", length($NBP_info{$allele}->[5]->[3]), "\n";
        print DELETE "-D PCR_product \"$allele"."_external\"\n";

        print DELETE "\nPCR_product : $allele"."_internal\n";
        print DELETE "-D Oligo $allele"."_internal_f\n";
        print DELETE "-D Oligo $allele"."_internal_b\n";
        print DELETE "-D Variation \"$allele\"\n";

        print DELETE "\nOligo : \"$allele"."_internal_f\"\n";
	print DELETE "-D Sequence \"$NBP_info{$allele}->[5]->[1]\"\n";
	print DELETE "-D Length ", length($NBP_info{$allele}->[5]->[1]), "\n";
        print DELETE "-D PCR_product \"$allele"."_internal\"\n";

      	print DELETE "\nOligo : \"$allele"."_internal_b\"\n";
	print DELETE "-D Sequence \"$NBP_info{$allele}->[5]->[2]\"\n";
	print DELETE "-D Length ", length($NBP_info{$allele}->[5]->[2]), "\n";
        print DELETE "-D PCR_product \"$allele"."_internal\"\n";

      }
    }
  }
}

__END__

=head2 NAME - parse_NBP_alleles.pl
                                                                                                                                                 
B<USAGE> 
    
    parse_NBP_alleles.pl -input <input_file>


B<Interval of update>

    Normally every 2-3 wks. Shohei Mitani sends out a full list of all tm* alleles in a tab 
    separated value (tsv) file.  This file will contain those alleles previously processed by 
    us (but some of their details may have changed) plus any new alleles and associated PCR 
    product information.

B<What the script does>	

    The deletion alleles in the input file are described using clone coordinates and so they 
    do not supply the 30 bp flanking sequences that we need to be able to curate the alleles 
    in our WormBase allele pipeline.  Hence, the script uses these coordinates to retrieve 
    30 bp flanking sequences (as well as processing other basic information from the file).
    The script does the following two steps to get flanking sequence info:

    (1) grep clone info (chrom, start, end, clone_name) from GFF file: 
        /wormsrv2/autoace/GFF_SPLITS/WSXXX/CHROMOSOME_*.clone_acc.gff`;
        (via get_clone_chrom_coords() of Geneace.pm).

    (2) fetch flanks from ~wormpub/DATABASES/current_DB/CHROMOSOMES/CHROMOSOME_*.dna file 
        using substr command


B<What you need to do>

     First you MUST remove the Mac carriage return in their input file and convert this to a
     UNIX style newline character:

     1) tr '\r' '\n' < input_file > NBP_alleles.yymmdd.txt

     This processed file should be saved in /wormsrv1/geneace/ALLELE_DATA/JAPANESE_KNOCKOUTS 
     (need to be wormpub for write access).  Then run the script, specifying the input file 
     that you have just generated:

     2) parse_NBP_alleles.pl -input NBP_alleles.yymmdd.txt
	
     The script will create three output files:

	i)   NBP_alleles.yymmdd.ace
	ii)  NBP_alleles.yymmdd.delete.ace
        iii) NBP_alleles_to_check_by_hand.yymmdd

     The first two files complement each other.  The first contains a suitable ace file of all 
     allele information from their input file (including flanking sequences).  The second is a 
     delete ace file which will remove only that information which is contained in the first file.

     This means that if there is a delete ace file from the *last* time this script was run then
     we can use this to first remove any tm* allele data that will be in geneace before updating 
     with the new information (plus any old information which has been changed).  So each time we 
     run this script we will make these two ace files, but we will always need to first load the 
     delete ace file from the last time the script was run.  E.g.

     If today is 041230 and we just run this script, we will make:

     NBP_alleles.041230.ace
     NBP_alleles.041230.delete.ace

     If we then run the script again on 050122 then we will make:
   
     NBP_alleles.050122.ace
     NBP_alleles.050122.delete.ace

     So we would first load NBP_alleles.041230.delete.ace to geneace to clear any old information
     and then load the new NBP_alleles.050122.ace file to add all of the new information.

     You should clear up the directory from time to time to remove very old files.


B<Fix by hand checks>

     The third file the script makes is called NBP_alleles_to_check_by_hand.yymmdd.txt.  These are
     things to be investigated by hand.  Typically, the fix-by-hand cases are due to 
       (1) the free-styled strings embedded in between two lesion coordinates of an allele.
           E.g. the line "tm325 -> 14441/14442- 33 bp addition-14971/14972  (530 bp deletion)" 
           indicates the script has no knowledge about the "33 bp addition" and what the 33 bp 
           additions are: so the script only knows it is a deletion allele.   Basically, you 
           need to look at tm325 allele in Geneace and see what info is missing from the output 
           line: i.e. you should change Deletion tag to Deletion_with_insertion tag and change the 
           Method to Deletion_and_insertion.


B<Contact person>

	 Shohei Mitani (mitani1@research.twmu.ac.jp)


