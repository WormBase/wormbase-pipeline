#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-04-30 11:09:01 $ 

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use GENEACE::Geneace;
use Getopt::Long;

# ----- command line options
my ($input, $debug);
GetOptions ("i|input=s"  => \$input,
	    "d|debug"    => \$debug,
           );

# ----- warn
if (!$input){
  print "\n\nUSAGE: -i path/filename OR -input path/filename\n\n";
  exit(0);
}

# ----- check user previledge
my $user = `whoami`; chomp $user;
if ($user ne "wormpub"){
  print "\n\nYou need to be wormpub to do this!\n\n";
  exit(0);
}

# ----- global variables
my $ga = init Geneace();
my $database = $ga->curr_db;
my $tace = &tace;
my $allele_dir = "/wormsrv1/geneace/ALLELE_DATA";
my $rundate = &rundate;

# ----- prepare a locus to gene_id conversion hash

my %Gene_info = $ga->gene_info();

# ----- prepare clone coords from GFF file
my %clone_info = $ga->get_clone_chrom_coords(); # key is clone, values: chrom, start_coord, end_coord

# ----- parse NBP flatfile
my (%NBP_info, %chrom_NBP_allele);

my @NBP = `cut -f 1-5,8,12 $input`;

open(CHECK, ">$allele_dir/NBP_alleles_to_check_by_hand") || die $!;

foreach(@NBP){
  chomp;
  my($allele, $locus, $clone, $indel, $pheno, $primers, $mapper)=split(/\t+/,$_);

  if ($allele =~ /^tm0.+/){
    $allele =~ /(tm)0(\d+)/;  # remove leading 0 in allele name
    $allele = $1.$2;	
  }

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
  $clone =~ s/\..+// if $clone ne "cosmis";

  #----- process primers info in $primer (4 primers)
  my ($ext_f, $int_b, $ext_b, $int_f) = split(/,/, $primers);

  # ----- get rid of hidden/strange space and double quotes
  $ext_f =~ s/\s|\"//g;
  $ext_b =~ s/\s|\"//g;
  $int_f =~ s/\s|\"//g;
  $int_b =~ s/\s|\"//g;

  my @primers;
  if ($ext_f && $int_f && $int_b && $ext_b){
    push(@primers, $ext_f, $int_f, $int_b, $ext_b);
  }
  else {
    @primers = "NA"
  }
  push(@indel_info, $site_L, $site_R, $remark_1, $insert, $L_clone, $R_clone, $indel);
  push(@{$NBP_info{$allele}}, $locus, $clone, \@indel_info, $pheno, \@primers, $mapper);
  push(@{$chrom_NBP_allele{$clone_info{$clone}->[0]}}, $allele); # key: chromosome letter
}

# ----- acefile name and location

my $acefile = $input;
$acefile =~ /.+\/(.+)\.txt/;
$acefile = $1.".ace";

open(ACE, ">$allele_dir/$acefile") || die $!;

my $log = "/wormsrv2/logs/$acefile.$rundate";
`rm -f $log` if -e $log;

open(LOG, ">$log") || die $!;
print LOG "\n\nLoaded $acefile to Geneace . . .\n" if !$debug;
print LOG "\n\nLoaded $acefile to CK1_testDB . . .\n" if $debug;
print LOG "--------------------------------------------------\n\n";

get_30_bp_flanks($database);

# ----- upload data to Geneace

my $command=<<END;
pparse $allele_dir/$acefile
save
quit
END

print LOG "\n\n";
$ga->upload_database($ga->geneace, $command, "NBP_allele", $log) if !$debug;
$ga->upload_database($ga->test_geneace, $command, "NBP_allele", $log) if $debug;

# ----- mail notice

my $recipients = "ck1\@sanger.ac.uk, krb\@sanger.ac.uk";
$recipients = "ck1\@sanger.ac.uk" if $debug;
mail_maintainer("Loading NBP_allele", $recipients, $log);


# ----- The End
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

      print ACE "\nAllele : \"$allele\"\n";

      my @indel_info = @{$NBP_info{$allele}->[2]};
      my $L_clone = $indel_info[4];
      my $R_clone = $indel_info[5];

      my $DNA_L =substr($line, $clone_info{$NBP_info{$allele}->[1]}->[1] + $indel_info[0]-1 -31, 30) if $L_clone eq "NA";
         $DNA_L =substr($line, $clone_info{$L_clone}->[1] + $indel_info[0]-1 -31, 30)                if $L_clone ne "NA";

      my $DNA_R =substr($line, $clone_info{$NBP_info{$allele}->[1]}->[1] + $indel_info[1]-1,     30) if $R_clone eq "NA";
         $DNA_R =substr($line, $clone_info{$R_clone}->[1] + $indel_info[1]-1,     30)                if $R_clone ne "NA";

      print ACE "Sequence \"$NBP_info{$allele}->[1]\"\n";

      if ( $NBP_info{$allele}->[0] =~ /\w{3,3}-.+/ ){
	my $locus = lc($NBP_info{$allele}->[0]);  # NBP data often use capitalized locus name
	print ACE "Gene \"$Gene_info{$locus}{'Gene'}\"  \/\/$NBP_info{$allele}->[0]\n" if exists $Gene_info{$locus}{'Gene'};
	print LOG "$allele is linked to non-existent locus ($locus)\n" if !exists $Gene_info{$locus}{'Gene'};
      }
      else {
	print ACE "Predicted_gene \"$NBP_info{$allele}->[0]\"\n";
      }

      # ----- remove tm and use the rest in as allele id to link back to NBP allele webpage
      my $allele_id = $allele;
      $allele_id =~ s/tm//;

      # ----- acefile
      print ACE "Flanking_PCR_product     \"$allele"."_external\"\n";
      print ACE "Flanking_PCR_product     \"$allele"."_internal\"\n";
      print ACE "Author \"$NBP_info{$allele}->[5]\"\n" if $NBP_info{$allele}->[5] ne "NA";
      print ACE "Flanking_sequences \"$DNA_L\" \"$DNA_R\"\n";
      print ACE "Deletion\n" if $indel_info[3] eq "NA";
      if ($indel_info[3] ne "NA"){
	print ACE "Deletion_with_insertion \"", lc($indel_info[3]),"\"\n";
	print ACE "Method \"Deletion_and_insertion_allele\"\n";
      }
      else {
	print ACE "Method \"Deletion_allele\"\n";
      }
      print ACE "NBP_allele\n";
      print ACE "Species \"Caenorhabditis elegans\"\n";
      print ACE "Phenotype \"$NBP_info{$allele}->[3]\"\n";
      print ACE "Mutagen \"TMP\/UV\"\n";
      print ACE "Location \"FX\"\n";
      print ACE "MAP \"$chrom\"\n"; 
      print ACE "Remark \"Mutations at cosmid coordinates: $indel_info[6]\"\n";
      print ACE "Remark \"<A href='http:\\/\/www.shigen.nig.ac.jp\/c.elegans\/mutants/DetailsSearch?lang=english&seq=$allele_id' target=_new> more on $allele...<\/A>\"\n";

      # ----- only for those tm allele having primer information
      if (scalar @{$NBP_info{$allele}->[4]} != 1){

	print ACE "\nPCR_product : $allele"."_external\n";
	print ACE "Oligo $allele"."_external_f\n";
	print ACE "Oligo $allele"."_external_b\n";
	print ACE "Flanks_deletion \"$allele\"\n";
	
	print ACE "\nOligo : \"$allele"."_external_f\"\n";
	print ACE "Sequence \"$NBP_info{$allele}->[4]->[0]\"\n"; 
	print ACE "Length ", length($NBP_info{$allele}->[4]->[0]), "\n";
        print ACE "PCR_product \"$allele"."_external\"\n";

      	print ACE "\nOligo : \"$allele"."_external_b\"\n";
	print ACE "Sequence \"$NBP_info{$allele}->[4]->[3]\"\n"; 
	print ACE "Length ", length($NBP_info{$allele}->[4]->[3]), "\n";
        print ACE "PCR_product \"$allele"."_external\"\n";

        print ACE "\nPCR_product : $allele"."_internal\n";
        print ACE "Oligo $allele"."_internal_f\n";
        print ACE "Oligo $allele"."_internal_b\n";
        print ACE "Flanks_deletion \"$allele\"\n";

        print ACE "\nOligo : \"$allele"."_internal_f\"\n";
	print ACE "Sequence \"$NBP_info{$allele}->[4]->[1]\"\n"; 
	print ACE "Length ", length($NBP_info{$allele}->[4]->[1]), "\n";
        print ACE "PCR_product \"$allele"."_internal\"\n";

      	print ACE "\nOligo : \"$allele"."_internal_b\"\n";
	print ACE "Sequence \"$NBP_info{$allele}->[4]->[2]\"\n"; 
	print ACE "Length ", length($NBP_info{$allele}->[4]->[2]), "\n";
        print ACE "PCR_product \"$allele"."_internal\"\n";
      }
    }
  }
}

__END__
