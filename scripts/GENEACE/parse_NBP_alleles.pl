#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: krb $
# Last updated on: $Date: 2004-09-14 09:32:30 $ 

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use GENEACE::Geneace;
use Getopt::Long;

# ----- command line options
my ($input, $debug);
GetOptions ("i|input=s"  => \$input);

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

open(CHECK, ">$allele_dir/NBP_alleles_to_check_by_hand.$rundate") || die $!;

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

  #print "(1) $allele\n(2)$locus\n(3)$clone\n(4)$indel\n(5)$pheno\n(6)$primers\n(7)$mapper\n";# if $allele eq "tm1223";
}

# ----- acefile names, log file and location

my $acefile = $input;
$acefile =~ /(.+)\.txt/;
$acefile = $1.".ace";

open(ACE, ">$acefile") || die $!;

# last update has the -D lines to delete everything from last update as the NBP sends full updates, not incremental
my @last_updates = glob("$allele_dir/NBP_last_update.ace.*");
my $last_update = $last_updates[-1];

open(ACE_del, ">$allele_dir/NBP_last_update.ace.$rundate") || die $!;

my $acelog = $acefile;
$acelog =~ s/^.+\///;

my $log = "/wormsrv2/logs/$acelog.$rundate";
`rm -f $log` if -e $log;

open(LOG, ">$log") || die $!;
print LOG "\n\nLoaded $acefile to Geneace . . .\n";
print LOG "--------------------------------------------------\n\n";

get_30_bp_flanks($database);

# ----- upload data to Geneace

# parse in the $last_update file first to remove everthing from last update and then upload new updates
#my $command="pparse $last_update\npparse $acefile\nsave\nquit\n";

print LOG "\n\n";
#$ga->upload_database($ga->geneace, $command, "NBP_allele", $log);

# ----- mail notice

my $recipients = "krb\@sanger.ac.uk";
mail_maintainer("NBP allele update", $recipients, $log);

# ----- moving old files to ARCHIVE dir

my @old_txt   = glob("$allele_dir/NBP_*txt"); pop @old_txt;
my @old_ace   = glob("$allele_dir/NBP_*ace"); pop @old_ace;	
my @old_check = glob("$allele_dir/NBP_alleles_to_check_by_hand*"); pop @old_check;
my @old_last  = glob("NBP_last_update.ace*"); pop @old_last;

`mv @old_txt @old_ace @old_check @old_last $allele_dir/ARCHIVE/`;


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
      print ACE_del "\nAllele : \"$allele\"\n";
      my @indel_info = @{$NBP_info{$allele}->[2]};
      my $L_clone = $indel_info[4];
      my $R_clone = $indel_info[5];

      my $DNA_L =substr($line, $clone_info{$NBP_info{$allele}->[1]}->[1] + $indel_info[0]-1 -31, 30) if $L_clone eq "NA";
         $DNA_L =substr($line, $clone_info{$L_clone}->[1] + $indel_info[0]-1 -31, 30)                if $L_clone ne "NA";

      my $DNA_R =substr($line, $clone_info{$NBP_info{$allele}->[1]}->[1] + $indel_info[1]-1,     30) if $R_clone eq "NA";
         $DNA_R =substr($line, $clone_info{$R_clone}->[1] + $indel_info[1]-1,     30)                if $R_clone ne "NA";

      print ACE "Sequence \"$NBP_info{$allele}->[1]\"\n";
      print ACE_del "-D Sequence \"$NBP_info{$allele}->[1]\"\n";

      if ( $NBP_info{$allele}->[0] =~ /\w{3,4}-.+/ ){
	my $locus = lc($NBP_info{$allele}->[0]);  # NBP data often use capitalized locus name
	print ACE "Gene \"$Gene_info{$locus}{'Gene'}\"  \/\/$NBP_info{$allele}->[0]\n" if exists $Gene_info{$locus}{'Gene'};
	print ACE_del "-D Gene \"$Gene_info{$locus}{'Gene'}\"  \/\/$NBP_info{$allele}->[0]\n" if exists $Gene_info{$locus}{'Gene'};
	print LOG "$allele is linked to a locus ($locus) which does not associate with a Gene id\n" if !exists $Gene_info{$locus}{'Gene'};
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
      print ACE_del "-D Flanking_PCR_product     \"$allele"."_external\"\n";
      print ACE_del "-D Flanking_PCR_product     \"$allele"."_internal\"\n";
      print ACE_del "-D Author \"$NBP_info{$allele}->[5]\"\n" if $NBP_info{$allele}->[5] ne "NA";
      print ACE_del "-D Flanking_sequences \"$DNA_L\" \"$DNA_R\"\n";
      print ACE_del "-D Deletion\n" if $indel_info[3] eq "NA";

      if ($indel_info[3] ne "NA"){
	print ACE "Deletion_with_insertion \"", lc($indel_info[3]),"\"\n";
	print ACE "Method \"Deletion_and_insertion_allele\"\n";
	print ACE_del "-D Deletion_with_insertion \"", lc($indel_info[3]),"\"\n";
	print ACE_del "-D Method \"Deletion_and_insertion_allele\"\n";
      }
      else {
	print ACE "Method \"Deletion_allele\"\n";
	print ACE_del "-D Method \"Deletion_allele\"\n";
      }
      print ACE "NBP_allele\n";
      print ACE "Species \"Caenorhabditis elegans\"\n";
      print ACE "Phenotype \"$NBP_info{$allele}->[3]\"\n";
      print ACE "Mutagen \"TMP\/UV\"\n";
      print ACE "Location \"FX\"\n";
      print ACE "MAP \"$chrom\"\n"; 
      print ACE "Remark \"Mutations at cosmid coordinates: $indel_info[6]\"\n";
      print ACE "Remark \"<A href='http:\\/\/www.shigen.nig.ac.jp\/c.elegans\/mutants/DetailsSearch?lang=english&seq=$allele_id' target=_new> more on $allele...<\/A>\"\n";
      print ACE_del "-D NBP_allele\n";
      print ACE_del "-D Species \"Caenorhabditis elegans\"\n";
      print ACE_del "-D Phenotype \"$NBP_info{$allele}->[3]\"\n";
      print ACE_del "-D Mutagen \"TMP\/UV\"\n";
      print ACE_del "-D Location \"FX\"\n";
      print ACE_del "-D MAP \"$chrom\"\n"; 
      print ACE_del "-D Remark \"Mutations at cosmid coordinates: $indel_info[6]\"\n";
      print ACE_del "-D Remark \"<A href='http:\\/\/www.shigen.nig.ac.jp\/c.elegans\/mutants/DetailsSearch?lang=english&seq=$allele_id' target=_new> more on $allele...<\/A>\"\n";


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

	print ACE_del "\nPCR_product : $allele"."_external\n";
	print ACE_del "-D Oligo $allele"."_external_f\n";
	print ACE_del "-D Oligo $allele"."_external_b\n";
	print ACE_del "-D Flanks_deletion \"$allele\"\n";
	
	print ACE_del "\nOligo : \"$allele"."_external_f\"\n";
	print ACE_del "-D Sequence \"$NBP_info{$allele}->[4]->[0]\"\n";
	print ACE_del "-D Length ", length($NBP_info{$allele}->[4]->[0]), "\n";
        print ACE_del "-D PCR_product \"$allele"."_external\"\n";

      	print ACE_del "\nOligo : \"$allele"."_external_b\"\n";
	print ACE_del "-D Sequence \"$NBP_info{$allele}->[4]->[3]\"\n";
	print ACE_del "-D Length ", length($NBP_info{$allele}->[4]->[3]), "\n";
        print ACE_del "-D PCR_product \"$allele"."_external\"\n";

        print ACE_del "\nPCR_product : $allele"."_internal\n";
        print ACE_del "-D Oligo $allele"."_internal_f\n";
        print ACE_del "-D Oligo $allele"."_internal_b\n";
        print ACE_del "-D Flanks_deletion \"$allele\"\n";

        print ACE_del "\nOligo : \"$allele"."_internal_f\"\n";
	print ACE_del "-D Sequence \"$NBP_info{$allele}->[4]->[1]\"\n";
	print ACE_del "-D Length ", length($NBP_info{$allele}->[4]->[1]), "\n";
        print ACE_del "-D PCR_product \"$allele"."_internal\"\n";

      	print ACE_del "\nOligo : \"$allele"."_internal_b\"\n";
	print ACE_del "-D Sequence \"$NBP_info{$allele}->[4]->[2]\"\n";
	print ACE_del "-D Length ", length($NBP_info{$allele}->[4]->[2]), "\n";
        print ACE_del "-D PCR_product \"$allele"."_internal\"\n";

      }
    }
  }
}

__END__

=head2 NAME - parse_NBP_alleles.pl
                                                                                                                                                 
B<USAGE>

	-i input_file

B<Interval of update>

        Normally every 2-3 wks. Shohei Mitani sends out update file (full update) in xls format.

B<get lesion flanks>	

	NBP allele lesions are marked up as clone coordinates and do not supply flanking sequences in the flatfile.
	Hence, the script takes these coordinates and retrieve 30 bp flanks.
	How:
	(1) grep clone info (chrom, start, end, clone_name) from GFF file: /wormsrv2/autoace/GFF_SPLITS/WSXXX/CHROMOSOME_*.clone_acc.gff`;
            (via get_clone_chrom_coords() of Geneace.pm).
	(2) fetch flanks from /nfs/disk100/wormpub/DATABASES/current_DB/CHROMOSOMES/CHROMOSOME_*.dna" file using substr()                                                                                                                                  

B<preparing files>

	Remove Mac carriage return in the flatfle from NBP: tr '\r' '\n' < NBP_mutants_file > NBP_alleles.yymmdd, 
	The processed file should be saved in /wormsrv1/geneace/ALLELE_DATA/ (need to be wormpub for write acess)
	
B<script generated files>

	(1) NBP_alleles.yymmdd.ace
	(2) NBP_last_update.ace.yymmdd

	The script uploads (1) and (2) from last update to geneace. The (2) file contains -D lines which deletes only last update 
	info from NBP for each tm alleles. Some alleles may have parsing problems due to their format. The script will flag them and 
	generate a file NBP_alleles_to_check_by_hand. 
	Typically, the fix-by-hand cases are due to 
	(1) the free-styled strings embedded in between two lesion coordinates of an 
	    allele. Eg, the line "tm325 -> 14441/14442- 33 bp addition-14971/14972  (530 bp deletion)" indicates the script has no 
	    knowledge about the "33 bp addition" and what the 33 bp additions are: so the script only knows it is a deletion allele. 
	    Basically, you basically need to look at tm325 allele in Geneace and see what info is missing from the output line: ie, 
	    you should change Deletion tag to Deletion_with_insertion tag and change the Method to Deletion_and_insertion.

B<Contact person>

	 Shohei Mitani (mitani1@research.twmu.ac.jp)


