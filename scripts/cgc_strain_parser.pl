#!/usr/local/bin/perl5.6.0 -w
# 
# cgc_strain_parser.pl
#
# by Keith Bradnam
#
# Script to convert cgc strain file into ace file for geneace
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2002-07-22 11:18:03 $

use strict;
use Getopt::Std;

##############################
# command-line options       #
##############################
our $opt_i = "";      # input file
our $opt_h = "";      # help
getopts ('i:h');
&usage if ($opt_h);

# set up two output streams, one to add ace file information, and one to delete information
# should you need to remove it in future.

my $rundate    = `date +%y%m%d`; chomp $rundate;
open(STRAIN,">cgc_strain_info_$rundate.ace") || die "cant create output file1\n";
open(DELETE_STRAIN,">cgc_strain_info_$rundate.delete.ace") || die "cant create output file2\n";
open(INPUT, "<$opt_i") || die "Couldn't open input file\n";

# Count how many strains to loop through
my $strain_count = `grep Strain: $opt_i | wc -l`;

# loop through input file using following record separator
$/ = "--------------------";

my $big_counter=0;

while(<INPUT>){
  # drop out of loop before you reach last line of file
  last if $big_counter == $strain_count;
  $big_counter++;

  my $ace_object;
  my $delete_ace_object;

  #remove carriage returns and new lines
  s/\015//g;
  s/\012//g;

  my $strain;
  m/\s+Strain: (.*?)Genotype:/;
  $strain = $1;
  $strain =~ s/\s+//g;
  $ace_object = "Strain : \"$strain\"\n";
  $delete_ace_object = "Strain : \"$strain\"\n";
 
  my $genotype;
  m/Genotype: (.*?)Description:/;
  $genotype = $1;
  $genotype =~ s/\s{2,}/ /g; # condense whitespace to single gap between words
  $genotype =~ s/\s+$//g; # remove trailing whitespace
  $ace_object .= "Genotype \"$genotype\"\n" unless ($genotype eq "");
  $delete_ace_object .= "-D Genotype \"$genotype\"\n" unless ($genotype eq "");

  my $clone;
  my @loci;
  my @loci2;
  my @loci3;
  my @alleles;
  my @alleles2;
  my @alleles3;
  my @rearrangements;
  my @transgenes;

  my $counter =0;
#  print "BEFORE: $genotype\n";

  # find simple locus allele combinations e.g. spt-3(hc184)
  while($genotype =~ m/([Ca-z\-]{3,6}\-\d+)\(([a-z]{1,2}\d+)\)/){
    $loci[$counter] = $1;
    $alleles[$counter] = $2;
    $genotype =~ s/[Ca-z\-]{3,6}\-\d+\([a-z]{1,2}\d+\)//;
    $counter++;
  }
  
  # find chromosomal aberrations e.g. szT1
  $counter = 0;
  while($genotype =~ m/([a-z]{1,2}(Dp|Df|In|T|C)\d+)/){
    $rearrangements[$counter] = $1;
    $genotype =~ s/[a-z]{1,2}(Dp|Df|In|T|C)\d+//;
    $counter++;
  }
  # find transgenes e.g. zhEx11
  $counter = 0;
  while($genotype =~ m/([a-z]{1,2}(Ex|Is)\d+)/){
    $transgenes[$counter] = $1;
    $genotype =~ s/[a-z]{1,2}(Ex|Is)\d+//;
    $counter++;
  }
  
  # find double barrelled alleles (revertants) e.g. daf-12(rh61rh412) 
  $counter = 0;
  while($genotype =~ m/([Ca-z\-]{3,6}\-\d+)\(([a-z]{1,2}\d+)([a-z]{1,2}\d+)\)/){
    $loci2[$counter] = $1;
    $alleles2[$counter] = $2;
    $counter++;
    $alleles2[$counter] = $3;
    $genotype =~ s/[Ca-z\-]{3,6}\-\d+\([a-z]{1,2}\d+[a-z]{1,2}\d+\)//;
    $counter++;
  }


  # find alleles attached to bogus loci e.g. let-?(h661)
  $counter =0;
  while($genotype =~ m/\(([a-z]{1,2}\d+)\)/){
    $alleles3[$counter] = $1;
    $genotype =~ s/\([a-z]{1,2}\d+\)//;
    $counter++;
  }

  # find any skulking loci missed by steps above
  $counter =0;
  while($genotype =~ m/[Ca-z\-]{3,6}\-\d+/){
    $loci3[$counter] = $1;
    $genotype =~ s/[Ca-z\-]{3,6}\-\d+//;
    $counter++;
  }
  



  foreach my $i (@loci) {$ace_object .= "Gene $i\n";}
  foreach my $i (@loci2) {$ace_object .= "Gene $i\n";}
  foreach my $i (@alleles){$ace_object .= "Allele $i\n";}
  foreach my $i (@alleles2){$ace_object .= "Allele $i\n";}
  foreach my $i (@alleles3){$ace_object .= "Allele $i\n";}
  foreach my $i (@rearrangements){$ace_object .= "Rearrangement $i\n";}
  foreach my $i (@transgenes){$ace_object .= "Transgene $i\n";}

  foreach my $i (@loci) {$delete_ace_object .= "-D Gene $i\n";}
  foreach my $i (@loci2) {$delete_ace_object .= "-D Gene $i\n";}
  foreach my $i (@alleles){$delete_ace_object .= "-D Allele $i\n";}
  foreach my $i (@alleles2){$delete_ace_object .= "-D Allele $i\n";}
  foreach my $i (@alleles3){$delete_ace_object .= "-D Allele $i\n";}
  foreach my $i (@rearrangements){$delete_ace_object .= "-D Rearrangement $i\n";}
  foreach my $i (@transgenes){$delete_ace_object .= "-D Transgene $i\n";}
#  print "AFTER:  $genotype\n\n";



  my $description;
  m/Description: (.*?)Mutagen:/;
  $description = $1;
  $description =~ s/\s{2,}/ /g;
  $description =~ s/\s+$//g;
  # get rid of any quotation marks
  $description =~ s/\"//g; 
  # change any URLs present else the double back slash will be treated as a comment
  $description =~ s/http:\/\//URL: /g;
  $ace_object .= "Remark \"$description\"\n" unless ($description eq "");
  $delete_ace_object .= "-D Remark \"$description\"\n" unless ($description eq "");

  my $mutagen;
  m/Mutagen: (.*?)Outcrossed:/;
  $mutagen = $1;
  $mutagen =~ s/\s{2,}/ /g;
  $mutagen =~ s/\s+$//g;
  $ace_object .= "Mutagen \"$mutagen\"\n" unless ($mutagen eq "");
  $delete_ace_object .= "-D Mutagen \"$mutagen\"\n" unless ($mutagen eq "");
  
  my $outcrossed;
  m/Outcrossed: (.*?)Reference:/;
  $outcrossed = $1;
  $outcrossed =~ s/\s{2,}/ /g;
  $outcrossed =~ s/\s+$//g;
  $ace_object .= "Outcrossed\n";
  $delete_ace_object .= "-D Outcrossed\n";

  my $reference;
  if(m/Reference: CGC \#(\d{1,4})\s+/){
    $reference = $1;
    $ace_object .= "Reference \"[cgc$reference]\"\n" unless ($reference eq "");
    $delete_ace_object .= "-D Reference \"[cgc$reference]\"\n" unless ($reference eq "");
  }

  my $made_by;
  m/Made by: (.*?)Received:/;
  $made_by = $1;
  $made_by =~ s/\s{2,}/ /g;
  $made_by =~ s/\s+$//g;
  $ace_object .= "Made_by \"$made_by\"\n" unless ($made_by eq "");
  $delete_ace_object .= "-D Made_by \"$made_by\"\n" unless ($made_by eq "");


  if (m/Received: (\d+\/\d+\/\d+)/){
    my @dates = split(/\//,$1);
    if($dates[2] > 60){
      $ace_object .= "CGC_received \"19$dates[2]-$dates[0]-$dates[1]\"\n";
      $delete_ace_object .= "-D CGC_received \"19$dates[2]-$dates[0]-$dates[1]\"\n";
    }
    else{
      $ace_object .= "CGC_received \"20$dates[2]-$dates[0]-$dates[1]\"\n";
      $delete_ace_object .= "-D CGC_received \"20$dates[2]-$dates[0]-$dates[1]\"\n";
    }
  }

  # always add CGC lab details
  $ace_object .= "Location \"CGC\"\n";
  $delete_ace_object .= "-D Location \"CGC\"\n";

  # print final output to two files
  print STRAIN "$ace_object\n";
  print DELETE_STRAIN "$delete_ace_object\n";

}
close(INPUT);
close(STRAIN);
close(DELETE_STRAIN);

#################################################################

sub usage {
    system ('perldoc',$0);
    exit;       
}


__END__

=pod

=head1 NAME - cgc_strain_parser.pl

=back


=head1 USAGE

=over 4

=item cgc_strain_parser.pl -i <input file>

=back

This script will convert the CGC file of strain information into ace format.
It should be run against the file available at

http://www.cbs.umn.edu/CGC/Strains/gophstrn

The script will write two ace files to your current directory, one to be loaded 
into geneace, and a second to be archived in /wormsrv1/geneace which will have 
delete instructions for removing all the data you have just added.

=over 4

=item MANDATORY arguments:

-i Specify input file (CGC strain file)

=back

=over 4

=item OPTIONAL arguments:

-h (this help page)

=head1 AUTHOR - Keith Bradnam

Email krb@sanger.ac.uk

=cut
