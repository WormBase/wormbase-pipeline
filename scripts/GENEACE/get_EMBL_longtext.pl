#!/usr/local/bin/perl5.8.0 -w

# get_EMBL_longtext.pl

# Usage: perl5.8.0 this_script -i input -o output


# Author: Chao-Kung Chen Jan 24 2003

# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-05-25 12:49:00 $ 

use strict;
use Getopt::Long;

#################################################
# variables and command line options with aliases
#################################################

my ($input, $output, $locus, $ac, @ACs, %AC_info, @all_longtext, %ac_longtext, $num, $start, $end, $lt);

$start = `date +%H:%M:%S`; chomp $start;

GetOptions ("i|input=s"       => \$input,
            "o|output=s"      => \$output,
            "l|locu=s"        => \$locus,
            "ac|accession=s"  => \$ac);

############################
# read in EMBL ACs from file
############################

if ($input){
  open(IN, $input) || die "Can't read file $input";
  while (<IN>){
    chomp;
    # $1 = AC
    # $2 = gene id
    # $3 = ID
    # $4 = SV
    # $5 = Prot_id (prefix.suffix)
    # $6 = DE
    # $7 = species

    if ($_ =~ /^(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)/){
      push(@ACs, $1);
      push(@{$AC_info{$1}}, $2, $3, $4, $5, $6, $7);
      #print "$1\t$2\t$3\t$4\t$5\t$6====\n";
    }
  }

  ############################
  # fetch EMBL entry in one go
  ############################

  print "Fetching all EMBL accessions for longtext in one go via getz.....\n\n";
  foreach (@ACs){
  #  push(@all_longtext, `getz -e "[embl-acc:$_]"`);
    push(@all_longtext, `pfetch -F $_`);
  }
}
else {
  push(@all_longtext, `getz -e "[embl-acc:$ac]"`);
}


############################################
# processing EMBL entry without DNA sequence
############################################

my ($sv, $version, $id, $pid);

$num=0;
foreach (@all_longtext){

  if ($_ =~ /^no match\n/){
    push(@{$ac_longtext{$ACs[$num]}}, "NA"); 
    print "Can't pfetch $ACs[$num]\n"; $num++; next
  }

  if ($_ =~ /^ID\s+(\w+)\s+.+/ ){$id = $1};
  if ($_ =~ /^SV\s+(\w+\.(\d+))/ ){$sv = $1; $version = $2}

  if ($_ !~ /^SQ.+|^\s+/ && $_ !~ /^\/\// && $_ !~ /^no match\n/){
    push(@{$ac_longtext{$ACs[$num]}}, $_) if !$ac;
    push(@{$ac_longtext{$ac}}, $_) if $ac;
    next;
  }

  if ($_ =~ /^\/\//){$num++}
}

#############################
# output acefile for longtext
#############################

print "\nWriting longtext ace file .....\n\n";

open(FH, ">$output") || die $!;

foreach (keys %ac_longtext){

  if ($ac){
    my $de = `pfetch -D $ac`;
    $de =~ s/$id|$ac|\.\d+|\n|\"//g; $de =~ s/^\s+//g;

    print FH "\n\nGene : \"$locus\"\n";
    print FH "Other_sequence \"$ac\"\n";
    print FH "\nSequence : \"$ac\"\n";
    print FH "Database  \"EMBL\" \"NDB_ID\" \"$id\"\n";
    print FH "Database  \"EMBL\" \"NDB_AC\" \"$ac\"\n";
    print FH "Database  \"EMBL\" \"NDB_SV\" \"$ac.$version\"\n";
    print FH "DB_annotation  \"EMBL\" \"$ac\"\n";
    print FH "Species  \"Caenorhabditis elegans\"\n";
    print FH "Title \"$de\"\n";
    print FH "Gene \"$locus\"\n";
  }
  else {
   # print $pid = $AC_info{$_}->[3], "\n";
    $pid = $AC_info{$_}->[3];
    $pid =~ /(.+)\.(\d+)/;
    $pid = $1; $version = $2;
    my $species = $AC_info{$_}->[5];

    print FH "\n\nGene : \"$AC_info{$_}->[0]\"\n";
    print FH "Other_sequence \"$_\"\n";
    print FH "\nSequence : \"$_\"\n";
    print FH "Database  \"EMBL\" \"NDB_ID\" \"$AC_info{$_}->[1]\"\n";
    print FH "Database  \"EMBL\" \"NDB_AC\" \"$_\"\n";
    print FH "Database  \"EMBL\" \"NDB_SV\" \"$_.$AC_info{$_}->[2]\"\n";
    print FH "Protein_id  \"$_\" \"$pid\" $version\n";
    print FH "DB_annotation  \"EMBL\" \"$_\"\n";
    print FH "Species  \"$species\"\n";
    print FH "Title \"$AC_info{$_}->[4]\"\n";
    print FH "Gene \"$AC_info{$_}->[0]\"\n";
  }

  $lt = join('',@{$ac_longtext{$_}});
  if($lt ne "NA"){	
    print FH "\nLongText : \"$_\"\n";
    print FH "$lt\/\/\n";
    print FH "***LongTextEnd***\n";
  }
}


close FH;

$end = `date +%H:%M:%S`; chomp $end;
print "\n$0 started at $start, finished at $end\n";

print "\n########## You need to manually enter Protein_id info ##########\n\n" if $ac;

__END__



