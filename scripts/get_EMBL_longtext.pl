#!/usr/local/bin/perl5.8.0 -w

# get_EMBL_longtext.pl 

# Usage: perl5.8.0 this_script -i input -o output


# Author: Chao-Kung Chen Jan 24 2003

# Last updated by $Author: ck1 $
# Last updated on: $Date: 2003-07-18 13:00:56 $ 

use strict;
use Getopt::Long;

#################################################
# variables and command line options with aliases
#################################################

my ($input, $output, $ac, @ACs, %AC_locus, @all_longtext, %ac_longtext, $num, $start, $end, $lt);

$start = `date +%H:%M:%S`; chomp $start;

GetOptions ("i|input=s"     => \$input,
            "o|output=s"    => \$output);

############################
# read in EMBL ACs from file
############################

print "\n\nRetrieving EMBL AC numbers .....\n\n";

open(IN, $input) || die "Can't read file $input";
while (<IN>){
  if ($_ =~ /^AC\s+(.+);/){
    push(@ACs, $1);
  }
} 

print "\nGot ", scalar @ACs, " EMBL AC numbers\n\n";

############################
# fetch EMBL entry in one go
############################

print "Fetching all EMBL entries for longtext in one go via getz.....\n\n";
foreach (@ACs){
  push(@all_longtext, `getz -e "[embl-acc:$_]"`); 
} 

############################################
# processing EMBL entry without DNA sequence
############################################

$num=0;
foreach (@all_longtext){
  if ($_ =~ /^no match\n/){
    push(@{$ac_longtext{$ACs[$num]}}, "NA"); 
    print "Can't pfetch $ACs[$num]\n"; $num++; next
  }
  if ($_ !~ /^SQ.+|^\s+/ && $_ !~ /^\/\// && $_ !~ /^no match\n/){
    push(@{$ac_longtext{$ACs[$num]}}, $_); next
  } 
  if ($_ =~ /^\/\//){$num++}
}   

#############################
# output acefile for longtext
#############################

print "\nWriting ace file .....\n\n";

open(FH, ">$output") || die $!;

foreach (keys %ac_longtext){
  $lt = join('',@{$ac_longtext{$_}});  	
  if($lt ne "NA"){	
    print FH "\n\nSequence :\t\"$_\"\n";
    print FH "DB_annotation\tEMBL\t\"$_\"\n";
    print FH "\nLongText : \"$_\"\n";
    print FH "$lt\/\/\n";
    print FH "***LongTextEnd***\n";
  }
}

close FH;

$end = `date +%H:%M:%S`; chomp $end;
print "\nJob started at $start\n";
print "Job finished at $end\n";

__END__



