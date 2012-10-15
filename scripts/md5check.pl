#!/usr/bin/env perl

use Ace;
use Digest::MD5 qw(md5_hex);
use Getopt::Long;
use strict;

my ($species,$database,$write,$species,$create,$outfile,$verbose,$autofix,);

GetOptions('species=s'  => \$species,
	   'database=s' => \$database,
           'overwrite'  => \$write,
	   'species:s'  => \$species,
	   'create'     => \$create,
	   'output:s'   => \$outfile,
	   'verbose'    => \$verbose,
);


if (!defined $species) {$species = "Caenorhabditis brenneri";}
elsif ($species !~ /\S+\s+\S+/){die "Need full species name not just $species\n";}

unless (defined $outfile){$outfile = "md5sum.".$species.".ace";}
open (OUTPUT, ">$outfile") if (($write) || ($create));

if (!defined $create) {&check_md5sum;}
else {&create_md5sum;}
close OUTPUT;
exit(0);


sub check_md5sum{
  print "Checking md5sum of genomic DNA of $species\n";
  if ($write) {"print - $outfile";}
  else {print "\n\n";}
  my $good = 0;
  my $bad = 0;
  my $db = Ace->connect(-path => $database) or die "Failed to connect to $database\n";
  my $query = "Find Sequence WHERE Species=\"$species\" AND Method=\"Genomic_canonical\"";
  print "$query\n";
  my $seqIt = $db->fetch_many(-query => $query);
  while(my $seq=$seqIt->next){
    my $ace_md5=$seq->MD5;
    my $s=$seq->asDNA;
    $s=~s/>\S+\n//;
    $s=~s/[\s\n]//g;
    my $seq_md5=md5_hex(uc($s));
    if ($seq_md5 eq $ace_md5){
      $good++;
      print "\nSequence $seq\nStored: $ace_md5\nDNAmd5: $seq_md5\n" if ($verbose);
    }
    else{
      print "\nERROR: $seq has a wrong MD5 checksum\nStored: $ace_md5\nDNAmd5: $seq_md5\n";
      $bad++;
	#$db->parse("Sequence : $seq\nMD5 $seq_md5\n\n")||print Ace->error;
      print OUTPUT "Sequence : $seq\nMD5 $seq_md5\n\n" if $write;
    }
  }
  print "\n$good Sequences passed\n$bad Sequences failed\n";
}

sub create_md5sum{
  print "Creating md5sums for genomic DNA of $species - $outfile\n\n";
  my $db = Ace->connect(-path => $database) or die "Failed to connect to $database\n";
  my $query = "Find Sequence WHERE Species=\"$species\" AND Method=\"Genomic_canonical\"";
  print "$query\n";
  my $seqIt = $db->fetch_many(-query => $query);
  while(my $seq=$seqIt->next){
    my $s=$seq->asDNA;
    $s=~s/>\S+\n//;
    $s=~s/[\s\n]//g;
    my $seq_md5=md5_hex(uc($s));
    print OUTPUT "Sequence : $seq\nMD5 $seq_md5\n\n";
  }
}
