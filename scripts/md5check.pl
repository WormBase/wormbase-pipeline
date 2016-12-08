#!/usr/bin/env perl

use lib $ENV{'CVS_DIR'};
use Ace;
use Digest::MD5 qw(md5_hex);
use Getopt::Long;
use strict;
use Species;
use Wormbase;
use Storable;

my ($database,$write,$species,$create,$outfile,$verbose,$autofix,$store,$wormbase,);

GetOptions('species=s'  => \$species,
	   'database=s' => \$database,
           'overwrite'  => \$write,
	   'create'     => \$create,
	   'output:s'   => \$outfile,
	   'verbose'    => \$verbose,
	   'store:s'    => \$store,
);

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} 
else {
  $wormbase = Wormbase->new( -organism  => $species,
                           );
}


my $full_species_name;
if (!defined $species) {die "Need to specify -species\n";}
else {$full_species_name = $wormbase->full_name;}

unless (defined $outfile){$outfile = "md5sum.".$species.".ace";}
open (OUTPUT, ">$outfile") if (($write) || ($create));

if (!defined $create) {&check_md5sum;}
else {&create_md5sum;}
close OUTPUT;
exit(0);


sub check_md5sum{
  print "Checking md5sum of genomic DNA of $full_species_name\n";
  if ($write) {"print - $outfile";}
  else {print "\n\n";}
  my $good = 0;
  my $bad = 0;
  my $db = Ace->connect(-path => $database) or die "Failed to connect to $database\n";
  my $query = "Find Sequence WHERE Species=\"$full_species_name\" AND Method=\"Genomic_canonical\"";
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
  print "Creating md5sums for genomic DNA of $full_species_name - $outfile\n\n";
  my $db = Ace->connect(-path => $database) or die "Failed to connect to $database\n";
  my $query = "Find Sequence WHERE Species=\"$full_species_name\" AND Method=\"Genomic_canonical\"";
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
