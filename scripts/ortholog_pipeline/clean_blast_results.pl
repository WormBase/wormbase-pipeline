#!/usr/local/bin/perl -w


use lib "/nfs/farm/Worms/Ensembl/bioperl-live";
use Bio::Seq;
use Bio::SeqIO;
use IO::File;
use strict;
use Getopt::Long;

my ($wormpep, $elegans, $brig);

GetOptions (
	    "wormpep:s" => \$wormpep,
	    "elegans:s" => \$elegans,
	    "briggsae:s"=> \$brig
	   );

my $IDS = {};
my %LONGEST;

get_lengths();
find_longest();
foreach (($elegans,$brig)) {
  clean_files($_);
}

sub get_lengths {
  my $flag = shift;
  my $seqin = Bio::SeqIO->new(-format=>'Fasta',-file=>"$wormpep");
  while (my $seqobj = $seqin->next_seq()) {
    
    # What is the absolute length of the current sequence
    my $seq_length = $seqobj->length();  
    
    # Create a BioSeqIO out object
    my $id = $seqobj->id;
    
    my $clean_id = ($id =~ /(.*\.\d+).*/) ? $1 : $id;
    $IDS->{$clean_id}->{$id} = $seq_length;
    #print $id,"\t",$seq_length,"\n";
  }
}


sub find_longest {
  foreach my $clean_id (keys %$IDS) {
    my ($prev,$longest_id);
    $prev = 0;
    foreach my $id (keys %{$IDS->{$clean_id}}) {
      my $length = $IDS->{$clean_id}->{$id};
      if ($length > $prev) {
	$longest_id = $id;
      } else {
	next;
      }
    }
    $LONGEST{$longest_id}++;
  }
}


sub clean_files {
  my $file =shift;
  my $outfile = $file . "_clean";

  open(IN, "<$file") or die "cant open $file\n";
  open(OUT, ">$outfile") or die "cant open $outfile\n";
  while ( <IN> ) {
    my @data = split;
    my $id;
    if( $data[7] == 3 and $data[1] =~ /[a..z]/ ) {
      $id = $data[1];
    }
    elsif  ( $data[7] == 2 and $data[6] =~ /[a..z]/ ) {
      $id = $data[6];
    }
    if ( $id ) {
      next unless (defined $LONGEST{$id});
    }

    print OUT;
  }
  close IN;
  close OUT;
}
