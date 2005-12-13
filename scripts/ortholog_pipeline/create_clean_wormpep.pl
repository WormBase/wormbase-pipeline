#!/usr/bin/perl

# Create a clean version of wormpep using the longest
# splice variant as the canonical one.

use Bio::Seq;
use Bio::SeqIO;
use IO::File;
use strict;

my $input  = shift;
my $output = shift;

chomp $output;
chomp $input;


my $IDS = {};
my %LONGEST;

# Create a SeqIO output
my $seqout = Bio::SeqIO->new( -format=>'Fasta',
			      -file  =>">$output");

get_lengths();
find_longest();
dump_seqs();

sub get_lengths {
  my $flag = shift;
  my $seqin = Bio::SeqIO->new(-format=>'Fasta',-file=>"$input");
  while (my $seqobj = $seqin->next_seq()) {
    
    # What is the absolute length of the current sequence
    my $seq_length = $seqobj->length();  
    
    # Create a BioSeqIO out object
    my $id = $seqobj->id;
    
    my $clean_id = ($id =~ /(.*\.\d+).*/) ? $1 : $id;
    $IDS->{$clean_id}->{$id} = $seq_length;
    #    print $id,"\t",$seq_length,"\n";
  }
}


sub find_longest {
  foreach my $clean_id (keys %$IDS) {
    my ($prev,$longest_id);
    foreach my $id (keys %{$IDS->{$clean_id}}) {
      my $length = $IDS->{$clean_id}->{$id};
      if ($length > $prev) {
	$prev = $length;
	$longest_id = $id;
      } else {
	next;
      }
    }
    $LONGEST{$longest_id}++;
  }
}


sub dump_seqs {
  print "LONGEST TRANSCRIPTS FOUND: ",scalar keys %LONGEST,"\n";
  my $seqin = Bio::SeqIO->new(-format=>'Fasta',-file=>"$input");
  while (my $seqobj = $seqin->next_seq()) {
    # Create a BioSeqIO out object
    my $id = $seqobj->id;

    next unless (defined $LONGEST{$id});

    print $id,"\n";

    # write the sequence out
    $seqout->write_seq($seqobj);
  }
}
