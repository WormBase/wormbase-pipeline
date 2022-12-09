#!/usr/bin/perl
use warnings;
use strict;

# Point to where your input fasta file is located
my $fasta_file = $ARGV[0];

# Declare that you want to read this file, and print the error message
# in case that the input fasta file does not exist at the specified path
open(FASTA, "<", $fasta_file) or die("Probably wrong path: $fasta_file\n");
print("$fasta_file has been read");

# We will linearize the sequences into this hash
my %singleLineSequences;

# Initialize a variable to store sequence ids
my $sequence_id;

# Read the fasta file line by line
while(<FASTA>){
    my $line = $_; chomp($line);

    # if this is a new sequence with the header at the beginnning
    # Extract sequence id using regular expression (\S+) and
    # store it into the variable $sequence_id

    if ($line =~ m/^>(\S+)/){

        $sequence_id = $1; # e.g., YKR054C

        # Reserve an entry in your hash
        # $sequence_id = YKR054C
        # $singleLineSequences{'YKR054C'} = ""
        # No sequence has been added yet, hence the empty value
        $singleLineSequences{$sequence_id} = "";

        }

    # if the line is not a header but part of a sequence,
    # append this part to the corresponding sequence id in
    # your hash entry

    else {

        # paste the current line to the end of the sequence
        $singleLineSequences{$sequence_id} = $singleLineSequences{$sequence_id} . $line;

        }
   }

# Now that your hash contains single line sequences, you can simply
# loop over each sequence in your hash, determine the sequence length,
# and print it out

foreach my $sequence_entry (keys %singleLineSequences){

    # grab a hash entry and store it into a variable
    my $currentSequence = $singleLineSequences{$sequence_entry};

    # determine length of the sequence
    my $lengthSequence = length($currentSequence);

    # print the result: id,length
    print $sequence_entry . "," . $lengthSequence . "\n";
}
