#! /usr/bin/perl 
# change mRNAs to pseudogenic_transcript
# remove CDSs these had as a parent
# what we're claiming biologically when we do this:
# - these transcripts belong to genes that have no function
# - they look like protein-coding genes but are not (because they have stops)
# - they are probably pseudogenes arising from retrotansposition, i.e. copying an re-integration back into somewhere else in the genome, followed by degradation
@ARGV or die "Usage: $0 in.gff transcripts.txt";

my $f = shift;
my $fh;
if ($f eq "-" ) {
  $fh = STDIN;
} else {
  open ($fh, '<', $f) or die;
}
my $g = shift;
my $gh;
if ($g eq "-") {
  $gh = STDIN;
} else {
  open ($gh, '<', $g) or die;
}
my %ts;
while (<$gh>) {
  chomp;
  $ts{$_}++;
}
my %gs;
sub extract_attribute {
   my ($attr,$str) = @_;
   $str =~ /$attr=(.*?);/;
   $str =~ /$attr=(.*)$/ unless $1;
   chomp $1;
   return $1;
}
while (<$fh>) {
    print and next if /^#/;
    my @c = split "\t";
    if ( @c[2] eq "mRNA" ) {
        my $id = &extract_attribute("ID", @c[8]);
        my $parent =  &extract_attribute("Parent", @c[8]);
        warn "@c[8] without an ID? " unless $id;
        @c[2] = "pseudogenic_transcript" if $ts{$id};
        warn "@c[8] pseudogene with multiple transcripts? " if $ts{$id} and ++$gs{$parent} > 1 ;
    } else {
        @c = qw// if $ts{&extract_attribute("Parent", @c[8])} && @c[2] eq "CDS" ;
    }
    print join "\t", @c;
}
