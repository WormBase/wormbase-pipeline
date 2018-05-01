#! /usr/bin/perl 
# change mRNAs to noncoding
# remove everything these had as a parent
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

while (<$fh>) {
    my @c = split "\t";
    if ( @c[2] eq "mRNA" ) {
        @c[8] =~ /ID=(.*?);/;
        @c[8] =~ /ID=(.*)$/ unless $1;
        warn "@c[8] without an ID? " unless $1;
        chomp $1;
        @c[2] = "ncRNA" if $ts{$1};
    } else {
        @c[8] =~ /Parent=(.*?);/;
        @c[8] =~ /Parent=(.*)$/ unless $1;
        chomp $1;
        @c = qw// if $ts{$1};
    }
    print join "\t", @c;
}
