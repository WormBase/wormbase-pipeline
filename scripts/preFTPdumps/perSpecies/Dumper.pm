#!/usr/bin/env perl
#
# Utility functions:
#

package Dumper;

use Exporter qw(import);
our @ISA =   qw(Exporter);
our @EXPORT = qw(get_date rewrap);

use Time::localtime;

# yyyy-mm-dd
sub get_date{
 my $tm = localtime;
 my ($day,$month,$year)=($tm->mday,$tm->mon,$tm->year+1900);
 return "$year-$month-$day";
}


# wrap a paragraph to 80 characters
sub rewrap {
    my $text = shift;
    $text =~ s/^\n+//gs;
    $text =~ s/\n+$//gs;
    my @words = split(/\s/,$text);
    my ($para,$line);
    foreach (@words) {
            $line .= "$_ ";
            next if length ($line) < 80;
            $para .= "$line\n";
            $line = undef;
    }
    $para .= "$line\n" if ($line);
    return $para;
}

1;
