#!/usr/local/bin/perl

my %seq = read_fasta (\*STDIN);
my $coor;
my $start = 1;
my $count = 0;
foreach my $name (sort {$a cmp $b} keys %seq) {
    $count++;
    my $length = length $seq{$name};
    my $end = $start+$length-1;
    print "999\t$start\t$end\t$count\tF\t$name\t1\t$length\t+\n";
    $start = $end+1;
}


sub read_fasta {
    # input: Filehandle of fasta file
    # returns: hash with key=name, value=sequence
    local (*FILE) = @_;
    my ($id , $seq , %name2seq);
    while (<FILE>) {
        chomp;
#        if (/^>\S+\s+(\S+)/) {
        if (/^>(\S+)/) {
            my $new_id = $1;
            if ($id) {
                $seq =~ tr/A-Z/a-z/;
                $name2seq{$id} = $seq;
            }
            $id = $new_id ; $seq = "" ;
        } 
        elsif (eof) {
            if ($id) {
                $seq .= $_ ;
                $seq =~ tr/A-Z/a-z/;
                $name2seq{$id} = $seq;
            }
        }
        else {
            $seq .= $_ ;
        }
    }
    return %name2seq;
}
