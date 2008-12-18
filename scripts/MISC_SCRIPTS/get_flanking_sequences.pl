#!/usr/bin/perl -w

use lib $ENV{'CVS_DIR'};
use Sequence_extract;
use Wormbase;
use Storable;
use Getopt::Long;

my ($help,$database, $input_file);
my ($seq,$start, $end, $length);
GetOptions ("help"       => \&help,
	    "database=s" => \$database,
	    "file:s"	 => \$input_file,
	    "seq:s"      => \$seq,
	    "start:i"    => \$start,
	    "end:i"      => \$end,
	    "length:i"   => \$length
	   );

$database = '/nfs/disk100/wormpub/DATABASES/current_DB' unless $database;

my $wb = retrieve("$database/Elegans.store");
my $seq_obj      =  Sequence_extract->invoke($database,undef,$wb);
$length = 30 unless $length;

if(defined $input_file){
    open(IF,"<$input_file") or die "cant open $input_file: $!\n";
    while(<IF>){
	($seq,$start, $end) = split(/\s+/,$_);
	&get_flanks;
    }
    close IF;
}
else {
    &get_flanks;
} 



sub get_flanks {
    &check_data($seq, $start, $end) or do{print "FAILED :\n";return};
    my ($c_seq, $c_start, $c_end) = $seq_obj->LocateSpan($seq, $start, $end);
    
    my ($flank1, $flank2);
    if($start < $end) {
	$flank1 = $seq_obj->Sub_sequence($c_seq, $c_start-($length+1), $length);
	$flank2 = $seq_obj->Sub_sequence($c_seq, $c_end,$length);
    }
    else {
	$flank1 = $wb->DNA_string_reverse($seq_obj->Sub_sequence($c_seq, $c_end+1, $length));
	$flank2 = $wb->DNA_string_reverse($seq_obj->Sub_sequence($c_seq, $c_start-($length+1),$length));
    }
    print "$seq\t$start\t$end\t$flank1\t$flank2\n";
}

sub check_data {
    my ($seq, $start, $end) = @_;
    return unless ($seq =~ /[a-zA-Z]/ and  $start =~ /\d/ and $end =~ /\d/);
    return 1;
}

sub help {
    system ("perldoc $0");
    exit;
}

__END__

=pod

=head1 get_flanking_sequences.pl

This is to get flanking sequences for a given piece of the genome.

=head1 Options

=over 4

    -database   Specificy which database to use. Defaults to current_DB
    -seq        Genomic sequence
    -start      start coordinate 
    -end        end coordinate
    -length     length of flanking sequences required
    -file       pass a file of coordinates to find flanks for.  File needs to be 3 columns (seq, start, end)
    -help       show this pod

=back

    e.g. perl  get_flanking_sequences.pl -seq V -start 39000 -end 40000 -length 35

=cut
