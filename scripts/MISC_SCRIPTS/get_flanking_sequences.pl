#!/usr/bin/perl -w

use lib $ENV{'CVS_DIR'};
use Sequence_extract;
use Wormbase;
use Storable;

#my $database = glob("~wormpub/DATABASES/");
my $database = '/nfs/disk100/wormpub/BUILD/autoace';

#my $wb = Wormbase->new( '-test' => 1,
#			'-debug' => 'ar2'
#);

my $wb = retrieve("$database/Elegans.store");

my $seq_obj      =  Sequence_extract->invoke($database,undef,$wb);

my ($seq,$start, $end) = @ARGV;
my ($c_seq, $c_start, $c_end) = $seq_obj->LocateSpan($seq, $start, $end);

my ($flank1, $flank2);
if($start < $end) {
    $flank1 = $seq_obj->Sub_sequence($c_seq, $c_start-31, 30);
    $flank2 = $seq_obj->Sub_sequence($c_seq, $c_end,30);
}
else {
    $flank1 = $wb->DNA_string_reverse($seq_obj->Sub_sequence($c_seq, $c_end+1, 30));
    $flank2 = $wb->DNA_string_reverse($seq_obj->Sub_sequence($c_seq, $c_start-31,30));
}
print "$seq\t$start\t$flank1\t$flank2\n";


	
