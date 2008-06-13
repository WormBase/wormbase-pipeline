#!/software/bin/perl -w
# convert chromosome base Homol_Data to clone based one

use lib $ENV{'CVS_DIR'};
use Getopt::Long;
use Wormbase;
use Coords_converter;
use Storable;

my ($store,$wb,$debug,$test);
GetOptions(
	'store=s' => \$store,
	'debug=s' => \$debug,
	'test'    => \$test,
) || die("bas commandline parameter\n");

if ($store) { 
	$wb = retrieve($store) or croak("Can't restore wormbase from $store\n")
} else {
	$wb = Wormbase->new( -debug    => $debug, 
		             -test     => $test, 
		     )
}

#my $wb = retrieve('/nfs/disk100/wormpub/BUILD/autoace/Elegans.store');
my $cc = Coords_converter->invoke($wb->orgdb,undef,$wb);

my $seq;
my $blastdb;
my %clones;
while (<>) {

#Homol_data : "CHROMOSOME_I:wublastx_fly"
#Pep_homol       "FLYBASE:CG16858" "wublastx_fly"  11.398 4571310 4571236 1492 1516

	if(/^Homol_data.*_(\w+):(wublastx_\w+)/) {
		$seq = $1;
		$blastdb = $2;
		#print "\nHomol_data : \"$seq:$blastdb\"\n";
		next;
	}
	elsif(/^Pep_homol/) {
		my @data = split;
		my($clone,$x, $y) = $cc->LocateSpan($seq, $data[4], $data[5]);
		die unless ($clone and $x and $y);
		$data[4] = $x;
		$data[5] = $y;
		if($data[9]) {
			@z = $cc->LocateSpan($seq, $data[9], $data[9]);
			$data[9] = $z[1];
		}
		#print join("\t",@data)."\n";
		push(@{$clones{$clone}}, \@data);
	}
}

foreach my $clone (keys %clones) {
	print "\nSequence : $clone\n";
	my $size = $cc->Superlink_length($clone);
	print "Homol_data 1 $size\n";
	print "\nHomol_data : \"$clone:$blastdb\"\n";
	
	foreach my $hit (@{$clones{$clone}}) {
		print join("\t",@$hit)."\n";
	}
}
