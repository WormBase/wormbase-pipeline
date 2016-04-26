#!/software/bin/perl -w
# convert chromosome base Homol_Data to clone based one

use lib $ENV{'CVS_DIR'};
use Getopt::Long;
use Wormbase;
use Coords_converter;
use Storable;

my ($store,$wb,$debug,$test,$species);
GetOptions(
	'store=s'  => \$store,
	'debug=s'  => \$debug,
	'test'     => \$test,
        'species=s'=> \$species,
) || die("bas commandline parameter\n");

if ($store) { 
	$wb = retrieve($store) or croak("Can't restore wormbase from $store\n")
} else {
	$wb = Wormbase->new( -debug    => $debug, 
		             -test     => $test,
                             -organism => $species,
		     )
}

my $cc = Coords_converter->invoke($wb->orgdb,undef,$wb);

my $seq;
my $blastdb;
my %clones;
while (<>) {
  
#Homol_data : "CHROMOSOME_I:wublastx_fly"
#Pep_homol       "FLYBASE:CG16858" "wublastx_fly"  11.398 4571310 4571236 1492 1516
  
  if(/^Homol_data\s+:\s+\"([\w\.]+):(wublastx_\w+)/) { # matches both chromosome names and briggsae contig names like "cb25.fpc4228b"
    $seq = $1;
    $blastdb = $2;
    #print "\nHomol_data : \"$seq:$blastdb\"\n";
    next;
  }
  elsif(/^Pep_homol/) {
    my @data = split;

    my ($orig_start, $orig_end) = ($data[4], $data[5]);

    my($clone,$x, $y) = $cc->LocateSpan($seq, $orig_start, $orig_end);
    die unless ($clone and $x and $y);

    ($data[4], $data[5]) = ($x, $y);
    $data[5] = $y;
    if($data[9]) {
      # for features that were reverse strand originally, 
      # col 10 marks the *end* of the of ungapped block
      if ($data[5] - $data[4] == $orig_end - $orig_start) {
        # simple case; transformation has preserved orientation,
        # so column9 is simple offset 
        $data[9] = $data[4] + ($data[9] - $orig_start);
      } else {
        if ($orig_end > $orig_start) {
          # feature was forward, is now reverse; col 9 aligned with start of block,
          # should now align with end
          my $offset = $data[9] - $orig_start;
          $data[9] = $data[4] - $offset;
        } else {
          # feature was reverse, is now forward; col 9 aligned with end of block,
          # should now align with start
          my $offset = $orig_start - $data[9];
          $data[9] = $data[4] + $offset;
        }
      }
    }
    push(@{$clones{$clone}}, \@data);
  }
}

foreach my $clone (keys %clones) {
  print "\nSequence : $clone\n";
  my $size = $cc->Superlink_length($clone);
  print "Homol_data \"$clone:$blastdb\" 1 $size\n";
  
  print "\nHomol_data : \"$clone:$blastdb\"\n";	
  foreach my $hit (@{$clones{$clone}}) {
    print join("\t",@$hit)."\n";
  }
}
