#!/software/bin/perl -w

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use Ace;
use File::Copy;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species);
my( $dbpath, $chrom, $allele);
GetOptions (
            "help"      => \$help,
            "debug=s"	=> \$debug,
	    "test"	=> \$test,
	    "store:s"	=> \$store,
	    "chrom:s"	=> \$chrom,
	    "allele:s"	=> \$allele,
	    "database:s"=> \$dbpath
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

my %alleles;
$dbpath = ($dbpath or $wormbase->orgdb);
my $db = Ace->connect('-path' => $dbpath) or $log->log_and_die("cant open Ace connection to $dbpath\n".Ace->error."\n");
my $query = "find Variation ";
$query .= defined $allele ? "\"$allele\"" : "where Natural_variant";

my $vars = $db->fetch_many('-query' => $query);
while (my $cgh = $vars->next){
	if($cgh->ThreePrimeGap and $cgh->FivePrimeGap) {
		$alleles{$cgh->name}->{'five'}  = $cgh->FivePrimeGap->name;
		$alleles{$cgh->name}->{'three'} = $cgh->ThreePrimeGap->name;
	}
	else {
		$log->write_to($cgh->name." missing 3/5 prime gap\n");
	}
}

undef $db;

#read GFF files
my @chroms = $chrom or $wormbase->get_binned_chroms;
foreach (@chroms) {
    my $gff_file = $wormbase->gff."/$_.gff";
    open (GFF,"<$gff_file") or $log->log_and_die("cant open GFF $gff_file: $!\n");
    open (NEW,">$gff_file.new") or $log->log_and_die("cant write new GFF file: $!\n");
    while(<GFF>){
	my @f = split(/\t/,$_);
	if(defined $f[8] and ($f[1] eq 'Allele') and ($f[2] eq 'deletion') ) {
	    my ($allele) = $f[8] =~ /Variation\s+\"([\w\(\)]+)\"/;
	    if($alleles{$allele}->{'five'} and $alleles{$allele}->{'three'}) {
		$f[1] = 'CGH_allele';
		$f[2] = 'deletion';
		$f[3] -= $alleles{$allele}->{'five'};
		$f[4] += $alleles{$allele}->{'three'};

		print NEW join("\t",@f);
	    }
	    else {
		#retain existing but also report missing data.
		$log->write_to("missing gaps for $allele\n") if (/\)/);
	    }
	}
	    print NEW;
    }
    close GFF;
    close NEW;
    move("$gff_file.new","$gff_file") or $log->error("cant move $gff_file.new :$!\n") ;
}


$log->mail();
exit(0);
