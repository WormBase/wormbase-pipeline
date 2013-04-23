#!/usr/local/ensembl/bin/perl -w
# ms2@sanger.ac.uk
# make a gsi index (using S. Eddy's GSI.pm) of one or several fasta files

use strict;
use lib "$ENV{'CVS_DIR'}/BLAST_scripts";
use lib "$ENV{'CVS_DIR'}";
use GSI;
use Getopt::Long;
use Wormbase;
use Storable;
use Log_files;

my( $opt_l, $opt_o, $opt_f, $directory);
my($species, $debug, $store, $test);

my $usage .= "fasta2gsi.pl\n";
$usage .= "-f [fasta file]\n";
$usage .= "-l [list of fasta files]\n";
$usage .= "-o [outfile]    DEFAULT: <-f or -l>.gsi\n";

GetOptions (
	    "l" => \$opt_l,
	    "o" => \$opt_o,
	    "file:s" => \$opt_f,
	    "directory:s"  => \$directory,
	    "test"   => \$test,
	    "store:s"=> \$store,
	    "debug:s"=> \$debug,
	    "species:s"=> \$species,
	    );

my $wormbase;
if( $store ) {
    $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
    $wormbase = Wormbase->new( -debug   => $debug,
			       -test     => $test,
			       -organism => $species
			       );
}
my $log = Log_files->make_build_log($wormbase);
$log->write_to("writing gsi files for $opt_f\n");


############################
# go to the directory that the fasta file lives in

if ($directory) {
  chdir $directory;
}

############################
# -f option:
# read the fasta file, determine the number of entries,
# and store the key and offset value in a hash

if ($opt_f) {

    my $offset = 0;
    my $nfiles = 1;
    my $nfile = 1;
    my $nkeys = 0;
    my %key2offset;
    my %seen;

    open (F , "$opt_f") || die "cannot read $opt_f\n";
    while (<F>) {
	chomp;
	if (/\>(\S+)/) {
	    my $id = $1;
	    # take the suffix away (like TR:Q17521)
	    if ($id =~ /^\w\w:(\S+)/) {
		$id = $1;
	    }
	    if (exists $seen{$id}) {
		die "duplicated key: $1\n";
	    }
	    else {
		$seen{$id}=1;
	    }
	    $nkeys++;
	    $key2offset{$id} = $offset;
	}
	# gives the current byte offset
	$offset = tell;

    }
    close F;
    undef %seen;

# write the GSI file
    if ($opt_o) {
	open (GSI , ">$opt_o") || die "cannot create $opt_o\n";
    }
    else {
	open (GSI , ">$opt_f.gsi") || die "cannot create $opt_f.gsi\n";
    }
    GSI::writeHeaderRecord (*GSI , $nfiles , $nkeys);
    GSI::writeFileRecord (*GSI , $opt_f , $nfile , $GSI::fmt_fasta);

    foreach my $key (sort {$a cmp $b} keys %key2offset) {
	GSI::writeKeyRecord (*GSI , $key , $nfile , $key2offset{$key});
      }

    close GSI;

}

#############################
# -l option:
#############################

if ($opt_l) {

    my $offset;
    my $nfiles = 0;
    my %nfile;
    my $nkeys = 0;
    my %key2offset;
    my %key2file;
    my %seen;

# get the list of databases to index
    my $string;
    open (LIST , "$opt_l") || die "cannot read $opt_l\n";
    while (<LIST>) {
	$string .= $_;
    }
    close LIST;
    $string =~ s/\n/ /g;
    my @files = split /\s+/ , $string;

# read the database files
    foreach my $file (@files) {
	warn "processing $file\n";
	$nfile{$file} = ++$nfiles;
	open (F , "$file") || die "cannot read $file\n";
	while (<F>) {
	    chomp;
	    if (/\>(\S+)/) {
		my $id = $1;
		# take the suffix away (like TR:Q17521)
		if ($id =~ /^\w\w:(\S+)/) {
		    $id = $1;
		}
		if (exists $seen{$id}) {
		    die "duplicated key in $file: $1\n";
		}
		else {
		    $seen{$id}=1;
		}
		$nkeys++;
		$key2offset{$id} = $offset;
		$key2file{$id} = $nfiles;
	    }
	    # gives the current byte offset
	    $offset = tell;
	}
	close F;
    }
    %seen = "";

# write the GSI file
    open (GSI , ">$opt_l.gsi") || die "cannot create $opt_f.gsi\n";
    GSI::writeHeaderRecord (*GSI , $nfiles , $nkeys);

    foreach my $file (sort {$nfile{$a} <=> $nfile{$b}} keys %nfile) {
	GSI::writeFileRecord (*GSI , $file , $nfile{$file} , $GSI::fmt_fasta);
      }

    foreach my $key (sort {$a cmp $b} keys %key2offset) {
	GSI::writeKeyRecord (*GSI , $key , $key2file{$key} , $key2offset{$key}); 
      }

    close GSI;
}

$log->mail;







