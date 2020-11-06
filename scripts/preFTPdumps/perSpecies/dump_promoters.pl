#!/usr/bin/env perl

use strict;
use Getopt::Long;
use IO::File;
use Storable;

use lib $ENV{CVS_DIR};

use Wormbase;
use Log_files;
use Ace::Sequence;

use constant UPSTREAM => 2_500;  # how many bases upstream to do

my ($species,$store,$debug,$test,$database,$outfile,$chunk);
GetOptions(
     'species=s'  => \$species,
     'store=s'    => \$store,
     'debug=s'    => \$debug,
     'test'       => \$test,
     'database=s' => \$database,
     'outfile=s'  => \$outfile,
     'chunk=s'    => \$chunk,
)||die(@!);


my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")
}else {
  $wormbase = Wormbase->new( -debug => $debug, 
                             -test => $test,
                             -organism => $species);
}

my $log = Log_files->make_build_log($wormbase);
$database = $wormbase->autoace if not defined $database;

$log->write_to("connecting to $database\n");
my $db = Ace->connect(-path => $database) || $log->log_and_die("Could not connect to $database (".Ace->error.")\n");

$outfile = $wormbase->reports . "/" . "potential_promoters.fa"
    if not defined $outfile;
my $of = IO::File->new($outfile,'w');
$log->write_to("writing to $outfile\n");

warn "finding predicted genes...\n";
my @genes = $db->fetch(-query => "find Gene Species=\"${\$wormbase->long_name}\";Live;Sequence");

$log->write_to("found ",scalar(@genes)," predicted genes\n");

# create a sequence from each one and find the closest upstream transcript
for my $g (@genes) {

    # chunks hook
    if ($chunk){
	    $g->name =~/(\d)$/;
	    next unless "$1" == $chunk;
    }

    my $s = Ace::Sequence->new(-seq    => $g,
			       -offset => -UPSTREAM(),
			       -length => UPSTREAM() 
	) 
	or die "Can't open sequence segment $g: ",Ace->error,"\n";
    
    my $dna = $s->dna;
    
    # find nearest upstream transcript
    if (my @foreign_sequences = 
	grep { $_->info =~ /\.t?\d+[a-z]?$/ and $_->info ne $g } $s->features('Sequence')) {
	my ($rightmost) = sort { $b->end <=> $a->end } @foreign_sequences;
	substr($dna,0,$rightmost->end) = '';  #truncate
    }
    
    print $of ">$g ${\$g->Public_name}\n";
    $dna =~ s/(.{1,60})/$1\n/g;
    print $of $dna;
}

$of->close;
$db->close();
$log->mail;
exit(0);
