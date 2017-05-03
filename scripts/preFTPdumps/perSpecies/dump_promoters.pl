#!/usr/bin/env perl

use Getopt::Long;
use IO::File;
use Storable;

use Wormbase;
use Log_files;

use strict;
use Ace::Sequence;

use constant UPSTREAM => 2_500;  # how many bases upstream to do

my ($species,$format,$store,$debug,$test,$database);
GetOptions(
     'species=s'  => \$species,
     'format=s'   => \$format,
     'store=s'    => \$store,
     'debug=s'    => \$debug,
     'test'       => \$test,
     'database=s' => \$database,
)||die(@!);


my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")
}else {
  $wormbase = Wormbase->new( -debug => $debug, -test => $test,-autoace=> $database,-organism => $species)
}

my $log = Log_files->make_build_log($wormbase);

# Establish a connection to the database.
$log->write_to("connecting to ${\$wormbase->autoace}\n");
my $dbh = Ace->connect(-path => $wormbase->autoace)||$log->log_and_die(Ace->error);

my $file = $wormbase->reports . '/'.
   join('.',$wormbase->gspecies_name,$wormbase->ncbi_bioproject,$wormbase->get_wormbase_version,'potential_promotors.fa');
my $of = IO::File->new($file,'w');
$log->write_to("writing to $file\n");

warn "finding predicted genes...\n";
my @genes = $db->fetch(-query => "find Gene Species=\"${\$wormbase->long_name}\";Live;Sequence");

$log->write_to("found ",scalar(@genes)," predicted genes\n");

# create a sequence from each one and find the closest upstream transcript
for my $g (@genes) {
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

$log->mail;
$of->close;
