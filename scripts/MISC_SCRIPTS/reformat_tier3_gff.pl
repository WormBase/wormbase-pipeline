#!/usr/bin/perl
#
# reformat_tier3_gff.pl
#
# Takes the GFF files dumped from EnsEMBL databases using ENSEMBL/scripts/dump_gff3.pl
# and changes the attributes column to match the format of the core species GFF3 files
#
########################################################################################

use strict;
use Wormbase;
use Getopt::Long;
use Log_files;
use Path::Class;

my ($in_gff_file, $out_gff_file, $log_file, $debug);

GetOptions(
    "in|i=s" => \$in_gff_file,
    "out|o=s" => \$out_gff_file,
    "log|l=s" => \$log_file,
    "debug|d:s" => \$debug
    );

my $log = Log_files->make_log($log_file, $debug);
my $wormbase = Wormbase->new(-debug => $debug);

my $target_mappings = get_target_mappings();
my $in_fh = file($in_gff_file)->openr;
my $out_fh = file($out_gff_file)->openw;
my %unknown_species;
while (my $line = $in_fh->getline()) {
    chomp $line;
    if ($line =~ /\s(\S+)\-BLASTX/) {
	my $species = $1;
	$species =~ s/_proteins//;
	my @cols = split("\t", $line);
	my %attr = split(/[;=]/, $cols[8]);
	$log->log_and_die("No name for $species BLASTX line $line\n") unless exists $attr{'Name'};

	if (!exists $target_mappings->{$species}) {
	    $log->log_and_die("No target for $species BLASTX line $line\n") unless exists $attr{'Target'};
	    if ($species eq 'scerevisiae') {
		$cols[8] = 'Target=SGD:' . $attr{'Target'};
	    } elsif ($species eq 'hsapiens' && ($attr{'Target'} =~ /^ENS/ || $attr{'Target'} =~ /^HIT/ || $attr{'Target'} =~ /^OTTHUMP/)){
		$cols[8] = 'Target=ENSEMBL:' . $attr{'Target'};
	    } elsif ($species eq 'dmelanogaster') {
		$cols[8] = 'Target=FLYBASE:' . $attr{'Target'};
	    } elsif ($species eq 'UniProt') {
		$cols[8] = 'Target=SW:' . $attr{'Target'};
	    } else {
		if (!exists $unknown_species{$species}) {
		    $log->write_to("WARNING - no accessions for $species\n");
		    $unknown_species{$species}++;
		}
	    }
	}
	else {
	    my ($name, $position) = $attr{'Target'} =~ /^(\S+)(\s.+)$/;
	    if (exists $target_mappings->{$species}{$name}) {
		$cols[8] = 'Target=' . $target_mappings->{$species}{$name} . $position;
	    }
	    else {
		$log->write_to("WARNING - no mapping found for $name in $species\n");
		$cols[8] = 'Target=' . $attr{'Target'};
	    }
	}
	$line = join("\t", @cols);
    }
    $out_fh->print($line . "\n");
}

$log->mail;
exit(0);

sub get_target_mappings {
    my %target_mappings;

    my %accessors = ($wormbase->species_accessors);
    my $wormpep_dir = dir($wormbase->basedir, 'WORMPEP');
    my @wormpep_children = $wormpep_dir->children();
    for my $wb (values %accessors, $wormbase) {
	$log->write_to('Getting ID to accession mapping for ' . $wb->short_name . "\n");
	my $prefix = $wb->pepdir_prefix;
	for my $subdir(@wormpep_children) {
	    next unless $subdir->is_dir;
	    next unless index($subdir->basename, $prefix . 'pep') == 0;
	    my ($version) = $subdir->basename =~ /pep(\d+)$/;
	    my $accessions_file = $wb->basedir . "/WORMPEP/${prefix}pep${version}/${prefix}pep.accession${version}";
	    next unless -e $accessions_file;
	
	    my $in_fh = file($accessions_file)->openr;
	    while (my $line = $in_fh->getline()) {
		chomp $line;
		my @col = split("\t", $line);
		my $acc = shift @col;
		for my $id (@col) {
		    my $species = $wb->gspecies_name;
		    $species =~ s/_//g;
		    $target_mappings{$species}{$id} = $acc;
		}
	    }
	}
    }

    return \%target_mappings;
}

1;
