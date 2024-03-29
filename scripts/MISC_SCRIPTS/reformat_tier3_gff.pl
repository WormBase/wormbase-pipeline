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
use File::Copy qw(move);
use Path::Class;

my ($species, $log_file, $debug);

GetOptions(
    "species:s" => \$species,
    "log|l=s" => \$log_file,
    "debug|d:s" => \$debug
    );

$species = 'all' unless $species;
my $log = Log_files->make_log($log_file, $debug);
my $wormbase = Wormbase->new(-debug => $debug);

my %unknown_species;
my $target_mappings = get_target_mappings();

my %tier3_accessors = $wormbase->tier3_species_accessors;
for my $t3 (values %tier3_accessors) {
    next unless $species eq 'all' || lc $species eq $t3->species;
    $log->write_to('Processing ' . $t3->species . " GFF3 file\n");
    my $gff3_file = $t3->sequences . '/' . $t3->species . '.gff3';
    my $backup_file = $gff3_file . '.reformat_backup';
    move $gff3_file, $backup_file;
    
    my $in_fh = file($backup_file)->openr;
    my $out_fh = file($gff3_file)->openw;
    while (my $line = $in_fh->getline()) {
	chomp $line;
	if ($line =~ /\s(\S+)\-BLASTX/) {
	    my $target_species = $1;
	    $target_species =~ s/_proteins//;
	    my @cols = split("\t", $line);
	    my %attr = split(/[;=]/, $cols[8]);
	    $log->log_and_die("No name for $species BLASTX line $line\n") unless exists $attr{'Name'};

	    if (!exists $target_mappings->{$target_species}) {
		$log->log_and_die("No target for $target_species BLASTX line $line\n") unless exists $attr{'Target'};
		if ($target_species eq 'scerevisiae') {
		    $cols[8] = 'Target=SGD:' . $attr{'Target'};
		} elsif ($target_species eq 'hsapiens' && ($attr{'Target'} =~ /^ENS/ || $attr{'Target'} =~ /^HIT/ || $attr{'Target'} =~ /^OTTHUMP/)){
		    $cols[8] = 'Target=ENSEMBL:' . $attr{'Target'};
		} elsif ($target_species eq 'dmelanogaster') {
		    $cols[8] = 'Target=FLYBASE:' . $attr{'Target'};
		} elsif ($target_species eq 'UniProt') {
		    $cols[8] = 'Target=SW:' . $attr{'Target'};
		} else {
		    if (!exists $unknown_species{$target_species}) {
			$log->write_to("WARNING - no accessions for $target_species\n");
			$unknown_species{$target_species}++;
		    }
		}
	    }
	    else {
		my ($name, $position) = $attr{'Target'} =~ /^(\S+)(\s.+)$/;
		if (exists $target_mappings->{$target_species}{$name}) {
		    $cols[8] = 'Target=' . $target_mappings->{$target_species}{$name} . $position;
		}
		else {
		    $log->write_to("WARNING - no mapping found for $name in $target_species\n");
		    $cols[8] = 'Target=' . $attr{'Target'};
		}
	    }
	    $line = join("\t", @cols);
	}
	$out_fh->print($line . "\n");
    }
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
