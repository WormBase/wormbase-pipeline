package Wormtest;

use strict;
use Wormbase;

use Const::Fast;
use Path::Class;

const my $MAX_DAYS_SINCE_CITACE_DUMP => 21;
const my $MAX_DAYS_SINCE_GENEACE_COPY => 7;

sub make_build_tester {
    my ($class, $wormbase, $log, $prev_release) = @_;

    # Previous release can be specified as WS release no (with/without WS prefix), or full database path
    # If not passed in the arguments, the environmental variable PREVREL wil be used if set
    if (defined $prev_release) {
	$prev_release = 'WS' . $prev_release if $prev_release =~ /^\d+$/;
	$prev_release = '/nfs/production/panda/ensemblgenomes/wormbase/DATABASES/' . $prev_release 
	    if $prev_release =~ /^WS\d+$/;
    }
    else {
	$prev_release = $ENV{'PREVREL'};
    }

    my $self = {};
    $self->{'wormbase'} = $wormbase;
    $self->{'previous_wormbase'} = Wormbase->new(
	-autoace => $prev_release,
	-species => $wormbase->species
	);
    
    $self->{'log'} = $log;

    bless($self, $class);
    
    return $self;
}


sub build_folder_contents_present {
    my $self = shift;
    
    my ($filenames, $subdirnames) = $self->_folder_contents($self->{'wormbase'}->autoace);
    my %files_present = map {$_ => 1} @$filenames;
    my %subdirs_present = map {$_ => 1} @$subdirnames;
    for my $build_file ('runlog', 'Elegans.store') {
	if (!exists $files_present{$build_file}) {
	    $self->{'log'}->log_and_die("$build_file not present in " . $self->{'wormbase'}->autoace . "\n");
	}
    }
    for my $build_subdir ('acefiles', 'BLAT', 'CHECKS', 'CHROMOSOMES', 'COMMON_DATA', 'database', 'GFF_SPLITS',
			  'logs', 'MISC_OUTPUT', 'ONTOLOGY', 'release', 'REPORTS', 'SEQUENCES', 'SPELL', 'TMP',
			  'TRANSCRIPTS', 'wgf', 'wquery', 'wspec') {  
	if (!exists $subdirs_present{$build_subdir}) {
	    $self->{'log'}->log_and_die("$build_subdir directory not present in " .
					$self->{'wormbase'}->autoace . "\n");
	}
    }
    $self->{'log'}->write_to('All expected files and subdirs present in ' . $self->{'wormbase'}->autoace . "\n");

    return;
}


sub dna_files_have_headers {
    my $self = shift;

    my $chr_dir = dir($self->{'wormbase'}->chromosomes);
    my $chr_count = 0;
    while (my $dna_file = $dir->next) {
	next unless $dna_file->stringify =~ /\.dna$/;
	$chr_count++;
	my $first_line;
	my $dna_fh = $dna_file->openr;
	while ($dna_fh->getline) {
	    chomp;
	    $first_line = $_;
	    last;
	}
	$dna_fh->close;
	$self->{'log'}->log_and_die("Header line not present for $dna_file\n")
	    unless $first_line =~ /^>.+/;
    }
    
    $self->{'log'}->write_to("$chr_count chromosome files found, all with headers present\n");

    return;
}


sub elegans_loaded_first {
    my $self = shift;

    return if $self->{'wormbase'}->species eq 'elegans';
    
    if (-e $self->{'wormbase'}->primary('camace')) {
	$self->write_to('Elegans primary database loaded before ' . $self->{'wormbase'}->species . "\n");
    }
    else {
	$self->log_and_die($self->{'wormbase'}->species . " primary database must be loaded after C. elegans\n");
    }
    
    return;
}


sub primary_seq_dumps_present {
    my $self = shift;

    my $filenames = $self->_folder_contains_files($self->{'wormbase'}->cdna_dir);

    return;
}


sub recent_citace_dump {
    my $self = shift;

    my $dump_dir = dir($self->{'wormbase'}->ftp_upload . '/citace');
    my $most_recent_file;
    my $most_recent_time = 0;
    for my $file ($dump_dir->children) {
	if ((stat($file->stringify))[9] > $most_recent_time) {
	    $most_recent_time = (stat($file))[9];
	    $most_recent_file = $file->basename;
	}
    }

    if ($most_recent_time == 0) {
	$self->{'log'}->log_and_die("No citace dump file found in $CITACE_DUMP_DIR\n");
    } 

    if ((time() - $most_recent_time) < ($MAX_DAYS_SINCE_CITACE_DUMP * 86400)) {
	$self->{'log'}->write_to('Latest citace dump ' . $most_recent_file . 
				 ' last modified @ ' . localtime($most_recent_time) .
				 "\n");				 
    }
    else {
	$self->{'log'}->log_and_die('Latest citace dump ' . $most_recent_file . 
				    " last modified more than $MAX_DAYS_SINCE_CITACE_DUMP days ago (" .
				    localtime($most_recent_time) . ")\n");
    }

    return;
}


sub recent_genace_dump {
    my $self = shift;
    
    my $db_file = $self->{'wormbase'}->primary('geneace') . '/database/ACEDB.wrm';
    my $geneace_last_mod = (stat($db_file))[9];
    if ((time() - $geneace_last_mod) < ($MAX_DAYS_SINCE_GENEACE_COPY * 86400)) {
	$self->{'log'}->write_to("Geneace dumped within last $MAX_DAYS_SINCE_GENEACE_DUMP days: " .
				 localtime($geneace_last_mod) . "\n");
    }
    else {
	$self->{'log'}->log_and_die("Geneace dumped more than $MAX_DAYS_SINCE_GENEACE_DUMP days ago: " .
				 localtime($geneace_last_mod) . "\n");
    }

    return;
}


sub tier2_contigs_dumped {
    my $self = shift;

    return if $self->{'wormbase'}->species eq 'elegans';

    my $contigs_in_db = $self->_nr_contigs;

    my $gff_dir = dir($self->{'wormbase'}->gff_splits);
    while (my $gff_file = $gff_dir->next) {
	next unless $gff_file->stringify =~ /\.gff$/;
	my $contig_count = 0;
	my $gff_fh = $gff_file->openr;
	while ($gff_fh->getline) {
	    $contig_count++ if $_ =~ /^##sequence-region/;
	}
	$gff_fh->close;

	if ($contig_count != $contigs_in_db) {
	    $self->{'log'}->log_and_die("Was expecting $contigs_in_db contigs in " .
					$gff_file->basename . ", but found $contig_count\n");
	}
    }

    $self->{'log'}->write_to("The expected number of contigs were found in all GFF files\n");

    return;
}


sub _folder_contents {
    my ($self, $dirpath) = @_;

    $self->{'log'}->log_and_die("$dirpath does not exist") unless -d $dirpath;
    
    my $dir = dir($dirpath);
    my (@filenames, @subdirnames);
    while (my $file_or_subdir = $dir->next) {
	if (-f $file_or_subdir) {
	    push @filenames, $file_or_subdir->basename;
	}
	else {
	    push @subdirnames, $file_or_subdir->basename;
	}
    }

    return (\@filenames, \@subdirnames);
}


sub _folder_contains_files {
    my ($self, $dirpath) = @_;

    my ($filenames, $subdirnames) = $self->_folder_contents($dirpath);
    if (@$filenames) {
	$self->{'log'}->write_to("$dirpath contains files: " . 
				 join(', ', @$filenames) . "\n");
    }
    else {
	$self->{'log'}->log_and_die("$dirpath is empty\n");
    }

    return $filenames;
}


sub _folder_contains_subdirs {
    my ($self, $dirpath) = @_;

    my ($filenames, $subdirnames) = $self->_folder_contents($dirpath);
    if (@$subdirnames) {
	$self->{'log'}->write_to("$dirpath contains files: " . 
				 join(', ', @$subdirnames) . "\n");
    }
    else {
	$self->{'log'}->log_and_die("$dirpath is empty\n");
    }

    return $subdir_names;
}

sub _nr_contigs {
    my $self = shift;

    my $db = Ace->connect(-path => $self->{'wormbase'}->autoace, -program => $self->{'wormbase'}->tace) or 
	die ('Connection failure: ' . Ace->error);

    my $it = $db->fetch_many(-query => 'find Sequence_collection WHERE Live');
    my $collection_count = 0;
    my $nr_sequences;
    while (my $analysis = $it->next){
	$collection_count++;
	my @sequences = $analysis->at('Sequences');
	$nr_sequences = scalar @sequences;
    }

    $self->{'log'}->log_and_die("More than one Live Sequence_collection object\n") if $collection_count > 1;
    
    return $nr_sequences;
}


1;
