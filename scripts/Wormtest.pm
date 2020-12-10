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


sub primary_seq_dumps_present {
    my $self;

    $self->_folder_contains_files($self->{'wormbase'}->cdna_dir);

    return;
}


sub recent_citace_dump {
    my $self;

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
    my $self;
    
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


sub _folder_contains_files {
    my ($self, $dirpath) = @_;

    $self->{'log'}->log_and_die("$dirpath does not exist") unless -d $dirpath;
    
    my $dir = dir($dirpath);
    my @filenames;
    for my $file ($dir->children) {
	push @filenames, $file->basename;
    }

    if (@filenames) {
	$self->{'log'}->write_to("$dirpath contains " . 
				 join(', ', @filenames) . "\n");
    }
    else {
	$self->{'log'}->log_and_die("$dirpath is empty\n");
    }

    return;
}



1;
