=head1 NAME                                                                                                                       

Wormtest - A package for WormBase build tests

=cut                                                                                                                                                                                                                                                                                                    
=head1 SYNOPSIS

    use Wormtest;
    my $wt = Wormtest->new($wb, $log, 'WS279');
    $wt->recent_citace_dump();
                
=head1 DESCRIPTION                                                                                                                                                                                                                                                                                      Various tests for different stages of the WormBase build process.                                                                                                                 

=cut 

package Wormtest;

use strict;
use Wormbase;

use Const::Fast;
use Path::Class;

const my $MAX_DAYS_SINCE_CITACE_DUMP => 21;
const my $MAX_DAYS_SINCE_GENEACE_COPY => 7;
const my $MIN_PERCENT_SEQS_WITH_HOMOLOGY_DATA => 80;
const my $MIN_PERCENT_NON_ELEGANS_GENES => 50;
const my $MIN_BRIGGSAE_CACHE2_SIZE => 550000;
const my $NR_DBXREF_REPORT_COLUMNS => 9;
const my %MOL_TYPES = ( 'elegans'          => [qw(EST mRNA ncRNA OST tc1 RST Trinity Nanopore)],
			'briggsae'         => [qw(mRNA EST Trinity)],
			'remanei'          => [qw(mRNA EST)],
			'brenneri'         => [qw(mRNA EST)],
			'japonica'         => [qw(mRNA EST Trinity)],
			'brugia'           => [qw(mRNA EST Trinity IsoSeq)],
			'pristionchus'     => [qw(mRNA EST)],
			'ovolvulus'        => [qw(mRNA EST Trinity)],
			'sratti'           => [qw(mRNA EST)],
			'tmuris'           => [qw(mRNA EST Trinity IsoSeq)],
			'nematode'         => [qw(EST)],
			'nembase'          => [qw(EST)],
			'washu'            => [qw(EST)],
    );
const my @NON_ELEGANS_CORE_SPECIES => ('Caenorhabditis briggsae', 'Caenorhabditis remanei',
'Caenorhabditis brenneri', 'Caenorhabditis japonica', 'Pristionchus pacificus', 'Brugia malayi',
'Brugia pahangi', 'Onchocerca volvulus', 'Strongyloides ratti', 'Trichuris muris');


=head2 make_build_tester()

    Function: creates Wormtest object for running build tests
    Args:     Wormbase object
              Log_files object
              Previous release string, e.g. WS278 (optional)
    Returns:  Wormtest object

=cut
  
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
	-debug   => $wormbase->debug,
	-test    => $wormbase->test
	);
    
    $self->{'log'} = $log;

    bless($self, $class);
    
    return $self;
}


=head2 blat_files_present()

    Function: checks for presence of of BLAT results files corresponding to all shattered masked sequence files
    Args:     n/a
    Returns:  n/a

=cut
  
sub blat_files_present {
    my $self = shift;

    my $blat_filenames = $self->_folder_contains_files($self->{'wormbase'}->blat);
    my %blat_files_present = map {$_ => 1} @$blat_filenames;

    my $seq_filenames = $self->folder_contains_files($self->{'wormbase'}->cdna);

    my %masked_files_found;
    for my $seq_filename (@$seq_filenames) {
	next unless $seq_filename =~ /^(.+)\.masked_(\d+)/;
	$masked_files_found{$1}{$2} = 1;
    }

    for my $species (keys %MOL_TYPES) {
	next if $species eq $self->{'wormbase'}->species;
	for my $mol_type (keys %masked_files_found) {
	    for my $shattered_file_nr (keys %{$masked_files_found{$mol_type}}) {
		my $expected_blat_file = "${species}_${mol_type}_${shattered_file_nr}.psl";
		$self->{'log'}->log_and_die("$expected_blat_file not found\n")
		    unless exists $blat_files_present{$expected_blat_file};
	    }
	}
    }
    $self->{'log'}->write_to("All expected BLAT output files found\n");

    return;
}


=head2 build_folder_contents_present() 

    Function: checks for presence of expected files and subdirectories within the build directory
    Args:     n/a
    Returns:  n/a

=cut
  
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


=head2 cache_size_sufficient()

    Function: ensures that cache2 size is sufficient for briggsae
    Args:     n/a
    Returns:  n/a

=cut

sub cache_size_sufficient {
    my $self = shift;

    return unless $self->{'wormbase'}->species eq 'briggsae';

    my $cache_def_fh = file($self->{'wormbase'}->autoace . '/wspec/cachesize.wrm')->openr;
    while ($cache_def_fh->getline) {
	next unless $_ =~ /^CACHE2\s=\s(\d+)\s/;
	if ($1 < $MIN_BRIGGSAE_CACHE2_SIZE) {
	    $self->{'log'}->log_and_die("Cache2 size is set as $1 - the minimum recommended size is $MIN_BRIGGSAE_CACHE2_SIZE\n");
	}
	else {
	    $self->{'log'}->write_to("Cache2 size of $1 is sufficient\n");
	}
	last;
    }

    return;
}


=head2 dbxref_report_correctly_formatted

    Function: checks that the DB Xrefs report is present and has the correct number of columns
    Args:     n/a
    Returns:  n/a

=cut

sub dbxref_report_correctly_formatted {
    my $self = shift;

    my $report_file = file($self->{'wormbase'}->reports . '/' . $self->{'wormbase'}->species .
	'.dbxrefs.txt');

    $self->_file_exists($report_file->stringify);

    my $report_fh = $report_file->openr;
    while ($report_fh->getline) {
	next if $_ =~ /^\/\//;
	my @columns = split("\t", $_);
	$self->{'log'}->log_and_die('The following line in ' . $report_file->basename . 
				    " does not have the required number of columns:\n$_\n")
	    unless $NR_DBXREF_REPORT_COLUMNS == @columns;
    }

    $self->{'log'}->write_to('All lines in ' . $report_file->basename .
			     " have the required number of columns\n");

    return;
}


=head2 dna_files_have_headers()

    Function: checks that header lines are present for all sequence files in the chromosome directory
    Args:     n/a
    Returns:  n/a

=cut
  
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


=head2 elegans_loaded_first() 

    Function: checks that primary database is loaded for elegans before any other species
    Args:     n/a
    Returns:  n/a

=cut
  
sub elegans_loaded_first {
    my $self = shift;

    return if $self->{'wormbase'}->species eq 'elegans';
    
    if (-e $self->{'wormbase'}->primary('camace')) {
	$self->{'log'}->write_to('Elegans primary database loaded before ' . $self->{'wormbase'}->species . "\n");
    }
    else {
	$self->{'log'}->log_and_die($self->{'wormbase'}->species . " primary database must be loaded after C. elegans\n");
    }
    
    return;
}


=head2 final_gff_dumps_present

    Function: checks that the final GFF2 and GFF3 dumps are present
    Args:     n/a
    Returns:  n/a

=cut

sub final_gff_dumps_present {
    my $self = shift;


    my $gff_dir = $self->{'wormbase'}->species eq 'elegans' ? $self->{'wormbase'}->chromosomes :
	$self->{'wormbase'}->sequences;

    my @gff_filestems;
    my $contigs_in_db;
    if ($self->{'wormbase'}->species eq 'elegans') {
	for my $chr ('I', 'II', 'III', 'IV', 'V', 'X', 'MtDNA') {
	    push @gff_filestems, $gff_dir . '/CHROMOSOME_' . $chr;
	}
    }
    else {
	$contigs_in_db = $self->_nr_contigs;
	push @gff_filestems, $gff_dir . '/' . $self->{'wormbase'}->species;
    }

    for my $suffix ('.gff', '.gff3') {
	for my $gff_filestem (@gff_filestems) {
	    $self->_file_exists($gff_filestem . $suffix);
	    $self->_expected_seq_region_count
	}
    }

    return;
}


=head2 homology_data_loaded()

    Function: checks that homology data has been loaded onto Sequence objects in AceDB
    Args:     n/a
    Returns:  n/a

=cut

sub homology_data_loaded {
    my $self = shift;

    my $db = Ace->connect(-path => $self->{'wormbase'}->autoace, -program => $self->{'wormbase'}->tace) or 
	die ('Connection failure: ' . Ace->error);

    my $seq_count = 0;
    my $seq_with_homol_data = 0;
    my $i = $db->fetch_many(-query => 'find Sequence_collection WHERE Live');
    while (my $analysis = $i->next){
	my @sequences = $analysis->at('Sequences');
	for my $sequence(@sequences) {
	    $seq_count++;
	    my $seq = $db->fetch(Sequence => $sequence);
	    my @homol = $seq->at('SMap.S_child.Homol_data');
	    $seq_with_homol_data++ if scalar @homol;
	}
    }
    
    my $percent_with_homol = ($seq_with_homol_data / $seq_count) * 100;
    my $min_percent_with_homol = $self->{'wormbase'}->species eq 'elegans' ? 100 : $MIN_PERCENT_SEQS_WITH_HOMOLOGY_DATA;
    if ($percent_with_homol < $min_percent_with_homol) {
	$self->{'log'}->log_and_die("Only $seq_with_homol_data of $seq_count sequence collection sequence objects have " .
				    "homology data, which is less than the cutoff of ${min_percent_with_homol}\%\n");
    print "${seq_with_homol_data}/${seq_count}\n";

}


=head2 masked_files_present()

    Function: checks for presence of unshattered masked sequence files
    Args:     n/a
    Returns:  n/a

=cut
  
sub masked_files_present {
    my $self = shift;

    my $filenames = $self->_folder_contains_files($self->{'wormbase'}->cdna_dir);
    my %files_present = map {$_ => 1} @$filenames;
    for my $mol_type (@{$MOL_TYPES{$self->{'wormbase'}->species}}) {
	$self->{'log'}->log_and_die("Masked $mol_type files not found\n" unless exists $files_present{$mol_type . '.masked'});
    }
    $self->{'log'}->write_to("Masked files for all expected molecule types present\n");

    return;
}


=head2 primary_seq_dumps_present() 

    Function: checks for presence of cDNA dumps for all expected molecule types
    Args:     n/a
    Returns:  n/a

=cut
  
sub primary_seq_dumps_present {
    my $self = shift;

    my $filenames = $self->_folder_contains_files($self->{'wormbase'}->cdna_dir);
    my %files_present = map {$_ => 1} @$filenames;
    for my $mol_type (@{$MOL_TYPES{$self->{'wormbase'}->species}}) {
	$self->{'log'}->log_and_die("$mol_type cDNA dump not present\n" unless exists $files_present{$mol_type});
    }
    $self->{'log'}->write_to("cDNA dumps for all expected molecule types present\n");

    return;
}


=head2 recent_citace_dump()

    Function: checks for presence of recent citace dump in the FTP directory
    Args:     n/a
    Returns:  n/a

=cut
  
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


=head2 recent_geneace_dump()

    Function: checks for recent dump of geneace
    Args:     n/a
    Returns:  n/a

=cut
  
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


=head2 species_merge_successful()

    Function: checks non-elegans core species have been successfully merged by checking for
              the presence of a selection of genes from the previous release
    Args:     n/a
    Returns:  n/a

=cut

sub species_merge_successful {
    my $self = shift;

    my $non_elegans_genes = $self->_previous_release_non_elegans_genes;

    my $db = Ace->connect(-path => $self->{'wormbase'}->autoace,
			  -program => $self->{'wormbase'}->tace) or 
			      die 'Connection failure: ' . Ace->error;
    
    for my $species (keys %$non_elegans_genes) {
	my $species_gene_count = 0;
	for my $gene_id (@{$non_elegans_genes->{$species}}) {
	    my $gene = $db->fetch(Gene => $gene_id);
	    $species_gene_count++ if defined $gene;
	}
	my $percent_found = ($species_gene_count / scalar @{$non_elegans_genes->{$species}}) * 100;
	$self->{'log'}->log_and_die(
	    "Less than ${MIN_PERCENT_NON_ELEGANS_GENES}\% of $species genes selected from " .
	    $self->{'previous_wormbase'}->get_wormbase_version_name . ' (' .
	    join('|', @{$non_elegans_genes->{$species}}) . ") found in merged database\n"
	    ) if $percent_found < $MIN_PERCENT_NON_ELEGANS_GENES;
    }
    $self->{'log'}->write_to("Presence of non-elegans genes suggests successful database merge\n");

    return;
}


=head2 split_gffs_present()
    
    Function: checks that the expected number of GFF files are present
    Args:     expected number of GFFs per chromosome/species
    Returns:  n/a

=cut

sub split_gffs_present {
    my ($self, $expected_nr) = @_;

    my $filenames = _folder_contains_files($self->{'wormbase'}->gff_splits);

    my @prefixes;
    if ($self->{'wormbase'}->species eq 'elegans') {
	for my $chr ('I', 'II', 'III', 'IV', 'V', 'X', 'MtDNA') {
	    push @prefixes, "CHROMOSOME_${chr}_";
	}
    }
    else {
	push @prefixes, undef;
    }

    for my $prefix (@prefixes) {
	my $nr_gff_files = $self->_nr_files_with_prefix_or_suffix($filenames, $prefix, '.gff');
	if ($nr_gff_files != $expected_nr) {
	    my $msg = "Was expecting $expected_nr GFF files ";
	    $msg .= "with prefix $prefix " if defined $prefix;
	    $msg .= "but found $nr_gff_files"; 
	    $self->{'log'}->log_and_die("$msg\n");
	}
    }
    $self->{'log'}->write_to("Found the expected number of GFF files\n");
 
    return;
}


=head2 tier2_contigs_dumped()

    Function: checks that contig dumps are present of all contigs present in the current
              Sequence_collection AceDB object for the species
    Args:     Partial filename string (optional) - only filenames containing this substring
              will be checked
    Returns:  n/a

=cut
  
sub tier2_contigs_dumped {
    my ($self, $substring) = @_;

    return if $self->{'wormbase'}->species eq 'elegans';

    my $contigs_in_db = $self->_nr_contigs;

    my $gff_dir = dir($self->{'wormbase'}->gff_splits);
    while (my $gff_file = $gff_dir->next) {
	next unless $gff_file->stringify =~ /\.gff$/;
	next if defined $substring and index($gff_file->basename, $substring) == -1;	
	my $contig_count = $self->_expected_seq_region_count($gff_file, $contigs_in_db);
    }

    $self->{'log'}->write_to("The expected number of contigs were found in all GFF files\n");

    return;
}


=head2 uniprot_ids_in_wormpep()

    Function: checks that wormpep file and table file contain Uniprot IDs
    Args:     n/a
    Returns:  n/a

=cut

sub uniprot_ids_in_wormpep {
    my $self = shift;

    my $wp_file = file($self->{'wormbase'}->wormpep . '/' . $self->{'wormbase'}->pepdir_prefix .
	'pep' . $self->{'wormbase'}->get_wormbase_version);
    my $table_file = file($self->{'wormbase'}->wormpep . '/' . $self->{'wormbase'}->pepdir_prefix .
	'pep.table' . $self->{'wormbase'}->get_wormbase_version);

    for my $file ($wp_file, $table_file) {
	my $fh = $file->openr;
	my $seq_count = 0;
	my $uniprot_acc_count = 0;
	while ($fh->getline) {
	    next unless $_ =~ /^>/;
	    $seq_count++;
	    $uniprot_acc_count++ if $_ =~ /uniprot=/;
	}
	if ($uniprot_acc_count == 0) {
	    $self->{'log'}->log_and_die($file->basename . " does not contain Uniprot IDs\n");
	}
	else {
	    $self->{'log'}->write_to($file->basename .
				     " has Uniprot IDs for $uniprot_acc_count of $seq_count sequences\n");
	}
    }

    return;
}


=head2 vep_output_present()

    Function: checks for presence of .ace file generated by VEP pipeline
    Args:     n/a
    Returns:  n/a

=cut

sub vep_output_present {
    my $self = shift;
    
    my $vep_ace_file = $self->{'wormbase'}->acefiles . '/mapped_alleles.' .
	$self->{'wormbase'}->get_wormbase_version_name . '.ace';

    $self->_file_exists($vep_ace_file);

    return;
}


sub _expected_seq_region_count {
    my ($self, $seq_regions_in_db) = @_;
    
    my $seq_region_count = 0;
    my $gff_fh = $gff_file->openr;
    while ($gff_fh->getline) {
	$seq_region_count++ if $_ =~ /^##sequence-region/;
    }
    $gff_fh->close;

	
    if ($seq_region_count != $seq_regions_in_db) {
	$self->{'log'}->log_and_die("Was expecting $seq_regions_in_db sequence regions in " .
				    $gff_file->stringify . ", but found $seq_region_count\n");
    }

    return;
}


sub _file_exists {
    my ($self, $filepath) = @_;

    if (-e $filepath) {
	$self->{'log'}->write_to("$filepath exists\n");
    }
    else {
	$self->{'log'}->log_and_die("$filepath does not exist\n");
    }

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

    my $db = Ace->connect(-path => $self->{'wormbase'}->autoace,
			  -program => $self->{'wormbase'}->tace) or 
			      die 'Connection failure: ' . Ace->error;

    my $it = $db->fetch_many(-query => 'find Sequence_collection WHERE Live');
    my $collection_count = 0;
    my $nr_sequences;
    while (my $collection = $it->next){
	$collection_count++;
	my @sequences = $collection->at('Sequences');
	$nr_sequences = scalar @sequences;
    }

    $self->{'log'}->log_and_die("More than one Live Sequence_collection object\n") if $collection_count > 1;
    
    return $nr_sequences;
}


sub _nr_files_with_prefix_or_suffix {
    my ($self, $filenames, $prefix, $suffix) = @_;

    my $nr_files = 0;
    for my $filename (@$filenames) {
	next if defined $prefix and index($filename, $prefix) != 0;
	next if defined $suffix and index($filename, $suffix) != length $filename - length $suffix;
	$nr_files++;
    }

    return $nr_files;
}


sub _previous_release_non_elegans_genes {
    my $self = shift;

    my $db = Ace->connect(-path => $self->{'previous_wormbase'}->autoace,
			  -program => $self->{'previous_wormbase'}->tace) or
			      die 'Connection failure: ' . Ace->error; 

    my %non_elegans_genes;
    for my $species (@NON_ELEGANS_CORE_SPECIES) {
	my $it = $db->fetch_many(-query => 'FIND Gene WHERE Live AND Species = "' . $species . '"')
	    or die Ace->error;
	my $obj_count = 0;
	while (my $obj = $it->next) {
	    $obj_count++;
	    push @{$non_elegans_genes{$species}}, $obj->name if $obj_count % 2500 == 0;
	}
    }

    return \%non_elegans_genes;
}


1;
