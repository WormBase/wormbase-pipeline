use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
use Modules::Curate;
use FileHandle;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species, $database, $curate_test);
my ($replace, $isoform, $newgene, $split, $merge, $to_cds, $to_pseudogene, $to_transcript, $make_transposon, $make_operon, $delete, $last_reviewed, $just_make_history, $check_gene_name, $show_next_sequence_name);
my ($noload, $force, $class, $oldseq, @newseq);


GetOptions ("help"              => \$help,
            "debug=s"           => \$debug,
            "test"              => \$test,
            "verbose"           => \$verbose,
            "store:s"           => \$store,
	    "species:s"         => \$species,           # this defaults to elegans
	    "database:s"        => \$database,          # the default curation database for the species is used if this is not specied
	    "curate_test"       => \$curate_test,       # test out the Curate package routines

	    "replace"           => \$replace,           # COMMAND: replace old structure by new structure
	    "isoform"           => \$isoform,           # COMMAND: add a new isoform to an existing sequence name
	    "create_gene"       => \$newgene,           # COMMAND: create a new gene
	    "split"             => \$split,             # COMMAND: split an existing locus to make a new gene(s)
	    "merge"             => \$merge,             # COMMAND: merge existing loci
	    "to_cds"            => \$to_cds,            # COMMAND: convert an object into a CDS class object
	    "to_pseudogene"     => \$to_pseudogene,     # COMMAND: convert an object into a Pseudogene class object
	    "to_transcript"     => \$to_transcript,     # COMMAND: convert an object into a Transcript class object (you have to edit it if you want a Type other than ncRNA)
	    "make_transposon"   => \$make_transposon,   # COMMAND: convert an object into a Transposon_CDS and new Transposon class object
	    "make_operon"       => \$make_operon,       # COMMAND: make an operon and add some genes to it
	    "delete"            => \$delete,            # COMMAND: delete an object (leaving a history object)
	    "last_reviewed"     => \$last_reviewed,     # COMMAND: add a Last_reviewed tag
	    "just_make_history"      => \$just_make_history,      # COMMAND: just create a history object
	    "check_gene_name"   => \$check_gene_name,   # COMMAND: report any Gene_name for the object
	    "show_next_sequence_name"   => \$show_next_sequence_name,   # COMMAND: show the next available Sequence-name on this clone or scaffold

	    "class:s"           => \$class,             # class of things being curated - defaults to CDS
	    "oldseq:s"          => \$oldseq,            # existing sequence name
	    "newseq:s"          => \@newseq,            # new sequence name(s)

	    "noload"            => \$noload,            # if this is set, do not load the ACE file
	    "force"             => \$force,             # force the command to succeed even if it would normally fail a warning or a sanity check
            ) || die(@!);


if (!defined $database) {
  if ($species eq 'elegans') {
    $database = "/nfs/panda/ensemblgenomes/wormbase/camace_$ENV{USER}";
  } else {
    $database = "/nfs/panda/ensemblgenomes/wormbase/${species}_curation";
  }
}
if (!-e $database || !-d $database) {die "Can't find database $database\n"}


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species,
			     -autoace  => $database,
                             );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);
}

# allow comma-separated lists of sequence names
@newseq = split(/,/,join(',',@newseq));

# establish log file.
my $log = Log_files->make_build_log($wormbase);

if (!defined $species) {$species = $wormbase->species}

my $outfile = "$database/$ENV{USER}_curate.ace";
## check to see if outfile exists already in which case it didn't get read in successfully in a previous session
#if (-e $outfile) {die "The output ACE file $outfile already exists\nWas there a problem parsing this existing file last time?\n";}
 
my $ace = Ace->connect (-path => $database) || die "cannot connect to database at $database\n";

my $out = FileHandle->new("> $outfile");
my $curate = Curate->new($out, $ace, $wormbase, $log);

if (defined $force) {$curate->set_force(1)}

if ($curate_test) {
  curate_test();

} elsif ($replace) {
  if (!defined $class) {$class = 'CDS'}
  $curate->replace_cmd($class, $oldseq, @newseq);
} elsif ($isoform) {
  if (!defined $class) {$class = 'CDS'}
  $curate->isoform_cmd($class, $oldseq, @newseq);
} elsif ($newgene) {
  if (!defined $class) {$class = 'CDS'}
  $curate->newgene_cmd($class, @newseq);
} elsif ($split) {
  if (!defined $class) {$class = 'CDS'}
  $curate->split_cmd($class, $oldseq, @newseq);
} elsif ($merge) {
  # allow comma-separated lists of existing sequence names to be merged
  my @oldseq = split(/,/,$oldseq);
  if (!defined $class) {$class = 'CDS'}
  $curate->merge_cmd($class, \@oldseq, @newseq);
} elsif ($to_cds) {
  if (!defined $class) {$class = 'Pseudogene'}
  $curate->to_cds_cmd($class, $oldseq, @newseq);
} elsif ($to_pseudogene) {
  if (!defined $class) {$class = 'CDS'}
  $curate->to_pseudogene_cmd($class, $oldseq, @newseq);
} elsif ($to_transcript) {
  if (!defined $class) {$class = 'CDS'}
  $curate->to_transcript_cmd($class, $oldseq, @newseq);
} elsif ($make_transposon) {
  if (!defined $class) {$class = 'CDS'}
  $curate->make_transposon_cmd($class, $oldseq, @newseq);
} elsif ($make_operon) {
  my @oldseq = split(/,/,$oldseq);
  $curate->make_operon_cmd(\@oldseq, @newseq);
  print "\nNow load $outfile into geneace\n\n";
} elsif ($delete) {
  if (!defined $class) {$class = 'CDS'}
  $curate->delete_cmd($class, $oldseq);
} elsif ($last_reviewed) {
  if (!defined $class) {$class = 'CDS'}
  $curate->last_reviewed_cmd($class, $oldseq);
} elsif ($just_make_history) {
  if (!defined $class) {$class = 'CDS'}
  $curate->just_make_history_cmd($class, $oldseq);
} elsif ($check_gene_name) {
  if (!defined $class) {$class = 'CDS'}
  $curate->check_gene_name_cmd($class, $oldseq);
} elsif ($show_next_sequence_name) {
  $curate->show_next_sequence_name_cmd($oldseq);
} else {
  die "command not recognised\n";
}

$out->close;

$curate->check_force_flag();

if (!$noload) {
  print "Loading ACE file ...\n";
  if ($curate->load_ace($outfile)) {print "ERROR The acefile $outfile was not parsed correctly\n";}
}

$ace->close;
$log->mail();
print "Finished.\n" if ($verbose);
exit(0);

##############################################################
#
# Subroutines
#
##############################################################

#################
# Test the routines
sub curate_test {
  if ($species eq 'elegans') {
    my $cgc = $curate->cgcname('CDS', 'F19G12.7');
    if ($cgc eq 'abu-2') {print "OK - cgcname()\n";} else {print "FAIL - cgcname() abu-2\n"}
    my $cgc = $curate->cgcname('CDS', 'AC3.12');
    if ($cgc eq '') {print "OK - cgcname()\n";} else {print "FAIL - cgcname() null\n"}

    my $gene = $curate->SeqName2Gene('AC3.12');
    if ($gene eq 'WBGene00077503') {print "OK - SeqName2Gene()\n"} else {print "FAIL - SeqName2Gene()\n"}
    my $seqname = $curate->Gene2SeqName('WBGene00077503');
    if ($seqname eq 'AC3.12') {print "OK - Gene2SeqName()\n"} else {print "FAIL - Gene2SeqName()\n"}

    my $next_seq_name = $curate->Next_CDS_ID('AC3');
    if ($next_seq_name eq 'AC3.17') {print "OK - Next_CDS_ID()\n"} else {print "FAIL - Next_CDS_ID()\n"}
    my $next_seq_name = $curate->Next_CDS_ID('AC3');
    if ($next_seq_name eq 'AC3.18') {print "OK - Next_CDS_ID()\n"} else {print "FAIL - Next_CDS_ID()\n"}


  } elsif ($species eq 'briggsae') {
    my $cgc = $curate->cgcname('CDS', 'CBG21774');
    if ($cgc eq 'Cbr-fem-3') {print "OK - cgcname()\n";} else {print "FAIL - cgcname() Cbr-fem-3\n"}
    my $cgc = $curate->cgcname('CDS', 'CBG22653');
    if ($cgc eq '') {print "OK - cgcname()\n";} else {print "FAIL - cgcname() null\n"}

    my $gene = $curate->SeqName2Gene('CBG16519');
    if ($gene eq 'WBGene00036432') {print "OK - SeqName2Gene()\n"} else {print "FAIL - SeqName2Gene() - this currently fails for non-elegans species because Gene's do not have Sequence_name data\n"}
    my $seqname = $curate->Gene2SeqName('WBGene00036432');
    if ($seqname eq 'CBG16519') {print "OK - Gene2SeqName()\n"} else {print "FAIL - Gene2SeqName() - this currently fails for non-elegans species because Gene's do not have Sequence_name data\n"}			    

    my $next_seq_name = $curate->Next_CDS_ID();
    if ($next_seq_name eq 'CBG30866') {print "OK - Next_CDS_ID()\n"} else {print "FAIL - Next_CDS_ID()\n"}
    my $next_seq_name = $curate->Next_CDS_ID();
    if ($next_seq_name eq 'CBG30867') {print "OK - Next_CDS_ID()\n"} else {print "FAIL - Next_CDS_ID()\n"}
  } else {
    die "Don't have tests for species: $species\n";
  }


}
##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################




# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - script_template.pl

=head1 USAGE

=over 4

=item script_template.pl  [-options]

=back

This script does...blah blah blah

script_template.pl MANDATORY arguments:

=over 4

=item None at present.

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item xxx (xxx@sanger.ac.uk)

=back

=cut
