#!/usr/local/bin/perl5.8.0 -w          
#
# map_features.pl
#
# by Dan Lawson
#
# This maps features to the genome based on their flanking sequence.
# Uses Ant's Feature_mapper.pm module
#
#
# Last updated by: $Author: gw3 $                      # These lines will get filled in by cvs and helps us
# Last updated on: $Date: 2009-12-04 10:46:17 $        # quickly see when script was last changed and by whom


$|=1;
use strict;
use lib $ENV{'CVS_DIR'};
use Feature_mapper;
use Wormbase;
use Ace;
use Getopt::Long;


my ($feature, $clone, $flanking_left, $flanking_right, $coords, $span,$store);

my $help;                    # Help menu
my $debug;                   # Debug mode 
my $verbose;                 # Verbose mode
my $all;                     # Do all the following features:
my $SL1;                     #  SL1 trans-splice leader acceptors
my $SL2;                     #  SL2 trans-splice leader acceptors
my $polyA_site;              #  polyA_site
my $polyA_signal;            #  polyA_signal
my $binding_site;            #  binding_site feature data.
my $binding_site_reg;        #  binding_site_region feature data
my $segmental_duplication;          #  segmental_duplication features
my $genome_sequence_error;   # genome_sequence_error
my $transcription_start_site;# transcription_start_site
my $transcription_end_site;  # transcription_end_site
my $promoter;                # promoter region
my $regulatory_region;       # regulatory region
my $adhoc;                   # Run against a file, output to screen
my $start;
my $stop;
my $test;

GetOptions (
	    "all"                   => \$all,
	    "SL1"                   => \$SL1,
	    "SL2"                   => \$SL2,
	    "polyA_site"            => \$polyA_site,
	    "polyA_signal"          => \$polyA_signal,
	    "binding_site"          => \$binding_site,
	    "binding_site_reg"      => \$binding_site_reg,
	    "segmental_duplication" => \$segmental_duplication,
	    "Genome_sequence_error" => \$genome_sequence_error,
	    "transcription_start_site"=> \$transcription_start_site,
	    "transcription_end_site"=> \$transcription_end_site,
	    "promoter"              => \$promoter,
	    "regulatory_region"     => \$regulatory_region,
	    "adhoc=s"               => \$adhoc,
            "debug=s"               => \$debug,
            "verbose"               => \$verbose,
	    "help"                  => \$help,
    	    'store=s'               => \$store,
	    'test'                  => \$test
		);

# Help pod if needed
exec ('perldoc',$0) if ($help);

# recreate configuration  
my $wb;
if ($store) { $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wb = Wormbase->new( -debug => $debug, 
			    -test  => $test 
			    ) }

my $log = Log_files->make_build_log($wb);

#######################
# ACEDB and databases #
#######################

my $tace   = $wb->tace;
my $dbdir  = $wb->autoace;
my $outdir = $wb->acefiles;
$log->write_to("writing to ".$wb->acefiles."\n\n");

# WS version for output files
our ($WS_version) = $wb->get_wormbase_version_name;

# coordinates for Feature_mapper.pm module

my $mapper      = Feature_mapper->new($dbdir,undef, $wb);

# sanity checks for the length of feature types
my %sanity = (
	      'SL1'          => 0,
	      'SL2'          => 0,
	      'polyA_site'   => 0,
	      'polyA_signal_sequence' => 6,
	      'binding_site' => -1,
	      'binding_site_region' => -1,
	      'segmental_duplication' => -1,
	      'Genome_sequence_error' => -1,
	      'transcription_end_site' => -1,
	      'transcription_start_site' => -1,
	      'promoter' => -1,
	      'regulatory_region' => -1,
	      );

# queue which Feature types you want to map
my @features2map;
push (@features2map, "SL1")           if (($SL1) || ($all));
push (@features2map, "SL2")           if (($SL2) || ($all));
push (@features2map, "polyA_signal_sequence")  if (($polyA_signal) || ($all));
push (@features2map, "polyA_site")    if (($polyA_site) || ($all));
push (@features2map, "binding_site")  if (($binding_site) || ($all));
push (@features2map, "binding_site_region")  if (($binding_site_reg) || ($all));
push (@features2map, "segmental_duplication")  if (($segmental_duplication) || ($all));
push (@features2map, "Genome_sequence_error")  if (($genome_sequence_error) || ($all));
push (@features2map, "transcription_end_site")  if (($transcription_end_site) || ($all));
push (@features2map, "transcription_start_site")  if (($transcription_start_site) || ($all));
push (@features2map, "promoter")  if (($promoter) || ($all));
push (@features2map, "regulatory_region")  if (($regulatory_region) || ($all));

#############
# main loop #
#############

foreach my $query (@features2map) {

  my $table_file = "/tmp/map_features_table_$query.def";
  $log->write_to("Mapping $query features\n");

  # open output files
  open (OUTPUT, ">$outdir/feature_${query}.ace") or die "Failed to open output file\n" unless ($adhoc);

  # start tace session for input data (or find file for adhoc run)
  if ($adhoc) {
    open (TACE, "<$adhoc") or die "Failed to open input file: $adhoc\n";
    print "// Opening a file for input: $adhoc\n" if ($verbose);
  }
  else {
    my $species_name = $wb->full_name;
    my $table = <<EOF;
Sortcolumn 1

Colonne 1
Width 12
Optional
Visible
Class
Class Feature
From 1
Condition Method = "$query" AND Species = "$species_name"

Colonne 2
Width 12
Optional
Visible
Class
Class Sequence
From 1
Tag Flanking_sequences

Colonne 3
Width 32
Optional
Visible
Text
Right_of 2
Tag  HERE

Colonne 4
Width 32
Optional
Visible
Text
Right_of 3
Tag  HERE

// End of these definitions
EOF
    open (TABLE, ">$table_file") || die "Can't open $table_file: $!";
    print TABLE $table;
    close(TABLE);

    open (TACE, "echo 'Table-maker -p $table_file\nquit\n' | $tace $dbdir | ");
  }
  while (<TACE>) {
    
    # when it finds a good line
    if (/^\"(\S+)\"\s+\"(\S+)\"\s+\"(\S+)\"\s+\"(\S+)\"/) {
      ($feature,$clone,$flanking_left,$flanking_right) = ($1,$2,$3,$4);
      print "NEXT FEATURE: $feature,$clone,$flanking_left,$flanking_right\n" if ($debug);
      print "Unintialised feature in line $_" if (!defined $feature);
      print "Unintialised $feature clone" if (!defined $clone);
      print "Unintialised $feature flanking_left" if (!defined $flanking_left);
      print "Unintialised $feature flanking_right" if (!defined $flanking_right);

      my @coords = $mapper->map_feature($clone,$flanking_left,$flanking_right);
      if (!defined $coords[2]) {
	$log->write_to("ERROR: Can't map feature $feature on clone $clone flanking sequences: $flanking_left $flanking_right\n");
	$log->error;
	next;
      }
      
      $log->write_to("Feature $feature maps to different clone than suggested $clone -> $coords[0]\n") if ($clone ne $coords[0]);
      $clone = $coords[0];
      $start = $coords[1];
      $stop  = $coords[2];
      
      # Deal with polyA_signal_sequence features which have the
      # flanking sequences outside of the feature, but
      # Feature_mapper::map_feature expects the flanking sequences to
      # overlap the feature by one base to enable correct mapping of
      # features less than 2bp

      # munge returned coordinates to get the span of the mapped feature
      
      if ($query eq "polyA_signal_sequence") {
	if ($start < $stop) {
	  $start++;
	  $stop--;
	  $span = $stop - $start + 1;
	}
	else {
	  $start--;
	  $stop++;
	  $span = $start - $stop + 1;
	}
      }
      # else deal with butt-ended features (e.g. SL1, SL2 & polyA_site)
      elsif ($start > $stop) {
	$span = $start - $stop - 1;
      }
      else {
	$span = $stop - $start - 1;
      }
      
      # check feature span is sane
      if ( exists $sanity{$query} && (($span == $sanity{$query}) || ($sanity{$query} < 0)) ) {
	
	if ($adhoc) {
	  print "$feature maps to $clone $start -> $stop, feature span is $span bp\n";
	}
	else {
	  print OUTPUT "//$feature maps to $clone $start -> $stop, feature span is $span bp\n";
	  print OUTPUT "\nSequence : \"$clone\"\n";
	  print OUTPUT "Feature_object $feature $start  $stop\n\n";
	}
      }
      else {
	$log->write_to("ERROR: $feature maps to $clone $start -> $stop, feature span is $span bp\n");
	$log->error;
      }
    } #_ if match line

    # lines that look like features but there is a problem eg. whitespace in flanks.
    elsif (/^\"(\S+)\"/) {
      $log->write_to("ERROR: $1 has a problem, please check flanking sequences!! (whitespace is one cause)\n");
      $log->error;
    }
  }				 
  close TACE;	
  unlink $table_file if (-e $table_file);

  $wb->load_to_database($wb->autoace, "$outdir/feature_${query}.ace", "feature_mapping", $log);
}


close OUTPUT unless ($adhoc);


###############
# hasta luego #
###############

$log->mail();
exit(0);

__END__

=pod

=head2 NAME - map_features.pl

=head1 USAGE

=over 4

=item map_features.pl [-options]

=back

map_features.pl mandatory arguments:

=over 4

=item none

=back

map_features.pl optional arguments:

=over 4

=item -all 

map all of the following feature types:

=item -SL1

map SL1 trans-splice leader acceptor sites (2 bp feature)

=item -SL2

map SL2 trans-splice leader acceptor sites (2 bp feature)

=item -polyA_site

map polyA attachement sites (2 bp feature)

=item -polyA_signal

map polyA signal sequence sites (6 bp feature)

=item -debug <user> 

Queries current_DB rather than autoace

=item -build

Assumes you are building, so therefore loads data into autoace and writes acefiles to
autoace/acefiles
=item -store <storefile>

specifies an storable file with saved options

=item -adhoc <file> 

Queries a flatfile for the Feature data. Flatfile format is
"<feature_name>"  "<clone>" "<flanking_sequence_left>" "<flanking_sequence_right>" 

=item -help      

This help page

=cut
