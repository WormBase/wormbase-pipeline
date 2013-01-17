#!/usr/local/bin/perl5.8.0 -w
#
# script_template.pl
#
# by Keith Bradnam
#
# This is a example of a good script template
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-01-17 09:33:47 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

######################################
# variables and command-line options #
######################################

my ( $help, $debug, $test, $verbose, $store, $wormbase, $species, $gff3, $infile, $outfile );

GetOptions(
    'help'      => \$help,
    'debug=s'   => \$debug,
    'test'      => \$test,
    'verbose'   => \$verbose,
    'store:s'   => \$store,
    'species:s' => \$species,
    'infile:s'  => \$infile,
    'outfile:s' => \$outfile,
    'gff3'      => \$gff3,
);

if ($store) {
    $wormbase = retrieve($store) or croak("Can't restore wormbase from $store\n");
} 
else {
    $wormbase = Wormbase->new(
        -debug => $debug,
        -test  => $test,
	-organism => $species,
    );
}

# Display help if required
&usage('Help') if ($help);

# in test mode?
print "In test mode\n" if ( $verbose && $test );

# establish log file.
my $log = Log_files->make_build_log($wormbase);

if ($infile and not $outfile or
    $outfile and not $infile) {
  $log->log_and_die("When using -infile or -outfile, you must use both\n");
}

#################################
# Set up some useful paths      #
#################################

# Set up top level base directories (these are different if in test mode)
my $gff_dir = $wormbase->gff;    # GFF
my $ace_dir = $wormbase->autoace;# Autoace

# other paths
my $tace = $wormbase->tace;      # TACE PATH

###################################
# get the species of the sequences
###################################

my $cmd1 = "Query Find Sequence Where Database = \"NEMATODE_NET\"\nshow -a Species\nquit";
my $cmd2 = "Query Find Sequence Where Database = \"NEMBASE\"\nshow -a Species\nquit";
my $cmd3 = "Query Find Sequence Where Database = \"EMBL\"\nshow -a Species\nquit";

my %species;
my ( $id, $db );

print "Finding BLAT_WASHU data\n";
open( TACE, "echo '$cmd1' | $tace $ace_dir |" );
while (<TACE>) {
    chomp;
    next if (/acedb\>/);
    next if (/\/\//);
    if (/Sequence\s+:\s+\"(\S+)\"/) {
        $id = $1;
    }
    elsif (/Species\s+\"(.+)\"/) {
        $species{'BLAT_WASHU'}->{$id} = $1;
    }
}
close TACE;

print "Finding BLAT_NEMBASE data\n";
open( TACE, "echo '$cmd2' | $tace $ace_dir |" );
while (<TACE>) {
    chomp;
    next if (/acedb\>/);
    next if (/\/\//);
    if (/Sequence\s+:\s+\"(\S+)\"/) {
        $id = $1;
    }
    elsif (/Species\s+\"(.+)\"/) {
        $species{'BLAT_NEMBASE'}->{$id} = $1;
    }
}
close TACE;

print "Finding BLAT_NEMATODE data\n";
open( TACE, "echo '$cmd3' | $tace $ace_dir |" );
while (<TACE>) {
    chomp;
    next if (/acedb\>/);
    next if (/\/\//);
    if (/Sequence\s+:\s+\"(\S+)\"/) {
        $id = $1;
    }
    elsif (/Species\s+\"(.+)\"/) {
        $species{'BLAT_NEMATODE'}->{$id} = $1;
    }
}
close TACE;

# Not all EMBL Sequences will be in the database when this script
# is run, because the objects for the other single-species sets
# (briggase, remanei etc) do not make their way into autoace
# until the merge. We therefore fill in the gaps by picking these
# up directly. 
my %accessors = $wormbase->all_species_accessors;
foreach my $owb (values %accessors) {
  my $cdnadir = $owb->maskedcdna;
  my $full_sp = $owb->full_name;

  foreach my $file (glob("$cdnadir/*.*")) {
    open(my $fh, $file);
    while(<$fh>) {
      /^\>(\S+)/ and do {
        if (not exists $species{BLAT_NEMATODE}->{$1}) {
          $species{BLAT_NEMATODE}->{$1} = $full_sp;
        }
      }
    }
  }
}



##########################
# MAIN BODY OF SCRIPT
##########################
my (@gff_files, $count, $gffout_fh, $gffin_fh, $gff_tmp_file);

# loop through the chromosomes
if ($infile) {
  @gff_files = ($infile);
  open($gffout_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");
} elsif ($wormbase->assembly_type eq 'contig') {
  @gff_files = $wormbase->GFF_file_name;
} else {
  my @chromosomes = $wormbase->get_chromosome_names( -mito => 1, -prefix => 1 );
  @gff_files = map { $wormbase->GFF_file_name($_) } @chromosomes;
}

foreach my $GFF_file_name (@gff_files) {
  print "Reading $GFF_file_name\n" if ($verbose);
  
  open($gffin_fh, "<$GFF_file_name") or $log->log_and_die("Can't open $GFF_file_name\n"); 
  if (not $outfile) {
    $gff_tmp_file = "${GFF_file_name}.new";
    open($gffout_fh, ">$gff_tmp_file") or $log->log_and_die("Failed to open gff file $gff_tmp_file\n");
  }  

  while ( my $line = <$gffin_fh> ) {
    chomp $line;
    if ( $line =~ /^#/ || $line !~ /\S/ ) {
      print $gffout_fh "$line\n";
      next;
    }
    my @f = split /\t/, $line;
    
    # It is possible that this script is being run on the same
    # input file multiple times (e.g. when sorting out problems)
    # in which case we do not want to add 'Species' multiple
    # times to the same line.
    if (defined $f[1] && $f[8] !~ /;\sSpecies/ && $f[8] !~ /Species\=/) {
      
      if (grep { $f[1] =~ /^$_/ } ('BLAT_WASHU', 'BLAT_NEMBASE', 'BLAT_NEMATODE', 'BLAT_Caen_EST_')) {
        my $id;
        
        if ($gff3) {
          ($id) = $f[8] =~ /Target\=(\S+)/;
        } else {
          ($id) = ( $f[8] =~ /Target \"Sequence:(\S+)\"/ );
        }
        
        my $hkey = $f[1];
        # the {'BLAT_NEMATODE'} hash holds the EMBL data which BLAT_Caen_EST_* uses as well
        if ($hkey =~ /BLAT_Caen_/) {
          $hkey = "BLAT_NEMATODE";
        }
        
        if ( exists $species{$hkey}->{$id} ) {
          my $suffix;
          if ($gff3) { 
            $suffix = ";Species=" . $species{$hkey}->{$id};
          } else {
            $suffix = " ; Species \"" . $species{$hkey}->{$id} . "\"";
          }           
 
          $line = $line . $suffix;
          $count++;
          print "$line\n" if ($verbose);
          
        } else {
          print "Cannot find species info for $id\n";
        }
      }
    }
    
    # write out the line
    print $gffout_fh "$line\n";
    
    # end of GFF loop
  }
  
  if (not $outfile) {
    close($gffout_fh);
    # copy new GFF files over
    system("mv -f $gff_tmp_file $GFF_file_name"); 
  }
}

if ($outfile and $gffout_fh) {
  close($gffout_fh);
}

# Close log files and exit
$log->write_to("\n\nStatistics\n");
$log->write_to("----------\n\n");
$log->write_to("Changed $count lines\n");

##################
# Check the files
##################
$wormbase->check_files($log) unless $infile; 

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);

##############################################################
#
# Subroutines
#
##############################################################

##########################################

sub usage {
  my $error = shift;
  
  if ( $error eq "Help" ) {
    
    # Normal help menu
    system( 'perldoc', $0 );
    exit(0);
  }
}

##########################################

# Add perl documentation in POD format
# This should expand on your brief description above and
# add details of any options that can be used with the program.
# Such documentation can be viewed using the perldoc command.

__END__

=pod

=head2 NAME - add_species_to_BLAT_GFF.pl

=head1 USAGE

=over 4

=item  add_species_to_BLAT_GFF.pl [-options]

=back

Todd wanted to have species names added to the GFF results for non-elegans EST BLAT hits.
This script gets the species for each EST sequence from the autoace database and adds this information to the GFF records.

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

=item Gary Williasm (gw3@sanger.ac.uk)

=back

=cut
