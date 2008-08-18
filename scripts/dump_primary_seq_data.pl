#!/software/bin/perl -w
#
# dump_primary_seq_data.pl
#
# Usage : dump_primary_seq_data.pl [-options]
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2008-08-18 13:28:45 $

my $script_dir = $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'};
use strict;
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Modules::Features;
use Species;

######################################
# variables and command-line options #
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $version, $organism, $output, $species);

GetOptions ("help"       => \$help, #
            "debug=s"    => \$debug, # debug option, turns on more printing and only email specified user.
            "test"       => \$test, # only genomic RNAs will be used for a quick test.
            "verbose"    => \$verbose, # additional printing
            "store:s"    => \$store, # wormbase storable object
	    "organism=s" => \$organism, # Specify an organism
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
			     -organism => $species,
			   );
}

#################################
# Set up some useful stuff      #
#################################
# establish log file.
my $log = Log_files->make_build_log($wormbase);
my $tace = $wormbase->tace; # TACE PATH

# main stuff goes here
my @species;
my $acefile;
my $Longtextfile;



##########################
# MAIN BODY OF SCRIPT
##########################
my %species = ($wormbase->species_accessors);
$species{$wormbase->species} = $wormbase;

if (!defined($organism)) {
@species = (keys %species);
}
elsif ($organism) {
push (@species, $organism);
}


#foreach my $organism(keys %species) {
foreach my $organism(@species) {
  $log->write_to("============================================\nDUMPING: $organism\n============================================\n");
# output dirs
  my   $output_dir = $wormbase->basedir."_DATA/cDNAace/$organism";
  $wormbase->run_command ("mkdir $output_dir", $log) if (!-e $output_dir);
  

  my $sourceDB;
  if ($organism eq "elegans") {
    $sourceDB = $wormbase->database('camace');
  }
  elsif ($organism eq "heterorhabditis") {
    $sourceDB = glob "~wormpub/DATABASES/heterorhabditis";
  }
  elsif ($organism eq "pristionchus"){
    $sourceDB = glob "~wormpub/DATABASES/pristionchus";
  }
  else {$sourceDB = $wormbase->database($organism);}
  # Fetch sequence data from primary database.
  $log->write_to("Fetching sequence data from $sourceDB:\n");
  

# not needed any more - Dump the data back to BUILD_DATA/cDNA/organism/
$log->write_to("Dumping BLAT sequences($organism)\n");
&dump_BLAT_data ($sourceDB, $organism, \%species, $output_dir);
  
  $log->write_to("============================================\nProcessed: $organism\n============================================\n\n");
}

$log->write_to("\nRefreshed........ ;) \n\n");
$log->write_to ("WORK DONE------------FINISHED\n\n");
#Close log/files and exit
$log->mail();
exit(0);




###################################################
#                 SUBROUTINES                     #
###################################################


sub dump_BLAT_data {
  my $dbdir = shift;
  my $subspecies = shift;
  my $speciesobj = shift;
  my $EST_dir =shift;
#  $EST_dir = $EST_dir.$subspecies;
    $log->write_to("Dumping $subspecies from $dbdir\n\n");

  # Remove stale data if it exists on disk.
  my @types = ('mRNA','ncRNA','EST','OST','tc1','RST');
  foreach my $type (@types) {
    $wormbase->run_command ("rm $EST_dir${type}", $log) if (-e $EST_dir."${type}");
    $log->write_to("Removed $EST_dir${type}\n\n")  if (-e $EST_dir."${type}" && $debug);
  }

  my $command=<<END;
query find Sequence where method = NDB & RNA AND NEXT = mRNA & !Ignore\n
Write $EST_dir/mRNA.ace\n
clear\n
query find Sequence where method = NDB & RNA AND NEXT != mRNA & !Ignore\n
Write $EST_dir/ncRNA.ace\n
clear\n
query find Sequence where method = EST_$subspecies & !OST* & !Ignore & !RST*\n
Write $EST_dir/EST.ace\n
clear\n
query find Sequence where method = EST_$subspecies & OST* & !Ignore\n
Write $EST_dir/OST.ace\n
clear\n
query find Sequence where method = EST_$subspecies & RST* & !Ignore\n
Write $EST_dir/RST.ace\n
clear\n
query find Sequence TC*\n
Write $EST_dir/tc1.ace\n
clear\n
query find longtext\n
Write $EST_dir/LongText.ace
quit\n
END
  print $command if ($debug);
  open (DB, "| $tace $dbdir") || die "Couldn't open $dbdir\n";
  print DB $command;
  close DB;
  $log->write_to("Finished $subspecies\n\n");
}


__END__

=pod

=head2 NAME - dump_primary_seq_data.pl

=head1 USAGE

=over 4

=item dump_primary_seq_data.pl [-options]

=back

This script dumps ace files for Transcript data for each species and stores them under BUILD_DATA/cDNAace/

dump_primary_seq_data.pl

MANDATORY arguments:

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

=item -test, Test mode.

=back

=over 4

=item -species, only performs dumping for that species only

=back

=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Paul Davis (pad@sanger.ac.uk)

=back

=cut
