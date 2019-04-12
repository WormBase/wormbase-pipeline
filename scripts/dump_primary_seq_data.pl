#!/software/bin/perl -w
#
# dump_primary_seq_data.pl
#
# Usage : dump_primary_seq_data.pl [-options]
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2008-10-08 10:52:20 $

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

my ($help, $debug, $test, $verbose, $store, $wormbase, $version, $organism, $species, $noseqdump, $noacedump);

GetOptions ("help"       => \$help, #
            "debug=s"    => \$debug, # debug option, turns on more printing and only email specified user.
            "test"       => \$test, # only genomic RNAs will be used for a quick test.
            "verbose"    => \$verbose, # additional printing
            "store:s"    => \$store, # wormbase storable object
	    "organism=s" => \$organism, # Specify an organism
	    "noace"      => \$noacedump, # don't update the ace data on disk.
            "noseq"      => \$noseqdump, # don't update the BLAT seq data on disk.
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
			     -organism => $organism,
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
  my   $ACEoutput_dir = $wormbase->basedir."_DATA/cDNAace/$organism";
  $wormbase->run_command ("rm -r $ACEoutput_dir", $log) if (-e $ACEoutput_dir);
  $log->write_to("Removing $ACEoutput_dir\n") if (-e $ACEoutput_dir);
  $wormbase->run_command ("mkdir $ACEoutput_dir", $log) if (!-e $ACEoutput_dir);
  $log->write_to("Creating $ACEoutput_dir\n") if (!-e $ACEoutput_dir);
  my $BLAToutput_dir = $wormbase->basedir."_DATA/cDNA/$organism";
  $wormbase->run_command ("rm -r $BLAToutput_dir", $log) if (-e $BLAToutput_dir);
  $log->write_to("Removing $BLAToutput_dir\n") if (-e $BLAToutput_dir);
  $wormbase->run_command ("mkdir $BLAToutput_dir", $log) if (!-e $BLAToutput_dir);
  $log->write_to("Creating $BLAToutput_dir\n") if (!-e $BLAToutput_dir);

  my $sourceDB = $wormbase->seq_db;
  
  # Fetch sequence data from primary database.
  $log->write_to("\n\nUsing: $sourceDB: for data retrieval\n");
  

# not needed any more - Dump the data back to BUILD_DATA/cDNA/organism/
$log->write_to("Dumping $organism Transcript .ace data.........\n") unless $noacedump;
&dump_BLAT_ace_data ($sourceDB, $organism, \%species, $ACEoutput_dir) unless $noacedump;
$log->write_to("Dumping $organism BLAT dna sequences...........\n") unless $noseqdump;
&dump_BLAT_data ($sourceDB, $organism, \%species, $BLAToutput_dir) unless $noseqdump; 
  $log->write_to("============================================\nPROCESSED: $organism\n============================================\n-\n-\n-\n-\n");
}

$log->write_to("\nRefreshed........ ;) \n\n");
#Close log/files and exit
$log->mail();
exit(0);




###################################################
#                 SUBROUTINES                     #
###################################################


sub dump_BLAT_ace_data {
  my $dbdir = shift;
  my $subspecies = shift;
  my $speciesobj = shift;
  my $EST_dir =shift;
#  $EST_dir = $EST_dir.$subspecies;
 
  # Remove stale data if it exists on disk.
  my @types = ('mRNA','ncRNA','EST','OST','tc1','RST','Trinity','IsoSeq','Nanopore');
  foreach my $type (@types) {
    my $file = "$EST_dir/${type}.ace";
    if (-e $file) {
      $wormbase->run_command ("rm $file", $log);
      $log->write_to("Removed $file\n\n");
    }
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
query find Sequence where method = RNASeq_trinity & !Ignore\n
Write $EST_dir/Trinity.ace\n
clear\n
query find Sequence where method = RNASeq_IsoSeq & !Ignore\n
Write $EST_dir/IsoSeq.ace\n
clear\n
query find Sequence where method = RNASeq_Nanopore & !Ignore\n
Write $EST_dir/Nanopore.ace\n
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

  # special case for Trinity
  if (not -e "$EST_dir/Trinity.ace" || not -e "$EST_dir/IsoSeq.ace") {
    $wormbase->run_command("touch $EST_dir/Trinity.ace", $log);
  }
  $log->write_to("$subspecies Transcripts dumped\n\n");
}

##########################
# dump data for BLATing  #
##########################
sub dump_BLAT_data {
  my $dbdir = shift;
  my $subspecies = shift;
  my $speciesobj = shift;
  my $EST_dir =shift;
#  $EST_dir = $EST_dir.$subspecies;

  # Remove stale data if it exists on disk.
  my @types = ('mRNA','ncRNA','EST','OST','tc1','RST','Trinity','IsoSeq','Nanopore');
  foreach my $type (@types) {
    if (-e "$EST_dir/$type") {
      $wormbase->run_command ("rm $EST_dir/$type", $log);
      $log->write_to("Removed $EST_dir/$type\n\n");
    }
  }

  my $command=<<END;
query find Sequence where method = NDB & RNA AND NEXT = mRNA & !Ignore\n
Dna -mismatch $EST_dir/mRNA\n
clear\n
query find Sequence where method = NDB & RNA AND NEXT != mRNA & !Ignore\n
Dna -mismatch $EST_dir/ncRNA\n
clear\n
query find Sequence where method = EST_$subspecies & !OST* & !Ignore & !RST*\n
Dna -mismatch $EST_dir/EST\n
clear\n
query find Sequence where method = EST_$subspecies & OST* & !Ignore\n
Dna -mismatch $EST_dir/OST\n
clear\n
query find Sequence where method = EST_$subspecies & RST* & !Ignore\n
Dna -mismatch $EST_dir/RST\n
clear\n
query find Sequence where method = RNASeq_trinity & !Ignore\n
Dna -mismatch $EST_dir/Trinity\n
clear\n
query find Sequence where method = RNASeq_IsoSeq & !Ignore\n
Dna -mismatch $EST_dir/IsoSeq\n
clear\n
query find Sequence where method = RNASeq_Nanopore & !Ignore\n
Dna -mismatch $EST_dir/Nanopore\n
clear\n
query find Sequence TC*\n
Dna -mismatch $EST_dir/tc1\n
clear\n
quit\n
END
  print $command if ($debug);
  open (DB, "| $tace $dbdir") || die "Couldn't open $dbdir\n";
  print DB $command;
  close DB;

  # special case for Trinity
  if (not -e "$EST_dir/Trinity" || not -e "$EST_dir/IsoSeq") {
    $wormbase->run_command("touch $EST_dir/Trinity", $log);
  }
  
  $log->write_to("$subspecies .ace data dumped.\n\n");
}


__END__

=pod

=head2 NAME - dump_primary_seq_data.pl

=head1 USAGE

=over 4

=item dump_primary_seq_data.pl [-options]

=back

This script dumps ace and sequence files for Transcript data for each species and stores them under BUILD_DATA/cDNAace/ and BUILD_DATA/cDNA

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
