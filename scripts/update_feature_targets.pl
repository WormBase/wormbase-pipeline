#!/usr/local/bin/perl5.8.0 -w
#
# update_feature_targets.pl
#
# by
#
#Refreshes some data from the build back to the canonical camace
#
# Last updated on: $Date: 2015-03-24 10:22:46 $
# Last updated by: $Author: pad $

use strict;
use Getopt::Long;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Storable;
use Log_files;
use Carp;

my ($db,$camace,$test,$noload,$store,$wormbase,$species, $debug, $acefiles,);

&GetOptions('db=s'     => \$db,
            'camace=s' => \$camace,
	    'test'     => \$test,
	    'noload'   => \$noload,
	    'debug=s'  => \$debug,
	    'species:s' => \$species,
	    'store:s'  => \$store,
	    'acefiles:s' => \$acefiles,
) || die('cant parse the command line parameter');

######################
# ----- globals -----
######################


my $user = `whoami`; chomp $user;

if ($store) {
  $wormbase = retrieve($store) or croak ("Can't restore wormbase from $store\n");
} 
else {
  $wormbase = Wormbase->new( -debug => $debug,
                             -test => $test,
                             -organism => $species
                           );
}

  $db = $wormbase->autoace unless $db;


my $tace     = $wormbase->tace;          # tace executable path
my $release  = $wormbase->get_wormbase_version_name(); # only the digits
my $ftype;
if ($wormbase->species eq 'elegans') {
  $camace = $wormbase->database('camace');
}
else {
  $camace = $wormbase->database($species);
}

my $log = Log_files->make_build_log($wormbase);

foreach $ftype ("SL1", "SL2", "three_prime_UTR", "transcription_end_site", "segmental_duplication", "polyA_site polyA_signal_sequence", "Genome_sequence_error", "DNAseI_hypersensitive_site", "Corrected_genome_sequence_error", "binding_site_region") {

  my $acedir;
  if (defined $acefiles){ 
    $acedir =  $acefiles;
  }
  else { 
    $acedir = $wormbase->acefiles;
  }
  
  my $ffile = $acedir."/feature_${ftype}.ace";
  $log->write_to ("Processing $ffile\n\n");
  if (-e $ffile) {
    my $tmp_file = &parse_out_parent_sequences($ffile, $ftype, $wormbase->species);
    unless ($noload){
      $wormbase->load_to_database($camace, $tmp_file, "feature_${ftype}_parents",$log);
      unlink $tmp_file;
    }
  }
  else {
    $log->write_to ("Failed to find $ffile\n\n");
  }
}

&remove_bogus_features unless ($noload);
print "I'm here\n";
$log->mail();

exit(0);

sub parse_out_parent_sequences {
  my $infile = shift;
  my $intype = shift;
  my $subspecies = shift;

  my $tmp_file = "/tmp/${subspecies}_load_file.$intype.ace";
  print "$tmp_file\n" if ($debug);
  open my $tmpfh, ">$tmp_file" or die("Could not open $tmp_file for writing\n");
  open my $infh, $infile or die("Could not open $infile for reading\n");

  my ($class, $obj_id);

  while(<$infh>) {
    if (/^(\S+)\s*:\s*\"(\S+)\"/) {
      ($class, $obj_id) = ($1, $2);
      next;
    } elsif (/^\s*$/) {
      ($class, $obj_id) = (undef, undef);
      next;
    }

    if (/^Feature_object\s+(\S+)/ and defined $obj_id and $class eq 'Sequence') {
      print $tmpfh "\nFeature : \"$1\"\nMapping_target $obj_id\n";
    }
  }
  return $tmp_file;
}

sub remove_bogus_features {
$log->write_to("Removing bogus Features from camace (No Flanks and No Method)\n");
my $command;
  $command  = "query find Feature where !method AND !Flanking_sequences\n";
  $command  .= "kill\n";
  $command  .= "y\n";
  $command  .= "save\n";
  $command  .= "quit\n";

  my $tsuser = "bogus_features";
  open (TACE, "| $tace $camace -tsuser $tsuser") || die "Couldn't open $camace\n";
  print TACE $command;
  close TACE;
  $log->write_to ("Removed bogus Features\n");
}


__END__
