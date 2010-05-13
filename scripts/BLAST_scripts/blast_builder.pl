#!/software/bin/perl -w
#
# ABC.pl                           
# 
# by Gary Williams                       
#
# This is a script to automate the sections A, B and C of the BLAST Build
#
# Last updated by: $Author: pad $     
# Last updated on: $Date: 2010-05-13 09:58:59 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Expect;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($species_list, $all, $dump, $blastp_dump, $motif_dump, $repeat_dump, $check, );

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "test"        => \$test,
	    "verbose"     => \$verbose,
	    "store:s"     => \$store,
	    "species:s"   => \$species_list, # the comma-delimted list of species to process
	    "all"         => \$all,          # run all the parts of the analysis
	    "dump"        => \$dump,         # run the dump part of the script
	    "blastp_dump" => \$blastp_dump,  # run the blastp_dump part of the script
	    "motif_dump"  => \$motif_dump,   # run the motif_dump part of the script
	    "repeat_dump" => \$repeat_dump,  # run the repeat_dump part of the script
	    "check"       => \$check,        # run the check part of the script
	    )||die(@!);


#############################################################################
#
# The user must be logged in as themself.
#
# We get complaints from System if we run jobs in LSF as wormpub
# because they don't known which one of us to contact when there is a
# problem.
#
#############################################################################
my $USER = $ENV{USER};
if ($USER eq "wormpub") {
  die "Sorry, you are logged in as 'wormpub'\nYou must run this as yourself.\n";
}

#############################################################################
#
# Check that the user has group set as 'worm' and file permission mask 002
# so that all others in this group can access and change the files
#
#############################################################################
my ($username, $pass, $uid, $gid, $quota, $comment, $gcos, $dir, $shell, $expire) = getpwnam($USER);
my ($groupname, $passwd, $groupid, $members) = getgrgid($gid);
if ($groupname ne 'worm') {
  die "Sorry, your group is set to be '$groupname'\n You must run this as group 'worm' (Run: 'newgrp worm').\n";
}
umask 002; # make all files group-read/writeable


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);
}

my @species = split /[,\s]+/, $species_list;
my $spDB;

# establish log file.
my $log = Log_files->make_build_log($wormbase);

##########################
# MAIN BODY OF SCRIPT
##########################

my $version = $wormbase->get_wormbase_version;
my %accessors = $wormbase->species_accessors;

foreach my $species (@species) {
  $log->write_to("Processing $species\n");
  print "Processing $species\n";
# get the species object from the list of accessors
  if (!exists $accessors{$species}) {
    if ($species eq 'elegans') {
      $spDB = $wormbase;
    } else {
      $log->error("The supplied species name '$species' is not a valid TierII species\n");
      next;
    }
  } else {
    $spDB = $accessors{$species};
  }


  &dump($species)          if ($dump || $all);
  &blastp_dump($species)   if ($blastp_dump || $all);
  &motif_dump($species)    if ($motif_dump || $all);
  &repeat_dump($species)   if ($repeat_dump || $all);
  &check($species)         if ($check || $all);
    
}

$log->mail();
print "\nFinished.\n";
exit(0);


##############################################################
#
# Subroutines
#
##############################################################
    
sub dump {

  my ($species) = @_;

  $log->write_to("  Running dump.pl . . .\n");

  my $sort_dir = "/lustre/scratch101/ensembl/wormpipe/sort_dump";
  chdir $sort_dir || die "Can't cd to $sort_dir: $!";
  $spDB->run_command('rm -f $sort_dir/junk*', $log);

  # check that we have enough free disk-space to work in
  open (DF, "df . |") || die "Can't run 'df .'\n";
  <DF>;<DF>;			# skip the first two lines
  my $df = <DF>;
  my ($available) = ($df =~ /\d+\s+(\d+)/);
  close(DF);
  if ($available < 30000000) { # NB df returns the available space in 1Kb blocks
    $log->error("\nThere is less than 30Gb of free disk space available on $sort_dir\nAborting - not running $species\n");
    next;
  }

  $ENV{PERL5LIB} =  ".:/software/worm/lib/site_perl:/software/worm/lib/bioperl-live:/software/worm/ensembl/ensembl-pipeline/modules:/software/worm/ensembl/old_ensembl/modules";
  
  #$spDB->run_command("/software/bin/perl ../script/dump.pl -db worm_ensembl_$species", $log);
  $spDB->run_script("BLAST_scripts/dump.pl -db worm_ensembl_$species", $log);

  $log->write_to("  Sorting $species.srt . . .\n");
  $ENV{SORT_OPTS} = "-k2,2 -k8,8n -k10,10nr";
  chdir "/lustre/scratch101/ensembl/wormpipe/sort_dump";
  $spDB->run_command("time sort -m -S 2G -T tmp \$SORT_OPTS junk*.srt -o $species.srt", $log); 
    
}

sub blastp_dump {

  my ($species) = @_;

  $log->write_to("  Running dump_blastp_from_file.pl . . .\n");
  $spDB->run_script("BLAST_scripts/dump_blastp_from_file.pl $species.srt -version $version -matches -database worm_$species", $log);  
  $spDB->run_command('rm -f /lustre/scratch101/ensembl/wormpipe/sort_dump/junk*.*', $log);
    
}

sub motif_dump {

  my ($species) = @_;

  $log->write_to("  Running Motif data . . .\n");
  $spDB->run_script("BLAST_scripts/dump_motif.pl -database worm_ensembl_$species", $log);  
  $spDB->run_script("BLAST_scripts/dump_interpro_motif.pl -database worm_ensembl_$species", $log);
    
}

sub repeat_dump {

  my ($species) = @_;

  $log->write_to("  Running Repeat data . . .\n");
  $spDB->run_script("BLAST_scripts/dump_repeats.pl -database worm_ensembl_$species", $log);
}


##################
# Check the files
##################
sub check {

  my ($species) = @_;

  my $dump_dir = "/lustre/scratch101/ensembl/wormpipe/dumps/";

  $spDB->check_file($dump_dir."/${species}_blastx.ace", 
		    $log,
		    minsize => 1000000,
		   );

  $spDB->check_file($dump_dir."/${species}_blastp.ace", 
		    $log,
		    minsize => 1000000,
		   );
  
  $spDB->check_file($dump_dir."/worm_ensembl_${species}_motif_info.ace", 
		    $log,
		    minsize => 1000000,
		   );
  
  $spDB->check_file($dump_dir."/worm_ensembl_${species}_interpro_motif_info.ace", 
		    $log,
		    minsize => 1000000,
		   );
  
  $spDB->check_file($spDB->acefiles."/repeat_homologies.ace", 
		    $log,
		    minsize => 1000000,
		   );

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
