#!/software/bin/perl -w
#
# geneace_patch.pl
# 
# A script to make patch files between the last upload and the current geneace.
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2011-04-06 09:25:23 $
#
#====================

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
my $wormbase;

##############################
# command-line options       #
##############################

my $debug;                 # Debug option
my $store;                 # Storable not needed as this is not a build script!
my $test;
my $email;                 # Option for child scripts that can take a user email option.
my $ldate;                  # Option for including a date to do an incremental diff

  GetOptions (
	      "debug:s"    => \$debug,
	      "store"      => \$store,
	      "test"       => \$test,
	      "email:s"    => \$email,
	      "ldate:s"    => \$ldate,
	     );

if ($store) {
  $wormbase = retrieve($store) or croak ("Can't restore wormbase from $store\n");
} 
else {
  $wormbase = Wormbase->new( -debug => $debug,
			     -test => $test,
			   );
}

my $log = Log_files->make_build_log($wormbase);
my $tace = $wormbase->tace;
my $acediff = "/nfs/users/nfs_a/acedb/RELEASE.2011_01_13.BUILD/bin.LINUX_64/acediff";
my $reformat = "/nfs/users/nfs_w/wormpub/wormbase/scripts/reformat_acediff";
my $WS_version = $wormbase->get_wormbase_version;
my $next_build = $WS_version + 1;

#my $date = `date +%y-%m%d-%H:%M`;
my $date = `date +%y%m%d`;
chomp $date;
#debugging
#$date++;


# directory paths
my $wormpub = $wormbase->wormpub;
our $curationdb = $wormbase->database('geneace');
our $builddb = $wormpub."/BUILD/PRIMARIES/geneace";
our $directory   = $curationdb."/Patch_data_WS${WS_version}";

#setup directory path for patch files
if ((-e $directory) && (!-e "$directory/COMPLETE")) {
  $log->log_and_die("ERROR: Previous run failed to complete as $directory/COMPLETE Lock file missing\n");
  $log->write_to ("Using existing directory : '$directory'\n\n");
}
elsif ((-e $directory) && (-e "$directory/COMPLETE")) {
  $wormbase->run_command("rm $directory/COMPLETE*\n");
  $log->write_to ("Using existing directory : '$directory'\n\n");
  #  $log->write_to ("Removing old data\n\n");
  #  $wormbase->run_command("rm $directory/*\n") && $log->log_and_die("Failed to remove old data.\n");
}
else {
  mkdir ($directory) or die "Failed to create ${directory}\n";
}


# what classes of data do we want to dump?
my @classes = ('Gene_name','Variation_name','Strain');
my @databases;
# don't re-dump reference data as there is already a dump on disk for comparison.
if ($ldate) {
  @databases = ('curationdb')
}
# dumps geneace data from the curation database and the one used in the current build.
else {
  @databases = ('curationdb', 'builddb');
}

my $database;
foreach $database(@databases) {
  if ($database =~ "curationdb") {
    $ENV{'ACEDB'} = $curationdb;
  }
  if ($database =~ "builddb") {
    $ENV{'ACEDB'} = $builddb;
  }
  open (TACE,"| $tace") or die "Failed to open database connection\n";
  foreach my $class (@classes) {
    $log->write_to ("dumping $class class from $database\n");
    my $path = "$directory/" . "${date}_${class}_${database}.ace";
    my $command = "nosave\nquery find $class\nshow -a -f $path\n";
    print TACE $command;
    $log->write_to ("dumped $class class from camace_${database}\n\n");
  }
  close TACE;
}

foreach my $class (@classes) {
  my $path_new = $directory."/${date}_${class}_curationdb.ace";
  my $path_ref;
  # diff to last dump.
  if ($ldate) {
    $path_ref = $directory."/${ldate}_${class}_curationdb.ace";
  }
  # diff to staging 
  else {
    $path_ref = $directory."/${date}_${class}_builddb.ace";
  }
  
  print "running.... acediff $path_ref $path_new > $directory/${date}_${class}_diff.ace\n\n";
  $wormbase->run_command("csh -c \"$acediff $path_ref $path_new >! $directory/${date}_${class}_diff.ace\"", $log) && die "acediff Failed for ${path_new}\n";
  print "running.... reformat_acediff -file $directory/${date}_${class}_diff.ace -fileout $directory/${date}_update_${class}.ace\n\n";
  $wormbase->run_command("/software/bin/perl $reformat -file $directory/${date}_${class}_diff.ace -fileout $directory/${date}_update_${class}.ace -sclass $class", $log) && die "reformat Failed for ${date}_${class}_diff.ace\n";
  $log->write_to ("update file ready for $class class\n\n");
}

$wormbase->run_command("touch $directory/COMPLETE");

$log->mail();
exit(0);
__END__

=pod

=head2 NAME - geneace_patch.pl

=head1 USAGE:

=over 4

=item  geneace_patch.pl [-options]

=back

geneace_patch.pl is a wrapper with options to automate the production of the 
incremental patch files for the web team to update the WS release mid cycle.

geneace_patch.pl mandatory arguments:

=over 4

=item none

=back

geneace_patch.pl optional arguments:

=over 4

=item -debug, Verbose/Debug mode

=back

=over 4

=item -ldate, compares the current annotation against the previous patch files giving an incremental update.

=back

=head1 RUN REQUIREMENTS:

Access to geneace and permission to write under the geneace directory

=back

=head1 EXAMPLES:

=over 4

=item geneace_patch.pl

Does a comparison between the current geneace instance and the one used in the 
current/last build

=back

=over 4

=item geneace_patch.pl -ldate 110405

Does a comparison between the current geneace instance and the data dumped from
the patch data done on 110405

=back

=over 4

=head1 AUTHOR - Paul Davis

Email pad@sanger.ac.uk

=cut
