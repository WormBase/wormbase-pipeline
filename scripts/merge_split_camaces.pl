#!/usr/local/bin/perl5.8.0 -w
#
# merge_split_camaces.pl
# 
# A script to make multiple copies of camace for curation, and merge them back again
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2003-09-26 09:19:35 $

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Getopt::Long;


##############################
# command-line options       #
##############################

my $all;                   # All
my $dan;                   # Use Daniel's split
my $ant;                   # Use Anthony's split
my $paul;                  # Use Paul's split
my $merge;                 # Merging databases
my $split;                 # Splitting databases
my $update;                # Update current database
my $debug;                 # Debug option
my $help;                  # Help menu

GetOptions (
            "all"        => \$all,
	    "dan"        => \$dan,
	    "ant"        => \$ant,
	    "paul"       => \$paul,
	    "merge"      => \$merge,
	    "split"      => \$split,
	    "update"     => \$update,
	    "help"       => \$help,
	    "debug"      => \$debug
	    );

# Help pod if needed
&usage("Help") if ($help);

our $tace       = &tace;
our $WS_version = &get_wormbase_version;

my $WS_previous = $WS_version - 1;
print "WS_version : $WS_version\tWS_previous : $WS_previous\n" if ($debug);

my @databases; #array to store what splits are to be merged.

# directory paths for the split databases

our $current     = "/wormsrv1/camace";
our $directory   = "/nfs/disk100/wormpub/camace_orig/WS${WS_previous}-WS${WS_version}";
our $camace_orig = "/nfs/disk100/wormpub/camace_orig";
our $camace_dl1  = "/nfs/disk100/wormpub/camace_dl1";
our $camace_pad  = "/nfs/disk100/wormpub/camace_pad";
our $camace_ar2  = "/nfs/disk100/wormpub/camace_ar2";


## (1) Merge split databases #1 - do the diffs ##
if ($merge) {

  print "Make a new directory : '$directory'\n" if ($debug);
  mkdir $directory;

  # load @databases array with users to create split databases for
  push(@databases,"orig");
  push(@databases,"ar2") if ($ant || $all);
  push(@databases,"dl1") if ($dan || $all);
  push(@databases,"pad") if ($paul || $all);

  print "You are merging Data from @databases\n\n";

  # dumps the Pseudogene, Transcript, Feature and Sequence classes from the database
  &dump_camace;

  ###################################################
  #   All of the raw data is now dumped to files    #
  ###################################################

  # run acediff on the files

  # Sequence
  my $path_orig = $directory . '/Sequence_orig.ace';
  my $path_dl1  = $directory . '/Sequence_dl1.ace';
  my $path_pad  = $directory . '/Sequence_pad.ace';
  my $path_ar2  = $directory . '/Sequence_ar2.ace';

  system ("acediff $path_orig $path_dl1 > $directory/sequence_diff_dl1.ace") if ($dan  || $all);
  system ("acediff $path_orig $path_pad > $directory/sequence_diff_pad.ace") if ($paul || $all);
  system ("acediff $path_orig $path_ar2 > $directory/sequence_diff_ar2.ace") if ($ant  || $all);

  # Transcript
  $path_orig = $directory . '/Transcript_orig.ace';
  $path_dl1  = $directory . '/Transcript_dl1.ace';
  $path_pad  = $directory . '/Transcript_pad.ace';
  $path_ar2  = $directory . '/Transcript_ar2.ace';

  system ("acediff $path_orig $path_dl1 > $directory/transcript_diff_dl1.ace") if ($dan  || $all);
  system ("acediff $path_orig $path_pad > $directory/transcript_diff_pad.ace") if ($paul || $all);
  system ("acediff $path_orig $path_ar2 > $directory/transcript_diff_ar2.ace") if ($ant  || $all);

  # Feature
  $path_orig = $directory . '/Feature_orig.ace';
  $path_dl1  = $directory . '/Feature_dl1.ace';
  $path_pad  = $directory . '/Feature_pad.ace';
  $path_ar2  = $directory . '/Feature_ar2.ace';

  system ("acediff $path_orig $path_dl1 > $directory/feature_diff_dl1.ace") if ($dan  || $all);
  system ("acediff $path_orig $path_pad > $directory/feature_diff_pad.ace") if ($paul || $all);
  system ("acediff $path_orig $path_ar2 > $directory/feature_diff_ar2.ace") if ($ant  || $all);

  # Pseudogene
  $path_orig = $directory . '/Pseudogene_orig.ace';
  $path_dl1  = $directory . '/Pseudogene_dl1.ace';
  $path_pad  = $directory . '/Pseudogene_pad.ace';
  $path_ar2  = $directory . '/Pseudogene_ar2.ace';

  system ("acediff $path_orig $path_dl1 > $directory/pseudogene_diff_dl1.ace") if ($dan  || $all);
  system ("acediff $path_orig $path_pad > $directory/pseudogene_diff_pad.ace") if ($paul || $all);
  system ("acediff $path_orig $path_ar2 > $directory/pseudogene_diff_ar2.ace") if ($ant  || $all); 


  ###################################################
  # all of the acediffs are now complete            #
  ###################################################

  # tidy up and reformat the diff files ready to be loaded

  if ($dan || $all) {
    system ("reformat_acediff $directory/sequence_diff_dl1.ace   > $directory/update_sequence_dl1.ace");
    system ("reformat_acediff $directory/transcript_diff_dl1.ace > $directory/update_transcript_dl1.ace");
    system ("reformat_acediff $directory/feature_diff_dl1.ace    > $directory/update_feature_dl1.ace");
    system ("reformat_acediff $directory/pseudogene_diff_dl1.ace > $directory/update_pseudogene_dl1.ace");
  }
  if ($paul || $all) {
    system ("reformat_acediff $directory/sequence_diff_pad.ace   > $directory/update_sequence_pad.ace");
    system ("reformat_acediff $directory/transcript_diff_pad.ace > $directory/update_transcript_pad.ace");
    system ("reformat_acediff $directory/feature_diff_pad.ace    > $directory/update_feature_pad.ace");
    system ("reformat_acediff $directory/pseudogene_diff_pad.ace > $directory/update_pseudogene_pad.ace");
  }
  if ($ant || $all) {
    system ("reformat_acediff $directory/sequence_diff_ar2.ace   > $directory/update_sequence_ar2.ace");
    system ("reformat_acediff $directory/transcript_diff_ar2.ace > $directory/update_transcript_ar2.ace");
    system ("reformat_acediff $directory/feature_diff_ar2.ace    > $directory/update_feature_ar2.ace");
    system ("reformat_acediff $directory/pseudogene_diff_ar2.ace > $directory/update_pseudogene_ar2.ace");
  }

  print "acediff's done and files can be found in $directory\n";
}


## (2) synchronises /wormsrv1/camace with the split versions ##
if ($update) {
    &update_camace;
}


## (3) TransferDB calls to move /wormsrv1/camace to the split databases ##
if ($split) {
    print "Removing old split databases and Copying /wormsrv1/database to the split camaces\n";
    &split_databases;
    exit(0);
}

print "hasta luego\n";

exit(0);

###################################################################################################
###################################################################################################

sub dump_camace {
  #dumps out subset of classes from camace splits and processes the files to be loaded back to /wormsrv1/camace
  #array of classes to be dumped
  my @classes = ('Pseudogene', 'Transcript', 'Feature', 'Sequence');

  my $camace_path;
  my $path;

  foreach my $database (@databases) {

    $camace_path = "/nfs/disk100/wormpub/camace_${database}";
    $ENV{'ACEDB'} = $camace_path;

    foreach my $class (@classes) {

      print "dumping $class class from camace_${database}\n";
      $path = "$directory/" . "${class}_${database}.ace";
      &dumpace("$class",$path);
      print "dumped $class class from camace_${database}\n\n";
    }
  }
}

sub dumpace {
    my $class    = shift;
    my $filepath = shift;

    my $command = "nosave\nquery find $class\nshow -a -f $filepath\nquit\n";

    # dump out from ACEDB
    print "\nFilename: $filepath\n";
    open (TACE,"| $tace");
    print TACE $command;
    close TACE;
}

sub loadace {
    my $filepath = shift;
    my $tsuser   = shift;

    my $command = "pparse $filepath\nsave\nquit\n";

    # dump out from ACEDB
    print "\nFilename: $filepath\n";
    open (TACE,"| $tace -tsuser $tsuser");
    print TACE $command;
    close TACE;
}

sub update_camace {
    # upload processed diff files into /wormsrv1/camace
    print "Upload diff files to /wormsrv1/camace";
    $ENV{'ACEDB'} = $current;

    if ($dan || $all) {
	&loadace("$directory/update_sequence_dl1.ace",'dl1');
	&loadace("$directory/update_transcript_dl1.ace",'dl1');
	&loadace("$directory/update_feature_dl1.ace",'dl1');
	&loadace("$directory/update_pseudogene_dl1.ace",'dl1');
    }
    if ($ant || $all) {
	&loadace("$directory/update_sequence_ar2.ace",'ar2');
	&loadace("$directory/update_transcript_ar2.ace",'ar2');
	&loadace("$directory/update_feature_ar2.ace",'ar2');
	&loadace("$directory/update_pseudogene_ar2.ace",'ar2');
    }
    if ($paul || $all) {
	&loadace("$directory/update_sequence_pad.ace",'pad');
	&loadace("$directory/update_transcript_pad.ace",'pad');
	&loadace("$directory/update_feature_pad.ace",'pad');
	&loadace("$directory/update_pseudogene_pad.ace",'pad');
    }

    # uplaod new mRNAs into camace
    print "Upload new mRNAs in /wormsrv1/camace\n";
    &loadace("/nfs/disk100/wormpub/analysis/ESTs/elegans_mRNAs.ace",'NDB_data');

    # upload BLAT results to database
    print "Update BLAT results in /wormsrv1/camace\n";
    system ("load_blat2db.pl -dbdir $current");

    # synchronize the locus - sequence connections
    print  "Update locus2sequence connections in /wormsrv1/camace\n";
    system ("locus2seq.pl -camace -update");
}

sub split_databases {
    # it has been decided that it is better to remove the database directory to make transfer db more stable #
    # initialise/copy to camace_orig (always do this)
    system("rm -rf /nfs/disk100/wormpub/camace_orig/database") && die "Failed to remove camace_orig/database\n";
    system ("TransferDB.pl -start /wormsrv1/camace -end $camace_orig -database -wspec -name camace_orig_WS$WS_version");

    # initialise/copy to camace_ar2 (-ant or -all)
    system("rm -rf /nfs/disk100/wormpub/camace_ar2/database") && die "Failed to remove camace_ar2/database\n";
    system ("TransferDB.pl -start /wormsrv1/camace -end $camace_ar2 -database -wspec -name camace_ar2_WS$WS_version") if ($ant || $all);

    # initialise/copy to camace_dl1 (-dan or -all)
    system("rm -rf /nfs/disk100/wormpub/camace_dl1/database") && die "Failed to remove camace_dl1/database\n";
    system ("TransferDB.pl -start /wormsrv1/camace -end $camace_dl1 -database -wspec -name camace_dl1_WS$WS_version") if ($dan || $all);

    # initialise/copy to camace_pad (-paul or -all)
    system("rm -rf /nfs/disk100/wormpub/camace_pad/database") && die "Failed to remove camace_pad/database\n";
    system ("TransferDB.pl -start /wormsrv1/camace -end $camace_pad  -database -wspec -name camace_pad_WS$WS_version") if ($paul || $all);

    print "CAMACE SPLITS UPDATED";
}

sub usage {
    my $error = shift;

    if ($error eq "Help") {
        # Normal help menu
        system ('perldoc',$0);
        exit (0);
    }
}

__END__

=pod

=head2 NAME - merge_split_camaces

=head1 USAGE:

=over 4

=item merge_split_camaces [-options]

=back

merge_split_camaces is a wrapper with options to automate the merging of the
working split copies of camace, load the updates into the current version
of camace (with some additional housekeeping), and finally running the
TransferDB.pl jobs to (re)generate the working split copies of camace.

merge_split_camaces mandatory arguments:

=over 4

=item none

=back

merge_split_camaces optional arguments:

=over 4

=item -merge, Generate diff files from split camace databases
 
=item -update, Upload diff files to /wormsrv1/camace and add BLAT, Locus data

=item -split, Transfer /wormsrv1/camace into split camace databases

=item -all, Work on all splits, i.e. ar2,dl1,pad

=item -ant, Work on split ar2 only

=item -dan, Work on split dl1 only

=item -paul, Work on split pad only

=item -help, Help page

=item -debug, Verbose/Debug mode

=back

=head1 RUN REQUIREMENTS:

=back

merge_split_camaces must be able to see the /wormsrv1 disk

=head1 RUN OUTPUT:

=back

=head1 EXAMPLES:

=over 4

=item merge_split_camaces -merge -all 

=back

Dumps the relevant classes from camace_orig and the split databases, runs an
acediff to calculate the changes. This acediff file is then reformated to take
account of the acediff bug.

=head1 AUTHOR - Daniel Lawson

Email dl1@sanger.ac.uk

=cut
