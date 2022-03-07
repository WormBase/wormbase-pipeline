#!/usr/bin/env perl
#
# RNASeq_cleanup.pl
# 
# by Gary Williams                        
#
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2015-04-27 13:37:33 $      

#
# The script takes a species name and then runs through all of the
# species' Experiments directories on /nfs/backup/ looking for data
# that should be backed up to the NGS server at the Sanger.
#
# It makes sym-links to these files in the directory $RNASeq->RNASeqTransferDir
# you can then scp the files from RNASeqTransferDir to the Sanger NGS server.
#
# for this to work, there must be a tunnel from the current machine to the Sanger (ssh -f -N -C -X ssh.sanger.ac.uk)
# and the following Host definition must be in the .ssh/config file.
# LocalForward 2229 web-bfint.internal.sanger.ac.uk:22
# the files will appear on the NGS server FTP site within an hour of this transfer
#
# Then give the command:
# ssh -f -N -C -X ssh.sanger.ac.uk
# scp -P 2229 /nfs/nobackup/ensemblgenomes/wormbase/BUILD/RNASeq/$species/Transfer/* ${USER}@localhost:/data/production/parasites/wormbase/RNASeq_alignments/$species\
#
# See also:
# https://www.ebi.ac.uk/seqdb/confluence/pages/viewpage.action?spaceKey=WORMBASE&title=Sanger+NGS+server
# https://www.ebi.ac.uk/seqdb/confluence/display/~jane/RNASeq+processing+for+WormBase+ParaSite
#
# To look at the NGS data:
# ssh to Sanger, then:
# ssh web-bfint
# cd /data/production/parasites/wormbase/RNASeq_alignments
#
# You should check that there is sufficient free space on web-bfint before commencing the scp.


use strict;
use lib $ENV{'CVS_DIR'};
use Carp;
use Modules::RNASeq;
use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use List::Util qw(sum); # for mean()
use Time::localtime;
use File::stat;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species, $experiment, $force);
GetOptions ("help"         => \$help,
            "debug=s"      => \$debug,
            "test"         => \$test,
            "verbose"      => \$verbose,
            "store:s"      => \$store,
            "species:s"    => \$species, # the default is elegans
	    "experiment:s" => \$experiment, # specify just one experiment to cleanup, for testing
	    "force"        => \$force, # stage the files to be copied even if there are problems or it exists already in the NGS server - mainly for debugging or if you can't wait for the NGS server http public view to catch up with its filesystem.
	   );


$species = 'elegans' unless $species;

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species,
			   );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


######################################

my $check = 1; # don't want to remove anything until it is backed up.
my $new_genome = 0; # this is not a new genome

my $RNASeq = RNASeq->new($wormbase, $log, $new_genome, $check);
my $cufflinks_transcripts = $wormbase->build_data."/BUILD_DATA/MISC_DYNAMIC/SHORT_READS/CUFFLINKS_TRANSCRIPTS/$species";
my $cufflinks_done_file = "cufflinks/genes.fpkm_tracking.done";

# foreach SRA Experiment directory, check to see if there are files to be copied to the NGS, move cufflinks transcripts.gtf to backup
my $SRA_dir = $RNASeq->{RNASeqSRADir};
my $TransferDir = $RNASeq->{RNASeqTransferDir};

chdir $SRA_dir;
mkdir ($TransferDir, 0775);
system("rm -f $TransferDir/*");

opendir(DIR, '.' ) || $log->log_and_die("Can't open the SRA directory");

mkdir ($cufflinks_transcripts, 0775);

my $mode = 0664;   

while ( my $file = readdir(DIR) ) {
  next if ( $file eq "." or $file eq ".." );
  if (defined $experiment && $file ne $experiment) {next} # for debugging when we only want one Experiment done
  if (-d $file) {
    $log->write_to("Staging $file: ");
    chdir $file;
    my $bam_done_file = $RNASeq->{'alignmentDir'}."/accepted_hits.bam.done";
    if (!$force && !$RNASeq->check_NGS_files($file)) {$log->write_to("already in the NGS server\n"); chdir ".."; next} # the files are already stored in the NGS server so don't copy them again
    if (!$force && !-e $bam_done_file) {$log->write_to("ERROR - ALIGNMENT FAILED!\n"); chdir ".."; next} # the alignment did not succeed
    if (!$force && !-e $cufflinks_done_file) {$log->write_to("ERROR - CUFFLINKS FAILED!\n"); chdir ".."; next} # the cufflinks analysis did not succeed
    
    $RNASeq->stage_files_for_NGS($file);
    
    system("cp cufflinks/transcripts.gtf $cufflinks_transcripts/$file.gtf");
    chmod $mode, "$cufflinks_transcripts/$file.gtf";
    $log->write_to("Done\n");
    chdir ".."; # back out to top level
  }
}


$log->write_to("\n\nNow you should check that there is sufficient free space on web-bfint.\nThen give the command:\nssh -f -N -C -X ssh.sanger.ac.uk\nscp -P 2229 $TransferDir/* $ENV{'USER'}\@localhost:/data/production/parasites/wormbase/RNASeq_alignments/$species\n(Sanger account password)\n\n");

# remove the Experiment directories and tidy up the Transfer directory
$log->write_to("\nThen if all transfers OK, give the command:\nrm -rf $SRA_dir/* $TransferDir\n\n");



$log->mail();
print "Finished.\n";
exit(0);

