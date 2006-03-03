#!/usr/local/bin/perl5.8.0 -w
#
# script_template.pl                           
# 
# by Keith Bradnam                         
#
# This is a example of a good script template
#
# Last updated by: $Author: ar2 $     
# Last updated on: $Date: 2006-03-03 14:48:01 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
#use Ace;
#use Sequence_extract;
#use Coords_converter;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    );

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

# establish log file.
my $log = Log_files->make_build_log($wormbase);

#################################
# Set up some useful paths      #
#################################

# Set up top level base directories (these are different if in test mode)
my $dbpath         = $wormbase->autoace;     # AUTOACE DATABASE DIR


&buildrelease;


##########################
# MAIN BODY OF SCRIPT
##########################

# main stuff goes here

$log->mail();
exit(0);

sub buildrelease{	

  $log->write_to(" Starting to build release files\n\n");
  my $WS_version = $wormbase->get_wormbase_version;
  my $WS_previous   = $WS_version - 1;
  #remove old files - shouldn't be any now that autoace is refreshed each build
  if (-e "$dbpath/release/database.WS"."$WS_previous".".4-0.tar.gz"){
    $log->write_to("Older WS version files exist, removing them\n");
    $wormbase->delete_files_from("$dbpath/release","*WS"."$WS_previous"."*") or $log->write_to("ERROR: Problems removing files from $dbpath/release: $!\n");
  }
  

  $log->write_to("Making distribution files for $WS_version\n\n");
  
  $wormbase->run_command("/bin/touch $dbpath/release/files_in_tar", $log);
  $wormbase->run_command("/bin/touch $dbpath/release/md5sum.${WS_version}", $log);
  
  my @tarfiles;
  $tarfiles[0] = "wspec/cachesize.wrm  wspec/constraints.wrm wspec/copyright wspec/database.wrm wspec/displays.wrm wspec/help.wrm wspec/layout.wrm wspec/models.wrm wspec/options.wrm wspec/passwd.wrm wspec/psfonts.wrm wspec/subclasses.wrm wspec/xfonts.wrm wgf wquery wscripts  pictures database/log.wrm database/database.map database/ACEDB.wrm" ;
  
  for (my $i = 1 ; -e "$dbpath/database/block$i.wrm" ; ++$i) {
    $tarfiles[($i+4)/5] .= " database/block$i.wrm" ;
  }
  $log->write_to("* Makedistr: beginning tar ..\n");
  $wormbase->delete_files_from("$dbpath/release","database\.$WS_version\..*\.tar","-");

  for (my $i = 0; $i < @tarfiles; ++$i) {
    $wormbase->run_command("cd $dbpath; tar -hcf $dbpath/release/database.$WS_version.4-$i.tar $tarfiles[$i]\"");
    
    # list files in the tar archive
    $wormbase->run_command("tar -tf $dbpath/release/database.$WS_version.4-$i.tar >> $dbpath/release/files_in_tar");
    
    # gzip the tar archive
    $wormbase->run_command("/bin/gzip $dbpath/release/database.$WS_version.4-$i.tar"); 
    
    # check consistency of gzip file
    $wormbase->run_command("/bin/gzip -t $dbpath/release/database.$WS_version.4-$i.tar.gz >> $dbpath/release/files_in_tar");
    
    # calculate md5sum for the gzip file
    $wormbase->run_command("/nfs/disk100/wormpub/bin.ALPHA/md5sum $dbpath/release/database.$WS_version.4-$i.tar.gz >> $dbpath/release/md5sum.WS$WS_version");
  }
  #check md5sums
  $wormbase->run_command("cd ".$wormbase->autoace."/release; md5sum -c md5sum.WS$WS_version");
  # zip up the dna and gff files
  $wormbase->run_script("chromosome_dump.pl --zipdna --zipgff", $log);
  $wormbase->run_script( "release_letter.pl -l"               , $log);
}




##############################################################
#
# Subroutines
#
##############################################################



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

=item Keith Bradnam (krb@sanger.ac.uk)

=back

=cut
