#!/software/bin/perl -w
#
# ABC.pl                           
# 
# by Gary Williams                       
#
# This is a script to automate the sections A, B and C of the BLAST Build
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2008-11-24 13:57:15 $      

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
my ($password, $species_list);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "password:s" => \$password,	# the wormpub password
	    "species:s"  => \$species_list, # the comma-delimted list of species to process
	    );


if ($ENV{USER} eq "wormpub") {
  die "Sorry, you are logged in as 'wormpub'\nYou must run this as yourself.\n";
}

if (!defined $password) {
  print "wormpub password> ";
  $password = <STDIN>;
}

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
  my $spDB;
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
  my $store_file = $spDB->build_store; # create and get the store file to use in all commands for this species


  # I don't think setting this SORT option is necessary until the sort
  # command is run under user wormpub, but the Build Guide says do it now,
  # so ...
  $ENV{SORT_OPTS} = "-k2,2 -k8,8n -k10,10nr";

  $log->write_to("  Running dump.pl . . .\n");
  print "  Running dump.pl . . .\n" if ($verbose);
  chdir "/lustre/work1/ensembl/wormpipe/sort_dump" || die "Can't cd to /lustre/work1/ensembl/wormpipe/sort_dump: $!";

  # check that we have enough free disk-space to work in
  open (DF, "df . |") || die "Can't run 'df .'\n";
  <DF>;<DF>;			# skip the first two lines
  my $df = <DF>;
  my ($available) = ($df =~ /\d+\s+(\d+)/);
  close(DF);
  if ($available < 30000000) { # NB df returns the available space in 1Kb blocks
    $log->error("\nThere is less than 30Gb of free disk space available on /lustre/work1/ensembl/wormpipe/sort_dump\nAborting - not running $species\n");
    next;
  }

  $ENV{PERL5LIB} =  ".:/software/worm/lib/site_perl:/software/worm/lib/bioperl-live:/software/worm/ensembl/ensembl-pipeline/modules:/software/worm/ensembl/old_ensembl/modules";
  $wormbase->run_command('rm -f /lustre/work1/ensembl/wormpipe/sort_dump/junk*', $log);
  $wormbase->run_command("/software/bin/perl ../script/dump.pl -db worm_ensembl_$species", $log);


  $log->write_to("  Logging in as wormpub . . .\n");
  print "  Logging in as wormpub . . .\n" if ($verbose);
  my $exp = new Expect;
  $exp->raw_pty(1);  
  $exp->spawn('ssh wormpub@farm-login')
      or die "Cannot spawn ssh: $!\n";
  my $timeout = 10;
  $exp->expect($timeout,
	       [ qr'Password: $' => sub { my $exp = shift;
					  $exp->send("$password\n");
					  exp_continue; } ],
	       );
  

  #my $logfile = "~/BUILD/autoace/logs/ABC_$species.log";
  #$log->write_to("  Log file for the wormpub session: $logfile\n");
  #print "  Log file for the wormpub session: $logfile\n";
  #$exp->send("script $logfile\n");

  $log->write_to("  Sorting $species.srt . . .\n");
  print "  Sorting $species.srt . . .\n" if ($verbose);
  $exp->send("setenv SORT_OPTS \"-k2,2 -k8,8n -k10,10nr\"\n");
  &wait_for_prompt($exp);
  $exp->send("cd /lustre/work1/ensembl/wormpipe/sort_dump\n");
  &wait_for_prompt($exp);
  $exp->send("time sort -m -S 2G -T tmp \$SORT_OPTS junk*.srt -o $species.srt\n"); 
  &wait_for_prompt($exp);

  $log->write_to("  Running dump_blastp_from_file.pl . . .\n");
  print "  Running dump_blastp_from_file.pl . . .\n" if ($verbose);
  my $cmd = "/software/bin/perl \$CVS_DIR/BLAST_scripts/dump_blastp_from_file.pl $species.srt -version $version -matches -database worm_$species";
  $cmd .= " -store $store_file";
  #print "cmd = $cmd\n";
  $exp->send("$cmd\n");
  &wait_for_prompt($exp);
  $exp->send("rm -f /lustre/work1/ensembl/wormpipe/sort_dump/junk*.*\n");
  &wait_for_prompt($exp);

  $log->write_to("  Running Motif data . . .\n") if ($verbose);
  print "  Running Motif data . . .\n";
  $cmd = "/software/bin/perl \$CVS_DIR/BLAST_scripts/dump_motif.pl -database worm_ensembl_$species";
  $cmd .= " -store $store_file";
  #print "cmd = $cmd\n";
  $exp->send("$cmd\n");
  &wait_for_prompt($exp);

  $cmd = "/software/bin/perl \$CVS_DIR/BLAST_scripts/dump_interpro_motif.pl -database worm_ensembl_$species";
  $cmd .= " -store $store_file";
  #print "cmd = $cmd\n";
  $exp->send("$cmd\n");
  &wait_for_prompt($exp);

  $log->write_to("  Running Repeat data . . .\n");
  print "  Running Repeat data . . .\n" if ($verbose);
  $cmd = "/software/bin/perl \$CVS_DIR/BLAST_scripts/dump_repeats.pl -database worm_ensembl_$species";
  $cmd .= " -store $store_file";
  #print "cmd = $cmd\n";
  $exp->send("$cmd\n");
  &wait_for_prompt($exp);

  $log->write_to("  Logging out of wormpub . . .\n");
  print "  Logging out of wormpub . . .\n" if ($verbose);

  #$exp->send("sleep 10\n"); # rest a while and survey the fruit of our labours

  #$exp->send("exit\n");		# from the 'script' command
  $exp->send("exit\n");		# from the ssh session
  # do a soft_close to nicely shut down the command
  $exp->soft_close();


}

$log->mail();
print "\nFinished.\n";
exit(0);






##############################################################
#
# Subroutines
#
##############################################################

sub wait_for_prompt {
  my $exp = shift @_;
  my $timeout = 36000;
  $exp->expect($timeout, [ qr'\[wormpub]\d+: $' ] ); # wait for the prompt

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
