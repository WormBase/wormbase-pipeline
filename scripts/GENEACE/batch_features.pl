#!/software/bin/perl -w
#
# features.pl
# 
# by Gary Williams                  
#
# This is for managing Feature IDs in the new (datomic) NameServer system
#
# Last updated by: $Author: pad $     
# Last updated on: $Date: 2013-08-14 12:19:59 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use lib '/software/worm/lib/perl';
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use FileHandle;          # Or `IO::Handle' or `IO::'-anything-else used for autoflush

use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
use NameDB_handler;


=pod

=head batch_pname_update.pl

=item Options:

  -action    one of "new", "kill", "resurrect", "help"
  -file      TAB or comma delimited file containing input IDs and old/new names  <Mandatory>
  -test      use the test nameserver

    for action "new":
    column 1 - number of IDs to create

    example:
    10100
    

    for action "kill":
    column 1 - Feature ID / name to kill
    
    example:
    WBsf1026946
    WBsf1026947
    WBsf1026948


    for action "resurrect":
    column 1 - Feature ID to resurrect
    
    example:
    WBsf1026946
    WBsf1026947
    WBsf1026948


    for action "help" - no input is required, instructions on how to set up authentication to use the Nameserver are output


  -output    output file holding ace results
  -why       optional string describing the reason for performing the action
  -debug     limits to specified user <Optional>
  -species   can be used to specify non elegans
  

e.g. perl batch_features.pl -species elegans -action new -file feature_name_data -output results_file.ace

=cut




######################################
# variables and command-line options # 
######################################

my ($test, $help, $debug, $verbose, $store, $wormbase);
my ($species, $file, $output, $action, $why);
my $BATCH_SIZE = 500; # maximum entries to put into any one batch API call

GetOptions (
	    "test"       => \$test,
	    "help"       => \$help,
            "debug=s"    => \$debug,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "file:s"     => \$file,
	    "output:s"   => \$output,
	    "action:s"   => \$action,
	    "why:s"      => \$why,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -organism => $species,
			     -test => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# establish log file.
my $log = Log_files->make_build_log($wormbase);



if (!defined $file) {die "-file file not specified\n"}
if (!-e $file) {die "-file file doesn't exist\n"}

if (!defined $output) {die "-output file not specified\n"}

if (!defined $action) {die "-action not specified\n"}


open (IN, "<$file") || $log->log_and_die("Can't open file $file");
open (OUT, ">$output") || $log->log_and_die("Can't open file $output");
OUT->autoflush(1);       # empty the buffer after every line when writing to OUT

my $db = NameDB_handler->new($wormbase, $test);


if ($action eq 'new') {
  new_ft();
} elsif ($action eq 'kill') {
  kill_ft();
} elsif ($action eq 'resurrect') {
  resurrect_ft();
} elsif ($action eq 'help') {
  help_authentication();
} else {
  die "-action action '$action' not recognised\n"
}



close(OUT);
close(IN);
$db->close;


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

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################

sub new_ft {
  my $number;
  
  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}
    print OUT "//\t'$line' features queued for creation\n";
    
    my ($new_ids, $batch) = $db->new_features($line);
    foreach my $id (@{$new_ids}) {
      print OUT "\nFeature : $id\n";
    }
    
    print OUT "// batch '$batch' created\n";
  }


}


##########################################

sub kill_ft {
  my @names;

  my $batch;
  my $count = 0;
  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}
    push @names, $line;
    print OUT "\n-D Feature $line\n";
    $count++;
    print OUT "//\tfeature '$line' queued for being killed\n";
    if ($count == $BATCH_SIZE) {
      my $info = $db->kill_features(\@names, $why);
      $count = 0;
      @names = ();
      $batch = $info->{'dead'}{'id'}; # was batch/id
      print OUT "// batch '$batch' killed\n";
    }
  }

  if ($count) {
    my $info = $db->kill_features(\@names, $why);
    $batch = $info->{'dead'}{'id'}; # was batch/id
    print OUT "// batch '$batch' killed\n";
  }

}

##########################################

sub resurrect_ft {
  my @names;

  my $batch;
  my $count = 0;
  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}
    push @names, $line;
    print OUT "\nFeature : $line\n";
    $count++;
    print OUT "//\tfeature '$line' queued for being resurrected\n";
    if ($count == $BATCH_SIZE) {
      my $info = $db->resurrect_features(\@names, $why);
      $count = 0;
      @names = ();
      $batch = $info->{live}{'id'}; # was batch/id
      print OUT "// batch '$batch' resurrected\n";
    }
  }

  if ($count) {
    my $info = $db->resurrect_features(\@names, $why);
    $batch = $info->{live}{'id'}; # was batch/id
    print OUT "// batch '$batch' resurrected\n";
  }

}


##########################################
sub help_authentication {
  $db->print_authentication_instructions();
}

##########################################
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
