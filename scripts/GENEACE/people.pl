#!/software/bin/perl -w
#
# batch_people.pl
# 
# by Gary Williams                  
#
# This is for managing Person IDs in the new (datomic) NameServer system
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

=head people.pl

=item Options:

  -test      use the test nameserver
  -action    one of "new", "update", "kill", "help"
  -file      TAB or comma delimited file containing input IDs and old/new names  <Mandatory>

    for action "new":
    column 1 - Person's name
    column 2 - Person's email
    column 3 - Person's WBPersonID

    example:
    Joe Bloggs	joe.bloggs@wormbase.org	WBPerson0029076    
    Jane Doe	jane.doe@wormbase.org	WBPerson0029231

    for action "update":
    column 1 - identifier (one of email or WBPerson)
    column 2 - new name, email or WBPerson to give this person - If it doesn't have a '@' or a 'WBPerson' in it then it is assumed to be a name

    example:
    Joe Bloggs	jb1@ebi.ac.uk
    WBPerson0029231	Jessy Doe


    for action "kill":
    column 1 - name, email or WBPerson of person to remove from the database
    
    example:
    WBPerson0029231
    joe.bloggs@wormbase.org


    for action "find":
    column 1 -  person ID to find
    person ID can be either of their email address or WBPersonID

    example name:
    Jane Doe

    for action "help" - no input is required, instructions on how to set up authentication to use the Nameserver are output


  -debug     limits to specified user <Optional>
  

e.g. perl batch_people.pl -action new -file new_name_data

=cut




######################################
# variables and command-line options # 
######################################

my ($test, $help, $debug, $verbose, $store, $wormbase);
my ($file, $action);

GetOptions (
	    "test"       => \$test,
	    "help"       => \$help,
            "debug=s"    => \$debug,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "file:s"     => \$file,
	    "action:s"   => \$action,
	    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# establish log file.
my $log = Log_files->make_build_log($wormbase);



if (!defined $file) {die "-file file not specified\n"}
if (!-e $file) {die "-file file doesn't exist\n"}

if (!defined $action) {die "-action not specified\n"}


open (IN, "<$file") || $log->log_and_die("Can't open file $file");

my $db = NameDB_handler->new($wormbase, $test);


if ($action eq 'new') {
  new_people();
} elsif ($action eq 'update') {
  update_people();
} elsif ($action eq 'kill') {
  kill_people();
} elsif ($action eq 'find') {
  info_person();
} elsif ($action eq 'help') {
  help_authentication();
} else {
  die "-action action '$action' not recognised\n"
}


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

#    column 1 - Person's name
#    column 2 - Person's email
#    column 3 - Person's WBPersonID

sub new_people {
  my @names;

  my $batch;
  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}
 
    my ($name, $email, $wbperson) = split /[\t,]/, $line;
    $name =~ s/ $//;
    $name =~ s/^ //;
    $email =~ s/ $//;
    $email =~ s/^ //;
    $wbperson =~ s/ $//;
    $wbperson =~ s/^ //;
    if ($name eq '') {die "line '$line' has no person name\n"}
    if ($email !~ /\@/) {die "line '$line' has an invalid email\n"}
    if ($wbperson !~ /^WBPerson/) {die "line '$line' has an invalid WBPerson\n"}
    
    print "//\tperson name '$line' queued for creation\n";

    my $info = $db->new_person($name, $email, $wbperson);

    print "// batch '$line' created\n";

  }
}
##########################################
#    column 1 - identifier (one of email or WBPerson)
#    column 2 - new name, email or WBPerson to give this person - If it doesn't have a '@' or a 'WBPerson' in it then it is assumed to be a name

sub update_people {

  while (my $line = <IN>) {
    my ($name, $email, $wbperson);
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}

    my ($person, $update) = split /[\t,]/, $line;
    $person =~ s/ $//;
    $person =~ s/^ //;
    $update =~ s/ $//;
    $update =~ s/^ //;

    if ($update =~ /\@/) {$email = $update}
    elsif ($update =~ /^WBPerson\d+/) {$wbperson = 'WBPerson'}
    else {$name = $update}

    my $info = $db->update_person($person, $name, $email, $wbperson);
    print "// batch '$line' updated\n";
    
  }

}

##########################################

sub kill_people {

  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}
    my ($person) = split /[\t,]/, $line;
    $person =~ s/ $//;
    $person =~ s/^ //;

      my $info = $db->kill_person($person);
      print "// batch '$person' killed\n";
    
  }

}


##########################################

sub info_person {

  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next}
    if ($line eq '') {next}

    my ($person) = split /[\t,]/, $line;
    $person =~ s/ $//;
    $person =~ s/^ //;

    my $info = $db->info_person($person);

    my $email = $info->{'email'};
    my $name = $info->{'name'};
    my $id = $info->{'id'};
    print "Person $person // Name: '$name' ID: '$id' Email: '$email'\n";
    
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
