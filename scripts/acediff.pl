#!/software/bin/perl -w
#
# acediff.pl : an improved acediff
#
# gw3
#
# Usage : acediff.pl [-options]
#
# Last edited by: $Author: gw3 $
# Last edited on: $Date: 2011-04-26 15:49:01 $
use lib $ENV{'CVS_DIR'};

use strict;
use Getopt::Long;
use Storable;
use Wormbase;
use Log_files;

my ($test, $debug, $store, $file1, $file2, $output, $logfile);

GetOptions (
	    "test"        => \$test,
	    "store:s"     => \$store,
	    "debug:s"     => \$debug,
	    "reference:s" => \$file1,
	    "new:s"       => \$file2,
	    "output:s"    => \$output,
            "logfile:s"   => \$logfile,
	    )
		    or die("invalid commandline option\n");
	    ;

my $wormbase;
# use the time and the process ID to make a unique file extension
my $time = time();
my $pid = "$$";
my $tmpDir= $ENV {'WB_SCRATCH'};
my $outfile1 = "/$tmpDir/acediff1.$pid.$time.new";
my $outfile2 = "/$tmpDir/acediff2.$pid.$time.new";
my $tmp = "/$tmpDir/acediff.$pid.$time.tmp";


if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}

my $log = $logfile ? Log_files->make_log($logfile, $debug) : Log_files->make_build_associated_log($wormbase);

# Change input separator to paragraph mode, but store old mode in $oldlinesep
my $oldlinesep = $/;
$/ = "";

open (FILE, "<$file1") or $log->log_and_die("cant open $file1 : $!\n");
open (PRE, ">$tmp") or $log->log_and_die("cant open $tmp : $!\n");
while (my $record = <FILE>) {
  while ($record =~ /^\s*\/\//) {$record =~ s/^\s*\/\/.*?\n//} # strip out any comments at the start of the record
  $record =~ s/\t/ /g;
  $record =~ s/\n/\t/g;
  print PRE $record, "\n";
}
close (PRE);
close (FILE);
  
# sort the file - sometimes the sorted file doesn't appear - this is very odd - try a few times to make it
my $tries = 5;
while ($tries-- && ! -e "$outfile1") {
  $wormbase->run_command("sort -S 4G $tmp -o $outfile1", $log);
  system('sleep 5'); # wait a few seconds for NFS to realise that there really is a file there
}

open (FILE, "<$file2") or $log->log_and_die("cant open $file2 : $!\n");
open (PRE, ">$tmp") or $log->log_and_die("cant open $tmp : $!\n");
while (my $record = <FILE>) {
  while ($record =~ /^\s*\/\//) {$record =~ s/^\s*\/\/.*?\n//} # strip out any comments at the start of the record
  $record =~ s/\t/ /g;
  $record =~ s/\n/\t/g;
  print PRE $record, "\n";
}
close (PRE);
close (FILE);
  
# sort the file - sometimes the sorted file doesn't appear - this is very odd - try a few times to make it
$tries = 5;
while ($tries-- && ! -e "$outfile2") {
  $wormbase->run_command("sort -S 4G $tmp -o $outfile2", $log);
  system('sleep 5'); # wait a few seconds for NFS to realise that there really is a file there
}

# reset input line separator
$/= $oldlinesep;

# get diff output - use 'system' rather than $wormbase->run_command
# because diff returns '1' if there is a difference, 2 means error
my $status = system("/usr/bin/diff $outfile1 $outfile2 > $tmp");
if ( ( $status >> 8 ) == 2 ) {
  $log->log_and_die("diff exited with an error\n");
}

# now read in $tmp
# when there is an object that is new (no previous instance): output object
# when there is an object that is deleted (no subsequent instance): output -D object
# when there is an object that is changed: read both instances in and out deleted tags as -D tags and add new tags

# '>' is a new object
# '<' is an old object

my %new;
my %reference;

open (ACE, "<$tmp") or $log->log_and_die("cant read from $tmp : $!\n");
while (my $line = <ACE>) {

  # put the data into the hashes
  my ($origin, $first, $rest) = ($line =~ /^([\>\<])\s+(.+?)\t(.+)/);
  if (!defined $origin) {next}
  if ($origin eq '>') {
    # put tags in %new, keyed by class and ID line
    $new{$first} = $rest;
  } elsif ($origin eq '<') {
    # put tags in %reference, keyed by class and ID line
    $reference{$first} = $rest;

  } else { # number or --- lines
    # do nothing
  }

}
close ACE;



# now go through %new and check in %reference to see if this is new or changed
# delete the %reference ones we find that are in %new

open (OUT, ">$output") || $log->log_and_die("cantopen $output : $!\n");
my $started;
foreach my $newkey (keys %new) {
  if (exists $reference{$newkey}) {
    # changed data - store the tags in hases for comparison
    my %n;
    my %r;
    my @new = split /\t/, $new{$newkey};
    my @ref = split /\t/, $reference{$newkey};
    foreach my $tag (@new) {$n{$tag}=1}
    foreach my $tag (@ref) {$r{$tag}=1}

    # removed tags
    $started = 0;
    foreach my $tag (keys %r) {
      if (! exists $n{$tag}) {
	if (!$started) {
	  print OUT "\n// Changed deleted data\n";    
	  print OUT "\n$newkey\n";
	  $started = 1;
	}
	print OUT "-D $tag\n";
      }
    }

    # added tags
    $started = 0;
    foreach my $tag (keys %n) {
      if (! exists $r{$tag}) {
	if (!$started) {
	  print OUT "\n// Changed added data\n";    
	  print OUT "\n$newkey\n";
	  $started = 1;
	}
	print OUT "$tag\n";
      }
    }

    delete $reference{$newkey};

  } else {
    # new data
    print OUT "\n// New data\n";
    print OUT "\n$newkey\n";
    my $rest = $new{$newkey};
    $rest =~ s/\t/\n/g;
    print OUT "$rest";
  }
}


# go through %reference to find the ones that are left - these are just deleted
foreach my $refkey (keys %reference) {
  if (exists $reference{$refkey}) {
    print OUT "\n// Deleted data\n";    
    print OUT "\n-D $refkey\n";
  }
}
close(OUT);

# tidy up
unlink $outfile1;
unlink $outfile2;
unlink $tmp;

$log->mail;
exit(0);




# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - acediff.pl

=head1 USAGE

=over 4

=item acediff.pl  -reference file1.ace -new file2.ace -out changed_file.ace

=back

This script is a rewrite of acediff in perl. It takes a reference ace fiekl and a new ace file and writes out an ace file that will patch the old database to reflect the new data.

script_template.pl MANDATORY arguments:

=over 4

=item -reference file, the old ace file

=back

=over 4

=item -new file, the new ace file

=back

=over 4

=item -output file, the difference ace file

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

=item Gary Williams (gw3@sanger.ac.uk)

=back

=cut
