#!/software/bin/perl -w
#
# changed_swissprot_proteins.pl
# 
# by Gary Williams
#
# This script finds SwissProt proteins that have been changed since the last Build
#
# Last updated by: $Author: pad $     
# Last updated on: $Date: 2013-08-14 12:19:59 $      

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

my ($help, $debug, $test, $verbose, $store, $wormbase, $output, $database, $version);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "output:s"   => \$output,   # optional output 
	    "database:s" => \$database, # optional database to use, default is autoace
	    "version:i"  => \$version,  # optional WS version to use, default is autoace version
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


##########################
# MAIN BODY OF SCRIPT
##########################

$output ||= $wormbase->reports . '/changed_swissprot_proteins.txt';
my $VER = (defined $version) ? $version : $wormbase->get_wormbase_version;
my $ace_dir =  (defined $database) ? $database : $wormbase->autoace;
my $tace = $wormbase->tace;
my $command = "query Find Wormpep; Database = SwissProt AND History = $VER\nShow -a\nquit\n";

my $title;
my $swissprot;
my $header_printed = 0;

open (OUT, ">$output") || $log->log_and_die("Can't open $output\n");
open (TACE, "echo '$command' | $tace $ace_dir |");
while (<TACE>) {
  next if (/^acedb\>/);
  next if (/^\/\//);
# some typical Histories:
# History 263 created W06H8.8g
# History 262 converted to isoform F13D12.2a F13D12.2
# History 262 reappeared coded by another gene B0218.1
# History 262 removed B0218.1a                            (in the same Protein object as the above line)
# History 239 reappeared as isoform T05H4.6b T05H4.6
# History 262 replaced by CE52351 T05H4.6a
# History 262 reappeared coded by another gene T07A9.9
# History 262 removed T07A9.9a                            (in the same Protein object as the above line)

  if ($_ =~ /^Protein/) {$title = $_}
  if ($_ =~ /^Database\s+\"SwissProt\"\s+\"UniProtAcc\"/) {$swissprot = $_}
  if ($_ =~ /^History\s+(\d+)\s/ && $1 == $VER) {
    if ($header_printed == 0) {
      print OUT $title;
      print OUT $swissprot;
      print "$title $swissprot $_";
      $header_printed = 1;
      $title = undef;
      $swissprot = undef;
    }
    print OUT $_;
  }
  if ($_ =~ /^\s+$/) {
    print OUT $_;
    $header_printed = 0;
  }
}
close TACE;
close OUT;




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
