#!/software/bin/perl -w
#
# remove_unwanted_GFF_lines.pl
# 
# by Gary Williams
#
# Last updated by: $Author: gw3 $
# Last updated on: $Date: 2011-04-28 15:22:41 $

# removes lines from GFF files that have no 'source' column and are of type 'intron'
# these come from the Confirmed_intron tags in the clone Sequence
# objects - they are created by curators but are not updated.

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species);
my ($chrom, $gffdir);

GetOptions (
            "help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"      => \$store,
	    "chrom:s"    => \$chrom,
	    "gff:s"      => \$gffdir,
	    "species:s"  => \$species	    
	    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species
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


$gffdir             = $wormbase->gff unless $gffdir;


# parse GFF lines
my @gff_files;

if (defined($chrom)){
    push(@gff_files,$chrom);
} else {
  @gff_files = $wormbase->get_chromosome_names('-prefix' => 1, '-mito' => 1);
  #for species in many contigs the gffs are all in a single file named after the species
  # e.g. remanei.gff
  if($wormbase->assembly_type eq 'contig') {
    @gff_files = ();
    push(@gff_files, lc($wormbase->species));
  }
}

foreach my $file (@gff_files) {
    
  # Check for existance of GFFfile
  $log->log_and_die("$file doesnt exist\n") unless (-e "$gffdir/$file.gff") ;
  
  open (OUT, ">$gffdir/${file}.gff.new");
  
  open (GFF, "<$gffdir/${file}.gff")  || die "Cannot open $file.gff\n";
  while (<GFF>) {
    chomp;
    
    #skip header lines of file
    
      unless  (/^\S+\s+(\.)\s+(intron)/) {
	print OUT "$_\n";
	next;
      }
	
    }
    close GFF; #_ end of input GFF file
    close OUT; #_ end of output GFF file

    # copy new GFF files over
    system("mv -f $gffdir/${file}.gff.new $gffdir/${file}.gff");

}


$log->mail();
print "Finished.\n" if ($verbose);
exit(0);

##############################################################
#
# Subroutines
#
##############################################################

sub usage {
 my $error = shift;
 
  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }

}

=pod

=head2 NAME - remove_unwanted_GFF_lines.pl

=head1 USAGE

=over 4

=item  remove_unwanted_GFF_lines.pl [-options]

=back

This script removes lines from GFF files that have no 'source' column
and are of type 'intron' these come from the Confirmed_intron tags in
the clone Sequence objects - they are created by curators but are not
updated.

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

=item Gary Williams

=back

=cut
