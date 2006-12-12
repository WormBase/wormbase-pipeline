#!/usr/local/bin/perl5.8.0 -w
#
# script_template.pl                           
# 
# by Keith Bradnam                         
#
# This is a example of a good script template
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2006-12-12 16:08:18 $      

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
	    "store:s"      => \$store,
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
my $gff_dir         = $wormbase->gff;         # AUTOACE GFF
###my $gff_dir = glob("~wormpub/DATABASES/current_DB/CHROMOSOMES/");

my $ace_dir = $wormbase->autoace;
###my $ace_dir = glob("~wormpub/DATABASES/current_DB/");

# other paths
my $tace            = $wormbase->tace;        # TACE PATH



###################################
# get the species of the sequences
###################################

my $cmd1 = "Query Find Sequence Where Database = \"NEMATODE_NET\"\nshow -a Species\nquit";
my $cmd2 = "Query Find Sequence Where Database = \"NEMBASE\"\nshow -a Species\nquit";
my $cmd3 = "Query Find Sequence Where Database = \"EMBL\"\nshow -a Species\nquit";

my %species;
my $id;
my $db;

print "Finding BLAT_WASHU data\n";
open (TACE, "echo '$cmd1' | $tace $ace_dir |");
while (<TACE>) {
  chomp;
  next if (/acedb\>/);
  next if (/\/\//);
  if (/Sequence\s+:\s+\"(\S+)\"/) {
    $id = $1;
  } elsif (/Species\s+\"(.+)\"/) {
    $species{'BLAT_WASHU'}->{$id} = $1;
  }
}
close TACE;

print "Finding BLAT_NEMBASE data\n";
open (TACE, "echo '$cmd2' | $tace $ace_dir |");
while (<TACE>) {
  chomp;
  next if (/acedb\>/);
  next if (/\/\//);
  if (/Sequence\s+:\s+\"(\S+)\"/) {
    $id = $1;
  } elsif (/Species\s+\"(.+)\"/) {
    $species{'BLAT_NEMBASE'}->{$id} = $1;
  }
}
close TACE;

print "Finding BLAT_NEMATODE data\n";
open (TACE, "echo '$cmd3' | $tace $ace_dir |");
while (<TACE>) {
  chomp;
  next if (/acedb\>/);
  next if (/\/\//);
  if (/Sequence\s+:\s+\"(\S+)\"/) {
    $id = $1;
  } elsif (/Species\s+\"(.+)\"/) {
    $species{'BLAT_NEMATODE'}->{$id} = $1;
  }
}
close TACE;



##########################
# MAIN BODY OF SCRIPT
##########################
my $count;

# loop through the chromosomes
  my @chromosomes =  $wormbase->get_chromosome_names(-mito => 1, -prefix => 0);
  foreach my $chromosome (@chromosomes) {
    print "Reading chromosome $chromosome\n" if ($verbose);

# loop through the GFF file
    my @f;
    open (GFF, "<$gff_dir/CHROMOSOME_${chromosome}.gff") || die "Failed to open gff file $gff_dir/CHROMOSOME_${chromosome}.gff\n";
    open (OUT, ">$gff_dir/CHROMOSOME_${chromosome}.gff.new") || die "Failed to open gff file $gff_dir/CHROMOSOME_${chromosome}.gff.new\n";
###    open (OUT, ">./CHROMOSOME_${chromosome}.gff.new") || die "Failed to open gff file ./CHROMOSOME_${chromosome}.gff.new\n";
    while (my $line = <GFF>) {
      chomp $line;
      if ($line =~ /^#/ || $line !~ /\S/) {
	print OUT "$line\n";
        next;
      }
      @f = split /\t/, $line;
      my $id;

# is this a BLAT_WASHU or BLAT_NEMBASE or BLAT_NEMATODE line?
      if ($f[1] eq 'BLAT_WASHU') {
	# get the ID name
	($id) = ($f[8] =~ /Target \"Sequence:(\S+)\"/);

	if (exists $species{'BLAT_WASHU'}->{$id}) {
	  $line = $line . " ; Species \"" . $species{'BLAT_WASHU'}->{$id} . "\"";
	  $count++;
	  #print "$line\n" if ($verbose);
	}
      } elsif ($f[1] eq 'BLAT_NEMBASE') {
	# get the ID name
	($id) = ($f[8] =~ /Target \"Sequence:(\S+)\"/);

	if (exists $species{'BLAT_NEMBASE'}->{$id}) {
	  $line = $line . " ; Species \"" . $species{'BLAT_NEMBASE'}->{$id} . "\"";
	  $count++;
	  #print "$line\n" if ($verbose);
	}
      } elsif ($f[1] eq 'BLAT_NEMATODE') {
	# get the ID name
	($id) = ($f[8] =~ /Target \"Sequence:(\S+)\"/);

	if (exists $species{'BLAT_NEMATODE'}->{$id}) {
	  $line = $line . " ; Species \"" . $species{'BLAT_NEMATODE'}->{$id} . "\"";
	  $count++;
	  #print "$line\n" if ($verbose);
	} else {
	  #print "BLAT_NEMATODE species doesn't exist for $id\n";
	}
      }


# write out the line
      print OUT "$line\n";

# end of GFF loop
    }

# close files
    close (GFF);
    close (OUT);

# end of chromosome loop
  }

# copy new GFF files over
  foreach my $chromosome (@chromosomes) {
    system("mv -f $gff_dir/CHROMOSOME_${chromosome}.gff.new $gff_dir/CHROMOSOME_${chromosome}.gff");
  }

# Close log files and exit
$log->write_to("\n\nStatistics\n");
$log->write_to("----------\n\n");
$log->write_to("Changed $count lines\n");


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

=head2 NAME - add_species_to_BLAT_GFF.pl

=head1 USAGE

=over 4

=item  add_species_to_BLAT_GFF.pl [-options]

=back

Todd wanted to have species names added to the GFF results for non-elegans EST BLAT hits.
This script gets the species for each EST sequence from the autoace database and adds this information to the GFF records.

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

=item Gary Williasm (gw3@sanger.ac.uk)

=back

=cut
