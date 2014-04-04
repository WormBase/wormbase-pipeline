#!/software/bin/perl -w
#
# get_non_core_EST.pl
# 
# by Gary Williams            
#

# This gets the EST sequences in fasta format from the ENA and formats
# them so that just their accession number of on the title line for
# use by the EST-STAR hive aligner

#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2014-04-04 08:46:47 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my (@only_species, %only_species);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "onlyspecies=s@" => \@only_species,
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

map { $only_species{$_} = 1 } @only_species;

my %tierIII = $wormbase->tier3_species_accessors;
foreach my $t3 (keys %tierIII){

  next if @only_species and not exists $only_species{$t3};

  my $wb = $tierIII{$t3};
  my $species  = $wb->species;
  my $full_species = $wb->full_name();
  my $taxon_id = $wb->ncbi_tax_id;

  print "$full_species (taxon ID: $taxon_id)\n";

  # the EST files go here
  my $target_dir = $wb->cdna_dir;
  if (!-d $target_dir) {mkdir $target_dir, 0777}
  my $outfile = "${target_dir}/EST";
  $outfile =~ s/\s/_/g;

  my $old_size = undef;
  if (-e $outfile) {
    $old_size = -s $outfile;
  }

  open (OUT, ">$outfile") || $log->log_and_die("Can't open the ouput file $outfile");

  # now get the EST/mRNA sequences
  my $query = 'http://www.ebi.ac.uk/ena/data/warehouse/search?query=%22mol_type=%22mRNA%22%20AND%20tax_tree(TAXON)%22&result=sequence_release&display=fasta';
  $query =~ s/TAXON/$taxon_id/;
  open (DATA, "wget -q -O - '$query' |") || $log->log_and_die("Can't get mRNA data for $full_species\n");

  my $entry_count=0;
  while(my $line = <DATA>) {
    if ($line =~ /^>ENA\|\w+\|(\S+)\s/) {
      $line = ">$1\n";
      $entry_count++;
    }
    print OUT $line;
  }

  $log->write_to("$entry_count EST/mRNA entries found\n");

  close(DATA);
  close (OUT);

  if (defined $old_size) {
    my $new_size = -s $outfile;
    if ($new_size < $old_size) {
      $log->write_to("WARNING: the new ${full_species} file is smaller than the old file (old size: $old_size bytes, new_size: $new_size bytes).\n");
      $log->error();
    } elsif ($new_size > $old_size) {
      $log->write_to("The ${full_species} file has been updated (old size: $old_size bytes, new_size: $new_size bytes).\n");
    } else {
      $log->write_to("The ${full_species} EST file has not changed - no updates to do.\n");
    }
  } else {
      $log->write_to("The ${full_species} EST file has been created.\n");
  }
  $log->write_to("\n");
}

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
