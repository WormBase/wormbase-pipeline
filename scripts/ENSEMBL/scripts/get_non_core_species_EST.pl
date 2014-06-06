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
# Last updated on: $Date: 2014-06-06 09:03:36 $      

use strict;                                      
use Getopt::Long;
use Carp;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


######################################
# variables and command-line options # 
######################################

my ($help, $project_dir);
my (@only_species, %only_species);

GetOptions ("help"           => \$help,
	    "project_dir:s"  => \$project_dir, # usually $wormbase or $parasite
	    "species=s@"     => \@only_species,
	    );

# Display help if required
&usage("Help") if ($help);

if (!defined $project_dir) {$project_dir = "/nfs/panda/ensemblgenomes/wormbase"}

##########################
# MAIN BODY OF SCRIPT
##########################

# Get the Registry
# This assumes the environment variable ENSEMBL_REGISTRY is set, this is
#     used as the name of the configuration file to read.
# Or, the file .ensembl_init exists in the home directory, it is
#     used as the configuration file.
Bio::EnsEMBL::Registry->load_all();

map { $only_species{$_} = 1 } @only_species;

foreach my $species (keys %only_species){

  my $meta_container = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'MetaContainer');
  my $taxon_id = $meta_container->get_taxonomy_id();

  print "$species taxon ID: $taxon_id\n";

  # the EST files go here
  my $target_dir = "$project_dir/BUILD_DATA/cDNA/$species";

  if (!-d $target_dir) {mkdir $target_dir, 0777}
  my $outfile = "${target_dir}/EST";
  $outfile =~ s/\s/_/g;

  my $old_size = undef;
  if (-e $outfile) {
    $old_size = -s $outfile;
  }

  open (OUT, ">$outfile") || die("Can't open the ouput file $outfile");

  # now get the EST/mRNA sequences
  my $query = 'http://www.ebi.ac.uk/ena/data/warehouse/search?query=%22mol_type=%22mRNA%22%20AND%20tax_tree(TAXON)%22&result=sequence_release&display=fasta';
  $query =~ s/TAXON/$taxon_id/;
  open (DATA, "wget -q -O - '$query' |") || die("Can't get mRNA data for $species\n");

  my $entry_count=0;
  while(my $line = <DATA>) {
    if ($line =~ /^>ENA\|\w+\|(\S+)\s/) {
      $line = ">$1\n";
      $entry_count++;
    }
    print OUT $line;
  }

  print "$entry_count EST/mRNA entries found\n";

  close(DATA);
  close (OUT);

  if (defined $old_size) {
    my $new_size = -s $outfile;
    if ($new_size < $old_size) {
      print "WARNING: the new ${species} file is smaller than the old file (old size: $old_size bytes, new_size: $new_size bytes).\n";
    } elsif ($new_size > $old_size) {
      print "The ${species} file has been updated (old size: $old_size bytes, new_size: $new_size bytes).\n";
    } else {
      print "The ${species} EST file has not changed - no updates to do.\n";
    }
  } else {
      print "The ${species} EST file has been created.\n";
  }
  print "\n";
}

print "Finished.\n";
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
