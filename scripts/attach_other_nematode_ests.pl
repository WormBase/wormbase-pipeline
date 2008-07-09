#!/software/bin/perl -w
#
# attach_other_nematode_ests.pl                   
# 
# by Gary Williams                        
#
# 
#
# Last updated by: $Author: mh6 $     
# Last updated on: $Date: 2008-07-09 12:29:15 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Modules::Overlap;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($species, $output, $database, $load);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "output:s"   => \$output,
	    "database:s" => \$database,
	    "load"       => \$load,
	    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species,
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


$database = $wormbase->autoace unless $database;
$output   = $wormbase->acefiles . "/attach_other_nematode_ests.ace" unless $output;

##########################
# MAIN BODY OF SCRIPT
##########################

open (ACE, "> $output") || die "cant open $output\n";

# get the Overlap object
my $ovlp = Overlap->new($database, $wormbase);


# We want the best match to a gene in the whole genome,
# not just the best in each chromosome, so have a hash of
# key: EST ID
# value: href with keys: {best score, gene}
# and update it it we find a better score
my %best_hits;			# hash of best-scoring matches to genes


# loop through the chromosomes
my @chromosomes = $wormbase->get_chromosome_names(-mito => 0, -prefix => 1);
foreach my $chromosome (@chromosomes) {

  print "reading GFF data for chromosome $chromosome\n";
  my @nematode_hsp = $ovlp->get_EST_NEMATODE($chromosome);
  my @nembase_hsp  = $ovlp->get_EST_NEMBASE($chromosome);
  my @washu_hsp    = $ovlp->get_EST_WASHU($chromosome);

  my @ests = (@nematode_hsp, @nembase_hsp, @washu_hsp);

# sort ests by chromosomal start position
  @ests = sort {$a->[1] <=> $b->[1]} @ests;

  my @genes = $ovlp->get_Genes($chromosome);

  &get_matching_nematode_est(\@genes, \@ests, \%best_hits, $chromosome);

}

print "writing output\n";
&write_output(\%best_hits);

close(ACE);


if ($load) {
  $wormbase->load_to_database($database, $output, 'attach_other_nematode_ests', $log);
}


$log->mail();
print "Finished.\n" if ($verbose);
exit(0);

##########################################


sub get_matching_nematode_est {

  my ($genes_aref, $ests_aref, $best_hits_href, $chromosome) = @_;

# set up the overlap compare objects for the secondary lists
  my $genes_obj = $ovlp->compare($genes_aref); # look for matches in either sense

  foreach my $est (@{$ests_aref}) { # $id, $chrom_start, $chrom_end, $chrom_strand, $hist_start, $hit_end, $score
    my $est_id = $est->[0];
    my $score = $est->[6];
    my @gene_matches = $genes_obj->match($est);
    my @genes = $genes_obj->matching_IDs;
    foreach my $gene (@genes) {
      if (exists $best_hits_href->{$est_id}{score} && $best_hits_href->{$est_id}{score} > $score) {next;}
      $best_hits_href->{$est_id}{score} = $score;
      $best_hits_href->{$est_id}{gene} = $gene;
    }
  }

}

##########################################
sub write_output {
  my ($best_hits_href) = @_;
  
  foreach my $est_id (keys %{$best_hits_href}) {
    my $gene_id = $best_hits_href->{$est_id}{gene};
    print "$gene_id\t$est_id\n" if ($verbose);
    print ACE "\n";
    print ACE "Sequence : \"$est_id\"\n";
    print ACE "Gene \"$gene_id\"\n";
  }
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

=head2 NAME - attach_other_nematode_ests.pl

=head1 USAGE

=over 4

=item attach_other_nematode_ests.pl  [-options]

=back

This script creates an ace file to attach other nematode ESTS (those with Method = EST_nematode) to the Gene that they align to best using BLAT.

attach_other_nematode_ests.pl MANDATORY arguments:

=over 4

=item none

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -output

This specifies an ace file to write to other than the default

The default is acefiles/attach_other_nematode_ests.ace

=over 4

=back

=item -load

This specifies that the output file should be loaded into the database.

=over 4

=back

=item -database

This specifies an ACeDB database to use to read GFF data etc.

The default is to use autoace.

=over 4

=back

=item -h, Help

=back

=over 4
 
=item -debug, Verbose/Debug mode
 
=back

=over 4

=item -test, Test mode, run the script using TEST_BUILD

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back
                                                                                             

=head1 REQUIREMENTS

=over 4

=item There must be the following GFF files in the GFF_SPLITS directory:
				*_BLAT_NEMATODE.gff
				*_BLAT_NEMBASE.gff
				*_BLAT_WASHU.gff
                                *_genes.gff

=back

=head1 AUTHOR

=over 4

=item Gary Williams (gw3@sanger.ac.uk)

=back

=cut
