#!/software/bin/perl -w
#
# script_template.pl                           
# 
# by xxx                         
#
# This is a example of a good script template
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

my ($help, $debug, $test, $verbose, $store, $wormbase, $outfile);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "outfile:s"  => \$outfile,
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

$outfile ||= $wormbase->misc_output . "/changed_CGC_names.dat";

my $RELEASE = $wormbase->get_wormbase_version;;
my $PREV_RELEASE = $RELEASE - 1;
my $prevdir = $wormbase->wormpub . "/DATABASES/WS${PREV_RELEASE}/COMMON_DATA/";
my %gene2cgc;
my %prev_gene2cgc;

$wormbase->FetchData('worm_gene2cgc_name', \%gene2cgc);
$wormbase->FetchData('worm_gene2cgc_name', \%prev_gene2cgc, $prevdir);

&remove_other_species(\%gene2cgc);
&remove_other_species(\%prev_gene2cgc);

my ($cgc2gene, $cgc2biotype, $sequence_name, $geneid2cgc) = &key_by_cgc(\%gene2cgc);
my ($prev_cgc2gene, $prev_cgc2biotype, $prev_sequence_name, $prev_geneid2cgc) = &key_by_cgc(\%prev_gene2cgc);


open (OUT, "> $outfile") || $log->log_and_die("Can't open $outfile\n");
&output_changed_cgc();
close OUT;


$log->mail();
print "Finished.\n" if ($verbose);
exit(0);






##############################################################
#
# Subroutines
#
##############################################################

sub remove_other_species {

  my ($hash_ref) = @_;

#          'F09E10.11b' => 'Transcript tts-1 WBGene00006650',
#          'F21F8.8' => 'Pseudogene str-137 WBGene00006185',
#          'W10G6.2a' => 'CDS sgk-1 WBGene00004789',
#          'CRE15400' => 'CDS Cre-epn-1 WBGene00058561',
#          'CBN05963' => 'CDS Cbn-srw-67 WBGene00144688',
#          'C05E4.1' => 'CDS srp-2 WBGene00005643',
 
  foreach my $key (keys %{$hash_ref}) {
    if ($hash_ref->{$key} =~ /^\S+\s+\w+\-\w+\-[\d\.]+\s+\S+$/) {
#      print "delete $hash_ref->{$key}\n";
      delete $hash_ref->{$key};
    }
  }


}


##########################################
#my %cgc2gene = &key_by_cgc(\%gene2cgc)
sub key_by_cgc {
  my ($hash_ref) = @_;
  my %gene;
  my %biotype;
  my %sequence_name;
  my %geneid2cgc; # key GeneID, value CGC_name


  foreach my $key (keys %{$hash_ref}) {
    my ($biotype, $cgc, $geneid) = ($hash_ref->{$key} =~ /^(\S+)+\s+(\S+)+\s+(\S+)+$/);
    $gene{$cgc} = $geneid;
    $biotype{$cgc} = $biotype;
    my $sequence_name = $key;
    $sequence_name =~ s/^(\S+?)[a-z]$/$1/;
    $sequence_name{$cgc} = $sequence_name;
    $geneid2cgc{$geneid} = $cgc;
  }

  return (\%gene, \%biotype, \%sequence_name, \%geneid2cgc);
}
##########################################

sub output_changed_cgc {

  # we want prev and new data output for:
  # cgc_name, GeneID, Sequence_name, Biotype

  # want to group by:
  # changed cgc_name on a gene
  # no previous cgc_name
  # no new cgc_name
  # GeneID different to previous one
  # Biotype changed


  my %all_cgcs; # unique set of all cgc_names
  foreach my $cgc (keys %{$cgc2gene}, keys %{$prev_cgc2gene}) {
    $all_cgcs{$cgc} = 1;
  }


  my @changed_cgc_name_on_a_gene;
  my @no_previous_cgc_name;
  my @no_new_cgc_name;
  my @geneid_different_to_prev;
  my @biotype_changed;
  my @other;

  # find the genes whose CGC name has changed
  foreach my $cgc (keys %all_cgcs) {
    if (!$all_cgcs{$cgc}) {next} # already dealt with this one
    my $prev = (exists $prev_cgc2gene->{$cgc}) ? "$cgc\t$prev_cgc2gene->{$cgc}\t$prev_sequence_name->{$cgc}\t$prev_cgc2biotype->{$cgc}" : ".\t.\t.\t.";
    my $new = (exists $cgc2gene->{$cgc}) ? "$cgc\t$cgc2gene->{$cgc}\t$sequence_name->{$cgc}\t$cgc2biotype->{$cgc}" : ".\t.\t.\t.";
    if ($prev ne $new) {
      if (exists $prev_cgc2gene->{$cgc} && !exists $cgc2gene->{$cgc}) { # CGC name no longer in use ...
	my $gene = $prev_cgc2gene->{$cgc};
	if (exists $geneid2cgc->{$gene}) { # ... but the Gene ID is still attached to a new CGC
	  my $new_cgc = $geneid2cgc->{$gene};
	  $all_cgcs{$cgc} = 0; # don't process this one again
	  $all_cgcs{$new_cgc} = 0; # don't process this one again
	  $new = "$new_cgc\t$cgc2gene->{$new_cgc}\t$sequence_name->{$new_cgc}\t$cgc2biotype->{$new_cgc}";
	  push @changed_cgc_name_on_a_gene, "$prev\t$new\n";
	} 
      }
    }
  }

  foreach my $cgc (keys %all_cgcs) {
    if (!$all_cgcs{$cgc}) {next} # already dealt with this one
    my $prev = (exists $prev_cgc2gene->{$cgc}) ? "$cgc\t$prev_cgc2gene->{$cgc}\t$prev_sequence_name->{$cgc}\t$prev_cgc2biotype->{$cgc}" : ".\t.\t.\t.";
    my $new = (exists $cgc2gene->{$cgc}) ? "$cgc\t$cgc2gene->{$cgc}\t$sequence_name->{$cgc}\t$cgc2biotype->{$cgc}" : ".\t.\t.\t.";
    if ($prev ne $new) {
      if (!exists $prev_cgc2gene->{$cgc}) {push @no_previous_cgc_name, "$prev\t$new\n"}
      elsif (!exists $cgc2gene->{$cgc}) {push @no_new_cgc_name, "$prev\t$new\n"}
      elsif ($prev_cgc2gene->{$cgc} ne $cgc2gene->{$cgc}) {push @geneid_different_to_prev, "$prev\t$new\n"}
      elsif ($prev_cgc2biotype->{$cgc} ne $cgc2biotype->{$cgc}) {push @biotype_changed, "$prev\t$new\n"}
      else {push @other, "$prev\t$new\n"} # can have a change in the Sequence name of a Gene
    }
  }

  print OUT "# Genes with changed CGC names\n";
  print OUT "# old CGC\told GeneID\told SequenceID\told Biotype\tnew CGC\tnew GeneID\tnew SequenceID\tnew Biotype\n";
  print OUT @changed_cgc_name_on_a_gene;
  print OUT "\n";

  print OUT "# new CGC names\n";
  print OUT "# old CGC\told GeneID\told SequenceID\told Biotype\tnew CGC\tnew GeneID\tnew SequenceID\tnew Biotype\n";
  print OUT @no_previous_cgc_name;
  print OUT "\n";
  
  print OUT "# removed CGC names\n";
  print OUT "# old CGC\told GeneID\told SequenceID\told Biotype\tnew CGC\tnew GeneID\tnew SequenceID\tnew Biotype\n";
  print OUT @no_new_cgc_name;
  print OUT "\n";
  
  print OUT "# assigned Gene changed\n";
  print OUT "# old CGC\told GeneID\told SequenceID\told Biotype\tnew CGC\tnew GeneID\tnew SequenceID\tnew Biotype\n";
  print OUT @geneid_different_to_prev;
  print OUT "\n";
  
  print OUT "# Biotype changed\n";
  print OUT "# old CGC\told GeneID\told SequenceID\told Biotype\tnew CGC\tnew GeneID\tnew SequenceID\tnew Biotype\n";
  print OUT @biotype_changed;
  print OUT "\n";
  
  print OUT "# other changes\n";
  print OUT "# old CGC\told GeneID\told SequenceID\told Biotype\tnew CGC\tnew GeneID\tnew SequenceID\tnew Biotype\n";
  print OUT @other;
  print OUT "\n";
  

    
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
