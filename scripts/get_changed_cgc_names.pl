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

my ($help, $debug, $test, $verbose, $store, $wormbase, $readable_outfile, $table_outfile, $database);
my ($RELEASE, $PREV_RELEASE);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "readable_outfile:s"  => \$readable_outfile, # human-readable format for people who like pretty outputs
	    "table_outfile:s"  => \$table_outfile, # table format used to present the data in the web site
	    "RELEASE:s"  => \$RELEASE, # optional - for use when looking at historical releases
	    "database:s" => \$database, # ditto
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

my $tace = $wormbase->tace;
$database ||= $wormbase->autoace;

$readable_outfile ||= $wormbase->misc_output . "/readable_changed_CGC_names.dat"; # human-readable format
$table_outfile ||= $wormbase->misc_output . "/changed_CGC_names.dat"; # table format

$RELEASE ||= $wormbase->get_wormbase_version;
$PREV_RELEASE = $RELEASE - 1;
my $releasedir = ($RELEASE == $wormbase->get_wormbase_version) ? $wormbase->common_data :  $wormbase->wormpub . "/DATABASES/WS${RELEASE}/COMMON_DATA/";
my $prevdir = $wormbase->wormpub . "/DATABASES/WS${PREV_RELEASE}/COMMON_DATA/";
my %gene2cgc;
my %prev_gene2cgc;

$wormbase->FetchData('worm_gene2cgc_name', \%gene2cgc, $releasedir);
$wormbase->FetchData('worm_gene2cgc_name', \%prev_gene2cgc, $prevdir);

&remove_other_species(\%gene2cgc);
&remove_other_species(\%prev_gene2cgc);

my ($cgc2gene, $cgc2biotype, $sequence_name, $geneid2cgc) = &key_by_cgc(\%gene2cgc);
my ($prev_cgc2gene, $prev_cgc2biotype, $prev_sequence_name, $prev_geneid2cgc) = &key_by_cgc(\%prev_gene2cgc);

$log->write_to("Connecting to database...\n");
my $db = Ace->connect (-path => "$database",
		       -program => $tace) || $log->log_and_die("Cannot connect to database at $database\n");


open (READABLE, "> $readable_outfile") || $log->log_and_die("Can't open $readable_outfile\n");
open (OUT, "> $table_outfile") || $log->log_and_die("Can't open $table_outfile\n");
&output_changed_cgc();
close OUT;
close READABLE;

$db->close();


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
    if ($hash_ref->{$key} =~ /^\S+\s+\S+\-\S+\-\S+\s+\S+$/) {
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


  # human-readable data
  my @changed_cgc_name_on_a_gene;
  my @no_previous_cgc_name;
  my @no_new_cgc_name;
  my @geneid_different_to_prev;
  my @biotype_changed;
  my @other;

  # table-format data
  my @table_changed_cgc_name_on_a_gene;
  my @table_no_previous_cgc_name;
  my @table_no_new_cgc_name;
  my @table_geneid_different_to_prev;
  my @table_biotype_changed;
  my @table_other;

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
	  my $details = get_cgc_name_changed_details($cgc); # see if the CGC_name is now an Other_name
	  push @changed_cgc_name_on_a_gene, "$prev\t$new\tCGC name changed on this gene$details\n"; # human-readable format
	  push @table_changed_cgc_name_on_a_gene, "$prev\t$new\n"; # table format
	} 
      }
    }
  }

  foreach my $cgc (keys %all_cgcs) {
    if (!$all_cgcs{$cgc}) {next} # already dealt with this one
    my $prev = (exists $prev_cgc2gene->{$cgc}) ? "$cgc\t$prev_cgc2gene->{$cgc}\t$prev_sequence_name->{$cgc}\t$prev_cgc2biotype->{$cgc}" : ".\t.\t.\t.";
    my $new = (exists $cgc2gene->{$cgc}) ? "$cgc\t$cgc2gene->{$cgc}\t$sequence_name->{$cgc}\t$cgc2biotype->{$cgc}" : ".\t.\t.\t.";
    if ($prev ne $new) {

      #######################
      # human readable format
      if (!exists $prev_cgc2gene->{$cgc}) { # new CGC-name - get details
	my $details = get_new_cgc_details($cgc);
	push @no_previous_cgc_name, "$new\t$details\n"
      }
      elsif (!exists $cgc2gene->{$cgc}) {
	my $details = get_retired_cgc_details($cgc); # try to discover why it has retired
	push @no_new_cgc_name, "$prev\t$new\tCGC_name retired$details\n"
      }
      elsif ($prev_cgc2gene->{$cgc} ne $cgc2gene->{$cgc}) {push @geneid_different_to_prev, "$prev\t$new\tCGC_name assigned to a new gene\n"}
      elsif ($prev_cgc2biotype->{$cgc} ne $cgc2biotype->{$cgc}) {push @biotype_changed, "$prev\t$new\tBiotype changed\n"}
      else {push @other, "$prev\t$new\tSequence ID changed\n"} # can have a change in the Sequence name of a Gene

      ##############
      # table format
      if (!exists $prev_cgc2gene->{$cgc}) {push @table_no_previous_cgc_name, "$prev\t$new\n"}
      elsif (!exists $cgc2gene->{$cgc}) {push @table_no_new_cgc_name, "$prev\t$new\n"}
      elsif ($prev_cgc2gene->{$cgc} ne $cgc2gene->{$cgc}) {push @table_geneid_different_to_prev, "$prev\t$new\n"}
      elsif ($prev_cgc2biotype->{$cgc} ne $cgc2biotype->{$cgc}) {push @table_biotype_changed, "$prev\t$new\n"}
      else {push @table_other, "$prev\t$new\n"} # can have a change in the Sequence name of a Gene

    }
  }

  ##############################
  # output human-readable format
  print OUT "# new CGC names\n";
  print OUT "# CGC_name\tGeneID\tSequenceID\tBiotype\tEvidence\n";
  print OUT sort @no_previous_cgc_name;
  print OUT "\n";
  
  print OUT "# Other changes\n";
  # Genes with changed CGC names
  print OUT "# old CGC\told GeneID\told SequenceID\told Biotype\tnew CGC\tnew GeneID\tnew SequenceID\tnew Biotype\n";
  print OUT sort @changed_cgc_name_on_a_gene;
  # removed CGC names
  print OUT sort @no_new_cgc_name;
  # assigned Gene changed
  print OUT sort @geneid_different_to_prev;
  # Biotype changed
  print OUT sort @biotype_changed;
  # other changes
  print OUT sort @other;
  print OUT "\n";

  #####################
  # output table format
  
  print READABLE "# Genes with changed CGC names\n";
  print READABLE "# old CGC\told GeneID\told SequenceID\told Biotype\tnew CGC\tnew GeneID\tnew SequenceID\tnew Biotype\n";
  print READABLE sort @table_changed_cgc_name_on_a_gene;
  print READABLE "\n";
  
  print READABLE "# new CGC names\n";
  print READABLE "# old CGC\told GeneID\told SequenceID\told Biotype\tnew CGC\tnew GeneID\tnew SequenceID\tnew Biotype\n";
  print READABLE sort @table_no_previous_cgc_name;
  print READABLE "\n";
  
  print READABLE "# removed CGC names\n";
  print READABLE "# old CGC\told GeneID\told SequenceID\told Biotype\tnew CGC\tnew GeneID\tnew SequenceID\tnew Biotype\n";
  print READABLE sort @table_no_new_cgc_name;
  print READABLE "\n";
  
  print READABLE "# assigned Gene changed\n";
  print READABLE "# old CGC\told GeneID\told SequenceID\told Biotype\tnew CGC\tnew GeneID\tnew SequenceID\tnew Biotype\n";
  print READABLE sort @table_geneid_different_to_prev;
  print READABLE "\n";
  
  print READABLE "# Biotype changed\n";
  print READABLE "# old CGC\told GeneID\told SequenceID\told Biotype\tnew CGC\tnew GeneID\tnew SequenceID\tnew Biotype\n";
  print READABLE sort @table_biotype_changed;
  print READABLE "\n";
  
  print READABLE "# other changes\n";
  print READABLE "# old CGC\told GeneID\told SequenceID\told Biotype\tnew CGC\tnew GeneID\tnew SequenceID\tnew Biotype\n";
  print READABLE sort @table_other;
  print READABLE "\n";
}

##########################################
# see if the old CGC_name is now an Other_name
sub get_cgc_name_changed_details {
  my ($cgc) = @_; 

  my $details;
  my $cgc_gene = $prev_cgc2gene->{$cgc};
  my $gene_obj = $db->fetch(Gene => $cgc_gene);

  print "CGC_name changed in gene\t$cgc\t$cgc_gene\n";
  my @other_names = $gene_obj->Other_name;
  if (grep /^$cgc$/, @other_names) {
    $details .= ", $cgc has been assigned as an Other_name of $cgc_gene";
  } else {
    $details .= ", $cgc is NOT an Other_name of $cgc_gene"; # something odd happening here
  }

  return $details;
}
##########################################
# get the evidence for the association of a CGC_name with the gene object
sub get_new_cgc_details {
  my ($cgc) = @_;

  my $details;
  my $new_cgc_gene = $cgc2gene->{$cgc};
  my $gene_obj = $db->fetch(Gene => $new_cgc_gene);

  my $cgc_escaped = $cgc;
  $cgc_escaped =~ s/\./\\./g;
  print "New CGC_name\t$cgc\t$new_cgc_gene\n";
  my @evidence = $gene_obj->at("Identity.Name.CGC_name.$cgc_escaped")->col(2);
  if (!scalar @evidence) {
    $details = "No evidence given for this!";
  } else {
    foreach my $evidence (@evidence) {
      if ($evidence =~ /^WBPaper/) {
	my $paper_obj = $db->fetch(Paper => $evidence);
	if (defined $details) {$details .= ', '}
	my $citation .= $paper_obj->Brief_citation;
	$citation =~ s/\s*".+//; # remove the quoted title
	$details .= $citation;
	$details .= ", ($evidence)";
      } elsif ($evidence =~ /^WBPerson/) {
	my $person_obj = $db->fetch(Person => $evidence);
	if (defined $details) {$details .= ', '}
	$details .= $person_obj->Standard_name;

      } else {
	if (defined $details) {$details .= ', '}
	$details .= "$evidence"; # not a paper or a person, just display the ID
      }
    }
  }
  return $details;
}

##########################################
# try to find out why it was retired
sub get_retired_cgc_details {
  my ($cgc) = @_;
  
  my $details = "";
  my $old_cgc_gene = $prev_cgc2gene->{$cgc};
  my $gene_obj = $db->fetch(Gene => $old_cgc_gene);

  # has the gene been merged and the CGC is now an Other_name?
  if ($gene_obj->Status->name eq 'Dead') {
    print "Retired CGC_name in Dead gene\t$cgc\t$old_cgc_gene\n";
    my @merged_genes = $gene_obj->at("Identity.History.Version_change")->col(6);
#     e.g. WBGene00006456 in WS252
    foreach my $merged (@merged_genes) {
      if ($merged =~ /^WBGene/) {
	my $merged_obj = $db->fetch(Gene => $merged);
	my @other_names = $merged_obj->Other_name;
	foreach my $other_name (@other_names) {
	  if ($other_name eq $cgc) {
	    $details .= ", old gene merged in $merged_obj and $cgc has been assigned as an Other_name of $merged_obj";
	  }
	}
      }
    }
  } else {
    $details .= ", $old_cgc_gene is still a Live gene";
    print "Retired CGC_name in Live gene\t$cgc\t$old_cgc_gene\n";
    my @other_names = $gene_obj->Other_name;
    if (grep /^$cgc$/, @other_names) {
      $details .= ", $cgc has been assigned as an Other_name of $old_cgc_gene";
    } else {
      $details .= ", $cgc is NOT an Other_name of $old_cgc_gene"; # something odd happening here
    }
  }

  return $details;
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
