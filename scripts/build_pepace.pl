#!/usr/local/bin/perl5.8.0 -w
#
# build_pepace.pl
#
# by Anthony Rogers
# 2019-01-04 re-written by Gary Williams
#
# This creates an acefile that can be loaded in to an empty
# database to completely recreate what was pepace. This is based
# solely in the wormpep.history file.
#
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-09-29 16:57:22 $

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

my ( $help, $debug, $test, $verbose, $store, $wormbase, $species, $one_protein );

GetOptions(
	   "help"      => \$help,
	   "debug=s"   => \$debug,
	   "test"      => \$test,
	   "verbose"   => \$verbose,
	   "store:s"   => \$store,
	   "species:s" => \$species,
	   "one_protein:s" => \$one_protein, # for specifying a single protein for testing purposes 
	  );

if ($store) {
  $wormbase = retrieve($store)
    or croak("Can't restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new(
			    -debug => $debug,
			    -test  => $test,
			    -organism => $species
			   );
}

# Display help if required
&usage("Help") if ($help);

# establish log file.
my $log = Log_files->make_build_log($wormbase);

##############
# variables  #
##############

my $ace_dir    = $wormbase->autoace;                  # AUTOACE DATABASE DIR
my $wormpepdir = $wormbase->wormpep;                  # CURRENT WORMPEP
my $ver        = $wormbase->get_wormbase_version();
my $PEP_PREFIX = $wormbase->pep_prefix;
my $PEPDIR     = $wormbase->pepdir_prefix;
my $acefile    = "$ace_dir/acefiles/pepace.ace";

my $Live = 9999999; # large number for the release number so it sorts to the end
my $regex = $wormbase->seq_name_regex;

# read history file
our ( $gene, $CE, $in, $out );


# get the sequence from the common data written earlier
my %CE_sequence = $wormbase->FetchData('cds2aa');
my %CDS_cgc = $wormbase->FetchData('cds2cgc');
if (scalar (keys %CDS_cgc) < 1) {$log->log_and_die("No CGC names found in COMMON_DATA/cds2cgc - this looks wrong!\n")}

#######################################################################################################

# read history file
my %protein_history = &read_history_file();

# construct histories of CDSs
my %cds_history = &get_cds_history(\%protein_history);

# process proteins
 # make History tags
 # add further details for Live genes
&process_proteins($acefile, $one_protein, \%protein_history, \%cds_history);

# sanity checks
check($acefile);

#load files in to autoace.
if ($debug) {$log->write_to("In debug mode so pepace will not be loaded into the database\n");}
$wormbase->load_to_database( $wormbase->autoace, "$ace_dir/acefiles/pepace.ace", 'pepace', $log ) if not $debug;

# update common data
$wormbase->run_script("update_Common_data.pl --cds2wormpep", $log) if not $debug;

$log->mail();
exit(0);


#######################################################################################################
# read the current history file into a hash
# the CDSs are held in a list because more than one CDS might be created coding for the same protein in the same range of releases
# $end is 0, $in is 1 so when they are sorted, we see 'remove' events in a release before we see 'add' events

sub read_history_file {

  my %protein_history;

  open( HISTORY, "$wormpepdir/${PEPDIR}pep.history$ver" ) or $log->log_and_die("Cant open wormpep.history$ver $!\n");
  while (<HISTORY>) {
    my @data = split( /\s+/, $_ );
    my ( $cds, $protein, $in, $out ) = @data;
    if (defined $out && $in >= $out) {$log->write_to("IN release is larger than OUT release in wormpep.history$ver\n@data\n");$log->error}
    if (!defined $out) {$out = $Live} # $Live is a large number used to make the 'Live' events sort to the end of the list of events
    push @{$protein_history{$protein}{$in}{1}}, $cds; #  1 = $in  Creation release number
    push @{$protein_history{$protein}{$out}{0}}, $cds; # 0 = $out Removal release number
  }

  return %protein_history;
}

#######################################################################################################
# invert the protein history hash to get the history of the CDS objects' proteins
sub get_cds_history {

  my ($protein_history) = @_;

  my %cds_history;

  foreach my $protein (keys %{$protein_history}) {
    foreach my $release (keys $protein_history->{$protein}) {
      foreach my $end (keys $protein_history->{$protein}{$release}) { # $end is: $in = 1, $out = 0
	my @cds = @{$protein_history->{$protein}{$release}{$end}};
	foreach my $cds (@cds) {
	  $cds_history{$cds}{$release}{$end} = $protein;
	}
      }
    }
  }

  return %cds_history;

}


#######################################################################################################
# process the protein histories to output the History tags
# where the protein is still Live, output the molecular weight, CGC name, Corresponding CDS and other details
sub process_proteins {

  my ($acefile, $one_protein, $protein_history, $cds_history) = @_;
  
  open( ACE, ">$acefile" ) || $log->log_and_die("cant write $acefile\n");
  
  # foreach $protein event
  foreach my $protein (keys %{$protein_history}) {
    
    if (defined $one_protein) {$protein = $one_protein} # for testing purposes

    my @History;
    my %live; # hash keyed by genes which are currently (as of the release being processed) Live
    
    
    # output Species, Wormpep tags
    &output_Species_Wormpep($protein);
    
    # foreach sort by $release, end
    foreach my $release (sort {$a <=> $b} keys $protein_history->{$protein}) {
      foreach my $end (sort {$a <=> $b} keys $protein_history->{$protein}{$release}) { # end is: $out is 0, $in is 1 so when they are sorted, we see 'remove' events in a release before we see 'add' events
	my @cds = @{$protein_history->{$protein}{$release}{$end}};
	foreach my $cds (@cds) {
	  
	  my $converting = 0;
	  
	  # add to History store:
	  # $release, text, $cds, state = 'Add', 'Dead', 'Convert', 'Live'
	  
	  
	  my $convert_type;
	  my $convert_cds;
	  ($convert_type, $convert_cds) = &check_for_convert($cds, $release, \@History);

	  my $state = 'unknown';
	  if ($end == 0 && $release != $Live) {
	    $state = 'Dead';
	  } elsif ($end == 0 && $release == $Live) {
	    $state = 'Live';
	  } elsif ($end == 1 && $convert_type) { # is there a Dead entry that matches the root name of $cds in the previous release?
	    $state = 'Convert';
	    $converting = 1;
	  } elsif ($end == 1) {
	    $state = 'Add';
	  } else {
	    $log->log_and_die("Unrecognised state in pepace - what should I do with Protein $protein, Release $release, End $end, CDS $cds ?\n");
	  }
	  
	  
	  if ($state eq 'Dead') {
	    push @History, [($release, "'removed' '$cds'", $cds, $state)];
	    delete $live{$cds};
	    
	  } elsif ($state eq 'Add') {
	    
	    my $old_cds;

	    if (&any_current_live_genes(\%live) || scalar @History == 0) { # if have live genes or this is the start of the History we "create"
	      push @History, [($release, "'created' '$cds'", $cds, $state)];
	    } else {
	      if (&matches_an_old_gene_name($cds, \@History)) {
		push @History, [($release, "'reappeared' '$cds'", $cds, $state)];
	      } elsif ($old_cds = &matches_an_old_gene_name_as_an_isoform($cds, \@History)) {
#		push @History, [($release, "'reappeared as isoform' '$cds' '$old_cds'", $cds, $state)];  #  the old format was like below, but this line looks better, I think
		push @History, [($release, "'reappeared as isoform $cds' '$old_cds'", $cds, $state)];
	      } elsif (&doesnt_match_an_old_gene_name($cds, \@History)) {
		push @History, [($release, "'reappeared coded by another gene' '$cds'",  $cds, $state)];
	      } else { # ??? do we ever get here ??? - I don't think so, but create it anyway as a last default
		push @History, [($release, "'created' '$cds'", $cds, $state)];
	      }
	    } 
	    $live{$cds} = 1;

	    
	  } elsif ($state eq 'Convert') {
	    if ($convert_type eq 'isoform') {
#	      push @History, [($release, "'converted to isoform' '$cds' '$convert_cds'", $cds, $state)]; #  the old format was like below, but this line looks better, I think
	      push @History, [($release, "'converted to isoform $cds' '$convert_cds'", $cds, $state)];
	    } else {
#	      push @History, [($release, "'converted to root cds name' '$cds' '$convert_cds'", $cds, $state)]; #  the old format was like below, but this line looks better, I think
	      push @History, [($release, "'converted to root cds name $cds' '$convert_cds'", $cds, $state)];
	    }
	    $live{$cds} = 1;
	    
	  } elsif ($state eq 'Live') {
	    push @History, [($release, "Live", $cds, $state)];
	    # the $protein is still coded for by $cds, so add $cds to %live - it should be there already, but better safe than sorry
	    if (!exists $live{$cds}) {$log->write_to("$cds is not Live at the end of ${protein}'s History - very odd - it is being set to Live anyway!\n")}
	    $live{$cds} = 1;
	    
	  } else {
	    $log->log_and_die("There is an unknown state: $state in processing $cds for $protein - very odd!\n");
	  }
	  
	  # amend previous state = 'Dead' entries in the History store where $release = previous release
	  if ($converting) {
	    #   change state 'Dead' to 'Ignore' if  "converted to isoform"
	    #   change state 'Dead' to 'Ignore' if  "converted to gene name"
	    &amend_history_convert($release, \@History, $convert_cds);
	  } else {
	    #   change text "remove" to "replaced by $protein $other_gene"
	    &check_if_replacing($cds, $release, \@History, $cds_history);
	  }
	  
	}
      }
    }	  # end $release loop
    
    # output History store for $protein
    # don't print anything for state='Live'
    &output_History($protein, \@History);
    
    # if protein is still Live, output further details: Corresponding_CDS for each %live, molecular weight, CGC name, Corresponding CDS and other details
    &output_further_details($protein, \%live);
    
    
    if ($one_protein) {last} # for testing puposes
  }    # end $protein loop


  close(ACE);

}

#######################################################################################################
# output all tags that are required for Live genes
sub output_further_details {
  my ($protein, $live) = @_;
  
  if (!any_current_live_genes($live)) {return}
  
  print ACE "Live\n";
  my %gnames;
  my $example_cds;
  foreach my $cds (keys %{$live}) {
    $example_cds = $cds; # we only need one example of a CDS name for this protein to get the sequence and so mol_weight

    
    print ACE "Corresponding_CDS \"$cds\"\n";
    
    my ($stem, $iso) = ($cds =~ /($regex)(\w*)/);
    
    my $gname = (exists $CDS_cgc{$cds}) ? uc($CDS_cgc{$cds}) : $stem;
    $gname=~s/([A-Z])([A-Z]{2})(-\S+-\S+)/$1\L$2\U$3/;
    $gname .= ", isoform $iso" if $iso;
    $gnames{$gname} = 1;
  }
  foreach my $gname (keys %gnames) {
    print ACE "Gene_name \"$gname\"\n";
  }
  
  my $seq = $CE_sequence{$example_cds};
  if (!defined $seq) {$log->log_and_die("No sequence found in COMMON_DATA/cds2aa for $example_cds\n");}
  my $mw = &get_mol_weight($seq);

  print ACE "Molecular_weight $mw Inferred_automatically \"build_pepace.pl\"\n";
  print ACE "\nPeptide : \"$protein\"\n";
  print ACE "$seq\n";
  
}

#######################################################################################################
# the History elements contain a list of: ($release, "text to output", $cds, $state)
sub output_History {
  my ($protein, $History) = @_;

  foreach my $data (@{$History}) {
    if ($data->[3] eq 'Live') {next} # still being Live is not an event
    if ($data->[3] eq 'Ignore') {next} # being Ignored (not Deleted) because it is in the process of being Converted into an isoform is not an event
    $data->[1] =~ s/\'/"/g; # change ' to " 
    print ACE "History \"" . $data->[0] . "\" " . $data->[1] . "\n";
  }
}

#######################################################################################################
# change a Dead entry to an Ignore entry because it is being converted to an isoform
# the History elements contain a list of: ($release, "text to output", $cds, $state)
sub amend_history_convert {
  my ($release, $History, $convert_cds) = @_;

  foreach my $data (@{$History}) {
    if ($data->[0] == $release && $data->[2] eq $convert_cds && $data->[3] eq 'Dead') {
      $data->[3] = 'Ignore';
    }
  }  
}

#######################################################################################################
# in a Dead entry, change text "remove" to "replaced by $protein $other_gene" if the new protein can be found
sub check_if_replacing {
  my ($cds, $release, $History, $cds_history) = @_;

  foreach my $data (@{$History}) {
    if ($data->[0] == $release && $data->[2] eq $cds && $data->[3] eq 'Dead') {
      if (exists $cds_history->{$cds}{$release}{1}) {
	my $new_protein = $cds_history->{$cds}{$release}{1};
	$data->[1] = "'replaced by $new_protein' '$cds'";
      }
    }
  }  

}
#######################################################################################################
# where $cds is currently Dead, check to see if its root sequence name has been seen before to code for this protein
# returns the matching name
# the History elements contain a list of: ($release, "text to output", $cds, $state)
sub matches_an_old_gene_name_as_an_isoform {
  my ($cds, $History) = @_;

  if (oldStyleName($cds)) {return ''}
  if ($cds !~ m/($regex)\w/) {return ''}
  my $root = $1;
	
  return &matches_an_old_gene_name($root, $History);
}

#######################################################################################################
# where $cds is currently Dead, check to see if its exact name has been seen before to code for this protein
# returns the matching name
# the History elements contain a list of: ($release, "text to output", $cds, $state)
sub matches_an_old_gene_name {
  my ($cds, $History) = @_;

  foreach my $data (@{$History}) {
    if ($data->[2] eq $cds) {return $data->[2]}
  }

  return '';
}

#######################################################################################################
# where $cds is currently Dead, check to see if its exact name has not been seen before to code for this protein
# returns true if it finds a match
# the History elements contain a list of: ($release, "text to output", $cds, $state)
sub doesnt_match_an_old_gene_name {
  my ($cds, $History) = @_;
  
  return ! matches_an_old_gene_name($cds, $History);
}

#######################################################################################################
# true if there are currently (as of the end of the previous release) any genes which are Live  
# we do not include genes which are Live in this current release because we sort @History by (Release, End) and we code $out releases as End='0', so 'remove' events are seen before 'add' events  so we don't see them now.
sub any_current_live_genes {
  my ($live) = @_;
  return scalar keys %{$live};
}
#######################################################################################################
# is there a Dead entry that matches the root name of $cds in the previous release?
# return "isoform" or "root" and the old name of the CDS - or null if not found.
sub check_for_convert {
  my ($cds, $release, $History) = @_;

  if (oldStyleName($cds)) {return ''}
  $cds =~ m/($regex)\w/; 
  my $root_cds = $1;
  foreach my $data (@{$History}) {
    if ($data->[0] == $release) {
      if ($data->[3] eq 'Dead') {
	$data->[2] =~ m/($regex)\w/;
	my $old_root_cds = $1;
	
	if (defined $root_cds && $root_cds eq $data->[2]) {
	  return ("isoform", $data->[2]); # the isoform converted back to the root matches the old cds
	} elsif (defined $old_root_cds && $cds eq $old_root_cds) {
	  return ("root", $data->[2]); # name converted back to the root of the old isoform matches
	}
      }
    }
  }
  return '';
}

#######################################################################################################
# the minimal required tags for a Protein object

sub output_Species_Wormpep {
  my ($protein, $species) = @_;
  print ACE "\n";
  print ACE "Protein : \"$protein\"\n";
  print ACE "Species \"".$wormbase->full_name."\"\n";
  if (!defined $species || $species eq 'elegans') {
    print ACE "Wormpep\n"; # this is a tag only used for indicating elegans proteins - used by the script changed_swissprot_proteins.pl
  }

}


#######################################################################################################
# do sanity checks
sub check {
  my ($acefile) = @_;

  my $live_peps  = `grep -c Live $acefile`;
  chomp $live_peps;
  my $table_peps = &countUniquePeptides("$wormpepdir/${PEPDIR}pep$ver.pep");
  
  $log->write_to("This file has $live_peps live peptides\n");
  $log->write_to("$wormpepdir/${PEPDIR}pep$ver suggests there should be $table_peps\n");
  
  if ( ($live_peps) == $table_peps ) {
    $log->write_to("\nso thats OK!\n\n");
  } else {
    $log->write_to("\nERROR: Investigate this!\n\n");
    $log->error;
  }
}

#######################################################################################################
sub countUniquePeptides {
  my ($fa) = @_;

  my (%peps, %unique_peps, $name);
  
  open(my $fah, $fa) or die "Could not open $fa for reading\n";
  while(<$fah>) {
    /^\>(\S+)/ and do {
      $name = $1;
      next
    };
    
    /^(\S+)/ and do {
      $peps{$name} .= $1;
    }
  }
  
  foreach my $pep (values %peps) {
    $unique_peps{$pep}++;
  }
  
  return scalar(keys %unique_peps);
}

#######################################################################################################



##########################################

sub usage {
  my $error = shift;
  
  if ( $error eq "Help" ) {
    # Normal help menu
    system( 'perldoc', $0 );
    exit(0);
  }
}

##########################################

sub get_mol_weight {
  my $pep = shift;



  my %mw = (
	    'A', '71.0788',  'R', '156.1876', 'D', '115.0886', 'N', '114.1039',
	    'C', '103.1448', 'E', '129.1155', 'Q', '128.1308', 'G', '57.0520',
	    'H', '137.1412', 'I', '113.1595', 'L', '113.1595', 'K', '128.1742',
	    'M', '131.1986', 'F', '147.1766', 'P', '97.1167',  'S', '87.0782',
	    'T', '101.1051', 'W', '186.2133', 'Y', '163.1760', 'V', '99.1326',
	    'U', '150.050'
	   );
  
  #Amino acids
  my $A = "";
  my $R = "";
  my $D = "";
  my $N = "";
  my $C = "";
  my $E = "";
  my $Q = "";
  my $G = "";
  my $H = "";
  my $I = "";
  my $L = "";
  my $K = "";
  my $M = "";
  my $F = "";
  my $P = "";
  my $S = "";
  my $T = "";
  my $W = "";
  my $Y = "";
  my $V = "";
  my $U = "";

  $A = $pep =~ tr/A/A/;    # count the number of each amino acids in peptide.
  $R = $pep =~ tr/R/R/;
  $D = $pep =~ tr/D/D/;
  $N = $pep =~ tr/N/N/;
  $C = $pep =~ tr/C/C/;
  $E = $pep =~ tr/E/E/;
  $Q = $pep =~ tr/Q/Q/;
  $G = $pep =~ tr/G/G/;
  $H = $pep =~ tr/H/H/;
  $I = $pep =~ tr/I/I/;
  $L = $pep =~ tr/L/L/;
  $K = $pep =~ tr/K/K/;
  $M = $pep =~ tr/M/M/;
  $F = $pep =~ tr/F/F/;
  $P = $pep =~ tr/P/P/;
  $S = $pep =~ tr/S/S/;
  $T = $pep =~ tr/T/T/;
  $W = $pep =~ tr/W/W/;
  $Y = $pep =~ tr/Y/Y/;
  $V = $pep =~ tr/V/V/;
  $U = $pep =~ tr/U/U/;
  
  #Calculate the Total Mw of the peptide by summing the subunits.
  my $sum =
    ( ( $A * $mw{A} ) +
      ( $R * $mw{R} ) +
      ( $D * $mw{D} ) +
      ( $N * $mw{N} ) +
      ( $C * $mw{C} ) +
      ( $E * $mw{E} ) +
      ( $Q * $mw{Q} ) +
      ( $G * $mw{G} ) +
      ( $H * $mw{H} ) +
      ( $I * $mw{I} ) +
      ( $L * $mw{L} ) +
      ( $K * $mw{K} ) +
      ( $M * $mw{M} ) +
      ( $F * $mw{F} ) +
      ( $P * $mw{P} ) +
      ( $S * $mw{S} ) +
      ( $T * $mw{T} ) +
      ( $W * $mw{W} ) +
      ( $Y * $mw{Y} ) +
      ( $V * $mw{V} ) +
      ( $U * $mw{U} ) ) / 1000;
  my $result = sprintf "%.1f", $sum;
  return $result;
}

#######################################################################################################
# sort helper routine
sub byRelease {
  # used by sort
  $a <=> $b;
}

#######################################################################################################
sub oldStyleName {
  my ($cds) = @_;
  if (( !defined $species || $species eq 'elegans') && $cds =~ m/\w+\.\p{IsAlpha}+/ ) {
    return 1;
  } else {
    return 0;
  }
}


#######################################################################################################
__END__


=pod

=head2 NAME - build_pepace.pl

=head1 USAGE 

=over 4

=item build_pepace.pl

=back

This creates an acefile that can be loaded in to an empty 
database to completely recreate what was pepace. This is based
solely in the wormpep.history file.

build_pepace.pl MANDATORY arguments:

=over 4

=item none

=back

build_pepace.pl  OPTIONAL arguments:

=over 4

=item none

=back

=head1 REQUIREMENTS

=over 4

=item None known.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
