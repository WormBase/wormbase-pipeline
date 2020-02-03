#!/software/bin/perl -w
#
# recent_changes.pl
# 
# by Gary Williams                  
#
# This is to read data about recent changes from the Nameserver and to
# output an ACE format file to update the geneace ACeDB database.
#
# Last updated by: $Author: pad $     
# Last updated on: $Date: 2013-08-14 12:19:59 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use lib '/software/worm/lib/perl';
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use FileHandle;          # Or `IO::Handle' or `IO::'-anything-else used for autoflush

use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
use NameDB_handler;


=pod

=head recent_changes.pl

=item Options:

  -load     load the output ACE file into geneace
  -from     start date to read data from the Nameserver
  -until    end date to read data from the Nameserver
  -outfile  name of the output ACE file - defaults to ./recent_changes_$from_$until.ace
  -debug    limits reports to specified user <Optional>
  -help     some help


e.g. perl recent_changes.pl -outfile results_file.ace -load

=cut




######################################
# variables and command-line options # 
######################################

my ($help, $debug, $verbose, $store, $wormbase, $geneace);
my ($outfile, $load, $from, $until);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "outfile:s"  => \$outfile,
	    "load"       => \$load,
	    "from:s"     => \$from,
	    "until:s"    => \$until,
	    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug);
}

# Display help if required
&usage("Help") if ($help);

# establish log file.
my $log = Log_files->make_build_log($wormbase);

# this stores the current versions of events in 
my %versions;

$from = get_from_date($from);
$until = get_until_date($until);
if ($debug) {
  print "In debug mode - not updating the '-from' date\n";
} else {
  write_from_date($until);
} 
print "Extracting data from dates from $from until $until, inclusive\n";

if (!defined $outfile) {$outfile = "recent_changes_${from}_to_${until}.ace";}

$geneace = Ace->connect('-path' => $wormbase->database('geneace')) or $log->log_and_die("Failed to connect to geneace\n");

open (OUT, ">$outfile") || $log->log_and_die("Can't open file $outfile");
print OUT "\n\n// Nameserver Recent changes\nFrom: $from\nUntil: $until\n\n\n";

my $db = NameDB_handler->new($wormbase);

# we wish to pull out recent activity on the *web* interface in order to read it into geneace in the correct order of events.
my $gene_data = $db->recent_gene($from, $until, 'web');
my $variation_data = $db->recent_variation($from, $until, 'web');
my $strain_data = $db->recent_strain($from, $until, 'web');

process_gene_data($gene_data);
process_variation_data($variation_data);
process_strain_data($strain_data);

close(OUT);
$db->close;

if ($load) {
  $wormbase->load_to_database($geneace, $outfile, "recent_changes.pl", $log)
} else {
  $log->write_to("Output file '$outfile' has not been loaded into geneace.\n");
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

sub get_from_date {
  my ($from) = @_;
  
  if (!defined $from) {

    # get the from date from the storage file
    my $dir = $wormbase->wormpub;
    my $storage = "$dir/.date_of_last_query_of_nameserver_for_recent_change.dat";
    open (DAT, "<$storage") || $log->log_and_die("Can't open $storage\n");
    $from = <DAT>;
    chomp $from;
    close(DAT);

  }
  if (check_date($from)) {$log->log_and_die("There is a problem with the format of the -from date: $from\nIt should be in ISO format e.g. '2019-06-07T15:04:15'\n")}
  return $from
}

##########################################
# if the date is not specified, then assume that we want today's date and we want to use this as the $from date next time
sub get_until_date {
  my ($until) = @_;

  if (!defined $until) {
    $until = `date --utc +%Y-%m-%dT%H:%M:%SZ`; # '2019-06-07T15:04:15Z'
    chomp $until;
  }

  if (check_date($until)) {$log->log_and_die("There is a problem with the format of the -until date: $until\n")}
  return $until;
}

##########################################

sub write_from_date {
  my ($until) = @_;
  
  my $dir = $wormbase->wormpub;
  my $storage = "$dir/.date_of_last_query_of_nameserver_for_recent_change.dat";
  open (DAT, ">$storage") || $log->log_and_die("Can't open $storage\n");
  print DAT $until;
  close(DAT);
  my $mode = 0664;   
  chmod $mode, $storage;
  
}

##########################################
# date should match ISO 8601 extended notation '2019-06-07' or '2019-06-07 15:04:15' or '2019-06-07T15:04:15Z'
sub check_date {
  my ($date) = @_;
  if ($date !~ /^[12]\d{3}\-[01]\d\-[0123]\d$/ && $date !~ /^[12]\d{3}\-[01]\d\-[0123]\d[T\s][012]\d\:[0-5]\d\:[0-5]\d$/) {return 1}
  return 0;
}

##########################################

sub help_authentication {
  $db->print_authentication_instructions();
}

##########################################

sub process_gene_data {
  my ($gene_data) = @_;
  
  my $from = $gene_data->{'from'};
  my $until = $gene_data->{'until'};
  foreach my $activity (@{$gene_data->{activities}}) {
    my $who = $activity->{'who'}{'id'}; # was {'provenance/who'}{'person/id'}
    my $when = $activity->{'when'}; # was {'provenance/when'}
    my $why = $activity->{'why'}; # was {'provenance/why'}
    if (!defined $why) {$why = ''}
    $why =~ s/\"//g; # remove quote marks
    $why =~ s/\n/ /g; # remove newlines
    my $what = $activity->{'what'}; # was {'provenance/what'}
    my $id = $activity->{'id'}; # was {'gene/id'}
    my %data;
    foreach my $thing (@{$activity->{'changes'}}) {
      my $attr = $thing->{'attr'};
      my $value = $thing->{value};
      # in events like an update both the new name and the old name
      # are provided with the same key - the first one in the array we
      # are iterating over appears to always be the new value
      if (!exists $data{$attr}) { 
	$data{$attr} = $value;
      }
    }
    print "gene what = $what ",%data,"\n";
    if ($what eq 'new-gene') { # was  'event/new-gene'
      new_gene($id, $when, $why, $who, \%data);
    } elsif ($what eq 'new-unnamed-gene') { # was 'event/new-unnamed-gene'
      new_unnamed_gene($id, $when, $why, $who, \%data);
    } elsif ($what eq 'kill-gene') { # was 'event/kill-gene'
      kill_gene($id, $when, $why, $who, \%data); # requires testing - do we have to allow for Pseudogenes/Transposons to be Suppressed instead of Dieing?
    } elsif ($what eq 'update-gene') { # was 'event/update-gene'
      update_gene($id, $when, $why, $who, \%data);
    } elsif ($what eq 'resurrect-gene') { # was 'event/resurrect-gene'
      resurrect_gene($id, $when, $why, $who, \%data);
    } elsif ($what eq 'merge-genes') { # was 'event/merge-genes'
      merge_genes($id, $when, $why, $who, \%data);
    } elsif ($what eq 'undo-merge-genes') { # was 'event/undo-merge-genes'
      undo_merge_genes($id, $when, $why, $who, \%data);
    } elsif ($what eq 'split-gene') { # was 'event/split-gene'
      split_gene($id, $when, $why, $who, \%data);
    } elsif ($what eq 'undo-split-gene') { # was 'event/undo-split-gene'
      undo_split_genes($id, $when, $why, $who, \%data);
    } elsif ($what eq 'suppress-gene') { # was 'event/suppress-gene'
      suppress_gene($id, $when, $why, $who, \%data);
    } elsif ($what eq 'remove-cgc-names') { # was 'event/remove-cgc-names'
      remove_cgc_names($id, $when, $why, $who, \%data);
    }
  }
  
  
}

##########################################

sub process_variation_data {
  my ($variation_data) = @_;
  
  my $from = $variation_data->{'from'};
  my $until = $variation_data->{'until'};
  foreach my $activity (@{$variation_data->{activities}}) {
    my $who = $activity->{'who'}{'id'}; # was {'provenance/who'}{'person/id'}
    my $when = $activity->{'when'}; # was {'provenance/when'}
    my $why = $activity->{'why'}; # was {'provenance/why'}
    if (!defined $why) {$why = ''}
    $why =~ s/\"//g; # remove quote marks
    $why =~ s/\n/ /g; # remove newlines
    my $what = $activity->{'what'}; # was {'provenance/what'}
    my $id = $activity->{'id'}; # was {'variation/id'}
    my %data;
    foreach my $thing (@{$activity->{'changes'}}) {
      my $attr = $thing->{'attr'};
      my $value = $thing->{value};
      # in events like an update both the new name and the old name
      # are provided with the same key - the first one in the array we
      # are iterating over appears to always be the new value
      if (!exists $data{$attr}) { 
	$data{$attr} = $value;
      }
    }
    #print "variation what = $what ",%data,"\n";
    if ($what eq 'new-variation') { # was 'event/new-variation'
      new_variation($id, \%data, $when, $why, $who);
    } elsif ($what eq 'kill-variation') { # was  'event/kill-variation'
      kill_variation($id, \%data, $when, $why, $who);
    } elsif ($what eq 'update-variation') { # was 'event/update-variation'
      update_variation($id, \%data, $when, $why, $who);
    } elsif ($what eq 'resurrect-variation') { # was 'event/resurrect-variation'
      resurrect_variation($id, \%data, $when, $why, $who);
    }    
    
  }
  
}

##########################################
sub process_strain_data {
  my ($strain_data) = @_;
  
  my $from = $strain_data->{'from'};
  my $until = $strain_data->{'until'};
  foreach my $activity (@{$strain_data->{activities}}) {
    my $who = $activity->{'who'}{'id'}; # was {'provenance/who'}{'person/id'}
    my $when = $activity->{'when'}; # was {'provenance/when'}
    my $why = $activity->{'why'}; # was {'provenance/why'}
    if (!defined $why) {$why = ''}
    $why =~ s/\"//g; # remove quote marks
    $why =~ s/\n/ /g; # remove newlines
    my $what = $activity->{'what'}; # was {'provenance/what'}
    my $id = $activity->{'id'}; # was {'strain/id'}
    my %data;
    foreach my $thing (@{$activity->{'changes'}}) {
      my $attr = $thing->{'attr'};
      my $value = $thing->{value};
      # in events like an update both the new name and the old name
      # are provided with the same key - the first one in the array we
      # are iterating over appears to always be the new value
      if (!exists $data{$attr}) { 
	$data{$attr} = $value;
      }
    }
    #print "strain what = $what ",%data,"\n";
    if ($what eq 'new-strain') { # was 'event/new-strain'
      new_strain($id, \%data, $when, $why, $who);
    } elsif ($what eq 'kill-strain') { # was  'event/kill-strain'
      kill_strain($id, \%data, $when, $why, $who);
    } elsif ($what eq 'update-strain') { # was 'event/update-strain'
      update_strain($id, \%data, $when, $why, $who);
    } elsif ($what eq 'resurrect-strain') { # was 'event/resurrect-strain'
      resurrect_strain($id, \%data, $when, $why, $who);
    }    
    
  }
  
}

##########################################

sub new_gene {
  my ($id, $when, $why, $who, $data) = @_;
  if ($debug) {print "New Gene: ID: $id, When: $when, Who: $who\n";}
  
  print OUT "\n";
  print OUT "// new_gene\n";
  print OUT "Gene : $id\n";
  print OUT "Remark \"[$when $who] Gene Created: $why\" Curator_confirmed $who\n";
  print OUT "CGC_name $data->{'cgc-name'}\n" if (exists $data->{'cgc-name'}); # was 'gene/cgc-name'
  print OUT "Sequence_name $data->{'sequence-name'}\n" if (exists $data->{'sequence-name'}); # was 'gene/sequence-name'
  print OUT "Public_name $data->{'cgc-name'}\n" if (exists $data->{'cgc-name'}); # was 'gene/cgc-name'
  print OUT "Public_name $data->{'sequence-name'}\n" if (exists $data->{'sequence-name'} && !exists $data->{'cgc-name'}); # was 'gene/sequence-name' 'gene/cgc-name'
  print OUT "Live\n";
  print OUT "Method Gene\n";
  $when =~ s/[T\s]/_/;
  $when =~ s/Z//;
  $when =~ s/\.\d+$//;
  print OUT "Version_change 1 \"$when\" $who Event Created\n";
  print OUT "Version 1\n";
  $versions{$id} = 1;
  print OUT "Species : \"$data->{'species'}\"\n"; # was 'gene/species'
  print OUT "\n";
}

sub new_unnamed_gene { # I think this event type is not used - both cloned and uncloned gene creation gives an event 'new-gene'.
  my ($id, $when, $why, $who, $data) = @_;
  if ($debug) {print "New Unnamed Gene: ID: $id, When: $when, Who: $who\n";}
  
  print OUT "\n";
  print OUT "// new_unnamed_gene\n";
  print OUT "Gene : $id\n";
  print OUT "Remark \"[$when $who] Gene Created: $why\" Curator_confirmed $who\n";
  print OUT "Public_name $data->{'name'}\n"; # was 'gene/name'
  print OUT "\n";
}

sub kill_gene { # We want to suppress Transposon and Pseudogene genes, not kill them - what does the Nameserver do, does it allow killing: Yes!
  my ($id, $when, $why, $who, $data) = @_;
  if ($debug) {print "Kill Gene: ID: $id, When: $when, Who: $who\n";}
  
  print OUT "\n";
  print OUT "// kill_gene\n";
  print OUT "Gene : $id\n";
  print OUT "Remark \"[$when $who] Gene Killed: $why\" Curator_confirmed $who\n";
  print OUT "Dead\n";
  update_event($id, $when, $who, 'Killed');
  print OUT "\n";
  
  my $gene_obj = $geneace->fetch(Gene => $id);
  my $biotype = $gene_obj->Biotype;
  # SO:0000336 pseudogene
  # SO:0000111 transposable_element_gene
  if ($biotype eq 'SO:0000336' || $biotype eq 'SO:0000111') { # Pseudogene or Transposon
    $log->error("Warning kill-genes: $id has been Killed, but it is a Pseudogene or Transposon so it should really have been Suppressed\n");
    return;
  }
}

sub update_gene {
  my ($id, $when, $why, $who, $data) = @_;
  if ($debug) {print "Update Gene: ID: $id, When: $when, Who: $who\n";}
  
  print OUT "\n";
  print OUT "// update_gene\n";
  if (exists $data->{'biotype'}) { # was 'gene/biotype'
    #my $biotype = $data->{'gene/biotype'};
    # We explicitly ignore the biotype change information.
    # The Nameserver Biotype data tends to be less accurate than the
    # camace and geneace data, so don't try to pprocess this.
    print OUT "// The nameserver records a change of Biotype in this gene, but we explicitly ignore the biotype change information in this script.\n";
  }
  print OUT "Gene : $id\n";
  print OUT "Remark \"[$when $who] Gene Updated: $why\" Curator_confirmed $who\n";
  if (exists $data->{'cgc-name'}) { # was 'gene/cgc-name'
    print OUT "CGC_name $data->{'cgc-name'}\n"; # was 'gene/cgc-name'
    print OUT "Public_name $data->{'cgc-name'}\n"; # was 'gene/cgc-name'
    my $class = gene_class($data->{'cgc-name'}); # was 'gene/cgc-name'
    print OUT "Gene_class $class\n";
    name_change($id, $when, $who, 'CGC_name', $data->{'cgc-name'}); # was {'gene/cgc-name'
  }
  if (exists $data->{'sequence-name'}) { # was 'gene/sequence-name'
    print OUT "Sequence_name $data->{'sequence-name'}\n"; # was 'gene/sequence-name'
    print OUT "Public_name $data->{'sequence-name'}\n" if (!exists $data->{'cgc-name'}); # was 'gene/sequence-name' 'gene/cgc-name'
    name_change($id, $when, $who, 'Sequence_name', $data->{'sequence-name'}); # was 'gene/sequence-name'
  }
  if (exists $data->{'species'}) { # was 'gene/species'
    my $species = $data->{'species'}; # was 'gene/species'
    print OUT "Species \"$species\"\n";
  }
  print OUT "\n";
}

sub resurrect_gene {
  my ($id, $when, $why, $who, $data) = @_;
  if ($debug) {print "Resurrect Gene: ID: $id, When: $when, Who: $who\n";}
  
  print OUT "\n";
  print OUT "// resurrect_gene\n";
  print OUT "Gene : $id\n";
  print OUT "Remark \"[$when $who] Gene Resurrected: $why\" Curator_confirmed $who\n";
  print OUT "Live\n";
  my $name;
  my $matches = $db->find_genes($id);
  foreach my $match (@{$matches}) {
    if ($match->{'id'} eq $id) { # was 'gene/id'
      $name = $match->{'sequence-name'} # was 'gene/sequence-name'
    }
  }
  if (defined $name) {
    print OUT "Public_name $name\n";
    print OUT "Sequence_name $name\n";
  }
  update_event($id, $when, $who, 'Resurrected');
  print OUT "\n";
}

sub merge_genes {
  my ($id, $when, $why, $who, $data) = @_;
  if ($debug) {print "Merge Gene: ID: $id, When: $when, Who: $who\n";}
  
  my $deadgene =  $data->{'merges'}; # was 'gene/merges'
  my $deadgene_obj = $geneace->fetch(Gene => $deadgene);
  my $livegene_obj = $geneace->fetch(Gene => $id);
  
  # check the livegene details
  if (!defined $livegene_obj) {$log->log_and_die("ERROR merge-genes: $id is not a valid gene id in geneace")}
  # is this a Live gene with no merge from the DEAD gene in the last Acquires_merge?
  my $status = $livegene_obj->Status->name;
  if ($debug) {print "LIVE GENE: $livegene_obj :  ".$livegene_obj->Status->name."\n";}
  if ($status ne 'Live') {
    $log->log_and_die("ERROR merge-genes: $id is not a Live gene\n");
  }
  # get the last acquires_merge
  my $acquires_merge;
  foreach my $acquires_merge_obj ($livegene_obj->at('Identity.History.Acquires_merge')) {
    if (defined $acquires_merge_obj) {
      ($acquires_merge) = $acquires_merge_obj->row;
    }
  }
  if (defined $acquires_merge && $acquires_merge eq $deadgene) {
    $log->error("Warning merge-genes: $id has a tag saying it has already merged with $deadgene - this merge will not be done again\n");
    return;
  }

  # check the deadgene details
  if (!defined $deadgene_obj) {$log->log_and_die("ERROR merge-genes: $deadgene is not a valid gene id in geneace")}
  # is this a Live gene with no merge into the LIVE gene in the last Merged_into?
  if ($debug) {print "DEAD GENE:$ deadgene :  ".$deadgene_obj->Status->name."\n";}
  $status = $deadgene_obj->Status->name;
  if ($status ne 'Live') {
    $log->log_and_die("ERROR merge-genes: $deadgene is not a Live gene\n");
  }
  # get the last Merged_into tag
  my $merged_into;
  foreach my $merged_into_obj ($deadgene_obj->at('Identity.History.Merged_into')) {
    if (defined $merged_into_obj) {
      ($merged_into) = $merged_into_obj->row;
    }
  }
  if (defined $merged_into && $merged_into eq $deadgene) {
    $log->error("Warning merge-genes: $deadgene has a tag saying it has already merged into $id - this merge will not be done again\n");
    return;
  }

  # get the CGC name
  my $dead_CGC_name = $deadgene_obj->CGC_name;
  if (defined $dead_CGC_name) {
    $log->error("Warning merge-genes: $deadgene has a CGC name ($dead_CGC_name) - this merge needs checking as cgc names are involved. Check this with Tim. and only then load the .ace file.\n");
    if (defined $load) {
      undef $load;
      $log->error("Warning merge-genes: Automatic loading of the ace file has been turned off.\n");
    }
  }

  # output the merge data
  print OUT "\n";
  print OUT "// merge_genes\n";
  print OUT "Gene : $id\n";
  print OUT "Remark \"[$when $who] Gene Merged: $why\" Curator_confirmed $who\n";
  update_event($id, $when, $who, "Acquires_merge $deadgene");
  print OUT "Acquires_merge $deadgene\n";
  print OUT "\n";

  print OUT "\n";
  print OUT "Gene : $deadgene\n";
  print OUT "Remark \"[$when $who] Gene Merged into $id: $why\" Curator_confirmed $who\n";
  update_event($deadgene, $when, $who, "Merged_into $id");
  print OUT "Merged_into $id\n";
  print OUT "Dead\n";
  print OUT "\n";
   
  #Stuff to be removed from the dead gene.
  print OUT "\n";
  print OUT "Gene : $deadgene\n";
  print OUT "-D Map_info\n";
  print OUT "-D Sequence_name\n";
  print OUT "-D Allele\n";
  print OUT "-D Reference\n";
  print OUT "-D Ortholog\n";
  print OUT "-D Paralog\n";
  print OUT "-D Map_info\n";
  print OUT "-D Method\n";
  print OUT "\n";

  # transfer operon connections.
  my $operon_connect = $deadgene_obj->at('Gene_info.Contained_in_operon');
  if (defined $operon_connect) {
    foreach my $operon_connect_list ($operon_connect->col) {
      print OUT "\n";
      print OUT "Gene : $id\n";
      print OUT "Contained_in_operon $operon_connect_list\n";
      print OUT "\n";
    }
  }

  # transfer the Other_names
  my $dead_Other_names_col = $deadgene_obj->at('Identity.Name.Other_name');
  if (defined $dead_Other_names_col) {
    foreach my $dead_Other_names ($dead_Other_names_col->col) {
      print OUT "\n";
      print OUT "Gene : $id\n";
      print OUT "Other_name $dead_Other_names\n";
      print OUT "\n";
    }
  } 


  # transfer Alleles
  foreach my $Alleles ($deadgene_obj->at('Gene_info.Allele')) {
    if (defined $Alleles) {
      print OUT "\n";
      print OUT "Gene : $id\n";
      print OUT "Allele $Alleles\n";
      print OUT "\n";
    }
  }

  # transfer references
  foreach my $references ($deadgene_obj->at('Reference')) {
    if (defined $references) {
      print OUT "\n";
      print OUT "Gene : $id\n";
      print OUT "Reference $references\n";
      print OUT "\n";
    }
  }
  
  # transfer Features
  foreach my $feature ($deadgene_obj->at('Associated_feature')) {
    if (defined $feature) {
      print OUT "\n";
      print OUT "Gene : $id";
      print OUT "Associated_feature $feature\n";
      print OUT "\n";
    }
  }

  # transfer the Ortholog tags
  foreach my $dead_Orthologs ($deadgene_obj->at('Gene_info.Ortholog')) {
    if (defined $dead_Orthologs) {
      my @row = $dead_Orthologs->row;
      my $row0=$row[0]->name;
      unless (defined $row[1]) {$log->log_and_die("\n\nERROR merge-genes: Missing ortholog species data in gene $deadgene\n\n");}
      my $row1=$row[1]->name;
      unless (defined $row[2]) {$log->log_and_die("\n\nERROR merge-genes: Missing ortholog Evidence data in gene $deadgene - $row0($row1)\n\n");}
      my $row2=$row[2]->name;
      my @col = $deadgene_obj->at("Gene_info.Ortholog.$row0.$row1.$row2")->col;
      $row1 = '"' . $row1 . '"'; # add quotes to the species name
      foreach my $col (@col) {
	print OUT "\n";
	print OUT "Gene : $id\n";
	print OUT "Ortholog $row0 $row1 $row2 $col\n";
	print OUT "\n";
      }
    }
  } 

}

sub undo_merge_genes {
  my ($id, $when, $why, $who, $data) = @_;
  if ($debug) {print "Undo Merge Gene: ID: $id, When: $when, Who: $who\n";}
  $log->error("Warning undo-merge-genes: $id should be unmerged - this might be best and most easily done by editing the merge out of the ACE file\n");
  if (defined $load) {
    undef $load;
    $log->error("Warning undo-merge-genes: Automatic loading of the ace file has been turned off.\n");
  }
}

sub split_gene {
  my ($id, $when, $why, $who, $data) = @_;
  if ($debug) {print "Split Gene: ID: $id, When: $when, Who: $who\n";}
  
  # two split-gene events are made for every split.
  # we only want the one with lots of information.
  if (!exists $data->{'id'}) {return} # was 'gene/id'

  my $oldgene = $data->{'splits'}; # was 'gene/splits'
  my $species = $data->{'species'}; # was 'gene/species'
  my $sequence_name = $data->{'sequence-name'}; # was 'gene/sequence-name'
  my $newgene = $data->{'id'}; # was 'gene/id'
  my $biotype = $data->{'cds'}; # was 'biotype/cds'

  # process LIVE
  my $oldgene_obj = $geneace->fetch('Gene', $oldgene);
  if ($oldgene_obj) {
    # is this a Live gene with no splits to the new gene in the last Split_into?
    my $status = $oldgene_obj->Status->name;
    if ($status ne 'Live') {
      $log->log_and_die("ERROR split-gene: $oldgene is not a Live gene\n");
    }
    # get the last acquires_merge
    my $split_into;
    foreach my $split_into_obj ($oldgene_obj->at('Identity.History.Split_into')) {
      if (defined $split_into_obj) {
	($split_into) = $split_into_obj->row;
      }
    }
    if (defined $split_into && $split_into eq $newgene) {
      $log->error("Warning split-gene: $oldgene has a tag saying it has already had $newgene split from it - this split will not be done again\n");
      return;
    }

    print OUT "\n";
    print OUT "// split_gene\n";
    print OUT "Gene : $oldgene\n";
    update_event($oldgene, $when, $who, "Split_into $newgene");
    print OUT "Split_into $newgene\n";
    print OUT "Remark \"[$when $who] Gene Split creating $sequence_name $newgene: $why\" Curator_confirmed $who\n";
    print OUT "\n";

    # process NEW gene
    print OUT "\n";
    print OUT "Gene : $newgene\n";
    print OUT "Sequence_name $sequence_name\n";
    print OUT "Public_name $sequence_name\n";
    print OUT "Species \"$species\"\n";
    update_event($newgene, $when, $who, "Split_from $oldgene");
    print OUT "Split_from $oldgene\n";
    print OUT "Live\n";
    print OUT "Method Gene\n";
    print OUT "Remark \"[$when $who] Gene Split from $oldgene: $why\" Curator_confirmed $who\n";
    print OUT "\n";
  } else {
    $log->log_and_die("ERROR split-gene: $oldgene is not a known gene in geneace\n");
  }
    
}

sub undo_split_gene {
  my ($id, $when, $why, $who, $data) = @_;
  if ($debug) {print "Undo Split Gene: ID: $id, When: $when, Who: $who\n";}

  $log->error("Warning undo-split-genes: $id should be unsplit - this might be best and most easily done by editing the split out of the ACE file\n");
  if (defined $load) {
    undef $load;
    $log->error("Warning undo-split-genes: Automatic loading of the ace file has been turned off.\n");
  }
}

sub suppress_gene {
  my ($id, $when, $why, $who, $data) = @_;
  if ($debug) {print "Supress Gene: ID: $id, When: $when, Who: $who\n";}

  print OUT "\n";
  print OUT "// suppress_gene\n";
  print OUT "Gene : $id\n";
  print OUT "Remark \"[$when $who] Gene Suppressed: $why\" Curator_confirmed $who\n";
  update_event($id, $when, $who, "Suppressed");
  print OUT "Suppressed\n";
  print OUT "\n";
}
sub remove_cgc_names {
  my ($id, $when, $why, $who, $data) = @_;
  if ($debug) {print "Removed CGC Names From Gene: ID: $id, When: $when, Who: $who\n";}

  my $gene_obj = $geneace->fetch(Gene => $id);
  my $sequence_name = $gene_obj->Sequence_name;

  print OUT "\n";
  print OUT "// remove_cgc_gene\n";
  print OUT "Gene : $id\n";
  print OUT "-D CGC_name\n";
  print OUT "-D Gene_class\n";
  print OUT "\n";

  print OUT "\n";
  print OUT "Gene : $id\n";
  print OUT "Public_name \"$sequence_name\"\n" if (defined $sequence_name);
  print OUT "Remark \"[$when $who] Gene CGC Name Removed: $why\" Curator_confirmed $who\n";
  name_change($id, $when, $who, 'CGC_name', ''); # blank CGC name to show it has been removed.
  print OUT "\n";
}

##########################################

sub new_variation {
  my ($id, $data, $when, $why, $who) = @_;
  if ($debug) {print "New Variation: ID: $id\n";}

  print OUT "\n";
  print OUT "// new_variation\n";
  print OUT "Variation : $id\n";
  print OUT "Public_name $data->{'name'}\n"; # was 'variation/name'
  print OUT "Remark \"[$when $who] New Variation: $why\" Curator_confirmed $who\n";
  print OUT "\n";
}

sub kill_variation {
  my ($id, $data, $when, $why, $who) = @_;
  if ($debug) {print "Kill Variation: ID: $id\n";}

  print OUT "\n";
  print OUT "// kill_variation\n";
  print OUT "Variation : $id\n";
  print OUT "Dead\n";
  print OUT "Remark \"[$when $who] Kill Variation: $why\" Curator_confirmed $who\n";
  print OUT "\n";
}

sub update_variation {
  my ($id, $data, $when, $why, $who) = @_;
  if ($debug) {print "Update Variation: ID: $id\n";}

  print OUT "\n";
  print OUT "// update_variation\n";
  print OUT "Variation : $id\n";
  print OUT "Public_name $data->{'name'}\n"; # was 'variation/name'
  print OUT "Remark \"[$when $who] Update Variation: $why\" Curator_confirmed $who\n";
  print OUT "\n";
}

sub resurrect_variation {
  my ($id, $data, $when, $why, $who) = @_;
  if ($debug) {print "Resurrect Variation: ID: $id\n";}

  print OUT "\n";
  print OUT "// ressurect_variation\n";
  print OUT "Variation : $id\n";
  my $name;
  my $matches = $db->find_variations($id);
  foreach my $match (@{$matches}) {
    if ($match->{'id'} eq $id) { # was 'variation/id'
      $name = $match->{'name'} # was 'variation/name'
    }
  }
  print OUT "Public_name $name\n";
  print OUT "Live\n";
  print OUT "Remark \"[$when $who] Ressurect Variation: $why\" Curator_confirmed $who\n";
  print OUT "\n";
}

##########################################

sub new_strain {
  my ($id, $data, $when, $why, $who) = @_;
  if ($debug) {print "New Strain: ID: $id\n";}

  print OUT "\n";
  print OUT "// new_strain\n";
  print OUT "Strain : $id\n";
  print OUT "Public_name $data->{'name'}\n"; # was 'strain/name'
  print OUT "Remark \"[$when $who] New Strain: $why\" Curator_confirmed $who\n";
  print OUT "\n";
}

sub kill_strain {
  my ($id, $data, $when, $why, $who) = @_;
  if ($debug) {print "Kill Strain: ID: $id\n";}

  print OUT "\n";
  print OUT "// kill_strain\n";
  print OUT "Strain : $id\n";
  print OUT "Dead\n";
  print OUT "Remark \"[$when $who] Kill Strain: $why\" Curator_confirmed $who\n";
  print OUT "\n";
}

sub update_strain {
  my ($id, $data, $when, $why, $who) = @_;
  if ($debug) {print "Update Strain: ID: $id\n";}

  print OUT "\n";
  print OUT "// update_strain\n";
  print OUT "Strain : $id\n";
  print OUT "Public_name $data->{'name'}\n"; # was 'strain/name'
  print OUT "Remark \"[$when $who] Update Strain: $why\" Curator_confirmed $who\n";
  print OUT "\n";
}

sub resurrect_strain {
  my ($id, $data, $when, $why, $who) = @_;
  if ($debug) {print "Resurrect Strain: ID: $id\n";}

  print OUT "\n";
  print OUT "// ressurect_strain\n";
  print OUT "Strain : $id\n";
  my $name;
  my $matches = $db->find_strains($id);
  foreach my $match (@{$matches}) {
    if ($match->{'id'} eq $id) { # was 'strain/id'
      $name = $match->{'name'} # was 'strain/name'
    }
  }
  print OUT "Public_name $name\n";
  print OUT "Live\n";
  print OUT "Remark \"[$when $who] Ressurect Strain: $why\" Curator_confirmed $who\n";
  print OUT "\n";
}

##########################################

sub gene_class {
  my ($cgc_name) = @_;
  return $cgc_name =~ /^(\S+)\-/;

}

##########################################

# update an event in a Gene object.  
# If we have seen this gene already in this run then use the current
# version held in %versions

sub update_event {
  my ($id, $when, $who, $event) = @_;

  $when =~ s/[T\s]/_/;
  $when =~ s/Z//;
  $when =~ s/\.\d+$//;
  my $version;
  my $gene_obj = $geneace->fetch(Gene => $id);
  if (defined $gene_obj) {
    if (exists $versions{$id}) {
      $version = $versions{$id};
    } else {
      $version = $gene_obj->Version->name;
    }
  } else {
    # no gene yet, so create the first Event tag entry
    $version = 0;
  }

  $version++;
  $versions{$id} = $version;
  print OUT "Version_change $version \"$when\" $who Event $event\n";
  print OUT "Version $version\n";

}

##########################################

# update a name change in a Gene object.  
# If we have seen this gene already in this run then use the current
# version held in %versions

sub name_change {
  my ($id, $when, $who, $name_type, $name) = @_;

  $when =~ s/[T\s]/_/;
  $when =~ s/Z//;
  $when =~ s/\.\d+$//;
  my $version;
  my $gene_obj = $geneace->fetch(Gene => $id);
  if (exists $versions{$id}) {
    $version = $versions{$id};
  } else {
    $version = $gene_obj->Version->name;
  }
  $version++;
  my $old_cgc_name = $gene_obj->CGC_name;
  if (defined $old_cgc_name && $name_type eq 'CGC_name') {
      print OUT "Version_change $version \"$when\" $who Other_name $old_cgc_name\n";
      print OUT "Other_name $old_cgc_name\n";
      $version++;
  }
  print OUT "Version_change $version \"$when\" $who $name_type $name\n";
  print OUT "Version $version\n";
  $versions{$id} = $version;

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
