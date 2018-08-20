#!/usr/bin/env perl
#
# load_panther_xrefs
# 
# A script to load xrefs to panther. The source file is expected to live here:
#
# BUILD_DATA/MISC_STATIC/orthologs/panther/[species]
# 
# C. elegans file needs to be downloaded from: 
#
# ftp://ftp.pantherdb.org/sequence_classifications/12.0/PANTHER_Sequence_Classification_files/PTHR12.0_nematode_worm
#
#
# and corresponds to a specific WormBase freeze. Since WormBase will have moved
# on since then, we take care not to add xrefs to objects that are no longer live
#
#==================

use strict;
use Carp;
use Getopt::Long;

use lib $ENV{CVS_DIR};
use Wormbase;
use Log_files;
use Storable;

##############################
# command-line options       #
##############################

my ($test,$debug,$wb,$store,$noload,$database,$pantherfile,$acefile,$species);

GetOptions (
  "debug:s"        => \$debug,
  "store=s"        => \$store,
  "test"           => \$test, 
  "database:s"     => \$database,
  "acefile=s"      => \$acefile,
  "pantherfile=s"  => \$pantherfile,
  "species=s"      => \$species,
  "noload"         => \$noload,
	   );

if ($store) {
  $wb = retrieve($store) or croak ("Can't restore wormbase from $store\n");
}
else {
  $wb = Wormbase->new( -debug => $debug,
                       -test => $test,
                       -organism => $species,
      );
}

my $log = Log_files->make_build_log($wb);

$database = $wb->autoace if not defined $database;
$acefile = $wb->acefiles . "/panther_xrefs.ace" if not defined $acefile;
$pantherfile = $wb->misc_static . "/ortholog_data/panther/sequence_classifications." . $wb->species . ".txt" if not defined $pantherfile;

$log->log_and_die("Could not find panther file $pantherfile\n") if not -e $pantherfile;

my (%genelookup, %wbgene2cds);

$wb->FetchData('wbgene_id2cds', \%wbgene2cds, $wb->common_data);
foreach my $wbgeneid (keys %wbgene2cds) {
  $genelookup{$wbgeneid} = $wbgeneid;
  foreach my $cds (@{$wbgene2cds{$wbgeneid}}) {
    my $seq_name = $cds;
    $seq_name =~ s/[a-z]$//; 

    $genelookup{$seq_name} = $wbgeneid;
    $genelookup{"CELE_${seq_name}"} = $wbgeneid;
  }
}


my $def = &get_tm_def($wb->full_name);
my $tm_query = $wb->table_maker_query($database,$def);
while(<$tm_query>) {
  chomp;
  next if (/>/ or /\/\//);
  s/\"//g;

  my ($gene_id,$uniprot_gcrp) = split(/\t/, $_);

  $genelookup{$uniprot_gcrp} = $gene_id;
}
unlink $def;

open(my $outfh, ">$acefile") or $log->log_and_die("Could not open $acefile for writing\n");

$log->write_to("Reading Panther file...\n");
open(my $fh, $pantherfile) or $log->log_and_die("Could not open $pantherfile for reading\n");
while(<$fh>) {
  my ($panther_gid, $panther_fid) = /^(\S+)\s+(\S+):\S+/; 
  my $wbgene;
  
  if ($panther_gid =~ /WormBase=(WBGene\d+)/) {
    # easy case; check
    $wbgene = $genelookup{$1} if exists $genelookup{$1};
  } elsif ($panther_gid =~ /WormBase=([^\|]+)/) {
    my $cds = $1;
    $cds =~ s/[a-z]$//; 
    $wbgene = $genelookup{$cds} if exists $genelookup{$cds};
  } elsif ($panther_gid =~ /Gene=([^\|]+)/) {
    $wbgene = $genelookup{$1} if exists $genelookup{$1};
  } elsif ($panther_gid =~ /Gene_ORFName=([^\|]+)/) {
    $wbgene = $genelookup{$1} if exists $genelookup{$1};
  } elsif ($panther_gid =~ /UniProtKB=(\S+)/) {
    $wbgene = $genelookup{$1} if exists $genelookup{$1};
  }

  if (not defined $wbgene) {
    $log->write_to("Could not find WB id for $panther_gid\n");
  } else {
    # write the cross-reference here
    print $outfh "\nGene : \"$wbgene\"\n";
    print $outfh "Database \"Panther\" \"gene\" \"$panther_gid\"\n";
    print $outfh "Database \"Panther\" \"family\" \"$panther_fid\"\n";
  }
}

close($outfh) or $log->log_and_die("Could not close $acefile\n");

if (not $noload) {
  $log->write_to("Loading results...\n");
  $wb->load_to_database($database,$acefile,'load_refseq_xrefs', $log);
}

$log->mail();
exit(0);


##########################################
sub get_tm_def {

  my ($spe) = @_;
  
  my $tm_def = "/tmp/load_refseq_xrefs.$$.tm_def";
  open(my $tfh, ">$tm_def") 
      or $log->log_and_die("Could not open $tm_def for writing\n");
  
  my $query =  <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Gene 
From 1 
Condition Species = "$spe"
 
Colonne 2 
Width 12 
Mandatory 
Hidden 
Class 
Class CDS 
From 1 
Tag Corresponding_CDS 
 
Colonne 3 
Width 12 
Mandatory 
Hidden 
Class 
Class Database 
From 2 
Tag Database 
Condition UniProt_GCRP
 
Colonne 4 
Width 12 
Mandatory 
Hidden 
Class 
Class Database_field 
Right_of 3 
Tag  HERE  
 
Colonne 5 
Width 12 
Optional 
Visible 
Class 
Class Text 
Right_of 4 
Tag  HERE  

EOF

  print $tfh $query;
  close($tfh);

  return $tm_def;
}


__END__
