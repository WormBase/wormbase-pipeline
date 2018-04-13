#!/usr/bin/env perl
#
# overload_gff_cds_wormpep.pl
#
# Overloads the CDS and Transcript lines with extra info (mostly wormpep)
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2014-08-27 21:50:10 $

use strict;                                      
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

my ($debug, $test, $store, $wormbase, $species, $database);
my ($infile, $outfile, $gff3, $changed_lines, %already_done_cds);

GetOptions (
  "debug=s"    => \$debug,
  "test"       => \$test,
  "store:s"    => \$store,
  "species:s"  => \$species,
  "infile:s"   => \$infile,
  "outfile:s"  => \$outfile,
  "gff3"       => \$gff3,
  "database:s" => \$database,
    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species
			     );
}

$database = $wormbase->autoace if not defined $database;
my $species_full = $wormbase->full_name;
my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}


my ($tran_status, $tran_wormpep, $tran_locus) = &get_data();

open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  

while (<$gff_in_fh>) {
  chomp;

  my ($chromosome,$source,$feature,$start,$stop,$score,$strand,$other,$attr) = split /\t/;

  unless ($source eq 'curated' and $feature eq 'CDS' or
          $source eq 'Coding_transcript' and ($feature eq 'mRNA' or $feature eq 'protein_coding_primary_transcript') or
          $source eq 'Non_coding_transcript' and $feature eq 'nc_primary_transcript' or
          $source eq 'miRNA_primary_transcript' and $feature eq 'miRNA_primary_transcript' or
          $source eq 'miRNA_mature' and $feature eq 'miRNA' or
          $source eq 'Pseudogene' and $feature eq 'pseudogenic_transcript' or
          $source =~ /RNA$/ and $source eq $feature) {
    print $gff_out_fh "$_\n";
    next;
  }
    
  print $gff_out_fh "$chromosome\t$source\t$feature\t$start\t$stop\t$score\t$strand\t$other\t";
  
  if ($gff3) {
    if ( $feature eq 'CDS' ) {
      my ($cds) = $attr =~ /ID=CDS:([^;]+);/;
      if (not defined $cds) {
        $log->log_and_die("Could not find CDS id in attribute field: $attr\n");
      } 
      # Note: in GFF3, CDS features are split across several lines. It is wasteful and unnecessary to 
      # decorate all of the segments, so only do the first
      if (not exists $already_done_cds{$cds}) {
        $attr .=  ";Name=CDS:$cds";
        $attr .=  ";prediction_status=$tran_status->{$cds}" if exists $tran_status->{$cds} and $tran_status->{$cds};
        $attr .=  ";wormpep=$tran_wormpep->{$cds}"           if exists $tran_wormpep->{$cds} and $tran_wormpep->{$cds};
        $attr .=  ";locus=$tran_locus->{$cds}"              if exists $tran_locus->{$cds} and $tran_locus->{$cds};
        
        $already_done_cds{$cds} = 1;
        $changed_lines++;
      }
    } else {
      my $tr;
      if (/ID=Transcript:([^;]+);/) {
        $tr = $1; 
      } elsif (/ID=Pseudogene:([^;]+);/) {
        $tr = $1;
      }

      if ($tr) {
        $attr .=  ";wormpep=$tran_wormpep->{$tr}"      if exists $tran_wormpep->{$tr} and $tran_wormpep->{$tr};
        $attr .=  ";locus=$tran_locus->{$tr}"          if exists $tran_locus->{$tr} and $tran_locus->{$tr};
        $changed_lines++;
      }
      elsif ($attr !~ /Name=Operon:/) {
        $log->log_and_die("Could not find Transcript id in attribute field: $attr\n");
      } 
    }
  } else {
    if( $feature eq 'CDS') {
      my ($cds) = $attr =~ (/CDS \"(\S+)\"/);

      $attr .=  " ; Prediction_status \"$tran_status->{$cds}\""  if exists $tran_status->{$cds} and $tran_status->{$cds};
      $attr .=  " ; WormPep \"$tran_wormpep->{$cds}\""           if exists $tran_wormpep->{$cds} and $tran_wormpep->{$cds};
      $attr .=  " ; Locus \"$tran_locus->{$cds}\""               if exists $tran_locus->{$cds} and $tran_locus->{$cds};

      $changed_lines++;
    } else {
      #non-coding genes
      my ($tr) = $attr =~ (/Transcript \"(\S+)\"/);
      if (not $tr) {
        ($tr) = $attr =~ (/Pseudogene \"(\S+)\"/);
      }

      if ($tr) {
        $attr .= " ; WormPep \"$tran_wormpep->{$tr}\"" if exists $tran_wormpep->{$tr} and $tran_wormpep->{$tr};
        $attr .= " ; Locus \"$tran_locus->{$tr}\""     if exists $tran_locus->{$tr} and $tran_locus->{$tr};
        $changed_lines++;
      } elsif ($attr !~ /Operon \"\S+\"/) {
        $log->log_and_die("Could not find transcript name in attr field : $attr\n");
      }
    }
  }
  
  print $gff_out_fh "$attr\n";
}
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n");
$log->write_to("Finished processing : $changed_lines lines modified\n");
$log->mail();
exit(0);

##############################################################
#
# Subroutines
#
##############################################################

sub get_data {
  my $db = Ace->connect('-path' => $database) 
      or $log->log_and_die("cant open Ace connection to db\n".Ace->error."\n");

  my %cds2cgc = $wormbase->FetchData('cds2cgc');
  my %cds2wormpep = $wormbase->FetchData('cds2wormpep');
  my %cds2status = $wormbase->FetchData('cds2status');

  my %rna2cgc = $wormbase->FetchData('rna2cgc');
  my %pseudo2cgc = $wormbase->FetchData('pseudo2cgc');


  $log->write_to("Fetching CDS info\n");

  my $query = "find CDS where Corresponding_protein AND Method = \"curated\" AND Species = \"$species_full\"";
  my $cds = $db->fetch_many('-query' => $query);
  while(my $cds = $cds->next) {
    my @coding_transcripts = $cds->Corresponding_transcript;

    foreach my $ct (@coding_transcripts) {
      if ($ct->name ne $cds->name) {
        my $wormpep = $cds2wormpep{$cds};
        $cds2wormpep{$ct->name} = $wormpep;

        if (exists $cds2cgc{$cds->name}) {
          $rna2cgc{$ct->name} = $cds2cgc{$cds->name};
        }
      }
    }
  }

  $db->close();
  $log->write_to("Closed connection to DB\n");


  my (%tran_status, %tran_wormpep, %tran_locus);

  # locus
  foreach my $cds (keys %cds2cgc) {
    $tran_locus{$cds} = $cds2cgc{$cds};
  }
  foreach my $tran (keys %rna2cgc) {
    $tran_locus{$tran} = $rna2cgc{$tran};
  }
  foreach my $pse (keys %pseudo2cgc) {
    $tran_locus{$pse} = $pseudo2cgc{$pse};
  }

  # status
  foreach my $cds (keys %cds2status) {
    $tran_status{$cds} = $cds2status{$cds};
  }

  # wormpep
  foreach my $cds (keys %cds2wormpep) {
    $tran_wormpep{$cds} = $wormbase->wormpep_prefix . ":" . $cds2wormpep{$cds};
  }
  
  # finally, need to get wormpep info for coding 

  return (\%tran_status, \%tran_wormpep, \%tran_locus);
}


1;
