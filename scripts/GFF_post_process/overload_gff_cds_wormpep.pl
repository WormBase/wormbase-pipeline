#!/usr/bin/env perl
#
# overload_gff_cds_wormpep.pl
#
# Overloads the CDS and Transcript lines with extra info (mostly wormpep)
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-07-21 11:07:59 $

#
#    1. Brief_identification
#    2. Locus names when available
#    3. WormPep IDs
#    4. Confirmed_by status
#    5. WBGene IDs
#    6. Genbank accesions (for region:Genomic_canonical)
#
# CDS "JC8.10a" ; Note "Inositol polyphosphate phosphatase, catalytic domain homologues" ; WormPep "WP:CE28239" ; Note "unc-26" ; Confirmed_by "cDNA" ; Gene "WBGene00006763"
#

use strict;                                      
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

my ($debug, $test, $store, $wormbase, $species);
my ($infile, $outfile, $gff3, $changed_lines, %already_done_cds);

GetOptions (
  "debug=s"    => \$debug,
  "test"       => \$test,
  "store:s"    => \$store,
  "species:s"  => \$species,
  "infile:s"   => \$infile,
  "outfile:s"  => \$outfile,
  "gff3"       => \$gff3,
    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species
			     );
}
my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}

# get data from wormpep release
my ($CDS, %wormpep, %geneID, %status, %locus, %briefID, %RNAgenes, %seqname2geneid);

$wormbase->FetchData("worm_gene2geneID_name",\%seqname2geneid);
&get_wormpep_info();

open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  

while (<$gff_in_fh>) {
  chomp;
  
  unless  (/^\S+\s+(WormBase|curated|ncRNA|snlRNA|tRNAscan\S+)\s+(\w+primary_transcript|CDS|tRNA)/ or
           /^\S+\s+(rRNA)\s+(\w+_primary_transcript)/ or
           /^\S+\s+(curated_miRNA)\s+(miRNA_primary_transcript)/ or
           /^\S+\s+(\w+_mature_transcript)\s+(snRNA|snoRNA|tRNA|scRNA|miRNA|stRNA)/) {
    print $gff_out_fh "$_\n";
    next;
  }
  
  my ($chromosome,$source,$feature,$start,$stop,$score,$strand,$other,$attr) = split /\t/;
  
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
        $attr .=  ";Note=$briefID{$cds}"                              if ($briefID{$cds} ne "");
        $attr .=  ";WormPep=".$wormbase->pep_prefix.":$wormpep{$cds}" if ($wormpep{$cds} ne "");
        $attr .=  ";Locus=$locus{$cds}"                               if ($locus{$cds} ne "");
        $attr .=  ";Status=$status{$cds}"                             if ($status{$cds} ne "");
        $attr .=  ";Gene=$geneID{$cds}"                               if ($geneID{$cds} ne "");
        
        $already_done_cds{$cds} = 1;
        $changed_lines++;
      }
    } else {
      my ($transcript) = $attr =~ /ID=Transcript:([^;]+);/;
      if (not defined $transcript) {
        $log->log_and_die("Could not find Transcript id in attribute field: $attr\n");
      } 
      $attr .= ";Note=".$RNAgenes{$transcript}->{remark} if $RNAgenes{$transcript}->{remark} ;
      $attr .= ";Locus=".$RNAgenes{$transcript}->{locus} if $RNAgenes{$transcript}->{locus} ;;
      $attr .= ";Gene=".$seqname2geneid{$transcript}     if $seqname2geneid{$transcript} ;
      $changed_lines++;
    }
  } else {
    if( $feature eq 'CDS') {
      my ($i) = $attr =~ (/CDS \"(\S+)\"/);
      $attr .=  " ; Note \"$briefID{$i}\""                              if ($briefID{$i} ne "");
      $attr .=  " ; WormPep \"".$wormbase->pep_prefix.":$wormpep{$i}\"" if ($wormpep{$i} ne "");
      $attr .=  " ; Locus \"$locus{$i}\""                               if ($locus{$i} ne "");
      $attr .=  " ; Status \"$status{$i}\""                             if ($status{$i} ne "");
      $attr .=  " ; Gene \"$geneID{$i}\""                               if ($geneID{$i} ne "");
      $changed_lines++;
    } else {
      #non-coding genes
      my ($transcript) = $attr =~ (/Transcript \"(\S+)\"/);
      $attr .= " ; Note \"".$RNAgenes{$transcript}->{remark}."\"" if $RNAgenes{$transcript}->{remark};
      $attr .= " ; Locus \"".$RNAgenes{$transcript}->{locus}."\"" if $RNAgenes{$transcript}->{locus};
      $attr .= " ; Gene \"".$seqname2geneid{$transcript}."\""     if $seqname2geneid{$transcript} ;
      $changed_lines++;
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

sub get_RNA_info {
  my $rna_file = $wormbase->wormrna."/wormrna".$wormbase->get_wormbase_version.".rna";
  open (RNA,"<$rna_file") or $log->log_and_die("cant open $rna_file\t$!\n");
  while(<RNA>) {
    #I couldnt think of a way to do this in one regex!
    my ($locus, $remark, $cds);
    if(/locus:(\S+)/){
      $locus = $1;
      />(\S+)\s+(.*)\s+locus/;
      $cds = $1;
      $remark = $2;
    }
    elsif(/>(\S+)\s+(.*)$/) {
      $cds = $1;
      $remark = $2;
    }
    $RNAgenes{$cds}->{remark} = $remark if $remark;
    $RNAgenes{$cds}->{locus} = $locus if $locus;
  }
  close RNA;
}


sub get_wormpep_info {
  my $file = $wormbase->wormpep."/".$wormbase->pepdir_prefix."pep".$wormbase->get_wormbase_version;
  open (my $wpfh, "<$file") or $log->log_and_die("cant open $file $!\n");
  while (<$wpfh>) {
    #>4R79.2 CE19650 WBGene00007067  Ras family      status:Partially_confirmed      UniProt:Q9XXA4_CAEEL    protein_id:CAA20282.1
    #>4R79.1b        CE39659 WBGene00003525  locus:nas-6     status:Partially_confirmed      UniProt:Q2HQL9_CAEEL    protein_id:CAJ76926.1
    if (/^>(\S+)\s+(\S+)\s+(WBGene\d+)(\s+locus:(\S+))*\s+([^\t]*?)\s*status:(\S+)/) {
      $CDS           = $1;
      $wormpep{$CDS} = $2;
      $geneID{$CDS}  = $3;
      $locus{$CDS}   = $5;
      $briefID{$CDS} = $6;
      $status{$CDS}  = $7;
      
    }
  }
  close($wpfh);
}

1;
