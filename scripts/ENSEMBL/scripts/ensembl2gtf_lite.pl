#!/usr/local/bin/perl -w

use strict;

use Getopt::Long;
use Carp;
use FileHandle;
use FindBin qw($Bin);

# Load libraries needed for reading config -----------------------------------

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SeqFeature;

my $host;#   = 'ens-livemirror';
my $user;#   = 'ensro';
my $pass;#   = '';
my $port;#   = 3306;
my $dbname;# = 'homo_sapiens_core_42_36d';

my @seq_regions;
my $gtf_file = undef;
my $include_codons = 1;

my $start = undef;
my $end = undef;
my $seqname = undef;
my $coordsystem = 'toplevel';
my $coordsysversion;

$| = 1;

&GetOptions(
            'host:s'        => \$host,
            'user:s'        => \$user,
            'dbname:s'      => \$dbname,
            'pass:s'        => \$pass,
            'port:n'        => \$port,
            'codons'        => \$include_codons,
            'seqregion=s@'  => \@seq_regions,
            'gtffile:s'     => \$gtf_file,
            'coordsystem:s' => \$coordsystem,
            'coordsysversion:s' => \$coordsysversion,
            
);



my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                            -host   => $host,
                                            -user   => $user,
                                            -port   => $port,
                                            -pass   => $pass,
                                            -dbname => $dbname
					    );

my $sa  = $db->get_SliceAdaptor();

my @slices;
if (@seq_regions) {
  foreach my $sr (@seq_regions) {
    push @slices, $sa->fetch_by_region($coordsystem,
                                       $sr,
                                       undef,
                                       undef,
                                       undef,
                                       $coordsysversion);
  }
} else {
  @slices = @{$sa->fetch_all(
                             $coordsystem, 
                             $coordsysversion,
                             )};
}

my $outfile = $gtf_file;
my $gtffp = new FileHandle;
$gtffp->open(">$outfile") or croak "Unable to open $outfile for write";
$gtffp->autoflush(1);


foreach my $slice (@slices) {
  my $genes = $slice->get_all_Genes(undef,undef,1);
  foreach my $gene (@$genes) {
    foreach my $trans (@{$gene->get_all_Transcripts}) {
      write_transcript_gtf($gtffp,$slice,$gene,$trans,$include_codons,$seqname);
    }
  }
}
 # end foreach species



###########################################################################
sub make_start_codon_features {
  my ($trans,$id) = @_;


  if (!$trans->translation) {
    return (());
  }

  my $cdna_start_pos = $trans->cdna_coding_start;

  my @pepgencoords = $trans->cdna2genomic($cdna_start_pos,$cdna_start_pos + 2);

  if(scalar(@pepgencoords) > 3) {
    croak "pep start does not map cleanly\n";
  }

  unless($pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
    croak "pep start maps to gap\n";
  }
  unless($pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
    croak "pep start (end of) maps to gap\n";
  }

  my @translateable = @{$trans->get_all_translateable_Exons};
  my @startc_feat;
  my $phase = 0;
  foreach my $pepgencoord (@pepgencoords) {
    push @startc_feat, new Bio::EnsEMBL::SeqFeature(
                             -seqname => $id,
                             -source_tag => 'starttrans',
                             -primary_tag => 'similarity',
                             -start => $pepgencoord->start,
                             -end   => $pepgencoord->end,
                             -phase => $phase,
                             -strand => $translateable[0]->strand);
    $phase = 3 - ($pepgencoord->end - $pepgencoord->start + 1);
  }
  if ($translateable[0]->strand == 1) {
    @startc_feat = sort {$a->start <=> $b->start } @startc_feat;
  } else {
    @startc_feat = sort {$b->start <=> $a->start } @startc_feat;
  }
  return @startc_feat;

}

sub write_transcript_gtf {
  my ($fh,$slice,$gene,$transcript,$include_codons,$seqname) = @_;
  my $sliceoffset = $slice->start-1;

  my ($hasstart, $hasend, @startcs, @endcs);
  
  if ($include_codons) {
    ($hasstart,$hasend) = check_start_and_stop($slice,$transcript);

    @startcs =  make_start_codon_features($transcript,$transcript->stable_id) if $hasstart;
    @endcs   =  make_stop_codon_features($transcript,$transcript->stable_id) if $hasend;
  }


  my $chrname;
  $chrname = $slice->seq_region_name;

  my $idstr;

  if (defined($seqname)) {
    $idstr = $seqname;
  } else {
    $idstr = $chrname;
  }

  my @translateable_exons = @{$transcript->get_all_translateable_Exons} if $transcript->translation;

  my $count=1;
  my $intrans = 0;
  my $instop = 0;

  foreach my $exon (@{$transcript->get_all_Exons}) {
    my $strand = $exon->strand;

    if ($exon->strand == -1) {
        $strand = "-";
    } elsif ($exon->strand == 1) {
        $strand = "+";
    } elsif ($exon->strand == 0) {
        $strand = ".";
    }

    if ($transcript->translation && $exon == $transcript->translation->start_Exon) {
      $intrans = 1;
    }

    print $fh $idstr . "\t" . 
              $gene->biotype . "\t" . 
              'exon' . "\t" . 
              ($exon->start+$sliceoffset) . "\t". 
              ($exon->end+$sliceoffset) . "\t". 
              "." . "\t". 
              $strand . "\t". 
              "." . "\t";
    print_attribs($fh,$gene,$transcript,$count,'exon');
    print $fh "\n";

    if ($intrans) {
      my $cdsexon = shift @translateable_exons;
      my $phase = $cdsexon->phase;
      if ($cdsexon->phase == 1) {
        $phase = 2;
      } elsif ($cdsexon->phase == 2) {
        $phase = 1;
      } elsif ($cdsexon->phase == -1) {
        $phase = 0;
      }

      my $exon_start = $cdsexon->start;
      my $exon_end   = $cdsexon->end;
      if ($transcript->translation && 
          $hasend && 
          ($exon->end >= $endcs[0]->start && $exon->start <= $endcs[0]->end)) {

        if ($cdsexon->strand == 1) {
          $exon_end = $cdsexon->end - $endcs[0]->length;
        } else {
          $exon_start = $cdsexon->start + $endcs[0]->length;
        }
      }

      if ($exon_start <= $cdsexon->end &&
          $exon_end >= $cdsexon->start &&
          !$instop) {
        print $fh $idstr . "\t" . 
                  $gene->biotype . "\t" . 
                  'CDS' . "\t" . 
                  ($exon_start+$sliceoffset) . "\t". 
                  ($exon_end+$sliceoffset) . "\t". 
                  "." . "\t". 
                  $strand . "\t". 
                  $phase . "\t";
        print_attribs($fh,$gene,$transcript,$count,'CDS');
        print $fh "\n";
      }
    }
    if ($transcript->translation && 
        $exon == $transcript->translation->start_Exon && $hasstart) {
      my $tmpcnt = $count;
      foreach my $startc (@startcs) {
        print $fh $idstr . "\t" . 
                  $gene->biotype . "\t" . 
                  'start_codon' . "\t" . 
                  ($startc->start+$sliceoffset) . "\t". 
                  ($startc->end+$sliceoffset) . "\t". 
                  "." . "\t". 
                  $strand . "\t". 
                  $startc->phase . "\t";
        print_attribs($fh,$gene,$transcript,$tmpcnt++,'start_codon');
        print $fh "\n";
      }
    }
    if ($transcript->translation && 
        ($exon == $transcript->translation->end_Exon)) {
      if ($hasend) {
        my $tmpcnt = $count - $#endcs;
        foreach my $endc (@endcs) {
          print $fh $idstr . "\t" . 
                    $gene->biotype . "\t" . 
                    'stop_codon' . "\t" . 
                    ($endc->start+$sliceoffset) . "\t". 
                    ($endc->end+$sliceoffset) . "\t". 
                    "." . "\t". 
                    $strand . "\t". 
                    $endc->phase . "\t";
          print_attribs($fh,$gene,$transcript,$tmpcnt++,'stop_codon');
          print $fh "\n";
        }
      }
      $intrans = 0;
    }

    if (scalar(@endcs) && 
        ($exon->end >= $endcs[0]->start && $exon->start <= $endcs[0]->end)) {
      $instop = 1;
    }

    $count++;
  }
}

sub make_stop_codon_features {
  my ($trans,$id) = @_;

  if (!$trans->translation) {
    return (());
  }
  my @translateable = @{$trans->get_all_translateable_Exons};

  my $cdna_endpos = $trans->cdna_coding_end;

  my @pepgencoords = $trans->cdna2genomic($cdna_endpos-2,$cdna_endpos);

  if(scalar(@pepgencoords) > 3) {
    croak "pep end does not map cleanly\n";
  }

  unless($pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
    croak "pep end maps to gap\n";
  }
  unless($pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
    croak "pep end (end of) maps to gap\n";
  }

  my @stopc_feat;
  my $phase = 0;
  foreach my $pepgencoord (@pepgencoords) {
    push @stopc_feat, new Bio::EnsEMBL::SeqFeature(
                             -seqname => $id,
                             -source_tag => 'endtrans',
                             -primary_tag => 'similarity',
                             -start => $pepgencoord->start,
                             -end   => $pepgencoord->end,
                             -phase => $phase,
                             -strand => $translateable[0]->strand);
    $phase = 3 - ($pepgencoord->end-$pepgencoord->start+1);
  }

  if ($translateable[0]->strand == 1) {
    @stopc_feat = sort {$a->start <=> $b->start } @stopc_feat;
  } else {
    @stopc_feat = sort {$b->start <=> $a->start } @stopc_feat;
  }
  return @stopc_feat;
}

sub print_attribs {
  my ($fh,$gene,$transcript,$count,$type) = @_;


  my $gene_name;
  $gene_name = $gene->external_name;

  my $trans_name;
  $trans_name = $transcript->external_name;

  print $fh " gene_id \"" .  get_gene_id($gene) . "\";" .
            " transcript_id \"" . get_transcript_id($transcript) . "\";";
  print $fh " exon_number \"$count\";";
  print $fh " gene_name \"" . $gene_name . "\";" if ($gene_name);
  print $fh " transcript_name \"" . $trans_name . "\";" if ($trans_name);
  if ($type eq 'CDS') {
    print $fh ' protein_id "' . get_translation_id($transcript->translation) . '";';
  }
}


sub get_gene_id {
  my $gene = shift;

  if (defined($gene->stable_id)) {
    return $gene->stable_id;
  }
  return $gene->dbID;
}

sub get_transcript_id {
  my $transcript = shift;

  if (defined($transcript->stable_id)) {
    return $transcript->stable_id;
  }
  return $transcript->dbID;
}
sub get_translation_id {
  my $translation = shift;

  if (defined($translation->stable_id)) {
    return $translation->stable_id;
  }
  return $translation->dbID;
}

sub check_start_and_stop {
  my ($slice,$trans) = @_;

  my ($has_start, $has_end);

  my $tln = $trans->translation;
  return (0,0) if not defined $tln;

  my $cdna_seq     = uc($trans->spliced_seq);

  if ($trans->get_all_Exons->[0]->phase > 0) {
    $has_start = 0;
  } else {
    my $coding_start = $trans->cdna_coding_start;
    my $startseq     = uc(substr($cdna_seq,$coding_start-1,3));

    $has_start = ($startseq eq 'ATG') ?  1 : 0;
  }

  if ($trans->get_all_Exons->[-1]->end_phase > 0) {
    $has_end = 0;
  } else {
    my $coding_end   = $trans->cdna_coding_end;
    my $endseq       = substr($cdna_seq,$coding_end-3,3);

    $has_end = ($endseq eq "TAG" or $endseq eq "TGA" or $endseq eq "TAA") ? 1 : 0;
  }

  return ($has_start, $has_end);
}

