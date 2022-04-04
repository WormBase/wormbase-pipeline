#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

use lib $ENV{CVS_DIR};

use Bio::SeqIO;
use Modules::Map_Helper;

use LSF RaiseError => 0, PrintError => 1, PrintOutput => 1;
use LSF::JobManager;


################################
# command-line options         #
################################

my (
  $outfile,
  $workdir,
  $in_genome,
  $in_gff,
  $mapping_script,
    );

GetOptions(
  "workdir=s"       => \$workdir,
  "outfile=s"       => \$outfile,
  "ingff=s"         => \$in_gff,
  "ingenome=s"      => \$in_genome,
  "mappingscript=s" => \$mapping_script,
    );

my ($gnome) = &read_genome( $in_genome );
my ($pcr_prod, $genes) = &read_gff( $in_gff );
my $gindex = &build_gene_segment_index( $genes );
my $rnai_objs = &get_rnai_sequences($gnome, $pcr_prod);
&map_rnai_segs_to_genome( $gnome, $rnai_objs );

open( my $out_fh, ">$outfile" ) or die "Could not open $outfile for writing\n";

foreach my $type (keys %$rnai_objs) {
  foreach my $rnai_id (sort keys %{$rnai_objs->{$type}}) {
    my %gene_mappings;
    
    if (exists $rnai_objs->{$type}->{$rnai_id}->{genome_mappings}) {
      foreach my $mapping_type (keys %{ $rnai_objs->{$type}->{$rnai_id}->{genome_mappings} }) {
        foreach my $seg (@{$rnai_objs->{$type}->{$rnai_id}->{genome_mappings}->{$mapping_type}}) {
          my @hits = @{$gindex->search_feature_segments( $seg->[0], $seg->[1], $seg->[2], undef, 50)};
          
          foreach my $hit (@hits) {
            my ($gene, $trans) = split(/\//, $hit);
            $gene_mappings{$mapping_type}->{$gene} = 1;
          }
        }
      }
    }
    
    my (@prim_genes, @scd_genes);
    
    if (exists $gene_mappings{RNAi_primary}) {
      foreach my $gid (sort keys %{$gene_mappings{RNAi_primary}}) {
        my $id = sprintf("%s(%s)(%s)",
                         $genes->{$gid}->{id},
                         ($genes->{$gid}->{locus}) ? $genes->{$gid}->{locus} : $genes->{$gid}->{sequence_name},
                         $genes->{$gid}->{biotype});
        push @prim_genes, $id;
      }
    }
    if (exists $gene_mappings{RNAi_secondary}) {
      foreach my $gid (sort keys %{$gene_mappings{RNAi_secondary}}) {
        if (not exists $gene_mappings{RNA_primary} or not exists $gene_mappings{RNA_primary}->{$gid}) {
          my $id = sprintf("%s(%s)(%s)",
                           $genes->{$gid}->{id},
                           ($genes->{$gid}->{locus}) ? $genes->{$gid}->{locus} : $genes->{$gid}->{sequence_name},
                           $genes->{$gid}->{biotype});
          push @scd_genes, $id;
        }
      }
    }
    
   printf $out_fh "%s\tprimary=%s\tsecondary=%s\n", $rnai_id, join(",", @prim_genes), join(",", @scd_genes);
  }
}

################################################################
sub map_rnai_segs_to_genome {
  my ($genome, $rnais) = @_;
  
  foreach my $type (sort keys %$rnais) {
    my @rnai_ids = sort keys %{$rnais->{$type}};
    
    my $chunk_size = 1000;      
    my $total_rnai = scalar(@rnai_ids);
    my $chunk_count = int(scalar(@rnai_ids) / $chunk_size);
    
    my (@rnai_chunks, @query);
    
    # randomise order for better distribution
    while(my $rnai_id = shift @rnai_ids) {
      my $chunk_id = int(rand($chunk_count));
      
      push @{$rnai_chunks[$chunk_id]}, $rnais->{$type}->{$rnai_id};
    }
    
    printf STDERR "Fetched %d objects in %d chunks for type %s\n", $total_rnai, scalar(@rnai_chunks), $type;
    
    for(my $i=0; $i < @rnai_chunks; $i++) {
      my $file = "$workdir/rnai_query.$type.$i.fa";
      
      my $seqout = Bio::SeqIO->new(-format => 'fasta', -file => ">$file");   
      foreach my $rnai (@{$rnai_chunks[$i]}) {
        my $sobj = Bio::PrimarySeq->new(-seq => $rnai->{sequence},
                                        -id  => $rnai->{name});
        
        $seqout->write_seq($sobj);
      }
      
      push @query, $file;
    }
    
    my $lsf = LSF::JobManager->new();
    
    my @outfiles;
    for(my $i=0; $i < @query; $i++) {
      my $query = $query[$i];
      my $jname    = "worm_rna2genome.$type.$i.$$";
      my $out_file = "$workdir/rnai_out.$type.$i.gff";
      my $lsf_out  = "$workdir/rnai_out.$type.$i.lsf_out";
      
      my $cmd = "perl $mapping_script -query $query -target $in_genome -outfile $out_file";
      
      my @bsub_options = (-J => $jname, 
                          -o => $lsf_out,
                          -E => 'test -w ' . $workdir,
                          -M => 2600,
                          -R => 'select[mem>=2600] rusage[mem=2600]'
          );
      if (not -e $out_file) {
        $lsf->submit(@bsub_options, $cmd);
      }
      push @outfiles, $out_file;
    }   
    
    $lsf->wait_all_children( history => (scalar(@query) == 1) ? 0 : 1 );
    # wait a few seconds here; the LSF job manager sometimes lags a little in
    # registering the job history
    sleep(5);
    
    for my $job ( $lsf->jobs ) {    
      print STDERR ("$job exited non zero\n") if $job->history->exit_status != 0;
    }
    $lsf->clear;

    foreach my $outfile (@outfiles) {
      open(my $rnaifh, $outfile) or die "Could not open $outfile for writing\n";
      while(<$rnaifh>) {
        /^\#/ and next;
        chomp;
        my @l = split(/\t/, $_);

        my ($hit) = $l[8] =~ /Target=(\S+)/;
        push @{$rnais->{$type}->{$hit}->{genome_mappings}->{$l[1]}}, [$l[0], $l[3], $l[4], $l[6]];
      }
    }
  }
}



################################################################
sub read_genome {
  my ($genome_file) = @_;

  my %genome;

  my $seqio = Bio::SeqIO->new(-format => 'fasta',
                              -file => $genome_file);
  while(my $seq = $seqio->next_seq) {
    $genome{$seq->id} = uc($seq->seq);
  }
  
  return \%genome;
}

################################################################
sub read_gff {
  my ($gff_file) = @_;

  my (%pcr_prod, %genes, %tran2gene);

  my $gff_fh;
  if ($gff_file =~ /\.gz$/) {
    open($gff_fh, "gunzip -c $gff_file |");
  } else {
    open($gff_fh, $gff_file);
  }

  while(<$gff_fh>) {
    /^\#/ and next;
    chomp;
    
    /^(\S+)\s+GenePair_STS\s+PCR_product\s+(\d+)\s+(\d+)\s+\S+\s+\S+\s+\S+\s+Name=(sjj\S+)/ and do {
      my ($chr, $st, $en, $pcr_prod) = ($1, $2, $3, $4);
      $pcr_prod =~ s/;vend\S+//; 

      my ($type) = $pcr_prod =~ /^([^_]+)_/;
      push @{$pcr_prod{$type}}, [$pcr_prod, $chr, $st, $en];

      next;
    };

    next unless /^\S+\s+WormBase\s+/; 
    
    my @l = split(/\t/, $_);
    my %props;
    foreach my $pair (split(/;/, $l[8])) {
      my ($k, $v) = split(/=/, $pair);
      $props{$k} = $v;
    }
    if ($l[2] eq 'gene') {
      my $gene = {
        gff_id => $props{ID},
        id => $props{Name},
        locus => $props{locus},
        biotype => $props{biotype},
        sequence_name => $props{sequence_name},
      };
      $genes{$props{ID}} = $gene;
    } elsif (exists $props{Parent} and $props{Parent} =~ /^Gene:\S+/) {
      # transcript-type feature
      my $tran = {
        gff_id => $props{ID},
        id => $props{Name},
        biotype => $props{biotype},
      };
      $genes{$props{Parent}}->{transcripts}->{$tran->{gff_id}} = $tran; 
      $tran2gene{$tran->{gff_id}} = $genes{$props{Parent}}->{gff_id} 
    } elsif ($l[2] eq 'exon') {
      push @{$genes{ $tran2gene{ $props{Parent} }}->{transcripts}->{ $props{Parent} }->{exons}}, [$l[0], $l[3], $l[4], $l[6]];
    }
  }

  return (\%pcr_prod, \%genes);
}


################################################################
sub build_gene_segment_index {
  my ($genesh) = @_;

  my $fm = Map_Helper->new();

  foreach my $gacc (keys %$genesh) {
    foreach my $tacc (keys %{$genesh->{$gacc}->{transcripts}}) {
      foreach my $ex ( @{$genesh->{$gacc}->{transcripts}->{$tacc}->{exons}} ) {
        $fm->register_segment( @$ex, join("/", $gacc, $tacc) );
      }
    }
  }
  
  $fm->build_index();
  
  return $fm;
}


################################################################
sub get_rnai_sequences {
  my ($genomeh, $pcr_prodh) = @_;

  my (%results);

  foreach my $type (keys %$pcr_prodh) {
    foreach my $prod (@{$pcr_prodh->{$type}}) {
      my ($nm, $chr, $st, $en) = @$prod;

      my $seq_string = substr( $genomeh->{$chr}, $st - 1, $en - $st + 1);

      $results{$type}->{$nm} = {
        name => $nm,
        sequence => $seq_string,
      };
    }
  }

  return \%results;
}
