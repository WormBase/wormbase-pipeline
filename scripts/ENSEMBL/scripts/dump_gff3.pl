#!/usr/bin/env perl

use strict;
use warnings;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use IO::File;
use Getopt::Long;
use List::MoreUtils qw/uniq/;

my ($dbname, $dbuser, $dbpass, $dbport, $dbhost, $debug,
    @dump_slice, $num_jobs, $out_file, $slim, $out_fh);

my $dumpdir = ".";

GetOptions(
  'host=s'     => \$dbhost,
  'port=s'     => \$dbport,
  'user=s'     => \$dbuser,
  'pass=s'     => \$dbpass,
  'dbname=s'     => \$dbname,
  'slice=s@'     => \@dump_slice,
  'submit=i'     => \$num_jobs,
  'dumpdir=s'    => \$dumpdir,
  'outfile:s'    => \$out_file,
  'slim'         => \$slim,
  'debug'        => \$debug,
          )or die ("Couldn't get options");

my $ensdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
	-dbname  => $dbname,
	-host    => $dbhost,
	-port    => $dbport,
	-user    => $dbuser,
	-pass    => $dbpass,
);

# Get all toplevel slices
my $sa = $ensdb->get_SliceAdaptor();
my $archive_adaptor = $ensdb->get_adaptor("ArchiveStableId");
my @slices;
@dump_slice = split(/,/,join(',',@dump_slice));

if (@dump_slice) {
  foreach (@dump_slice) {
    my $sl = $sa->fetch_by_region('toplevel',$_);
    if (not defined $sl) {
      die "Could not fetch slice for $_\n";
    }
    push @slices, $sl;
  }
}
else {
  @slices = sort { $b->length <=> $a->length } @{$sa->fetch_all('toplevel')};
}

if ($num_jobs) {
  &submit_and_collate();
  exit(0);
}


if (defined $out_file) {
  open($out_fh, ">$out_file") or die "Could not open $out_file for writing\n";
  print $out_fh "##gff-version 3\n";
}


while( my $slice = shift @slices) {
  my $slice_name = $slice->seq_region_name();
  my $slice_size = $slice->length;

  $debug and print STDERR " Fetching and processing data for $slice_name...\n";

  if (not defined $out_fh) {
    my $sl_file = "$dumpdir/${slice_name}.gff3";
    open($out_fh, ">$sl_file") or die "Could not open $sl_file for writing\n";  
  }

  print $out_fh "##sequence-region $slice_name 1 $slice_size\n";
  
  # Get all the genes on this slice
  $debug and print STDERR " Fetching and processing genes...\n";
  my $genes = $slice->get_all_Genes();
  while( my $gene=shift @$genes) {
    my $gene_gff_id = 'gene:'.$gene->stable_id();
    my %gene_to_dump = (
      gff_id    => $gene_gff_id,
      seqname   => $slice_name,
      start     => $gene->seq_region_start(),
      end       => $gene->seq_region_end(),
      strand    => $gene->strand(),
      display => $gene->stable_id(), 
      gff_source  => (defined $gene->analysis->gff_source) ? $gene->analysis->gff_source : "WormBase",
      attribs   => {
           biotype => $gene->biotype,
      },
    );
    if ($gene->display_xref()){
      $gene_to_dump{attribs}{locus} = $gene->display_xref->display_id();
    }
    if (my ($description, $description_source, $description_source_acc) = ($gene->description // "") =~ /(.*?)\s*\[Source:(.*);Acc:(.*)\]/){
      $gene_to_dump{attribs}{description} = $description;
      $gene_to_dump{attribs}{description_source} = $description_source;
      $gene_to_dump{attribs}{description_source_acc} = $description_source_acc;
    }
    my $archive_stable_id = $archive_adaptor->fetch_by_stable_id($gene->stable_id);
    if ($archive_stable_id) {
# EnsEMBL code doesn't quite handle non-linear history, e.g.  Smp_120050+Smp_199230->Smp_335780+Smp_315690
       my @archived_ids =  grep {$_ ne $gene->stable_id} uniq  map {$_->{old_id} ? $_->{old_id}->stable_id : ()}  @{$archive_stable_id->get_history_tree->get_all_StableIdEvents };
       if (@archived_ids){
          $gene_to_dump{attribs}{previous_stable_id} = join ",", @archived_ids;
       }
    }


    my $all_transcripts = $gene->get_all_Transcripts();
    while( my $transcript) {
      
      my $transcript_gff_id = 'transcript:'.$transcript->stable_id();

      my $tr_obj =  {
        gff_id    => $transcript_gff_id,
        seqname   => $slice_name,
        start     => $transcript->start(),
        end       => $transcript->end(),
        strand    => $transcript->strand(),
        display     => $transcript->stable_id(),
        gff_source  => (defined $transcript->analysis->gff_source) ?  $transcript->analysis->gff_source : "WormBase",
        gff_type    => (defined $transcript->analysis->gff_feature ) ? $transcript->analysis->gff_feature : $transcript->biotype,
        attribs     => { info => get_info($transcript) },
      };

      my $translation = $transcript->translation;

      if (defined $translation) {
        my $translation_gff_id = 'cds:'.$transcript->translation->stable_id();

        $tr_obj->{cds_start} = $translation->genomic_start();
        $tr_obj->{cds_end}   = $translation->genomic_end();
        $tr_obj->{translation_stable_id} = $translation_gff_id;

        my $all_t_exons = $transcript->get_all_CDS();
        
        # Note that the get_all_CDS method returns exons with phases that
        # have already been converted to GFF3 style. No need for conversion then. 
        while (my $cds = shift @{$all_t_exons}) {

          push @{$tr_obj->{'cds'}}, {
            gff_id    => $translation_gff_id,
            seqname   => $slice_name,
            start     => $cds->seq_region_start(),
            end       => $cds->seq_region_end(),
            strand    => $cds->strand(),
            phase     => $cds->phase(),
          };
        }

        my $cdna_coding_start = $transcript->cdna_coding_start();
        if ($cdna_coding_start > 1) {
          my @coords = $transcript->cdna2genomic(1, $cdna_coding_start - 1);
          @coords = grep { $_->isa("Bio::EnsEMBL::Mapper::Coordinate") } @coords;
          foreach my $coord (@coords) {
            push @{$tr_obj->{'utr5'}}, {
              seqname => $slice_name,
              start   => $coord->start, 
              end     => $coord->end,
              strand  => $coord->strand,
            };
          }
        }

        my $cdna_coding_end = $transcript->cdna_coding_end();
        if ($cdna_coding_end < $transcript->length) {
          my @coords = $transcript->cdna2genomic($cdna_coding_end + 1, $transcript->length);
          @coords = grep { $_->isa("Bio::EnsEMBL::Mapper::Coordinate") } @coords;
          foreach my $coord (@coords) {
            push @{$tr_obj->{'utr3'}}, {
              seqname => $slice_name,
              start   => $coord->start, 
              end     => $coord->end,
              strand  => $coord->strand,
            };
          }
        }
      }
      
      # get all exons and translatable_exons of the transcript
      my $all_exons=$transcript->get_all_Exons();

      my $exon_count = 1;
      while( my $exon = shift @{$all_exons}) {
        # dont use exon_stable_id because it may or may not be shared
        # between exons shared between transcripts depending on how the database was
        # built, and we want consistency in the dumps
        my $exon_gff_id = sprintf("exon:%s.%d", $transcript->stable_id, $exon_count++);
        
        push @{$tr_obj->{exon}}, {
          gff_id    => $exon_gff_id,
          seqname   => $slice_name,
          start     => $exon->seq_region_start(),
          end       => $exon->seq_region_end(),
          strand    => $exon->strand(),
        };
      }
      

      push @{$gene_to_dump{transcript}}, $tr_obj;
    }
    print $out_fh dump_gene(\%gene_to_dump);
  }

    
  #
  # Prediction transcripts
  #
  $debug and print STDERR " Fetching and processing prediction transcripts...\n";
  
  my $pts = $slice->get_all_PredictionTranscripts;    

  while(my $pt = shift @$pts) {
    my $label =  $pt->display_label;
    my $gff_id = sprintf("GBG_%d.%s", $pt->dbID, $label);

    foreach my $pe (@{$pt->get_all_Exons}) {
      my $feat = {
        seqname      => $pt->slice->seq_region_name,
        strand       => ($pt->strand > 0?'+':'-'),
        start        => $pe->seq_region_start,
        end          => $pe->seq_region_end,
        score        => ($pe->score||'.'),
        phase        => (defined $pe->phase) ? (3 - $pe->phase) % 3 : ".",
        logic_name   => $pt->analysis->logic_name,
        gff_source   => (defined $pt->analysis->gff_source) ? $pt->analysis->gff_source : "WormBase_prediction",
        feature_type => "CDS",
        display      => $label,
        gff_id       => $gff_id,
      };
      print $out_fh dump_feature($feat);
    }
  }
  
  next if $slim;

  {
    $debug and print STDERR " Fetching and processing protein alignments...\n";

    my @logics = &get_feature_logics($ensdb, "protein_align_feature", $slice);

    foreach my $logic (@logics) {
      my $ana = $ensdb->get_AnalysisAdaptor->fetch_by_logic_name($logic);
      if (not defined $ana) {
        print "Skipping $logic because not found in the analysis table\n";
        next;
      }      
      if (not $ana->gff_source or not $ana->gff_feature) {
        print "Skipping $logic because it does not have a defined gff_source and/or gff_feature\n";
        next;
      }

      my @blastx_features;

      $debug and print STDERR "  Fetching protein alignments for $logic...\n";
      
      my $features = $slice->get_all_ProteinAlignFeatures($logic);  
      while(my $feat = shift @$features) {

        my $cigar_line = flipCigarReference($feat->cigar_string); # for Lincoln
        if ($feat->strand < 0) {
          $cigar_line = reverse_cigar($cigar_line);
        }
        $cigar_line = cigar_to_almost_cigar($cigar_line);
        
        my $stripped_feature = {
          seqname      => $slice->seq_region_name,
          hit_id       => $feat->hseqname, 
          hit_start    => $feat->hstart,
          hit_end      => $feat->hend,
          display      => $feat->hseqname,
          strand       => ($feat->strand > 0?'+':'-'),
          start        => $feat->start,
          end          => $feat->end,
          score        => $feat->score,
          phase        => ".",
          p_value      => $feat->p_value,
          gff_id       => $feat->analysis->logic_name . "." . $feat->dbID,
          logic_name   => $feat->analysis->logic_name,
          cigar        => $cigar_line,
          gff_source   => (defined $feat->analysis->gff_source) ? $feat->analysis->gff_source : $feat->analysis->logic_name,
          feature_type => (defined $feat->analysis->gff_feature)? $feat->analysis->gff_feature : 'protein_match',
        };
        push @blastx_features, $stripped_feature;
      }

      $debug and print STDERR "  Filtering protein alignments for $logic...\n";

      my @filtered_features=filter_features(\@blastx_features, $logic);
      map {print $out_fh dump_feature($_)} @filtered_features;
    }
  }   
  
  $debug and print STDERR " Fetching and processing DNA alignments...\n";
  
  my @logics = &get_feature_logics($ensdb, "dna_align_feature", $slice);
  
  foreach my $logic (@logics) {
    my $features = $slice->get_all_DnaAlignFeatures($logic);
    while(my $feat = shift @$features) {
      my $cigar_line = flipCigarReference($feat->cigar_string); # for Lincoln
      if ($feat->strand < 0) {
        $cigar_line = reverse_cigar($cigar_line);
      }
      $cigar_line = cigar_to_almost_cigar($cigar_line);
      
      my $stripped_feature = {
        seqname      => $feat->slice->seq_region_name,
        hit_id       => $feat->hseqname,
        hit_start    => $feat->hstart,
        hit_end      => $feat->hend,
        display      => $feat->hseqname,
        strand       => ($feat->strand > 0?'+':'-'),
        start        => $feat->start,
        end          => $feat->end,
        score        => $feat->score,
        phase        => ".",
        p_value      => $feat->p_value,
        gff_id       => $feat->analysis->logic_name . "." . $feat->dbID,
        logic_name   => $feat->analysis->logic_name,
        cigar        => $cigar_line,
        gff_source   => (defined $feat->analysis->gff_source) ? $feat->analysis->gff_source : $feat->analysis->logic_name,
        feature_type => (defined $feat->analysis->gff_feature) ? $feat->analysis->gff_feature : 'nucleotide_match',
      };
      print $out_fh dump_feature($stripped_feature);
    }
  }

  $debug and print STDERR " Fetching and processing Repeat features...\n";
  
  my $features = $slice->get_all_RepeatFeatures;
  while(my $feature = shift @$features) {
    
    my $stripped_feature = {
      seqname      => $feature->slice->seq_region_name,
      strand       => ($feature->strand > 0?'+':'-'),
      start        => $feature->seq_region_start,
      end          => $feature->seq_region_end,
      score        => ($feature->score||'.'),
      phase        => ".",
      logic_name   => $feature->analysis->logic_name,
      gff_source   => (defined $feature->analysis->gff_source) ? $feature->analysis->gff_source : "WormBase",
      feature_type => (defined $feature->analysis->gff_feature) ? $feature->analysis->gff_feature : 'repeat_region',
    };


    if ($feature->repeat_consensus->repeat_class and $feature->repeat_consensus->name ne $feature->analysis->logic_name){
      $stripped_feature->{attributes}->{repeat_class} = $feature->repeat_consensus->repeat_class;
      $stripped_feature->{display} = $feature->repeat_consensus->name;
    }

    print $out_fh dump_feature($stripped_feature);
  }

  $debug and print STDERR " Fetching and processing Simple features...\n";  

  # get all simple features stored in the database (Operons etc. etc.)
  $features = $slice->get_all_SimpleFeatures;
  while(my $simpfeature = shift @$features) {
    my $stripped_simpfeature = {
      seqname      => $simpfeature->slice->seq_region_name,
      strand       => ($simpfeature->strand > 0?'+':'-'),
      start        => $simpfeature->seq_region_start,
      end          => $simpfeature->seq_region_end,
      score        => ($simpfeature->score||'.'),
      phase        => ".",
      logic_name   => $simpfeature->analysis->logic_name,
      gff_source   => (defined $simpfeature->analysis->gff_source) ? $simpfeature->analysis->gff_source : "WormBase",
      feature_type => (defined $simpfeature->analysis->gff_feature) ? $simpfeature->analysis->gff_feature : $simpfeature->analysis->logic_name,
    };
    print $out_fh dump_feature($stripped_simpfeature);
  }
  

  print $out_fh "###\n";
  close($out_fh) unless defined $out_file;
}

# close the file handle of the primary gff as you are done
close($out_fh) if defined $out_file;


############################

sub submit_and_collate {

  eval "use LSF::JobManager";
  
  if (not defined $out_file) {
    $out_file = "$dumpdir/collated_dump.gff3";
  }
  open(my $out_fh, ">$out_file") or die "Could not open $out_file for writing\n";
  print $out_fh "##gff-version 3\n";

  my (@batches);
  my $idx = 0;
  while(@slices) {
    push @{$batches[$idx]}, shift @slices;

    $idx++; $idx = 0 if $idx == $num_jobs;
  }

  my $lsf_man = LSF::JobManager->new();

  my $base_cmd = "perl $0 -dbhost $dbhost -dbname $dbname -dbuser $dbuser -dbport $dbport -dumpdir $dumpdir";
  $base_cmd .= " -dbpass $dbpass" if $dbpass;
  $base_cmd .= " -slim" if $slim;

  my (@out_files, @err_files, @gff3_files);

  my @base_bsub_opts = (-M => 3000000,
                        -R => 'select[mem>=3000] rusage[mem=3000]');
                        
  for(my $batch_idx = 0; $batch_idx < @batches; $batch_idx++) {
    my @sl = map { $_->seq_region_name } @{$batches[$batch_idx]};

    my $gff3_file = "$dumpdir/wormbase_gff3_dump.batch_${batch_idx}.gff3";
    my $lsf_out = "$dumpdir/wormbase_gff3_dump.lsfreport.${batch_idx}.out";
    my $lsf_err = "$dumpdir/wormbase_gff3_dump.lsfreport.${batch_idx}.err";

    my @local_bsub_opts = (@base_bsub_opts, 
                           -J => "wormbase_gff3_dump.$$.$batch_idx",
                           -o => $lsf_out,
                           -e => $lsf_err);


    my $local_cmd = sprintf("%s -slice %s -outfile $gff3_file", $base_cmd, join(",", @sl), $gff3_file);

    $lsf_man->submit(@local_bsub_opts, $local_cmd);

    push @out_files, $lsf_out;
    push @err_files, $lsf_err;
    push @gff3_files, $gff3_file;
  }
                           
  $lsf_man->wait_all_children( history => 1 );
  sleep(5);

  my @errors;
  foreach my $job ($lsf_man->jobs) {
    if ($job->history->exit_status != 0) {
      push @errors, $job;
    }
  }
  
  if (@errors) {
    die sprintf("%d jobs failed! Not bothering to collate output. Check output files!\n", scalar(@errors));
  }

  foreach my $gff3 (@gff3_files) {
    open(my $in_gff3, $gff3) or die "Could not open $gff3 for reading\n";
    while(<$in_gff3>) {
      next if /^\#\#gff-version/;
      print $out_fh $_;
    }
  }
  close($out_fh);

  unlink @gff3_files, @out_files, @err_files;
}


# CIGAR to old GFF3 CIGAR format converter
sub cigar_to_almost_cigar {
    my $i=shift;
    $i=~s/(\d*)(\w)/$2$1 /g;
    return $i;
}

#convert refernece strand for Lincoln
sub flipCigarReference {
  my $i=shift;
  $i=~tr/ID/DI/;
  return $i;
}

# reverse cigar string
sub reverse_cigar {
  my $i=shift;
  my @pairs=$i=~/(\d*[MIDFR])/g;
  my $reversed_cigar = join('',reverse @pairs);
  return $i;
}

# print the feature using some funky template
sub dump_feature {
  my $i=shift;
  my %feature=%{$i};
  my $gff_line= join("\t", 
                     $feature{seqname},
                     $feature{gff_source},
                     $feature{feature_type},
                     $feature{start},
                     $feature{end},
                     $feature{score},
                     $feature{strand},
                     $feature{phase});
  my @group;
  push @group, "ID=$feature{gff_id}" if $feature{gff_id};
  push @group, "Name=$feature{display}" if $feature{display};
  push @group, "Target=$feature{hit_id} $feature{hit_start} $feature{hit_end}" if $feature{hit_id};
  push @group, "Gap=$feature{cigar}" if $feature{cigar};

  if (exists $feature{attributes}) {
    foreach my $k (sort keys %{$feature{attributes}}) {
      push @group, "$k=" . $feature{attributes}->{$k};
    }
  }

  $gff_line .= "\t" . join(";", @group);
  $gff_line .= "\n";

  return $gff_line;
}


# build the info tag including protein features and interpro
sub get_info {
  my $transcript= shift;
  my $info = '';
  my @info;
  
  if (defined $transcript->translation) {
    # get all protein_features on the transcript
    my $features=$transcript->translation->get_all_ProteinFeatures();
    # get logic_name and hit_id
    my %domains;
    foreach my $feat (@$features) {
      if ($feat->interpro_ac) {
        $domains{$feat->interpro_ac} = 1;
      }
    }

    foreach my $acc (sort keys %domains) {
      my $dbe = $ensdb->get_DBEntryAdaptor->fetch_by_db_accession('InterPro', $acc );
      if (defined $dbe and $dbe->description) {
        my $desc = $dbe->description;
        $desc =~ s/\&/\%26/g;
        $desc =~ s/\,/\%2C/g;
        $desc =~ s/\=/\%3D/g;
        $desc =~ s/\;/\%3B/g;
        push @info, "method:InterPro accession:$acc description:$desc";
      }
    }

    $info = join(" %0A", @info);
  }
  return $info;
}

# print the gene including transcripts and exons
sub dump_gene {
  my ($gene) = @_;
  
  my $output = '';
  
  # Dump gene
  $output .= "# Gene " . $gene->{gff_id} . "\n";
  $output .= gff_line($gene->{seqname}, 
                      $gene->{gff_source}, 
                      'gene', 
                      $gene->{start}, 
                      $gene->{end},
                      $gene->{strand}, 
                      $gene->{gff_id},
                      undef, 
                      $gene->{display}, 
                      $gene->{attribs});
  
  # Dump transcripts
  my $parent = $gene->{gff_id};
  my %exon_parent;
  foreach my $transcript (@{$gene->{transcript}}) {
    $output .= gff_line($transcript->{seqname}, 
                        $transcript->{gff_source}, 
                        (exists $transcript->{cds}) ? 'mRNA' : $transcript->{gff_type}, 
                        $transcript->{start}, 
                        $transcript->{end},
                        $transcript->{strand}, 
                        $transcript->{gff_id}, 
                        undef,
                        ($transcript->{display} || undef), 
                        $transcript->{attribs},
                        $parent);
    
    # Store the parent of this transcript's exons
    foreach my $exon (@{$transcript->{exon}}) {
      ${$exon_parent{$exon->{gff_id}}}{$transcript->{gff_id}} = 1;
    }
  }
  
  # Dump exons
  foreach my $transcript (@{$gene->{transcript}}) {
    foreach my $exon (@{$transcript->{exon}}) {
      next if !$exon_parent{$exon->{gff_id}}; # If there are no parents then we've already dumped this exon
      my @parents = keys %{$exon_parent{$exon->{gff_id}}};
      delete $exon_parent{$exon->{gff_id}};
      $output .= gff_line($exon->{seqname}, 
                          $transcript->{gff_source}, 
                          'exon', 
                          $exon->{start}, 
                          $exon->{end},
                          $exon->{strand}, 
                          $exon->{gff_id}, 
                          undef,
                          undef,
                          undef,
                          @parents);
    }
  }
  
  # Dump coding_exons
  foreach my $transcript (@{$gene->{transcript}}) {
    my $parent = $transcript->{gff_id};
    if (exists $transcript->{utr5}) {
      foreach my $utr_seg (@{$transcript->{utr5}}) {
        $output .= gff_line($utr_seg->{seqname},
                            $transcript->{gff_source}, 
                            'five_prime_UTR', 
                            $utr_seg->{start}, 
                            $utr_seg->{end}, 
                            $utr_seg->{strand}, 
                            undef, 
                            undef, 
                            undef, 
                            undef, 
                            $transcript->{gff_id});
      }
    }
    if (exists $transcript->{cds}) {
      foreach my $cds (@{$transcript->{cds}}) {
        $output .= gff_line($cds->{seqname}, 
                            $transcript->{gff_source}, 
                            'CDS', 
                            $cds->{start}, 
                            $cds->{end},
                            $cds->{strand}, 
                            $cds->{gff_id}, 
                            $cds->{phase},
                            undef,
                            undef,
                            $transcript->{gff_id});
      }
    }
    if (exists $transcript->{utr3}) {
      foreach my $utr_seg (@{$transcript->{utr3}}) {
        $output .= gff_line($utr_seg->{seqname}, 
                            $transcript->{gff_source}, 
                            'three_prime_UTR', 
                            $utr_seg->{start}, 
                            $utr_seg->{end}, 
                            $utr_seg->{strand}, 
                            undef, 
                            undef, 
                            undef, 
                            undef, 
                            $transcript->{gff_id});
      }
    }
  }
  
  $output .= "###\n";
  
  return $output;
}

# a template for a GFF line
sub gff_line {
  my ($seqid, $source, $type, $start, $end, $strand, $stable_id, $phase, $name, $attribs, @parents) = @_;
  
  my $output = '';
  
  $phase='.' unless defined $phase;
  $strand= $strand>0? '+' : '-';
 
  $output .= "$seqid\t$source\t$type\t$start\t$end\t.\t$strand\t$phase\t";
  my @tags;
  push @tags, "ID=$stable_id" if defined $stable_id;
  if (@parents) {
    my $parent = join(',', @parents);
    push @tags, "Parent=$parent";
  }
  push @tags, "Name=$name" if defined $name;
  if (defined $attribs) {
    foreach my $k (sort keys %$attribs) {
      if ($attribs->{$k}) {
        push @tags, "$k=" . $attribs->{$k};
      }
    }
  }

  $output .= join(';', @tags);
  $output .= "\n";
  
  return $output;
}

# remove < 75% of evalue features from 100bp windows
sub filter_features {
  my ($features, $logic)=@_;
  my %f_features;
  
  # 50bp bin size
  my $size=100;
  
  # should I bin them instead?
  my @bins;
  
  # put the feature in all bins between start and stop
  foreach my $f(@$features){
    my ($_start,$_end)=sort {$a <=> $b } ($f->{start},$f->{end});
    for (my $n=int($_start/$size);$n <=int($_end/$size);$n++){
      push @{$bins[$n]},$f;
    }
  }
  
  # get the best 5 within 25% of the best hsp and add them to a hash
  my $c=0;
  foreach my $bin(@bins) {
    $c++;
    next unless $bin; # skip empty bins
    my $best=0;
    my $max_hsp=0;
    my @sorted_bin = sort {$a->{p_value} <=> $b->{p_value} or $b->{score} <=> $a->{score} } @$bin;

    $best=&p_value($sorted_bin[0]->{p_value});
    map { $f_features{$_->{gff_id}}=$_ if (&p_value($_->{p_value}) > $best*0.75 && $max_hsp++ <5)} @sorted_bin; 
  }
  
  # attempt to group the hits by id so that we can "join them up"
  my (%hits_by_name, @final_list);
  foreach my $hit (values %f_features) {
    push @{$hits_by_name{$hit->{hit_id}}->{$hit->{strand}}}, $hit;
  }
  foreach my $hid (keys %hits_by_name) {
    foreach my $strand (keys %{$hits_by_name{$hid}}) {
      my @hits = sort { $a->{start} <=> $b->{start} } @{$hits_by_name{$hid}->{$strand}};

      my @grouped_hits;
      while(my $hit = shift @hits) {
        if (@grouped_hits and
            $hit->{start} > $grouped_hits[-1]->[-1]->{end} and
            (($strand eq '+' and $hit->{start} > $grouped_hits[-1]->[-1]->{end}) or
             ($strand eq '-' and $hit->{end} < $grouped_hits[-1]->[-1]->{start}))) {
          push @{$grouped_hits[-1]}, $hit;
        } else {
          push @grouped_hits, [$hit];
        }
      }
      foreach my $group (@grouped_hits) {
        my $id = $logic . "." . $group->[0]->{gff_id};
        foreach my $hit (@$group) {
          $hit->{id} = $id;
          push @final_list, $hit;
        }
      }
    }
  }

  return @final_list;
}

# p_value shorthand
sub p_value {
  my ($p)=@_;
  my $log = (eval($p) > 0 ? -(log(eval($p)))/log(10) : 999.9); # if e-value is zero set it to 999.9
  return $log;
}

##
sub get_feature_logics {
  my ($db, $table, $slice) = @_;

  my %logics;

  my $sql = "SELECT distinct logic_name from $table,analysis WHERE $table.analysis_id = analysis.analysis_id";
  if (defined $slice) {
    my $sid = $slice->get_seq_region_id();
    $sql .= " AND seq_region_id = $sid";
  }
  my $sth = $db->dbc->prepare($sql);
  $sth->execute();
  while(my ($row) = $sth->fetchrow_array) {
    $logics{$row} = 1;
  }

  return sort keys %logics;
}
