#!/usr/bin/perl
# based on the ZFish GFF3 VEGA dumper from Ian Sealy
# the Blast Filter is lifted from the WormBlast processing pipeline
#
# To lower the memory footprint and converts the rather large
# Ensembl objects into lighter hashes for later use


use strict;
use warnings;

use lib '../lib';
use lib '/software/worm/ensembl/ensembl/modules';
use lib $ENV{'CVS_DIR'}."/ENSEMBL/lib";

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use IO::File;
use Domain2Interpro;
use Getopt::Long;

use LSF RaiseError => 0, PrintError => 1, PrintOutput => 1;
use LSF::JobManager;


my $debug =1;

# Vega database
my $dbname = 'bmalayi2';
my $dbhost = 'eagle';
my $dbport =  3306;
my $dbuser = 'wormro'; 
my $dbpass = '';
my $dumpdir = ".";

my (@dump_slice, $num_jobs, $out_file, $slim, $out_fh);

GetOptions(
  'dbhost=s'     => \$dbhost,
  'dbname=s'     => \$dbname,
  'dbuser=s'     => \$dbuser,
  'dbpass=s'     => \$dbpass,
  'dbport=s'     => \$dbport,
  'slice=s@'     => \@dump_slice,
  'submit=i'     => \$num_jobs,
  'dumpdir=s'    => \$dumpdir,
  'outfile:s'    => \$out_file,
  'slim'         => \$slim,
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

my @slices;
@dump_slice = split(/,/,join(',',@dump_slice));

if (@dump_slice) {
  foreach (@dump_slice) {
    push @slices, $sa->fetch_by_region('toplevel',$_);
  }
}
else {
  @slices = sort { $b->length <=> $a->length } @{$sa->fetch_all('toplevel')};
}

if ($num_jobs) {
  &submit_and_collate();
  exit(0);
}

my $mapper = Domain2Interpro->new();

if (defined $out_file) {
  open($out_fh, ">$out_file") or die "Could not open $out_file for writing\n";
}


while( my $slice = shift @slices) {
  my $slice_name = $slice->seq_region_name();
  my $slice_size = $slice->length;

  print STDERR "Proccesing slice ", $slice->seq_region_name, "\n";

  if (not defined $out_fh) {
    my $sl_file = "$dumpdir/${slice_name}.gff3";
    open($out_fh, ">$sl_file") or die "Could not open $sl_file for writing\n";  
  }

	
  print $out_fh "##gff-version 3\n";
  print $out_fh "##sequence-region $slice_name 1 $slice_size\n";
  
  # Get all the genes on this slice
  my $genes = $slice->get_all_Genes();
  while( my $gene=shift @$genes) {
    my $gene_stable_id = 'gene:'.$gene->stable_id();
    
    my %gene_to_dump = (
      stable_id => $gene_stable_id,
      name      => $slice_name,
      start     => $gene->seq_region_start(),
      end       => $gene->seq_region_end(),
      strand    => $gene->strand(),
      note      => ($gene->status()||'PREDICTED' ). " " . $gene->biotype(),
      public_name => $gene->stable_id(),
      display     => $gene->stable_id(), # wrong but fixes db's without xref_mapping
      gff_source  => (defined $gene->analysis->gff_source) ? $gene->analysis->gff_source : "WormBase",
        );
    
    # get all transcripts of the gene
    my $all_transcripts = $gene->get_all_Transcripts();
    while( my $transcript = shift @{$all_transcripts}) {
      
      my $transcript_stable_id = 'transcript:'.$transcript->stable_id();

      my $tr_obj =  {
        stable_id => $transcript_stable_id,
        name      => $slice_name,
        start     => $transcript->start(),
        end       => $transcript->end(),
        strand    => $transcript->strand(),
        note      => $transcript->biotype(),
        info      => get_info($transcript),
        public_name => $transcript->stable_id(),
        display     => $transcript->stable_id(),
        gff_source  => (defined $transcript->analysis->gff_source) ?  $transcript->analysis->gff_source : "WormBase",
        gff_type    => (defined $transcript->analysis->gff_feature ) ? $transcript->analysis->gff_feature : $transcript->biotype . '_primary_transcript',
      };

      my $translation = $transcript->translation;

      if (defined $translation) {
        my $translation_id = 'cds:'.$transcript->translation->stable_id();

        $tr_obj->{cds_start} = $translation->genomic_start();
        $tr_obj->{cds_end}   = $translation->genomic_end();
        $tr_obj->{translation_stable_id} = $translation_id;

        my $all_t_exons = $transcript->get_all_translateable_Exons();
        
        while (my $cds = shift @{$all_t_exons}) {
          push @{$tr_obj->{'cds'}}, {
            stable_id => $translation_id,
            name      => $slice_name,
            start     => $cds->seq_region_start(),
            end       => $cds->seq_region_end(),
            strand    => $cds->strand(),
            phase     => (3-$cds->phase())%3, # phase/frame conversion to a sane system
          };
        }

        my $cdna_coding_start = $transcript->cdna_coding_start();
        if ($cdna_coding_start > 1) {
          my @coords = $transcript->cdna2genomic(1, $cdna_coding_start - 1);
          @coords = grep { $_->isa("Bio::EnsEMBL::Mapper::Coordinate") } @coords;
          foreach my $coord (@coords) {
            push @{$tr_obj->{'utr5'}}, {
              name => $slice_name,
              start => $coord->start, 
              end   => $coord->end,
              strand => $coord->strand,
            };
          }
        }

        my $cdna_coding_end = $transcript->cdna_coding_end();
        if ($cdna_coding_end < $transcript->length) {
          my @coords = $transcript->cdna2genomic($cdna_coding_end + 1, $transcript->length);
          @coords = grep { $_->isa("Bio::EnsEMBL::Mapper::Coordinate") } @coords;
          foreach my $coord (@coords) {
            push @{$tr_obj->{'utr3'}}, {
              name => $slice_name,
              start => $coord->start, 
              end   => $coord->end,
              strand => $coord->strand,
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
        my $exon_stable_id = sprintf("exon:%s.%d", $transcript->stable_id, $exon_count++);
        
        push @{$tr_obj->{'exon'}}, {
          stable_id => $exon_stable_id,
          name      => $slice_name,
          start     => $exon->seq_region_start(),
          end       => $exon->seq_region_end(),
          strand    => $exon->strand(),
        };
      }
      

      push @{$gene_to_dump{'transcript'}}, $tr_obj;
    }
    print $out_fh dump_gene(\%gene_to_dump);
  }
  
  next if $slim;
  # get all protein align features on the slice
  my %blastx_features;
  
  my $features = $slice->get_all_ProteinAlignFeatures;
  
  while(my $feat = shift @$features) {
    my $cigar_line = flipCigarReference($feat->cigar_string); # for Lincoln
    if ($feat->strand < 0) {
      $cigar_line = reverse_cigar($cigar_line);
    }
    $cigar_line = cigar_to_almost_cigar($cigar_line);

    my $stripped_feature = {
      hit_id      => $feat->hseqname, 
      target_id   => $slice->seq_region_name,
      target_start=> $feat->hstart,
      target_stop => $feat->hend,
      strand      => ($feat->strand > 0?'+':'-'),
      hit_start   => $feat->start,
      hit_stop    => $feat->end,
      score       => $feat->score,
      p_value     => $feat->p_value,
      dbid        => $feat->dbID,
      logic_name  => $feat->analysis->logic_name,
      gff_source  => ($feat->analysis->gff_source) ? $feat->analysis->gff_source : $feat->analysis->logic_name,
      cigar       => $cigar_line,
      feature_type=> 'protein_match',
    };
    push @{$blastx_features{$stripped_feature->{logic_name}}}, $stripped_feature;
  }
  
  while (my($k,$v)=each %blastx_features){
    my @filtered_features=filter_features($v, $k);
    map {print $out_fh dump_feature($_)} @filtered_features;
  }
  
  # get all dna align features on the slice
  
  $features=$slice->get_all_DnaAlignFeatures();

  while(my $feat = shift @$features) {
    my $cigar_line = flipCigarReference($feat->cigar_string); # for Lincoln
    if ($feat->strand < 0) {
      $cigar_line = reverse_cigar($cigar_line);
    }
    $cigar_line = cigar_to_almost_cigar($cigar_line);
    
    my $stripped_feature = {
      hit_id      => $feat->hseqname,
      target_id   => $feat->slice->seq_region_name,
      target_start=> $feat->hstart,
      target_stop => $feat->hend,
      strand      => ($feat->strand > 0?'+':'-'),
      hit_start   => $feat->start,
      hit_stop    => $feat->end,
      score       => $feat->score,
      p_value     => $feat->p_value,
      dbid        => $feat->dbID,
      logic_name  => $feat->analysis->logic_name,
      gff_source  => ($feat->analysis->gff_source) ? $feat->analysis->gff_source : $feat->analysis->logic_name,
      cigar       => $cigar_line,
      feature_type=> 'nucleotide_match',
    };
    print $out_fh dump_feature($stripped_feature);
  }
  
  
  # get all repeat features on the slice
  my $repeats = $slice->get_all_RepeatFeatures;
  foreach my $feature (@$repeats){
    my $stripped_feature = {
      target_id   => $feature->slice->seq_region_name,
      strand      => ($feature->strand > 0?'+':'-'),
      hit_start   => $feature->seq_region_start,
      hit_stop    => $feature->seq_region_end,
      score       => ($feature->score||'.'),
      dbid        => $feature->dbID,
      logic_name  => $feature->analysis->logic_name,
      gff_source  => (defined $feature->analysis->gff_source) ? $feature->analysis->gff_source : "WormBase",
      feature_type=> 'repeat_region',
    };
    print $out_fh dump_feature($stripped_feature);
  }
  

  # get all simple features stored in the database (Operons etc. etc.)
  my $simp_features = $slice->get_all_SimpleFeatures;
  foreach my $simpfeature (@$simp_features){
    my $stripped_simpfeature = {
      target_id   => $simpfeature->slice->seq_region_name,
      strand      => ($simpfeature->strand > 0?'+':'-'),
      hit_start   => $simpfeature->seq_region_start,
      hit_stop    => $simpfeature->seq_region_end,
      score       => ($simpfeature->score||'.'),
      dbid        => $simpfeature->dbID,
      logic_name  => $simpfeature->analysis->logic_name,
      gff_source  => (defined $simpfeature->analysis->gff_source) ? $simpfeature->analysis->gff_source : "WormBase",
      feature_type=> $simpfeature->analysis->logic_name,
    };
    print $out_fh dump_feature($stripped_simpfeature);
  }
  
  print $out_fh '#'x80;
  print $out_fh "\n";
  close($out_fh) unless defined $out_file;
}

# close the file handle of the primary gff as you are done
close($out_fh) if defined $out_file;


############################

sub submit_and_collate {

  if (not defined $out_file) {
    $out_file = "$dumpdir/collated_dump.gff3";
  }
  open(my $out_fh, ">$out_file") or die "Could not open $out_file for writing\n";

  my (@batches, $idx);
  while(@slices) {
    push @{$batches[$idx]}, shift @slices;

    $idx++; $idx = 0 if $idx == $num_jobs;
  }

  my $lsf_man = LSF::JobManager->new();

  my $base_cmd = "perl $0 -dbhost $dbhost -dbname $dbname -dbuser $dbuser -dbport $dbport -dumpdir $dumpdir";
  $base_cmd .= " -dbpass $dbpass" if $dbpass;
  $base_cmd .= " -slim" if $slim;

  my (@out_files, @err_files, @gff3_files);

  my @base_bsub_opts = (-M => 1000000,
                        -R => 'select[mem>=1000] rusage[mem=1000]');
                        
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
    my $gff_line=
	"$feature{target_id}\t$feature{gff_source}\t$feature{feature_type}\t$feature{hit_start}\t$feature{hit_stop}\t".
	"$feature{score}\t$feature{strand}\t.".
	($feature{id} ? "\tID=$feature{id}"  : "\tID=$feature{logic_name}.$feature{dbid}").
	($feature{cigar}?";Name=$feature{hit_id};Target=$feature{hit_id} $feature{target_start} $feature{target_stop};Gap=$feature{cigar}\n":"\n");
    return $gff_line;
}


# build the info tag including protein features and interpro
sub get_info {
  my $transcript= shift;
  my $info='';

  if (defined $transcript->translation) {
    # get all protein_features on the transcript
    my $features=$transcript->translation->get_all_ProteinFeatures();
    # get logic_name and hit_id
    
    my %plain_features;
    my $rest_features;
    map {
      if ($mapper->get_method2database($_->analysis->logic_name())) {
        push @{$plain_features{$_->analysis->logic_name()}},
        [$_->display_id(),$_->start(), $_->end(),$_->hstart(),$_->hend(),$_->score(),$_->p_value()]
      }
      else {
        push @$rest_features,$_
      }
    } @$features;
    my @interpros=$mapper->get_mapping(\%plain_features);
    map {$info.=sprintf( "position:%d-%d method:%s accession:%s description:%s %%0A", $_->[1], $_->[2],
                         'InterPro', $_->[0] , $_->[7]) if $_->[1]} @interpros;
    
    while ( my $pfeature = shift @$rest_features ) {
      my $logic_name = $pfeature->analysis()->logic_name();
      $info.=sprintf( "position:%d-%d %s method:%s accession:%s %%0A", $pfeature->start(), $pfeature->end(), 
                      $pfeature->p_value(),$logic_name, $pfeature->display_id());
    }
  }
  return $info;
}

# print the gene including transcripts and exons
sub dump_gene {
  my ($gene) = @_;
  
  my $output = '';
  
  # Dump gene
  $output .= "# Gene " . $gene->{'stable_id'} . "\n";
  $output .= gff_line(
    $gene->{'name'}, $gene->{'gff_source'}, 'gene', $gene->{'start'}, $gene->{'end'},
    $gene->{'strand'}, $gene->{'stable_id'},undef, $gene->{'display'}, $gene->{'note'},undef,$gene->{'public_name'}
      );
  
  # Dump transcripts
  my $parent = $gene->{'stable_id'};
  my %exon_parent;
  foreach my $transcript (@{$gene->{'transcript'}}) {
    $output .= gff_line(
      $transcript->{'name'}, 
      $transcript->{'gff_source'}, 
      (exists $transcript->{cds}) ? 'mRNA' : $transcript->{gff_type}, 
      $transcript->{'start'}, 
      $transcript->{'end'},
      $transcript->{'strand'}, 
      $transcript->{'stable_id'}, 
      undef,
      ($transcript->{'display'} || undef), 
      undef, 
      ($transcript->{info}||undef),
      $transcript->{'public_name'},$parent);
    
    # Store the parent of this transcript's exons
    foreach my $exon (@{$transcript->{'exon'}}) {
      ${$exon_parent{$exon->{'stable_id'}}}{$transcript->{'stable_id'}} = 1;
    }
  }
  
  # Dump exons
  foreach my $transcript (@{$gene->{'transcript'}}) {
    foreach my $exon (@{$transcript->{'exon'}}) {
      next if !$exon_parent{$exon->{'stable_id'}}; # If there are no parents then we've already dumped this exon
      my @parents = keys %{$exon_parent{$exon->{'stable_id'}}};
      delete $exon_parent{$exon->{'stable_id'}};
      $output .= gff_line(
        $exon->{'name'}, $transcript->{'gff_source'}, 'exon', $exon->{'start'}, $exon->{'end'},
        $exon->{'strand'}, $exon->{'stable_id'}, undef,undef, undef,undef,undef ,@parents);
    }
  }
  
  # Dump coding_exons
  foreach my $transcript (@{$gene->{'transcript'}}) {
    my $parent = $transcript->{'stable_id'};
    if (exists $transcript->{'utr5'}) {
      foreach my $utr_seg (@{$transcript->{'utr5'}}) {
        $output .= gff_line(
          $utr_seg->{'name'}, $transcript->{'gff_source'}, 'five_prime_UTR', $utr_seg->{'start'}, 
          $utr_seg->{'end'}, $utr_seg->{'strand'}, undef, undef, undef, 
          undef, undef, undef, $transcript->{'stable_id'});
      }
    }
    if (exists $transcript->{'cds'}) {
      foreach my $cds (@{$transcript->{'cds'}}) {
        $output .= gff_line(
          $cds->{'name'}, $transcript->{'gff_source'}, 'CDS', $cds->{'start'}, $cds->{'end'},
          $cds->{'strand'}, $cds->{'stable_id'}, $cds->{'phase'},undef, 
          undef,undef ,undef, $transcript->{'stable_id'});
      }
    }
    if (exists $transcript->{'utr3'}) {
      foreach my $utr_seg (@{$transcript->{'utr3'}}) {
        $output .= gff_line(
          $utr_seg->{'name'}, $transcript->{'gff_source'}, 'three_prime_UTR', $utr_seg->{'start'}, 
          $utr_seg->{'end'}, $utr_seg->{'strand'}, undef, undef, undef, 
          undef, undef, undef, $transcript->{'stable_id'});
      }
    }
  }
  
  $output .= "###\n";
  
  return $output;
}

# a template for a GFF line
sub gff_line {
  my ($seqid, $source, $type, $start, $end, $strand, $stable_id, $phase,$name, $note, $info,$public,@parents) = @_;
  
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
  push @tags, "Note=$note" if defined $note;
  push @tags, "info=$info" if defined $info;
  push @tags, "public_name=$public" if defined $public;
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
    my ($_start,$_end)=sort {$a <=> $b } ($f->{hit_start},$f->{hit_stop});
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
    map { $f_features{$_->{dbid}}=$_ if (&p_value($_->{p_value}) > $best*0.75 && $max_hsp++ <5)} @sorted_bin; 
  }
  
  # attempt to group the hits by id so that we can "join them up"
  my (%hits_by_name, @final_list);
  foreach my $hit (values %f_features) {
    push @{$hits_by_name{$hit->{hit_id}}->{$hit->{strand}}}, $hit;
  }
  foreach my $hid (keys %hits_by_name) {
    foreach my $strand (keys %{$hits_by_name{$hid}}) {
      my @hits = sort { $a->{hit_start} <=> $b->{hit_start} } @{$hits_by_name{$hid}->{$strand}};

      my @grouped_hits;
      while(my $hit = shift @hits) {
        if (@grouped_hits and
            $hit->{hit_start} > $grouped_hits[-1]->[-1]->{hit_stop} and
            (($strand eq '+' and $hit->{target_start} > $grouped_hits[-1]->[-1]->{target_stop}) or
             ($strand eq '-' and $hit->{target_stop} < $grouped_hits[-1]->[-1]->{target_start}))) {
          push @{$grouped_hits[-1]}, $hit;
        } else {
          push @grouped_hits, [$hit];
        }
      }
      foreach my $group (@grouped_hits) {
        my $id = $logic . "." . $group->[0]->{dbid};
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



