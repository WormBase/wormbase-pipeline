
package Physical_mapper;

use strict;
use warnings;
use YAML;

sub new {
  my ( $class, $yfile, $acedb, @files ) = @_;

  my $self = {};
  bless $self, $class;
  
  my ($map, $gmap, %sorted_map);

  if ($yfile && -e $yfile) { 
    $self->thaw($yfile);
  } elsif(grep {-e $_} @files)  {
    my ($pmap, $gmap) = Map_func::build($acedb, @files );    
    $self->gmap($gmap);
    $self->pmap($pmap);    
  } else{
    die("ERROR: cannot build mapper from @files\n");
  }

  foreach my $key ( keys %{$self->pmap} ) {
    @{ $sorted_map{$key} } = sort { $a <=> $b } keys %{ $self->pmap->{$key} };
  }
  $self->sorted_pmap(\%sorted_map);
  
  return $self;
}

sub map {
  my ( $self, $pos, $chr ) = @_;

  if ($pos < $self->sorted_pmap->{$chr}->[0]) { 
    # before first marker
    my $current   = $self->sorted_pmap->{$chr}->[1];
    my $fake_last = $self->sorted_pmap->{$chr}->[0];
    my $pdiff     = ($fake_last - $pos) / ( $current - $fake_last );
    my $mlast    = $self->pmap->{$chr}->{$fake_last}->[0];
    my $mcurrent = $self->pmap->{$chr}->{$current}->[0];
    my $mpos = $mlast - (( $mcurrent - $mlast ) * $pdiff);
    return $mpos;
  } 
  elsif ($pos > $self->sorted_pmap->{$chr}->[-1]){
    # after last marker
    my $current  = $self->sorted_pmap->{$chr}->[-2];
    my $fake_last = $self->sorted_pmap->{$chr}->[-1];
    my $pdiff    = ( $pos - $fake_last ) / ( $fake_last - $current );
    my $mlast    = $self->pmap->{$chr}->{$fake_last}->[0];
    my $mcurrent = $self->pmap->{$chr}->{$current}->[0];
    my $mpos = $mlast + (($mlast - $mcurrent) * $pdiff );
    return $mpos;
  } else {
    #.................................
    my $last = $self->sorted_pmap->{$chr}->[0];

    for ( my $i = 1 ; $i < scalar @{ $self->sorted_pmap->{$chr} } ; $i++ ) {
      my $next  = $self->sorted_pmap->{$chr}->[$i];

      if ( $pos <= $next && $pos >= $last ) {
        my $pdiff    = ( $pos - $last ) / ( $next - $last );
        my $mlast    = $self->pmap->{$chr}->{$last}->[0];
        my $mnext = $self->pmap->{$chr}->{$next}->[0];
        my $mpos = ( $mnext - $mlast ) * $pdiff + $mlast;
        return $mpos;
      }
      else {
        $last = $next;
      }
    }
  }
  # glorious correction routine for the ends
  # should be f(x)=dx/dy + x1
  # my $x1 = $last
  # my $y1 = $self->{'pmap'}->{$chr}->{$last};
  # my $x2 = $next
  # my $y2 = $self->{'pmap'}->{$chr}->{$next};
  #    my $f = sub { my $x = shift; return ( $x1 - $x2 ) / ( $y1 - $y2 ) + $x1 }
  
  return undef;
}

sub x_to_ace {
  my ( $self, $id, $map, $chr, $x ) = @_;

  my $mpos = $self->map( $map, $chr );
  if ($mpos) {
    return "$x : \"$id\"\nInterpolated_map_position \"$chr\" $mpos\n\n";
  }
}

sub thaw {
  my ($self, $yfile) = @_;

  my %obj = YAML::LoadFile($yfile);
  $self->pmap($obj{pmap});
  $self->gmap($obj{gmap});
}

sub freeze {
  my($self,$file)=@_;

  my $obj = {
    gmap => $self->gmap,
    pmap  => $self->pmap,
  };
  YAML::DumpFile( $file, %$obj);
}

sub check_and_fix_mapping {
  my ( $self, $database, $fixes_file, $fixes_log, $attempt_fix, $log ) = @_;
  
  open(my $revh, ">$fixes_log") or $log->log_and_die("Could not open $fixes_log for writing\n");
  
  my $pmap = $self->pmap;
  my $smap = $self->sorted_pmap;
  
  my $errors = 0;
  if ($attempt_fix) {        
    open(my $aceout, ">$fixes_file") or $log->log_and_die("Could not open $fixes_file for writing\n");

    foreach my $key ( sort keys %{$pmap} ) {
      my @genes;
      
      foreach my $i ( @{ $smap->{$key} } ) {
        push @genes, {
          chr       => $key, 
          ppos      => $i,
          orig_gpos => $pmap->{$key}->{$i}->[0],
          gpos      => $pmap->{$key}->{$i}->[0],
          gene      => $pmap->{$key}->{$i}->[1],
          landmark  => $pmap->{$key}->{$i}->[2],
        };
      }
      
      $log->write_to("Inconsistencies present for $key before attempted fix:\n");
      $self->report_inconsistencies($key, $log);
      $log->write_to("-----------------------------------\n");
      
      &fix_gmap( \@genes );       

      foreach my $col (@genes) {
        next if $col->{gpos} == $col->{orig_gpos};
        
        $pmap->{$key}->{ $col->{ppos} }->[0] = $col->{gpos}; # change gmap
        my $_chrom = $key;
        my $_pos   = $col->{gpos};
        my $_gene  = $col->{gene};
        
        # create acefile
        print $aceout "\n";
        print $aceout "Gene : $_gene\n";
        print $aceout "Map $_chrom Position $_pos\n";
      }

      $errors += $self->report_inconsistencies($key, $log, $revh);
    }
    close($aceout);
  } else {
    foreach my $key ( sort keys %{$pmap} ) {
      $errors += $self->report_inconsistencies($key, $log, $revh);
    }
  }    
  
  close($revh);

  return $errors;
}

sub report_inconsistencies {
  my ($self, $chr, $log, $revfh) = @_;

  my $pmap = $self->pmap;
  my $smap = $self->sorted_pmap;
    
  my $last;
  my $errors = 0;

  foreach my $i ( @{ $smap->{$chr} } ) {   # sorted pmap positions
    if ( $last and $pmap->{$chr}->{$i}->[0] < $pmap->{$chr}->{$last}->[0] ) {
      if (defined $revfh) {
        print $revfh "----------------------------------\n";
        printf($revfh "%s: %d\t%s\t%s ERROR (conflict with last line)", 
               $chr, 
               $i, 
               $pmap->{$chr}->{$i}->[0],
               $pmap->{$chr}->{$i}->[1]);
      }
      $log->write_to( "$chr: $i\t"
                      . $pmap->{$chr}->{$i}->[0] . "\t"
                      . $pmap->{$chr}->{$i}->[1]
                      . " ERROR (conflict with last line) \n" );
      if (defined $revfh) {
        print $revfh "----------------------------------\n";
      }
      $errors++;
    } else {
      if (defined $revfh) {
        printf($revfh "%s: %d\t%s\t%s\n",
               $chr, 
               $i,
               $pmap->{$chr}->{$i}->[0],
               $pmap->{$chr}->{$i}->[1]);
      }
    }
    
    $last = $i;
  }

  return $errors;
}

# class functions for Physical_mapper

sub fix_gmap {
  my ($genes) = @_;
  
  #
  # Find genes that have map positions inconsistent with the landmarks
  # These are the ones we want to fix in the first pass
  #
  for(my $i = 0; $i < @$genes; $i++) {
    $genes->[$i]->{needs_fix} = 0;

    next if $genes->[$i]->{landmark};
    my ($left_lm, $right_lm);
    for(my $j = $i-1; $j >= 0; $j--) {
      if ($genes->[$j]->{landmark}) {
        $left_lm = $j;
        last;
      }
    }
    for(my $j = $i+1; $j < @$genes; $j++) {
      if ($genes->[$j]->{landmark}) {
        $right_lm = $j;
        last;
      }
    }
    if (defined $left_lm and 
        $genes->[$i]->{orig_gpos} < $genes->[$left_lm]->{orig_gpos}) {
      $genes->[$i]->{needs_fix} = 1;
    }
    if (defined $right_lm and 
        $genes->[$i]->{orig_gpos} > $genes->[$right_lm]->{orig_gpos}) {
      $genes->[$i]->{needs_fix} = 1;
    }
  }
  
  my $changed_in_this_iteration;
  my $iteration = 0;
  do {
    $iteration++;
    $changed_in_this_iteration = 0;
    my ($prev_gpos, $next_gpos, $prev_ppos, $next_ppos);
    
    for ( my $i = 0 ; $i < @$genes ; $i++ ) {

      my $ppos  = $genes->[$i]->{ppos};
      my $gpos  = $genes->[$i]->{gpos};
      my $gname = $genes->[$i]->{gene};
      my $orig_gmap    = $genes->[$i]->{orig_gpos};
      my $is_landmark  = $genes->[$i]->{landmark};

      # get the next position
      # are we at the end of the array or end of the chromosome?
      next if $i == 0;
      last if $i == scalar(@$genes) - 1;

      if ($iteration == 1) {
        for(my $j=$i-1; $j >= 0; $j--) {
          $prev_ppos = $genes->[$j]->{ppos};
          $prev_gpos = $genes->[$j]->{gpos};
          last if not $genes->[$j]->{needs_fix};
        }
        
        for(my $j=$i+1; $j < @$genes; $j++) {
          $next_ppos = $genes->[$j]->{ppos};
          $next_gpos = $genes->[$j]->{gpos};
          last if not $genes->[$j]->{needs_fix};
        }
      } else {
        $prev_ppos = $genes->[$i-1]->{ppos};
        $prev_gpos = $genes->[$i-1]->{gpos};
        $next_ppos = $genes->[$i+1]->{ppos};
        $next_gpos = $genes->[$i+1]->{gpos};
      }

      if ( ($iteration == 1 and $genes->[$i]->{needs_fix}) or
           ($iteration > 1 and $gpos < $prev_gpos)) {
        # get the difference between the previous and next positions
        my $frac = ($ppos - $prev_ppos) / ($next_ppos - $prev_ppos);
        my $gdiff = $frac * ($next_gpos - $prev_gpos);
        $genes->[$i]->{gpos} = $prev_gpos + $gdiff;
        $changed_in_this_iteration = 1;
        $gpos = $genes->[$i]->{gpos};
      }

      #printf(STDERR "%-5d\t%-5d\t%s\t%2s\t%s\t%s\t%-10s\t%.6f\t%.6f\n", 
      #       $iteration,
      #       $i,
      #       $genes->[$i]->{gene},
      #       $genes->[$i]->{chr},
      #       $genes->[$i]->{landmark},
      #       $genes->[$i]->{needs_fix},
      #       $genes->[$i]->{ppos},
      #       $genes->[$i]->{gpos}, 
      #       $genes->[$i]->{orig_gpos},
      #    );

    }
  } while ($changed_in_this_iteration);
}


sub gmap {
  my ($self, $val) = @_;
  if (defined $val) {
    $self->{_gmap} = $val;
  }
  return $self->{_gmap};
}

sub pmap {
  my ($self, $val) = @_;
  if (defined $val) {
    $self->{_pmap} = $val;
  }
  return $self->{_pmap};
}

sub sorted_pmap {
  my ($self, $val) = @_;
  if (defined $val) {
    $self->{_sorted_pmap} = $val;
  }
  return $self->{_sorted_pmap};
}

######################
package Map_func;

use Ace;

#
#
#
sub get_phys {
  my ($database) = @_;
  my %map;

  my $db = Ace->connect( -path => $database );

  #push @genes, $db->find('find Gene * where Map AND SMap AND NOT Pseudo_map_position');
  my @genes = $db->find('find Gene * where Map AND SMap');

  while(my $gene = shift @genes) {
    if ( !defined $gene->Map(2) ) {
      print STDERR "cannot find genetic map position for $gene ...\n";
      next;
    }
    my $name = "$gene";
    my $pos  = $gene->Map(3)->name;
    if (not $pos) {
      $pos = "0.0000";
    }
    my $landmark = 0;
    if ($gene->at('Map_info.Landmark_gene')) {
      $landmark = 1;
    }
    $map{$name} = [$pos, $landmark];
  }
  $db->close();
  return \%map;
}

# 
#
#
sub build {
  my ( $acedb, @infiles ) = @_;
  my %genes;    # gene -> phys_map
  my %gen_map = %{ get_phys($acedb) };
  
  foreach my $file (@infiles) {
    open IN, $file;
    while (<IN>) {
      next if /^\#/;
      s/\"//g;
      my @a = split(/\t/, $_);
      next if !( $a[1] eq 'gene' && $a[2] eq 'gene' );
      my ($gene_id) = $a[8] =~ /^\S+\s+(\S+)/;
      next if not exists $gen_map{$gene_id};
      
      my $chrom   = $a[0];
      my $map_pos = ( $a[3] + $a[4] ) / 2;
      my $gen_pos = $gen_map{$gene_id}->[0];
      my $is_landmark = $gen_map{$gene_id}->[1];
      $chrom =~ s/CHROMOSOME_//;
      
      #should be chromosome->map_pos->gen_pos
      $genes{$chrom}->{$map_pos} = [$gen_pos,$gene_id,$is_landmark];
    }
    close IN;
  }
  
  return \%genes, \%gen_map;
}

1;
