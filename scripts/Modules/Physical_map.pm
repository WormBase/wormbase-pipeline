
package Physical_mapper;

use strict;
use warnings;
use YAML;

sub new {
  my ( $class, $yfile, $acedb, @files ) = @_;

  my $self = {};
  bless $self, $class;
  
  my ($map, $gmap, %sorted_map);

  if ( defined $yfile and -e $yfile ) { 
    $self->thaw($yfile);
  } else {
    my ($pmap, $gmap) = Map_func::build($acedb, @files );    
    $self->gmap($gmap);
    $self->pmap($pmap);    
  }

  foreach my $key ( keys %{$self->pmap} ) {
    @{ $sorted_map{$key} } = sort { $a <=> $b } keys %{ $self->pmap->{$key} };
  }
  $self->sorted_pmap(\%sorted_map);
  
  return $self;
}

sub map {
  my ( $self, $pos, $chr ) = @_;
  my $last = 0;
  my $next = 0;
  # leaves alot for refactoring ...
  if ($pos < @{ $self->sorted_pmap->{$chr} }[0]) { # before first marker
    my $current  = @{ $self->sorted_pmap->{$chr} }[1];
    my $fake_last = @{ $self->sorted_pmap->{$chr} }[0];
    my $pdiff    = ( $pos - $fake_last ) / ( $current - $fake_last );    # might  be wrong prefix ->better now?
    my $mlast    = $self->pmap->{$chr}->{$fake_last}->[0];
    my $mcurrent = $self->pmap->{$chr}->{$current}->[0];
    my $mpos = ( $mcurrent - $mlast ) * $pdiff + $mlast;
    return $mpos;
  } 
  elsif ($pos>@{ $self->sorted_pmap->{$chr} }[-1]){
    my $current  = @{ $self->sorted_pmap->{$chr} }[-1];
    my $fake_last = @{ $self->sorted_pmap->{$chr} }[-2];
    my $pdiff    = ( $pos - $fake_last ) / ( $current - $fake_last );    # might  be wrong prefix ->better now?
    my $mlast    = $self->pmap->{$chr}->{$fake_last}->[0];
    my $mcurrent = $self->pmap->{$chr}->{$current}->[0];
    my $mpos = ( $mcurrent - $mlast ) * $pdiff + $mlast;
    return $mpos;
  }
  else {
    #.................................
    for ( my $i = 0 ; $i < scalar @{ $self->sorted_pmap->{$chr} } ; $i++ ) {
      my $current  = @{ $self->sorted_pmap->{$chr} }[$i];
      my $mlast    = $self->pmap->{$chr}->{$last}->[0];
      my $mcurrent = $self->pmap->{$chr}->{$current}->[0];
      
      $next = $current;
      if ( $pos <= $current && $pos >= $last ) {
        my $pdiff    = ( $pos - $last ) / ( $current - $last );    # might  be wrong prefix ->better now?
        my $mlast    = $self->pmap->{$chr}->{$last}->[0];
        my $mcurrent = $self->pmap->{$chr}->{$current}->[0];
        my $mpos = ( $mcurrent - $mlast ) * $pdiff + $mlast;
        return $mpos;
      }
      else {
        $last = $current;
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

sub check_mapping {
  my ( $self, $database, $fixes_file, $fixes_log, $log ) = @_;
  
  my $revh = IO::File->new( $fixes_log, "w" );
  my $acefile = IO::File->new($fixes_file,'w');

  $log->write_to("writing genetic map fixes to $fixes_file\n");
  $log->write_to("have a look at $fixes_log to resolve\n");
  $log->write_to("-----------------------------------\n");
  
  my $pmap = $self->pmap;
  my $smap = $self->sorted_pmap;

  my $errors = 0;
  foreach my $key ( sort keys %{$pmap} ) {

    # need to build a @genes list [chromosome,ppos,gpos,geneid]
    
    my @genes;
    
    foreach my $i ( @{ $smap->{$key} } ) {
      push @genes, [
        $key,
        $i,
        $pmap->{$key}->{$i}->[0],
        $pmap->{$key}->{$i}->[1],
        $pmap->{$key}->{$i}->[2],
        $pmap->{$key}->{$i}->[0]
      ];
    }
    
    # call genetic fix function
    
    &fix_gmap( \@genes );       
    
    foreach my $col (@genes) {
      next if $col->[2] == $col->[5];
      $pmap->{$key}->{ $col->[1] }->[0] = $col->[2]; # change gmap
      my $_chrom = $key;
      my $_pos   = $col->[2];
      my $_gene  = $col->[3];
      
      # create acefile
      print $acefile "\n";
      print $acefile "Gene : $_gene\n";
      print $acefile "Map $_chrom Position $_pos\n";
    }
    
    ##
    
    my $last;
    foreach my $i ( @{ $smap->{$key} } ) {   # sorted pmap positions
      if ( $last and $pmap->{$key}->{$i}->[0] < $pmap->{$key}->{$last}->[0] ) {
        print $revh "----------------------------------\n";
        printf($revh "%s: %d\t%s\t%s ERROR (conflict with last line)", 
               $key, 
               $i, 
               $pmap->{$key}->{$i}->[0],
               $pmap->{$key}->{$i}->[1]);
        $log->write_to( "$key: $i\t"
                        . $pmap->{$key}->{$i}->[0] . "\t"
                        . $pmap->{$key}->{$i}->[1]
                        . " ERROR (conflict with last line) \n" );
        print $revh "----------------------------------\n";
        $errors++;
      } else {
        printf($revh "%s: %d\t%s\t%s\n",
               $key, 
               $i,
               $pmap->{$key}->{$i}->[0],
               $pmap->{$key}->{$i}->[1]);
      }
      
      $last = $i;
    }
  }
  return $errors;
}

# class functions for Physical_mapper

sub fix_gmap {
  my ($genes) = @_;
  
  my $changed_in_this_iteration;
  do {
    $changed_in_this_iteration = 0;
    my $prev_chrom = "";
    my $prev_pos;
    my $next_pos;
    for ( my $i = 0 ; $i < @$genes ; $i++ ) {

      my $chrom = $$genes[$i]->[0];
      my $pos   = $$genes[$i]->[2];
      my $gname = $$genes[$i]->[3];
      
      my $is_landmark = $$genes[$i]->[4];
      if ( $prev_chrom ne $chrom ) {
        $prev_pos   = $pos;
        $prev_chrom = $chrom;
        next;    # always skip the first gene in the chromosome
      }

      #printf "%-10d\t%s\t%2s\t%.6f\t\t%s\t%.6f\t%s\t%.6f\n", $i, $gname, $chrom, $pos, $prev_gname, $prev_pos, $$genes[$i+1]->[3], $$genes[$i+1]->[2];
      
      # get the next position
      # are we at the end of the array or end of the chromosome?
      if ( $i + 1 < @$genes && $$genes[ $i + 1 ]->[0] eq $chrom ) {
        $next_pos = $$genes[ $i + 1 ]->[2];
      }
      else {
        $next_pos = $pos + 0.5;
      }
      
      # should this position be changed? Test for this position less than previous.
      if ( $prev_pos > $pos ) {
        if (not $is_landmark) {
          # get the difference between the previous and next positions
          my $diff = $next_pos - $prev_pos;
          if ( $diff > 0.0005 ) {
            $$genes[$i]->[2] = $prev_pos + ( $diff / 2 );
          }
          else {
            $$genes[$i]->[2] = $prev_pos + 0.0005;
          }
          $changed_in_this_iteration = 1;
        }
        $pos = $$genes[$i]->[2];
      }
      
      # should this position be changed? Test for this position greater than next
      if ( $pos > $next_pos && $next_pos > $prev_pos ) {
        
        # get the difference between the previous and next positions
        if (not $is_landmark) {
          my $diff = $next_pos - $prev_pos;
          if ( $diff > 0.0005 ) {
            $$genes[$i]->[2] = $prev_pos + ( $diff / 2 );
          }
          else {
            $$genes[$i]->[2] = $prev_pos + 0.0005;
          }
          $changed_in_this_iteration = 1;
        }
        $pos = $$genes[$i]->[2];
      }
      $prev_pos = $pos;
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

  my @genes;
  my $db = Ace->connect( -path => $database );
  push @genes, $db->find('find Gene * where Map & SMap');
  foreach my $gene (@genes) {
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
      s/\"//g;
      my @a = split;
      next if !( $a[1] eq 'gene' && $a[2] eq 'gene' );
      my $gene_id = $a[9];
      next if not exists $gen_map{$gene_id};
      
      my $chrom   = $a[0];
      my $map_pos = ( $a[3] + $a[4] ) / 2;
      my $gen_pos = $gen_map{$gene_id}->[0];
      my $is_landmark = $gen_map{$gene_id}->[1];
      $chrom =~ s/CHROMOSOME_//;
      
      #should be chromosome->map_pos->gen_pos
      $genes{$chrom}->{$map_pos} = [$gen_pos,$gene_id, $is_landmark];
    }
    close IN;
  }
  
  return \%genes, \%gen_map;
}

1;
