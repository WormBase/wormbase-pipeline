#/software/bin/perl -w
#
# blat2ace.pl
# 
# by Kerstin Jekosch
#
# Exporter to map blat data to genome and to find the best match for each EST, mRNA, OST, etc.
#
# Last edited by: $Author: klh $
# Last edited on: $Date: 2011-03-21 11:01:49 $

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use File::Copy;

#########################
# Command line options  #
#########################

my ($help, $debug, $test, $verbose, $store, $wormbase, $species);
my ($database, $confirmed_introns, $virtualobjs, $qtype, $qspecies);
my ($acefile, $pslfile, $virtualobjsfile, $confirmedfile);

GetOptions (
	    "help"       => \$help,
	    "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "database:s" => \$database,
	    "intron"     => \$confirmed_introns,
            "cifile:s"   => \$confirmedfile,
	    "virtual"    => \$virtualobjs,
            "vfile:s"    => \$virtualobjsfile,
	    "type:s"     => \$qtype,
	    "qspecies:s" => \$qspecies, #query species
	    "ace:s"      => \$acefile,
	    "psl:s"      => \$pslfile
	   );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species
			     );
}
# establish log file.
my $log = Log_files->make_build_log($wormbase);

$log->log_and_die("no type specified\n") unless $qtype;

#############################
# variables and directories #
#############################

# set database paths, default to autoace unless -camace
my $ace_dir   = $wormbase->orgdb;
my $blat_dir = (defined $database) ? $database."/BLAT" : $wormbase->blat;
my @nematodes = qw(nematode washu nembase);

#############################
# CommonData hash retrieval #
#############################
my ($accessor, %estorientation);

if ($qspecies ne $wormbase->species) {
  my (%sa) = $wormbase->species_accessors;
  if (exists $sa{$qspecies}) {
    $accessor = $sa{$qspecies};
  }
} else {
  $accessor = $wormbase;
}

if (defined $accessor) {
  eval {
    %estorientation   = $accessor->FetchData('estorientation');
  };
  if ($@) {
    $log->write_to("Could not successfully retrieve orientation info for $qspecies\n");
    %estorientation = ();
  }
}


##########################################################################################
# map the blat hits to ace - i.e. process blat output (*.psl) file into set of ace files #
##########################################################################################
$log->write_to($wormbase->runtime.": Start mapping $qspecies $qtype\n\n");

# open input and output filehandles
$acefile = sprintf("%s/%s.blat.%s_%s.ace", $blat_dir, $wormbase->species, $qspecies, $qtype)
    unless $acefile;
$pslfile = sprintf("%s/%s_%s_out.psl", $blat_dir, $qspecies, $qtype) 
    unless $pslfile;
$virtualobjsfile = sprintf("%s/virtual_objects.%s.blat.%s.%s.ace", $blat_dir, $wormbase->species, $qtype, $qspecies) 
    unless defined $virtualobjsfile;


open(BLAT, "<$pslfile") or $log->log_and_die("Could open $pslfile $!\n");

my (%hits, %hits_by_tname, %target_lengths);

while (<BLAT>) {
  
  next unless (/^\d/);
  my @f            = split "\t";
  
  my $match        = $f[0];                    # number of bases matched by blat
  my $strand       = $f[8];                    # strand that match is on
  my $qname        = $f[9];                    # query sequence name
  my $qsize        = $f[10];                   # query sequence length
  my $qstart       = $f[11];
  my $qend         = $f[12];
  my $tname        = $f[13];                   # name of superlink that was used as blat target sequence
  my $tsize        = $f[14];                   # target seq size (used to be superlink hence sl)
  my $tstart       = $f[15];                   # target (superlink) start coordinate...
  my $tend         = $f[16];                   # ...and end coordinate
  my $block_count  = $f[17];                   # block count
  my @lengths      = split (/,/, $f[18]);      # sizes of each blat 'block' in any individual blat match
  my @q_starts     = split (/,/, $f[19]);      # start coordinates of each query block
  my @t_starts     = split (/,/, $f[20]);      # start coordinates of each target (superlink) block
  
  if (not exists $target_lengths{$tname}) {
    $target_lengths{$tname} = $tsize;
  }
  
  # -- hits are usually shadows of ++ hits (apparently)
  # next if ($strand eq '--');
  
  # calculate (acedb) score for each blat match
  # new way of calculating score, divide by query size rather than sum of matching blocks, 
  my $score = ($match/$qsize)*100;
  
  if ($strand =~ /^\S$/) {
    # single strand reported; BLAT always reports forward on the target in this case, 
    # flipping the query if necessary. 
    $strand .= "+";
  }
  
  my ($qstrand, $tstrand) = $strand =~ /^(\S)(\S)$/;
  
  my @segs = ();  
  
  for (my $x = 0;$x < $block_count; $x++) {
    
    my ($query_start,$query_end, $target_start, $target_end);
    
    $query_start = $q_starts[$x] + 1;
    $query_end   = $query_start + $lengths[$x] - 1;
    
    $target_start = $t_starts[$x] + 1;
    $target_end   = $target_start + $lengths[$x] - 1;
    
    if ($qstrand eq '-') {
      $query_end = $qsize - $q_starts[$x];
      $query_start   = $query_end - $lengths[$x] + 1;
    }
    
    if ($tstrand eq '-') {
      $target_end = $tsize - $t_starts[$x];
      $target_start = $target_end - $lengths[$x] + 1;
    }
    
    # 1 out of 8000,000 cDNAs ends up with 0 as a start and breaks the mapping, this hack is to stop it
    $target_start++ if $target_start == 0;
    $query_start++ if  $query_start == 0;
    
    push @segs, {
      tstart => $target_start,
      tend   => $target_end,
      qstart => $query_start,
      qend   => $query_end,
    }			
  }
  
  # we want to flip the strands such that the query is always forward (easier that way)
  if ($qstrand eq '-') {
    $qstrand = '+';
    $tstrand = ($tstrand eq '+') ? "-" : "+";
  }
  
  my $hit =  {
    tname    => $tname,
    tstart   => $tstart,
    tend     => $tend,
    tstrand  => $tstrand,
    qname    => $qname,
    qstart   => $qstart,
    qend     => $qend,
    qstrand  => $qstrand,
    score    => $score,
    #match    => $match,
    segments => \@segs,
    isbest   => 0,
  };
  
  push @{$hits{$qname}}, $hit;
 
}
close(BLAT);

foreach my $qname (keys %hits) {
  # BLAT is prone to produce "shadow" alignments (trivial variants of another
  # higher scoring alignment, usually with the query and target strand reversed)
  # Check for these here and filter them out

  my @nr_hits;
  foreach my $hit (sort { $b->{score} <=> $a->{score} } @{$hits{$qname}}) {
    my $overlaps = 0;

    KEPT: foreach my $kept (@nr_hits) {
      if ($hit->{tname} eq $kept->{tname} and 
          $hit->{tstrand} eq $kept->{tstrand} and 
          $hit->{tstart} <= $kept->{tend} and
          $hit->{tend} >= $kept->{tstart}) {
        # check for exon overlap
        foreach my $seg (@{$hit->{segments}}) {
          foreach my $kseg (@{$kept->{segments}}) {
            if ($kseg->{tstart} <= $seg->{tend} and
                $kseg->{tend} >= $seg->{tstart} and
                $kseg->{qstart} <= $seg->{qend} and
                $kseg->{qend} >= $seg->{qstart}) {
              $overlaps = 1;
              last KEPT;
            }
          }
        }

      }
    }

    if (not $overlaps) {
      push @nr_hits, $hit;
    }
  }

  # finally, mark the best-scoring hit as best iff it is the single best-scoring alignment

  if (scalar(@nr_hits) == 1 or
      $nr_hits[0]->{score} > $nr_hits[1]->{score}) {
    $nr_hits[0]->{isbest} = 1;
  }

  foreach my $hit (@nr_hits) {
    push @{$hits_by_tname{ $hit->{tname} }}, $hit;
  }
  delete $hits{$qname};
}


if ($confirmed_introns) {
  my $out_file = (defined $confirmedfile)
      ? $confirmedfile
      : sprintf("%s/%s.ci.%s_%s.ace", $blat_dir, $wormbase->species, $qspecies, $qtype);
  &confirm_introns($out_file, $qtype, \%hits_by_tname);
}


my ($bestm, $otherm);
 
if (grep(/$qspecies/, @nematodes)) {
  $bestm = $otherm = "BLAT_".uc($qspecies);

} else {
  my $mpre = $qspecies eq $wormbase->species 
    ? "BLAT_${qtype}" 
    : "BLAT_Caen_${qtype}";
  $bestm = $mpre . "_BEST";
  $otherm = $mpre . "_OTHER";
}

# need to write separate out hits into BEST and OTHER
my $virt_hash = &write_ace($acefile, $qtype, $bestm, $otherm, \%hits_by_tname);

if ($virtualobjs) {
  # write virtuals sequences here
  open my $vfh, ">$virtualobjsfile" or $log->log_and_die("Could not open $virtualobjsfile for writing\n");

  foreach my $tname (keys %$virt_hash) {
    print $vfh "\nSequence : \"$tname\"\n";
    foreach my $child (sort { my ($na) = ($a =~ /_(\d+)$/); my ($nb) = ($b =~ /_(\d+)$/); $na <=> $nb } keys %{$virt_hash->{$tname}}) {
      printf $vfh "S_Child Homol_data %s %d %d\n", $child, @{$virt_hash->{$tname}->{$child}};
    }
  }
  close($vfh);
}



$log->mail;
print "Finished.\n" if ($verbose);
exit(0);




#################################################################################
#                                                                               #
#                          Subroutines                                          #
#                                                                               #
#################################################################################


sub write_ace {
  my ($out_file, $type, $best_m, $other_m, $hits) = @_;

  # strategy; divide the parent sequence into 150000-sized bins, and place each alignment
  # into a bin. Also create an SChild for the whole parent sequence containing alignments that
  # do not fit completely inside a bin
  my $binsize = 100000;
  
  my %virtuals;

  foreach my $tname (keys %$hits) {
    foreach my $hit (@{$hits->{$tname}}) {
      my $bin = 1 +  int( $hit->{tstart} / $binsize );
      my $bin_start = ($bin - 1) * $binsize + 1;
      my $bin_end   = $bin_start + $binsize - 1;
      
      if ($bin_end > $target_lengths{$tname}) {
        $bin_end = $target_lengths{$tname};
      }
      my $bin_of_end = 1 +  int( $hit->{tend} / $binsize );

      if ($bin == $bin_of_end or
          ($bin == $bin_of_end - 1 and $hit->{tend} - $bin_end < ($binsize / 2))) {
        # both start and end lie in the same bin or adjacent bins - okay
        foreach my $seg (@{$hit->{segments}}) {
          $seg->{tstart} = $seg->{tstart} - $bin_start + 1;
          $seg->{tend}   = $seg->{tend} - $bin_start + 1;
        }
      } else {
        # hit has too great a span; skip it
        next;
      }

      # the code below places hits that span multiple bins onto
      # a special Homol_data that spans the whole parent sequence. 
      # However, this is not what happens in the old code, so I 
      # will not activiate it (for now) - klh 110317

      #if ($hit->{tend} > $bin_end) {
      #  $bin = 0;
      #  $bin_start = 1;
      #  $bin_end   = $target_lengths{$tname};
      #} else {
      #  foreach my $seg (@{$hit->{segments}}) {
      #    $seg->{tstart} = $seg->{tstart} - $bin_start + 1;
      #    $seg->{tend}   = $seg->{tend} - $bin_start + 1;
      #  }
      #}
    
      $hit->{bin} = $bin;
      my $parent = "BLAT_$type:$tname";
      if ($bin) {
        $parent .= "_$bin";
      }
      
      if (not exists $virtuals{$tname}->{$parent}) {
        $virtuals{$tname}->{$parent} = [$bin_start, $bin_end];
      }
    }
  }

 
  # finally, write out the results
  open(my $out_fh, ">$out_file") or $log->log_and_die("Could not open $out_file for writing\n");
  foreach my $tname (keys %$hits) {
    my @list =  grep { exists $_->{bin} } @{$hits->{$tname}};
    @list = sort { $a->{bin} <=> $b->{bin} or $a->{tstart} <=> $b->{tstart} } @list;

    my $prev_bin;
    foreach my $hit (@list) {

      if (not defined $prev_bin or $hit->{bin} != $prev_bin) {
        printf($out_fh 
               "\nHomol_data : \"BLAT_%s:%s%s\"\n", 
               $type,
               $tname,
               $hit->{bin} ? "_" . $hit->{bin} : "");
        $prev_bin = $hit->{bin};
      }
      
      
      foreach my $seg (@{$hit->{segments}}) {
        printf($out_fh "DNA_homol\t\"%s\"\t%s\t%.1f\t%d\t%d\t%d\t%d\n",
               $hit->{qname},
               ($hit->{isbest}) ? $best_m : $other_m,
               $hit->{score},
               $hit->{tstrand} eq "+" ? $seg->{tstart} : $seg->{tend},
               $hit->{tstrand} eq "+" ? $seg->{tend}   : $seg->{tstart},
               $seg->{qstart},
               $seg->{qend});
        
      }
    }
  }

  close($out_fh);

  return \%virtuals;
}


#########################
# confirm introns
########################
sub confirm_introns {
  my ($outfile, $type, $hits) = @_;

  my (%ci, $no_direction);

  foreach my $tname (keys %$hits) {
    
    my @alns = grep { $_->{isbest} } @{$hits->{$tname}};

    foreach my $aln (@alns) {
      my $qname = $aln->{qname};

      if ($type eq 'mRNA' and not exists $estorientation{$qname}) {
        $estorientation{$aln->{qname}} = 5;
      }

      my @segs = sort { $a->{tstart} <=> $b->{tstart} } @{$aln->{segments}};

      for(my $y=1; $y < @segs; $y++) {
        my $first = $segs[$y-1]->{tend} + 1;
        my $second = $segs[$y]->{tstart} - 1;

        if ($aln->{qstrand} eq $aln->{tstrand}) {
          if ($segs[$y-1]->{qend} + 1 == $segs[$y]->{qstart} and ($second - $first) > 2) {
            if (exists $estorientation{$qname} && $estorientation{$qname} eq '3') {
              push @{$ci{$tname}}, [$second,$first,$qname];
            } elsif (exists $estorientation{$qname} && $estorientation{$qname} eq '5') {
              push @{$ci{$tname}}, [$first,$second,$qname];
            } else {
              $no_direction++;
            }
          }
        } else {
          if ($segs[$y-1]->{qstart} - 1 == $segs[$y]->{qend} and (($second - $first) > 2)) {
            if (exists $estorientation{$qname} && $estorientation{$qname} eq '3') {
              push @{$ci{$tname}}, [$first,$second,$qname];
            } elsif (exists $estorientation{$qname} && $estorientation{$qname} eq '5') {
              push @{$ci{$tname}}, [$second,$first,$qname]; 
            } else {
              $no_direction++;
            }
          }
        }
      }
    }
  }
    
  $log->write_to("WARNING: Direction not found for $no_direction transcripts\n\n") if ($no_direction);

  ###################################
  # produce confirmed intron output #
  # Write them out as-is (i.e. not on virtuals)
  # because they are not loaded directly but processed 
  # another script
  ###################################

  open(my $ci_fh, ">$outfile");

  foreach my $tname (sort keys %ci) {
    my %double;
      
    print $ci_fh "\nSequence : \"$tname\"\n";
      
    foreach my $intr (@{$ci{$tname}}) {
      my $merge = $intr->[0].":".$intr->[1];
      if (not exists $double{$merge}) {
        # If RST or OST modify output.
        if (($type eq "RST") or ($type eq "OST")) {
          printf $ci_fh "Confirmed_intron %d %d EST %s\n", @$intr;
        } elsif ($type eq "mRNA") {
          printf $ci_fh "Confirmed_intron %d %d cDNA %s\n", @$intr;
        } else {
          printf $ci_fh "Confirmed_intron %d %d $type %s\n", @$intr;
        }
      }
      $double{$merge} = 1;
    }
  }
  close($ci_fh);
}

__END__

