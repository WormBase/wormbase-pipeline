#!/usr/local/bin/perl5.8.0 -w
#
# transcript_builder.pl
# 
# by Anthony Rogers
#
# Script to make ?Transcript objects
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2005-03-15 12:10:27 $

use strict;
use lib $ENV{'CVS_DIR'};
use Getopt::Long;
use Data::Dumper;
use Coords_converter;
use Wormbase;
use Modules::SequenceObj;
use Modules::Transcript;
use Modules::CDS;
use Modules::Strand_transformer;

my $tace = &tace;

my ($debug, $help, $verbose, $really_verbose, $est, $gff, $database, $build, $new_coords, $test, $UTR_range, @chromosomes, $gff_dir, $test_cds);

my $gap = 15;			# $gap is the gap allowed in an EST alignment before it is considered a "real" intron

# to track failings of system calls
my $errors = 0;

GetOptions ( "debug:s"          => \$debug,
	     "help"             => \$help,
	     "verbose"          => \$verbose,
	     "really_verbose"   => \$really_verbose,
	     "est:s"            => \$est,
	     "gap:s"            => \$gap,
	     "gff:s"            => \$gff,
	     "database:s"       => \$database,
	     "build"            => \$build,
	     "new_coords"       => \$new_coords,
	     "test"             => \$test,
	     "utr_size:s"       => \$UTR_range,
	     "chromosome:s"     => \@chromosomes,
	     "gff_dir:s"        => \$gff_dir,
	     "cds:s"            => \$test_cds
	   ) ;

@chromosomes = split(/,/,join(',',@chromosomes));

# Who will log be emailed to?
my $maintainers = "All";
if ($debug) {
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
  my $set_debug = SequenceObj->new();
  $set_debug->debug($debug);
}

*STDERR = *STDOUT;

&check_opts;
die "no database\n" unless $database;

my $log = Log_files->make_build_log($debug);

#setup directory for transcript
my $transcript_dir = "$database/TRANSCRIPTS";
&run_command("mkdir $transcript_dir") unless -e "$transcript_dir";

die "cant create log file\n\n" unless $log;

my $coords;
# write out the transcript objects
# get coords obj to return clone and coords from chromosomal coords
$coords = Coords_converter->invoke($database);


#Load in Feature_data : cDNA associations from COMMON_DATA
my %feature_data;
&load_features( \%feature_data );


# process chromosome at a time
@chromosomes = qw(I II III IV V X) unless @chromosomes;

foreach my $chrom ( @chromosomes ) {

  # links store start /end chrom coords
  my $link_start;
  my $link_end;

  my %genes_exons;
  my %genes_span;
  my %cDNA;
  my %cDNA_span;
  my %cDNA_index;
  my %features;

  my $transformer;

  my @cdna_objs;
  my @cds_objs;

  my $index = 0;

  my $gff_file;

  $gff_file = $gff ? $gff : "$gff_dir/CHROMOSOME_$chrom.gff";

  open (GFF,"<$gff_file") or $log->log_and_die("gff_file\n");
  $log->write_to("reading gff file $gff_file\n");

  # parse GFF file to get CDS, exon and cDNA info
  while (<GFF>) {
    my @data = split;
    
    next if( $data[1] eq "history" );
    #  GENE STRUCTURE
    if ( $data[1] eq "curated" ) {
      $data[9] =~ s/\"//g;
      next if( defined $test_cds and ($data[9] ne $test_cds)) ;
      if ( $data[2] eq "CDS" ) {
	# GENE SPAN
	$genes_span{$data[9]} = [($data[3], $data[4], $data[6])];

      } elsif ($data[2] eq "exon" ) {
	# EXON 
	$genes_exons{$data[9]}{$data[3]} = $data[4];
      }
    }
    
    # BLAT DATA
    elsif ( ($data[1] eq "BLAT_EST_BEST") or ($data[1] eq "BLAT_OST_BEST") or ($data[1] eq "BLAT_mRNA_BEST") ) {
      $data[9] =~ s/\"//g;
      $data[9] =~ s/Sequence:// ;
      $cDNA{$data[9]}{$data[3]} = $data[4];
      # keep min max span of cDNA
      if ( !(defined($cDNA_span{$data[9]}[0])) or ($cDNA_span{$data[9]}[0] > $data[3]) ) {
	$cDNA_span{$data[9]}[0] = $data[3];
	$cDNA_span{$data[9]}[2] = $data[6]; #store strand of cDNA
      } 
      if ( !(defined($cDNA_span{$data[9]}[1])) or ($cDNA_span{$data[9]}[1] < $data[4]) ) {
	$cDNA_span{$data[9]}[1] = $data[4];
      }
    }

    #create Strand_transformer for '-' strand coord reversal
    elsif ( $data[1] eq "Link" and $data[9] =~ "CHROMOSOME") {
      $transformer = Strand_transformer->new($data[3],$data[4]);
    }
    
    # add feature_data to cDNA
    #CHROMOSOME_I    SL1     SL1_acceptor_site       182772  182773  .       -       .       Feature "WBsf016344"
    elsif ( $data[9] and $data[9] =~ /(WBsf\d+)/) { #Feature "WBsf003597"
      my $feat_id = $1;
      my $dnas = $feature_data{$feat_id};
      if( $dnas ) {
	foreach my $dna ( @{$dnas} ) {
	  #	print "$dna\t$data[9]  --- $data[6] ---  ",$cDNA_span{"$dna"}[2],"\n";
	  next unless ( $cDNA_span{"$dna"}[2] and $data[6] eq $cDNA_span{"$dna"}[2] ); # ensure same strand
	  $cDNA_span{"$dna"}[3]{"$data[1]"} = [ $data[3], $data[4], $1 ]; # 182772  182773 WBsf01634
	}
      }
    }
  }

  # need to sort the cds's into ordered arrays + and - strand genes are in distinct coord space so they need to be kept apart
  my %fwd_cds;
  my %rev_cds;
  foreach ( keys %genes_span ) {
    if( $genes_span{$_}->[2] eq "+" ) {
      $fwd_cds{$_} = $genes_span{$_}->[0];
    }
    else {
      $rev_cds{$_} = $genes_span{$_}->[0];
    }
  }

  close GFF;
  
  &load_EST_data(\%cDNA_span, $chrom);
  # &checkData(\$gff,\$%cDNA_span, \%genes_span); # this just checks that there is some BLAT and gene data in the GFF file
  &eradicateSingleBaseDiff(\%cDNA);

  #create transcript obj for each CDS
  # fwd strand cds will be in block first then rev strand
  foreach (sort { $fwd_cds{$a} <=> $fwd_cds{$b} } keys  %fwd_cds ) {
    #next if $genes_span{$_}->[2] eq "-"; #only do fwd strand for now
    my $cds = CDS->new( $_, $genes_exons{$_}, $genes_span{$_}->[2], $chrom, $transformer );
    push( @cds_objs, $cds);
    $cds->array_index("$index");
    $index++;
  }
  foreach ( sort { $rev_cds{$b} <=> $rev_cds{$a} } keys  %rev_cds ) {
    #next if $genes_span{$_}->[2] eq "-"; #only do fwd strand for now
    my $cds = CDS->new( $_, $genes_exons{$_}, $genes_span{$_}->[2], $chrom, $transformer );
    push( @cds_objs, $cds);
    $cds->array_index("$index");
    $index++;
  }



  $index = 0;
  foreach ( keys %cDNA ) {

    my $cdna = SequenceObj->new( $_, $cDNA{$_}, $cDNA_span{$_}->[2] );
    $cdna->array_index($index);
    if ( $cDNA_span{"$_"}[3] ) {
      foreach my $feat ( keys %{$cDNA_span{"$_"}[3]} ) {
	$cdna->$feat( $cDNA_span{"$_"}[3]->{"$feat"} );
      }
    }
    # add paired read info
    if( $cDNA_span{"$_"}[4] ) {
      $cdna->paired_read( $cDNA_span{"$_"}[4] );
    }

    # index info
    $cDNA_index{"$_"} = $index;
    $index++;
    $cdna->transform_strand($transformer,"transform") if ( $cdna->strand eq "-" );

    #check for and remove ESTs with internal SL's 
    if( $cdna->SL ) {
      if( $cdna->start < $cdna->SL->[0] ) {
	my $gap = $cdna->SL->[0] - $cdna->start;
	$log->write_to($cdna->name." has internal SL ".$cdna->SL->[2]."gap = $gap\n");
	next;
      }
    }

    push(@cdna_objs,$cdna);
  }


  # these are no longer needed so free memory !
  %genes_exons = ();
  %genes_span= ();
  %cDNA = ();
  %cDNA_span = ();

  ################################################
  #            DATA LOADED           #
  ################################################


 CDNA:
  foreach my $CDNA ( @cdna_objs) {
    next if ( defined $est and $CDNA->name ne "$est"); #debug line

    #sanity check features on cDNA ie SLs are at start
    next if ( &sanity_check_features( $CDNA ) == 0 );

    foreach my $cds ( @cds_objs ) {
      if ( $cds->map_cDNA($CDNA) ) {
	$CDNA->mapped($cds );
	next CDNA;
      }
    }
    last if $est;
  }

  # tie up read pairs
  foreach my $cdna ( @cds_objs ) {
    my $name = $cdna->name;
    if($cDNA_index{$name}) {}
  }
      
  print "Second round additions\n";

  foreach my $CDNA ( @cdna_objs) {
    next if ( defined $est and $CDNA->name ne "$est"); #debug line
    next if ( defined($CDNA->mapped) );
    foreach my $cds ( @cds_objs ) {
      if ( $cds->map_cDNA($CDNA) ) {
	$CDNA->mapped($cds);
	print $CDNA->name," overlaps ",$cds->name,"on 2nd round \n" if $verbose;
      }
    }
  }
  # try and attach paired reads that dont overlap
 PAIR: foreach my $CDNA ( @cdna_objs) {
    next if ( defined $est and $CDNA->name ne "$est"); #debug line
    next if $CDNA->mapped;

    # get name of paired read
    next unless (my $mapped_pair_name = $CDNA->paired_read );

    # retreive object from array
    next unless ($cDNA_index{"$mapped_pair_name"} and (my $pair = $cdna_objs[ $cDNA_index{"$mapped_pair_name"} ] ) );

    # get cds that paired read maps to 
    if (my $cds = $pair->mapped ) {
      my $index = $cds->array_index;

      # find next downstream CDS - must be on same strand
      my $downstream_CDS;
    DOWN: while (! defined $downstream_CDS ) {
	$index++;
	if ( $downstream_CDS = $cds_objs[ $index ] ) {
	
	  unless ( $downstream_CDS ) {
	    last;
	    print "last gene in array\n";
	  }
	  # dont count isoforms
	  my $down_name = $downstream_CDS->name;
	  my $name = $cds->name;
	  $name =~ s/[a-z]//;
	  $down_name =~ s/[a-z]//;
	  if( $name eq $down_name ) {
	    undef $downstream_CDS;
	    next;
	  }
	  # @cds_objs is structured so that + strand genes are in a block at start, then - strand
	  last DOWN if( $downstream_CDS->strand ne $cds->strand );

	  # check unmapped cdna ( $CDNA ) lies within 1kb of CDS that paired read maps to ( $cds ) and before $downstream_CDS
	
	  print "trying ",$cds->name, " downstream is ", $downstream_CDS->name," with ",$CDNA->name,"\n";
	  if( ($CDNA->start > $cds->gene_end) and ($CDNA->start - $cds->gene_end < 1000) and ($CDNA->end < $downstream_CDS->gene_start) ) {
	    print " adding 3' cDNA ",$CDNA->name," to ",$cds->name,"\n";
	    $cds->add_3_UTR($CDNA);
	    $CDNA->mapped($cds);
	    last;
	  }
	}
	else {
	  last DOWN;
	}
      }
    }
  }

  print "Third round additions\n";

  foreach my $CDNA ( @cdna_objs) {
    next if ( defined $est and $CDNA->name ne "$est"); #debug line
    next if ( defined($CDNA->mapped) );
    foreach my $cds ( @cds_objs ) {
      if ( $cds->map_cDNA($CDNA) ) {
	$CDNA->mapped($cds);
	print $CDNA->name," overlaps ",$cds->name,"on 2nd round \n" if $verbose;
      }
    }
  }


  my $out_file = "$transcript_dir/chromosome_$chrom.ace";
  print "writing output to $out_file\n";
  open (FH,">$out_file") or $log->log_and_die("cant open $out_file\n");
  foreach my $cds (@cds_objs ) {
    print "reporting : ",$cds->name,"\n" if $debug;
    $cds->report(*FH, $coords, $transformer);
  }

  last if $gff;			# if only doing a specified gff file exit after this is complete
}

$log->mail("$maintainers");
exit(0);


######################################################################################################
#
#
#                           T  H  E        S  U  B  R  O  U  T  I  N  E  S
#
#
#######################################################################################################


sub eradicateSingleBaseDiff
  {
    my $cDNA = shift;
    $log->write_to( "removing small cDNA mismatches\n\n\n");
    foreach my $cdna_hash (keys %{$cDNA} ) {
      my $last_key;
      my $check;
      foreach my $exons (sort { $$cDNA{$cdna_hash}->{$a} <=> $$cDNA{$cdna_hash}->{$b} } keys %{$$cDNA{$cdna_hash}}) {
	print "\n############### $cdna_hash #############\n" if $really_verbose;
	print "$exons -> $$cDNA{$cdna_hash}->{$exons}\n" if $really_verbose;
	my $new_last_key = $exons;
	if ( $last_key ) {
	  if ( $$cDNA{$cdna_hash}->{$last_key} >= $exons - $gap ) { #allows seq error gaps up to $gap bp
	    $$cDNA{$cdna_hash}->{$last_key} = $$cDNA{$cdna_hash}->{$exons};
	    delete $$cDNA{$cdna_hash}->{$exons};
	    $check = 1;
	    $new_last_key = $last_key;
	  }
	}
	$last_key = $new_last_key;
      }
      if ( $check ) {
	foreach my $exons (sort keys  %{$$cDNA{$cdna_hash}}) {
	  print "single base diffs removed from $cdna_hash\n" if $really_verbose;
	  print "$exons -> $$cDNA{$cdna_hash}->{$exons}\n" if $really_verbose;
	}
      }
    }
  }

sub check_opts {
  # sanity check options
  if ( $help ) {
    system("perldoc $0");
    exit(0);
  }

  $database = glob("~wormpub/TEST_BUILD/autoace") if $test;

}

sub checkData
  {
    my $file = shift;
    my $cDNA_span = shift;
    my $genes_span = shift;
    die "There's no BLAT data in the gff file $$file\n" if scalar keys %{$cDNA_span} == 0;
    die "There are no genes in the gff file $$file\n" if scalar keys %{$genes_span} == 0;
  }



###################################################################################

sub run_command
  {
    my $command = shift;
    $log->write_to( &runtime, ": started running $command\n");
    my $status = system($command);
    if ($status != 0) {
      $errors++;
      $log->write_to("ERROR: $command failed\n");
    }
    $log->write_to( &runtime, ": finished running $command\n");

    # for optional further testing by calling subroutine
    return($status);
  }

#####################################################################################


sub load_EST_data
  {
    my $cDNA_span = shift;
    my $chrom = shift;
    my %est_orient;
    &FetchData("estorientation",\%est_orient,glob("~wormpub/DATABASES/autoace/COMMON_DATA"));
    foreach my $EST ( keys %est_orient ) {
      if ( $$cDNA_span{$EST} ) {
	my $GFF_strand = $$cDNA_span{$EST}->[2];
	my $read_dir = $est_orient{$EST};
      CASE:{
	  ($GFF_strand eq "+" and $read_dir eq "5") && do {
	    $$cDNA_span{$EST}->[2] = "+";
	    last CASE;
	  };
	  ($GFF_strand eq "+" and $read_dir eq "3") && do {
	    $$cDNA_span{$EST}->[2] = "-";
	    last CASE;
	  };
	  ($GFF_strand eq "-" and $read_dir eq "5") && do {
	    $$cDNA_span{$EST}->[2] = "-";
	    last CASE;
	  };
	  ($GFF_strand eq "-" and $read_dir eq "3") && do {
	    $$cDNA_span{$EST}->[2] = "+";
	    last CASE;
	  };
	}
      }
    }

    # load paired read info
    print STDERR "Loading EST paired read info\n";
    my $pairs = "$database/EST_pairs.txt";
    
    if ( -e $pairs ) {
      open ( PAIRS, "<$pairs") or $log->log_and_die("cant open $pairs :\t$!\n");
      while ( <PAIRS> ) {
	chomp;
	s/\"//g;
	s/Sequence://g;
	next if( ( $_ =~ /acedb/) or ($_ =~ /\/\//) );
	my @data = split;
	$$cDNA_span{$data[0]}->[4] = $data[1];
      }
      close PAIRS;
    }
    else {
      my $cmd = "select cdna, pair from cdna in class cDNA_sequence where exists_tag cdna->paired_read, pair in cdna->paired_read";
      
      open (TACE, "echo '$cmd' | $tace $database |") or die "cant open tace to $database using $tace\n";
      open ( PAIRS, ">$pairs") or die "cant open $pairs :\t$!\n";
      while ( <TACE> ) {
	chomp;
	s/\"//g;
	my @data = split;
	$$cDNA_span{$data[0]}->[4] = $data[1];
	print PAIRS "$data[0]\t$data[1]\n";
      }
      close PAIRS;
    }
  }

sub load_features 
  {
    my $features = shift;
    my %tmp;
    &FetchData("Featurelist",\%tmp,glob("~wormpub/DATABASES/autoace/COMMON_DATA") );
    foreach my $seq ( keys %tmp ) {
      my @feature = split(/,/,$tmp{$seq});
      foreach ( @feature ) {
	push(@{$$features{$_}},$seq);
      }
    }
  }

sub sanity_check_features
  {
    my $cdna = shift;
    my $return = 1;
    if( my $SL = $cdna->SL ) {
      $return = 0 if( $SL->[0] != $cdna->start );
      print STDERR $SL->[2]," inside ",$cdna->name,"\n";
    }
    if( my $polyA = $cdna->polyA_site ){
      $return = 0 if( $polyA->[1] != $cdna->end );
      print STDERR $polyA->[2]," inside ",$cdna->name,"\n";
    }

    return $return;
  }

__END__

=pod

=head2 NAME - transcript_builder.pl

=head1 USAGE

=over 4

=item transcript_builder.pl  [-options]

=back

This script "does exactly what it says on the tin". ie it builds transcript objects for each gene in the database

To do this it ;

1) Determines matching_cDNA status for each gene. Goes through gff files and examines each cDNA to see if it matches any gene that it overlaps.

2) For each gene that has matching cDNAs it then confirms that every exon of the gene that lies within the region covered by the cDNA is covered correctly.  So if a gene has an extra exon that is in the intron of a cDNA, that cDNA will NOT be linked to that gene.


=back

=head2 transcript_builder.pl arguments:

=over 4

=item * verbose and really_verbose  -  levels of terminal output
  
=item *  est:s     - just do for single est 

=item * gap:s      - when building up cDNA exon structures from gff file there are often single / multiple base pair gaps in the alignment. This sets the gap size that is allowable [ defaults to 5 ]

=item * gff:s         - pass in a gff file to use

=item * database:s      - database directory to use. Expects gff files to be in CHROMSOMES subdirectory.

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
