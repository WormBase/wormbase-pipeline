#!/usr/local/bin/perl5.8.0 -w
#
# transcript_builder.pl
# 
# by Anthony Rogers
#
# Script to make ?Transcript objects
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2004-05-04 09:09:02 $

use strict;
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
use Getopt::Long;
use Data::Dumper;
use Coords_converter;
use Wormbase;
use Modules::Transcript;
use Modules::CDS;

my $tace = &tace;

my ($debug, $help, $verbose, $really_verbose, $est, $count, $report, $transcript, $gff, $show_matches, $database, $overlap_check, $load_matches, $load_transcripts, $build, $new_coords, $test, $UTR_range);

my $gap = 5; # $gap is the gap allowed in an EST alignment before it is considered a "real" intron

# to track failings of system calls
my $errors = 0;

GetOptions ( "debug:s"          => \$debug,
	     "help"             => \$help,
	     "verbose"          => \$verbose,
	     "really_verbose"   => \$really_verbose,
	     "est:s"            => \$est,
	     "count"            => \$count,
	     "report"           => \$report,
	     "gap:s"            => \$gap,
	     "transcript"       => \$transcript,
	     "gff:s"            => \$gff,
	     "show_matches"     => \$show_matches,
	     "database:s"       => \$database,
	     "overlap"          => \$overlap_check,
	     "load_transcripts" => \$load_transcripts,
	     "load_matches"     => \$load_matches,
	     "build"            => \$build,
	     "new_coords"       => \$new_coords,
	     "test"             => \$test,
	     "utr_size:s"       => \$UTR_range
	   ) ;

# Who will log be emailed to?
my $maintainers = "All";
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}


my $log = Log_files->make_build_log($debug);

&check_opts; # if -build set, this will set all relevant opts to works as if in build. Will NOT overwrite the others (eg -count)

$database = glob("~wormpub/DATABASES/current_DB") if $test;
die "no database\n" unless $database;


my %genes_exons;
my %genes_span;
my %cDNA;
my %cDNA_span;
my %transcript_span;
my $gff_file;
my @ordered_genes;
my %matches_3;

my %EST_pairs;
my @non_overlapping_ESTs;
my %cdna2gene;
my %gene2cdnas;
my @cds_objs;
my %cds_index;


#setup directory for +transcript
my $transcript_dir = "$database/TRANSCRIPTS";
&run_command("mkdir $transcript_dir") unless -e "$transcript_dir";
&run_command("rm -f $transcript_dir/*");


my $coords;
# write out the transcript objects
# get coords obj to return clone and coords from chromosomal coords
$coords = Coords_converter->invoke($database, $new_coords);

$gff_file = $gff if $gff;

my @chromosomes = qw(I II III IV V X);
foreach my $chrom ( @chromosomes ) {

  #Make sure these are empty
  %genes_exons = ();
  %genes_span= ();
  %cDNA = ();
  %cDNA_span = ();
  %transcript_span = ();
  @ordered_genes = ();
  %gene2cdnas = ();
  my @cdna_objs;

  my $index = 0;
  
  unless ($gff) {
    $gff_file = "$database/CHROMOSOMES/CHROMOSOME_$chrom.gff";
  }

  open (GFF,"<$gff_file") or die "gff_file\n";
  $log->write_to("reading gff file $gff_file\n");
  # C43G2	  curated	  exon         10841   10892   .       +       .       CDS "C43G2.4"
  # C43G2   curated         coding_exon          10841   10892   .       +       0       CDS "C43G2.4"
  # C43G2   curated         CDS     10841   13282   .       +       .       CDS "C43G2.4"
  # C43G2   BLAT_EST_BEST   similarity   8877    9029    100     +       .       Target "Sequence:yk1125e01.5" 37 189

  # parse GFF file to get CDS, exon and cDNA info
  while (<GFF>) {
    my @data = split;
    
    #  GENE STRUCTURE
    if( $data[1] eq "curated" ) {
      $data[9] =~ s/\"//g;
      if( $data[2] eq "CDS" ) {
	# GENE SPAN
	$genes_span{$data[9]} = [($data[3], $data[4], $data[6])];

      }
      elsif ($data[2] eq "exon" ) {
	# EXON 
	$genes_exons{$data[9]}{$data[3]} = $data[4];
      }
    }
    
    # BLAT DATA
    elsif( ($data[1] eq "BLAT_EST_BEST") or ($data[1] eq "BLAT_OST_BEST") or ($data[1] eq "BLAT_mRNA_BEST") ){
      $data[9] =~ s/\"//g;
      $data[9] =~ s/Sequence:// ;
      $cDNA{$data[9]}{$data[3]} = $data[4];
      # keep min max span of cDNA
      if( !(defined($cDNA_span{$data[9]}[0])) or ($cDNA_span{$data[9]}[0] > $data[3]) ) {
	$cDNA_span{$data[9]}[0] = $data[3];
	$cDNA_span{$data[9]}[2] = $data[6]; #store strand of cDNA
      } 
      if( !(defined($cDNA_span{$data[9]}[1])) or ($cDNA_span{$data[9]}[1] < $data[4]) ) {
	$cDNA_span{$data[9]}[1] = $data[4];
      }
    }
  }

  close GFF;

  &load_EST_dir(glob("~/TRANSCRIPTS/EST_dir"));
 # &checkData(\$gff); # this just checks that there is some BLAT and gene data in the GFF file
  &eradicateSingleBaseDiff;
 # &getESTpairs;

  #create transcript obj for each CDS
  foreach ( keys %genes_exons ) {
    next if $genes_span{$_}->[2] eq "-"; #only do fwd strand for now
    my $cds = CDS->new( $_, $genes_exons{$_}, $genes_span{$_}->[2], $chrom );
    push( @cds_objs, $cds);
    $cds_index{$_} = $index;
    $index++;
  }

  foreach ( keys %cDNA ) {
    my $cdna = SequenceObj->new( $_, $cDNA{$_}, $cDNA_span{$_}->[2] );
    push(@cdna_objs,$cdna);
  }



################################################
  #            DATA LOADED           #
#################################################


  foreach my $CDNA ( @cdna_objs) {
   next if ( defined $est and $CDNA->name ne "$est"); #debug line
    foreach my $cds ( @cds_objs ) {
      if( $cds->map_cDNA($CDNA) ) {
	$CDNA->mapped(1);
	#print $CDNA->name," overlaps ",$cds->name,"\n" if $really_verbose;
      }
    }
    last if $est;
  }

  print "Second round additions\n";

  foreach my $CDNA ( @cdna_objs) {
   next if ( defined $est and $CDNA->name ne "$est"); #debug line
    next if ( defined($CDNA->mapped) ); 
    foreach my $cds ( @cds_objs ) {
      if( $cds->map_cDNA($CDNA) ) {
	$CDNA->mapped(1);
	print $CDNA->name," overlaps ",$cds->name,"on 2nd round \n" if $really_verbose;
      }
    }
  }

  my $out_file = glob("~ar2/TRANSCRIPTS/out.ace");
  open (FH,">$out_file") or die "cant open $out_file\n";
  foreach my $cds (@cds_objs ) {
    $cds->report(*FH, $coords);
  }

  last if $gff; # if only doing a specified gff file exit after this is complete
}

if( $load_transcripts ) {
  &run_command("cat $transcript_dir/transcripts_*.ace > $transcript_dir/transcripts_all.ace");
  &run_command("echo \"pparse $transcript_dir/transcripts_all.ace\nsave\nquit\" | $tace -tsuser transcripts $database");
}

if( $load_matches ) {
  &run_command("cat $transcript_dir/chromosome*_matching_cDNA.ace > $transcript_dir/matching_cDNA_all.ace");
  &run_command("echo \"pparse $transcript_dir/matching_cDNA_all.ace\nsave\nquit\" | $tace -tsuser matching_cDNA $database");
}

#$log->mail("$maintainers");

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
    $log->write_to( "removing small cDNA mismatches\n\n\n");
    foreach my $cdna_hash (keys %cDNA ) {
      my $last_key;
      my $check;
      foreach my $exons (sort { $cDNA{$cdna_hash}->{$a} <=> $cDNA{$cdna_hash}->{$b} } keys %{$cDNA{$cdna_hash}}) {
	print "\n############### $cdna_hash #############\n" if $really_verbose;
	print "$exons -> $cDNA{$cdna_hash}->{$exons}\n" if $really_verbose;
	my $new_last_key = $exons;
	if( $last_key ) {
	  if( $cDNA{$cdna_hash}->{$last_key} >= $exons - $gap ) { #allows seq error gaps up to $gap bp
	    $cDNA{$cdna_hash}->{$last_key} = $cDNA{$cdna_hash}->{$exons};
	    delete $cDNA{$cdna_hash}->{$exons};
	    $check = 1;
	    $new_last_key = $last_key;
	  }
	}
	$last_key = $new_last_key;
      }
      if ( $check ) {
	foreach my $exons (sort keys  %{$cDNA{$cdna_hash}}) {
	  print "single base diffs removed from $cdna_hash\n" if $really_verbose;
	  print "$exons -> $cDNA{$cdna_hash}->{$exons}\n" if $really_verbose;
	}
      }
    }
  }




sub checkOverlappingTranscripts  {
  my $chrom = shift;
  $log->write_to( "checking overlapping transcripts for chromosome_$chrom\n");
  my @ordered_transcripts;
  foreach my $transcript ( sort { $transcript_span{$a}[0]<=>$transcript_span{$b}[0]  } keys %transcript_span) {
    push (@ordered_transcripts, $transcript);
  }

  my $trans_count = scalar (@ordered_transcripts);

  open (OLT,">$transcript_dir/overlapping_transcripts_$chrom") or warn "cant open $transcript_dir/overlapping_transcripts\n";
  TRANS:
    for( my $i = 0;$i < $trans_count; $i++) {
      for( my $j = 1; $j < 6; $j++ ) {
	next unless ($i-$j >= 0 and $transcript_span{ $ordered_transcripts[$i-$j] });
	if ( $transcript_span{ $ordered_transcripts[$i] }->[0] < $transcript_span{ $ordered_transcripts[$i-$j] }->[1] ) {
	  {
	    my $gene1 = substr ($ordered_transcripts[($i-$j)] ,0,-6);
	    my $gene2 = substr ($ordered_transcripts[$i] ,0,-6);
	    next TRANS if &checkIsos(\$gene1, \$gene2, *OLT) == 1;
	  }
	}
	next if (($i+$j)>$trans_count);
	next unless $transcript_span{ $ordered_transcripts[$i+$j] };
       if ( $transcript_span{ $ordered_transcripts[$i] }->[1] > $transcript_span{ $ordered_transcripts[$i+$j] }->[0] ) {
	 {
	   my $gene1 = substr ($ordered_transcripts[($i+$j)] ,0,-6);
	   my $gene2 = substr ($ordered_transcripts[$i] ,0,-6);
	   next TRANS if &checkIsos(\$gene1, \$gene2, *OLT) == 1;	 
	}
       }
      }
    }
  sub checkIsos
    {
      my $gene1 = shift;
      my $gene2 = shift;
      my $file  = shift;
      if( $genes_span{$$gene1}->[2] eq $genes_span{$$gene2}->[2] ){ # check same strand strand
	# dismiss isoforms
	my ($gene1_name) = $$gene1 =~ /(\w+\.\d+)\w?/;
	my ($gene2_name) = $$gene2 =~ /(\w+\.\d+)\w?/;
	
	unless( $gene1_name eq $gene2_name ) {
	  print $file "overlapping transcripts\t$$gene1\t$$gene2\n";
	  return 1;
	}
      }
      return 0;
    }
  close OLT;
}

sub check_opts {
  $log->write_to("checking options . . . \n\n\n");
  # sanity check options
  if( $help ) {
    system("perldoc $0");
    exit(0);
  }

  unless ($transcript ) {
    if( $load_transcripts ) {
      die print "\n\n\nyou have chosen to load_transcripts without generating them\n Try adding -transcript if you want to build them first\n\n\n\n";
    }
    if( $overlap_check ) {
      die print "\n\n\nyou have chosen to check for overlapping transcripts without generating them\n Try adding -transcript if you want to build them first\n\n\n\n";
    }
    if( $build ) {
      $database = "/wormsrv2/autoace";
      $transcript = 1;
      $show_matches = 1;
      $load_transcripts = 1;
      $load_matches = 1;
      $overlap_check = 1;
      $new_coords = 1;
    }
  }

  if( defined($load_matches) & !defined($show_matches ) ) {
    die print "\n\n\nyou have chosen to load_matches without generating them\n Try adding -show_matches if you want to generate them first\n\n\n\n";
  }
}

sub checkData
  {
    my $file = shift;
    die "There's no BLAT data in the gff file $$file\n" if scalar keys %cDNA_span == 0;
    die "There are no genes in the gff file $$file\n" if scalar keys %genes_span == 0;
  }


sub getESTpairs
  {
    unless ( -e "$database/wquery/cDNApairs.def") {
      warn "$database/wquery/cDNApairs.def does not exist - inclusion of 3'EST not overlapping CDS will fail  . . ."  ;
      return 1;
    }

    my $command="Table-maker -p $database/wquery/cDNApairs.def\nquit\n";    
    open (TACE, "echo '$command' | $tace $database |");
    while (<TACE>) {
      chomp;
      s/\"//g;
      my @data = split;
      $EST_pairs{$data[0]} = $data[1];
    }
    close TACE;
  }

sub check_3_ends 
  {
    foreach my $cdna ( @non_overlapping_ESTs ) {
      if ( $EST_pairs{$cdna} ) {
	# paired read matches a gene
	my $paired_read = $EST_pairs{$cdna};
	my $potential_gene = $cdna2gene{$paired_read};
	next unless $potential_gene;

	#orientation.
	next unless $genes_span{$potential_gene}[2] eq $cDNA_span{$cdna}[2];
	if ( $genes_span{$potential_gene}[2] eq "+" ) {
	  if ( $genes_span{$potential_gene}[1] < $cDNA_span{$cdna}[0] and
	       $cDNA_span{$cdna}[0] - $genes_span{$potential_gene}[1] < $UTR_range ) {
	    push( @{$gene2cdnas{$_}}, $cdna);
	  }
	}
	elsif ( $genes_span{$potential_gene}[0] > $cDNA_span{$cdna}[1] and
		$genes_span{$potential_gene}[0] - $cDNA_span{$cdna}[1]  < $UTR_range ) {
	  push( @{$gene2cdnas{$_}}, $cdna);
	}
      }
    }
  }
	


###################################################################################

sub run_command{
  my $command = shift;
  $log->write_to( &runtime, ": started running $command\n");
  my $status = system($command);
  if($status != 0){
    $errors++;
    $log->write_to("ERROR: $command failed\n");
  }
  $log->write_to( &runtime, ": finished running $command\n");

  # for optional further testing by calling subroutine
  return($status);
}

#####################################################################################


sub load_EST_dir
  {
    my $file = shift;
    open (EST, "<$file") or die "cant open EST file $file\t$!\n";

    while (<EST>) {
      s/\"//g;
      if( /(\S+)\s+(EST_\d)/ ) {
	my $EST = $1;
	my $read_dir = $2;
	
	if( $cDNA_span{$EST} ) {
	  my $GFF_strand = $cDNA_span{$EST}->[2];
	CASE:{
	    ($GFF_strand eq "+" and $read_dir eq "EST_5") && do {
	      $cDNA_span{$EST}->[2] = "+";
	      last CASE;
	    };
	    ($GFF_strand eq "+" and $read_dir eq "EST_3") && do {
	      $cDNA_span{$EST}->[2] = "-";
	      last CASE;
	    };
	    ($GFF_strand eq "-" and $read_dir eq "EST_5") && do {
	      $cDNA_span{$EST}->[2] = "-";
	      last CASE;
	    };
	    ($GFF_strand eq "-" and $read_dir eq "EST_3") && do {
	      $cDNA_span{$EST}->[2] = "+";
	      last CASE;
	    };
	  }
	}
      }
    }
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


The script outputs up to 3 files for each chromosome.

  chromosome*_transcripts.ace    -     transcript object ace file 

  chromosome*_matching_cDNA.ace  -     matching_cDNA assignments

  chromosome*_matching_cDNA.dat  -     matching_cDNA assignments as Data::Dumper hash   gene => [cdna, cdna]

=back

=head2 transcript_builder.pl arguments:

=over 4

=item * verbose and really_verbose  -  levels of terminal output
  
=item *  est:s     - just do for single est 

=item * count     - reports number of matching_cDNAs for each gene 

=item * report    - reports list of matching_cDNA for each gene

=item * gap:s      - when building up cDNA exon structures from gff file there are often single / multiple base pair gaps in the alignment. This sets the gap size that is allowable [ defaults to 5 ]

=item *  transcript   - write transcipt files

=item * gff:s         - pass in a gff file to use

=item * show_matches   - write the matching_cDNA ace files

=item * database:s      - database directory to use. Expects gff files to be in CHROMSOMES subdirectory.

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
