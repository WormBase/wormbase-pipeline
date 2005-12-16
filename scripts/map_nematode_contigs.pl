#!/usr/local/bin/perl5.8.0 -w
#
# map_nematode_contigs.pl
# 
# by Gary Williams
#
# Map WashU and Nembase EST contigs to genome and write out file to send to authors
# of these two data sets for their web sites to point back at us
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2005-12-16 11:18:55 $

use strict;
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
#use Wormbase;
use Getopt::Long;
use Coords_converter;
use Carp;
use Log_files;
#use Ace;
#use Sequence_extract;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose);
my ($all, $washu, $nembase);
my $maintainers = "All";

my $washu_email = "jmartin\@watson.wustl.edu";
my $nembase_email = "mark.blaxter\@ed.ac.uk";



GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "all"      	 => \$all,
	    "washu"      => \$washu,
	    "nembase"    => \$nembase,
	    );

my $log = Log_files->make_build_log($debug);

# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if ($debug) {
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

# in test mode?
if ($test) {
  print "In test mode\n";
  $washu_email = "gw3\@sanger.ac.uk";
  $nembase_email = "gw3\@sanger.ac.uk";
}



##########################
# MAIN BODY OF SCRIPT
##########################

my $type;
my @best_hits;
my $coords = Coords_converter->invoke;
#my $gffdir      = glob("~wormpub/DATABASES/current_DB/CHROMOSOMES/");   # chromosome GFF files
my $gffdir      = glob("/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS/");   # chromosome GFF files

if ($washu || $all) {
  print "Mapping WashU contigs\n" if ($verbose);
  $log->write_to("Mapping WashU contigs\n");
  @best_hits = &map_hits_to_genome("washu");
  &map_to_gene("washu", @best_hits);
  &mail_author($washu_email, "/tmp/washu_result.dat");
}

if ($nembase || $all) {
  print "Mapping Nembase contigs\n" if ($verbose);
  $log->write_to("Mapping Nembase contigs\n");
  @best_hits = &map_hits_to_genome("nembase");
  &map_to_gene("nembase", @best_hits);
  &mail_author($nembase_email, "/tmp/nembase_result.dat");
}




# Close log files and exit
$log->write_to("\n\nFinished.\n");


$log->mail();

print "Finished.\n" if ($verbose);
exit(0);






##############################################################
#
# Subroutines
#
##############################################################

# $type = "washu" (do mapping for WashU contigs)
# or "nembase" (do mapping for Nembase)

sub map_hits_to_genome {

  my ($type) = @_;

# set database paths
  my $blat_dir  = "/wormsrv2/autoace/BLAT";

  my %best;



##########################################################################################
# map the blat hits to ace - i.e. process blat output (*.psl) file into set of ace files #
##########################################################################################

 
  print "Reading BLAT data\n";

  open(BLAT, "<$blat_dir/${type}_out.psl")    or die "Cannot open $blat_dir/${type}_out.psl $!\n";

  # loop through each blat hit
  while (<BLAT>) {
    next unless (/^\d/);
    my @f            = split "\t";

    my $match        = $f[0];                    # number of bases matched by blat
    my $strand       = $f[8];                    # strand that match is on
    my $query        = $f[9];                    # query sequence name
    my $query_size   = $f[10];                   # query sequence length
    my $superlink    = $f[13];                   # name of superlink that was used as blat target sequence
    my $slsize       = $f[14];                   # superlink size
    my $matchstart   = $f[15];                   # target (superlink) start coordinate...
    my $matchend     = $f[16];                   # ...and end coordinate
    my $block_count  = $f[17];                   # block count
    my @lengths      = split (/,/, $f[18]);      # sizes of each blat 'block' in any individual blat match
    my @query_starts = split (/,/, $f[19]);      # start coordinates of each query block
    my @slink_starts = split (/,/, $f[20]);      # start coordinates of each target (superlink) block
    
    
    # calculate (acedb) score for each blat match
    # new way of calculating score, divide by query size rather than sum of matching blocks, 
    my $score = ($match/$query_size)*100;
    
    my @exons = ();
    my ($chrom, $chrom_start);
    my $chrom_end;
    
    
    for (my $x = 0;$x < $block_count; $x++) {
      my $start = $slink_starts[$x];
      
      my ($query_start,$query_end);
      my $gff_strand;
      
      ($chrom, $chrom_start) = $coords->Coords_2chrom_coords($superlink, $start);
      $chrom_end = $chrom_start + $lengths[$x] -1;
      if (($strand eq '--') || ($strand eq '+-')) {
	$gff_strand = '-';
      } else {
	$gff_strand = '+';
      }
      
      # blatx 6-frame translation v 6-frame translation
      
      my $temp;
      if (($strand eq '++') || ($strand eq '-+')) {
	$query_start = $query_starts[$x] +1;
	$query_end   = $query_start + $lengths[$x] -1;
	if ($strand eq '-+') {
	  $temp        = $query_end;
	  $query_end   = $query_start;
	  $query_start = $temp; 
	}
      }
      elsif (($strand eq '--') || ($strand eq '+-')) {
	$temp         = $chrom_start;
	$chrom_start = $chrom_end;
	$chrom_end   = $temp;
	
	$query_start  = $query_size  - $query_starts[$x];
	$query_end    = $query_start - $lengths[$x] +1;
	
	if ($strand eq '--') {
	  $temp        = $query_end;
	  $query_end   = $query_start;
	  $query_start = $temp; 
	}
      }			
      push @exons, [$chrom_start,$chrom_end,$query_start,$query_end];
    }
    
    # collect best hits for each query sequence 
    # Choose hit with highest score (% of query length which are matching bases) 
    # If multiple hits have same scores (also meaning that $match must be same) store 
    # details of extra hits against same primary key in %best
    if (exists $best{$query}) {
      if (($score > $best{$query}->{'score'})) { 
	# Add all new details if score is better...
	$best{$query}->{'score'} = $score;
	$best{$query}->{'match'} = $match;
	@{$best{$query}->{'entry'}} = ({'chrom' => $chrom,'exons' => \@exons});
      }
      elsif($score == $best{$query}->{'score'}){
	#...only add details (name and coordinates) of extra hits if scores are same
	push @{$best{$query}->{'entry'}}, {'chrom' => $chrom,'exons' => \@exons};
      }
    }
    else {
      $best{$query}->{'match'} = $match;
      $best{$query}->{'score'} = $score;
      @{$best{$query}->{'entry'}} = ({'chrom' => $chrom,'exons' => \@exons});
    }
  }
  close(BLAT);
  

  # produce list of best matches 
  
  my @best_hits;		# LoL of best-match details
  
  foreach my $found (sort keys %best) {
    if (exists $best{$found}) {
      foreach my $entry (@{$best{$found}->{'entry'}}) {
	if (@{$best{$found}->{'entry'}} < 2) { # only want best hits at a unique position
	  my $chrom = $entry->{'chrom'};
	  foreach my $ex (@{$entry->{'exons'}}) {
	    my $score        = $best{$found}->{'score'};
	    my $chrom_start  = $ex->[0];
	    my $chrom_end    = $ex->[1];
	    my $query_start  = $ex->[2];
	    my $query_end    = $ex->[3];
	    
	    # print output 
	    my $query = $found;
	    my $gff_strand;
	    if ($chrom_start < $chrom_end) {
	      $gff_strand = '+';
	    } else {
	      $gff_strand = '-';
	      my $temp = $chrom_start;
	      $chrom_start = $chrom_end;
	      $chrom_end = $temp;
	  }
	    push @best_hits, [ ($chrom, $chrom_start, $chrom_end, $gff_strand, $score, $query) ];
	  }
	  
	}
      }
    }
  }

  # sort the best hist by the chromosome and then the end position
  @best_hits = sort {$a->[0] cmp $b->[0] or $a->[2] <=> $b->[2]} @best_hits;
  
  return @best_hits;
}

##########################################

# maps the best hit to a gene
# writes the results to a file

sub map_to_gene() {
  my ($type, @best_hits) = @_;

  my @chromosomes = qw( I II III IV V X );                            # chromosomes
  my $resultsfile = "/tmp/${type}_result.dat";
  my %unique_results = ();
  
  if ($test) {
    @chromosomes = qw( III );
  }

  foreach my $chromosome (@chromosomes) {
 
    print "Reading chromosome $chromosome\n" if ($verbose);
    print "gffdir = $gffdir\n" if ($verbose);

    # loop through the GFF file
    my @f;
    my @chrom;

 #   open (GFF, "<$gffdir/CHROMOSOME_${chromosome}.gff") || die "Failed to open chromosome gff file $gffdir/CHROMOSOME_${chromosome}.gff\n";
   open (GFF, "<$gffdir/CHROMOSOME_${chromosome}.WBgene.gff") || die "Failed to open chromosome gff file $gffdir/CHROMOSOME_${chromosome}.WBgene.gff\n";
    while (<GFF>) {
      chomp;
      s/^\#.*//;
      next unless /\S/;
      @f = split /\t/;
      next unless ($f[1] eq "gene");

      my ($name) = ($f[8] =~ /Gene\s+\"(\S+)\"/); # get the gene name
      #print "gene name = $name\n";
      push @chrom, [ ($f[3], $f[4], $f[6], $name) ]; # start, end, strand, gene name
    }
    close(GFF);
                  
    print "Read chromosome $chromosome, sorting...\n" if ($verbose);
    # sort by start pos
    @chrom = sort {$a->[0] <=> $b->[0]} @chrom;


  ###############
  # map contigs #
  ###############
  
    print "Find overlaps for BLATed contigs\n" if ($verbose);
    my $result;
    
    my $previous_start=0;		# remember the first gene that matched the previous BLAT hit
    my $set_flag;		        # true if we have set $previous_start to the first gene match for this BLAT 

    # get the first gene
    my $current_gene_count = 0;
    my $current_gene = $chrom[$current_gene_count];
    my ($gene_start, $gene_end, $gene_strand, $gene_name) = @{$current_gene};
    print "First gene = $gene_start, $gene_end, $gene_strand, $gene_name\n" if ($verbose);
    print "Have ",scalar(@best_hits), " best hits to look at\n" if ($verbose);

    # loop through BLAT results
    # @best_hits is array-reference of [ ($chrom, $chrom_start, $chrom_end, $gff_strand, $score, $query) ];
    foreach my $hit_aref (@best_hits) {
    
      my @this_hit = @{$hit_aref};
      # now have @this_hit with the BLAT results in and @chrom with genes for this chromosome
      my ($blat_chrom, $blat_start, $blat_end, $blat_strand, $blatscore, $blat_query) = @this_hit;

      # only want to look at this chromosome
      if ($blat_chrom ne "CHROMOSOME_$chromosome") {
	next;			                          # next BLAT
      }

      print "Next BLAT hit = $blat_chrom, $blat_start, $blat_end, $blat_strand, $blatscore, $blat_query\n" if ($verbose);

      # step down here
LOOP:
      if ($blat_end < $gene_start) {                      # BLAT is before current gene - get next BLAT
	print "next BLAT (gene=$gene_name)\n" if ($verbose);
	$set_flag = 0;		                          # not yet found the first match for the next BLAT
	$current_gene_count = $previous_start;            # go back to the first gene to match the previous BLAT 
	$current_gene = $chrom[$current_gene_count];
	($gene_start, $gene_end, $gene_strand, $gene_name) = @{$current_gene};
	next;			                          # get next BLAT
      }
      if ($blat_start > $gene_end) {                      # BLAT is after current gene - get next gene
	print "next gene (BLAT=$blat_query)\n" if ($verbose);
	$current_gene_count++;	                          # next gene
	$current_gene = $chrom[$current_gene_count];
	($gene_start, $gene_end, $gene_strand, $gene_name) = @{$current_gene};
	goto LOOP;		                          # go to test this next gene
      }

      # we have a match
      # save results as key of hash to get unique ones
      print "Got a hit to a gene $blat_query = $gene_name\n" if ($verbose);
      $unique_results{"${blat_query}:$gene_name"} = 1; 
      if (!$set_flag) {
	$set_flag = 1;
	$previous_start = $current_gene_count;            # remember this first gene that matches
      }
      next;			                          # next BLAT

    }
  }
 
  # print out the unique contig-gene results
  # using a hash to get the unique results
  print "Writing results to $resultsfile\n" if ($verbose);
  open (OUT, ">$resultsfile") || die "Can't open $resultsfile\n";
  foreach my $key (keys %unique_results) {
    my ($c, $t) = split (/:/,$key);
    print OUT "$c\t$t\n";
  }
  close (OUT);
}

##########################################

# returns the index into @list of the gene that overlaps the contig
# else it returns -1
sub binary_chop() {
  # start and end and strand of contig to find overlaps to
  # list of details of genes to overlap with
  my ($contig_start, $contig_end, $contig_strand, @list) = @_;

  my $list_length = scalar(@list);

  my $chop_start = 0;
  my $chop_end = $list_length-1;
  my $i;

  while ($chop_start < $chop_end-1) {
    #print "$chop_start $chop_end\t";
    $i = $chop_start + int (($chop_end - $chop_start)/2);
    my ($gene_start, $gene_end, $gene_strand, $gene_name) =  @{$list[$i]};
    my $contig_middle = $gene_start + (($gene_end - $gene_start)/2);
    if (#$contig_strand eq $gene_strand &&                                     # no - don't check for strand
        (($contig_start >= $gene_start && $contig_start <= $gene_end) || # do we have an overlap?
        ($contig_end >= $gene_start && $contig_end <= $gene_end)) &&
        ($contig_middle >= $gene_start && $contig_middle <= $gene_end) # at least half of the contig match overlaps
        ) {
      return $i;                # found
    } elsif ($contig_start < $gene_start) {
      $chop_end = $i;
    } else {
      $chop_start = $i;
    }
  }

  return -1;                    # not found

}

##########################################

sub mail_author {
  my ($address, $file) = @_;

  open (OUTLOG,  "|/bin/mailx -r \"wormbase\@sanger.ac.uk\"  -s \"Your contigs mapped to Wormbase genes\" $address ");
  if ( $file ) {
    open (READLOG, "<$file") || die "Didn't send mail properly\n\n";
    while (<READLOG>) { 
      print OUTLOG "$_";
    }
    close READLOG;
  }
}


##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################




# Add perl documentation in POD format
# This should expand on your brief description above and add details of any options
# that can be used with the program.  Such documentation can be viewed using the perldoc
# command.


__END__

=pod

=head2 NAME - map_nematode_contigs.pl

=head1 USAGE

=over 4

=item map_nematode_contigs.pl [-options]

=back

This script maps WashU and Nembase EST contigs to genome and writes out files and sends them to authors
of these two data sets for their web sites to point back at us.


script_template.pl MANDATORY arguments:

=over 4

=item none

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -all, do both WashU and Nembase contigs

=back

=over 4

=item -nembase, do the Nembase contigs 

=back

=over 4

=item -washu, do the WashU contigs

=back

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Verbose/Debug mode
 
=back

=over 4

=item -test, Test mode, generate the acefile but do not upload themrun the script, but don't change anything

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back
                                                                                             

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Gary Williams

=back

=cut


