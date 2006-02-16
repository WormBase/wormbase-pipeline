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
# Last edited on: $Date: 2006-02-16 15:42:47 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($all, $washu, $nembase);


# the people to send the results to
my $washu_email; #= "jmartin\@watson.wustl.edu";
my $nembase_email; #= "mark.blaxter\@ed.ac.uk";

print "DEBUG addresses\n";
$washu_email = "gw3\@sanger.ac.uk";
$nembase_email = "gw3\@sanger.ac.uk";



GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "all"      	 => \$all,
	    "washu"      => \$washu,
	    "nembase"    => \$nembase,
	    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($wormbase->test) {
  print "In test mode\n" if ($verbose);

  # mail the results to this user
  my $me = `whoami`;
  $washu_email = "${me}\@sanger.ac.uk";
  $nembase_email = "${me}\@sanger.ac.uk";
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);



##########################
# MAIN BODY OF SCRIPT
##########################

my $type;
my $gffdir = $wormbase->gff_splits;         # AUTOACE GFF_SPLIT

if ($washu || $all) {
  $type = "WASHU";
  print "Mapping $type contigs\n" if ($verbose);
  print "Using temporary file: /tmp/${type}_result.dat\n";
  $log->write_to("Mapping $type contigs\n");
  &map_to_gene($type);
  &mail_author($washu_email, "/tmp/${type}_result.dat");
}

if ($nembase || $all) {
  $type = "NEMBASE";
  print "Mapping $type contigs\n" if ($verbose);
  print "Using temporary file: /tmp/${type}_result.dat\n";
  $log->write_to("Mapping $type contigs\n");
  &map_to_gene($type);
  &mail_author($nembase_email, "/tmp/${type}_result.dat");
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


##########################################

# maps the best hit to a gene
# writes the results to a file

sub map_to_gene() {
  my ($type) = @_;

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
    my @best_hits = ();

    print "Reading genes\n" if ($verbose);
    open (GFF, "<$gffdir/CHROMOSOME_${chromosome}.WBgene.gff") || die "Failed to open gene gff file $gffdir/CHROMOSOME_${chromosome}_WBgene.gff\n";
    while (<GFF>) {
      chomp;
      s/^\#.*//;
      next unless /\S/;
      @f = split /\t/;

      my ($name) = ($f[8] =~ /Gene\s+\"(\S+)\"/); # get the gene name
      #print "gene name = $name\n";
      push @chrom, [ ($f[3], $f[4], $f[6], $name) ]; # start, end, strand, gene name
    }
    close(GFF);
                  
    # sort by start pos
    @chrom = sort {$a->[0] <=> $b->[0]} @chrom;


    print "Reading BLAT hits\n" if ($verbose);
    open (GFF, "<$gffdir/CHROMOSOME_${chromosome}.BLAT_$type.gff") || die "Failed to open BLAT gff file $gffdir/CHROMOSOME_${chromosome}_BLAT_$type.gff\n";
    while (<GFF>) {
      chomp;
      s/^\#.*//;
      next unless /\S/;
      @f = split /\t/;

      my ($query) = ($f[8] =~ /Target\s+\"Sequence:(\S+)\"/); # get the EST contig name
      push @best_hits, [ ($f[3], $f[4], $f[6], $query) ];
    }
    close(GFF);

    # sort the best hits by the start position
    @best_hits = sort {$a->[0] <=> $b->[0]} @best_hits;

  ###############
  # map contigs #
  ###############
  
    print "Find overlaps for BLATed contigs\n" if ($verbose);
    my $result;
    
    my $previous_start=0;		# remember the first gene that matched the previous BLAT hit
    my $set_flag = 0;		        # true if we have set $previous_start to the first gene match for this BLAT 

    # get the first gene
    my $current_gene_count = -1;

    # loop through BLAT results
    # @best_hits is array-reference of [ ($chrom_start, $chrom_end, $gff_strand, $query) ];
    foreach my $hit_aref (@best_hits) {
    
      my @this_hit = @{$hit_aref};
      # now have @this_hit with the BLAT results in and @chrom with genes for this chromosome
      my ($blat_start, $blat_end, $blat_strand, $blat_query) = @this_hit;
      print "Next BLAT hit = $blat_start, $blat_end, $blat_strand, $blat_query\n" if ($verbose);
 
# exhaustive all vs all mapping - takes seevral hours
      foreach my $current_gene (@chrom) {
	my ($gene_start, $gene_end, $gene_strand, $gene_name) = @{$current_gene};
	if ($gene_start < $blat_start && $gene_end > $blat_end) {
	  print "Got a hit to a gene $blat_query = $gene_name\n" if ($verbose);
	  $unique_results{"${blat_query}:$gene_name"} = 1;
	}
      }

      # step down through genes mapping algorithm - much faster
#      while (++$current_gene_count <= scalar(@chrom)-1) { # point to next gene
#	#print "Current gene = $current_gene_count, max = ".scalar(@chrom)."\n";
#	my $current_gene = $chrom[$current_gene_count];
#	my ($gene_start, $gene_end, $gene_strand, $gene_name) = @{$current_gene};
#	print "Next gene = $gene_start, $gene_end, $gene_strand, $gene_name\n" if ($verbose);

#	if ($gene_end < $blat_end) {
#	  next;			# get next gene

#	} elsif ($gene_start > $blat_end) {
#	  $current_gene_count = $previous_start-1;
#	  $set_flag = 0;	# reset flag for 'got a first match'
#	  last;			# get next blat

#	} else {
#	  # we have a match
#	  # save results as key of hash to get unique ones
#	  print "Got a hit to a gene $blat_query = $gene_name\n" if ($verbose);
#	  $unique_results{"${blat_query}:$gene_name"} = 1;
#	  if (!$set_flag) {	# if not yet got a match to a gene for this blat
#	    $set_flag = 1;
#	    $previous_start = $current_gene_count;            # remember this first gene that matches
#	  }
#	}
#      }
    }				# blat hit loop
  }				# chromosome loop

  # prepare hash for getting cds IDs from WBgene IDs
  my %wbgene_id2cds;
  $wormbase->FetchData("wbgene_id2cds",\%wbgene_id2cds);

  # print out the unique contig-gene results
  # using a hash to get the unique results
  print "Writing results to $resultsfile\n" if ($verbose);
  open (OUT, ">$resultsfile") || die "Can't open $resultsfile\n";
  foreach my $key (keys %unique_results) {
    my ($contig, $wbgene) = split (/:/,$key);
    # want to also get the CDS name (this misses out pseudogenes and RNA genes etc, so there are some blanks output)
    if (exists $wbgene_id2cds{$wbgene} ) {
      my $cds = $wbgene_id2cds{$wbgene};
      print OUT "$contig\t$wbgene\t@{$cds}\n";
    } else {
      print OUT "$contig\t$wbgene\t\n";
    }
  }
  
  close (OUT);

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

This script maps WashU and Nembase EST contigs to genome and writes
out files and sends them to authors of these two data sets for their
web sites to point back at us.


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

=item -test, Test mode, run the script only for chromosome III and mail the person running the script, not the collaborators.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back
                                                                                             

=head1 REQUIREMENTS

=over 4

=item 

=back

=head1 AUTHOR

=over 4

=item Gary Williams

=back

=cut


