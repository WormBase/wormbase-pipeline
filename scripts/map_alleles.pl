#!/usr/local/bin/perl5.6.0 -w                    # perl5.6.0 and -w flag
#
# map_alleles.pl
#
# by Anthony Rogers
#
# This maps alleles to the genome based on their flanking sequence
#
# Last updated by: $Author: ar2 $                      # These lines will get filled in by cvs and helps us
# Last updated on: $Date: 2002-08-22 10:46:50 $        # quickly see when script was last changed and by whom


use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Ace;

use Getopt::Std;
#######################################
# command-line options                #
#######################################

use vars qw($opt_d $opt_c);
# $opt_d debug   -  all output goes to ar/allele_mapping

getopts ('dc');

##############
# variables  #
##############

# Most checking scripts should produce a log file that is a) emailed to us all
# and b) copied to /wormsrv2/logs

my $maintainers = "All";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
my $log;
my $ver = &get_wormbase_version;

my $allele_fa_file;
my $genome_fa_file;
my $scan_file;
my $database;

my %chromosomeI_clones;
my %chromosomeII_clones;
my %chromosomeIII_clones;
my %chromosomeIV_clones;
my %chromosomeV_clones;
my %chromosomeX_clones;

if( defined($opt_d) ) {
  $allele_fa_file = glob("~ar2/allele_mapping/alleles.fa");
  $genome_fa_file = glob("~ar2/allele_mapping/genome.fa");
  $scan_file = glob("~ar2/allele_mapping/alleles.scan");
  $log        = glob("~ar2/allele_mapping/map_alleles.$rundate");
  $database = glob("~ar2/testace/");
 # $database = "/wormsrv2/current_DB/";
  $ver--;
}
else { 
  $log        = "/wormsrv2/logs/map_alleles.$rundate";
  $scan_file = glob("~ar2/allele_mapping/alleles.scan");
  $allele_fa_file = "/wormsrv2/autoace/BLATS/alleles.fa";
  $database = "/wormsrv2/autoace/";
}

if( $opt_c ){
  #update clone position hashes
  &UpdateHashes(\%chromosomeI_clones, "CHROMOSOME_I.clone_ends.gff");
  &UpdateHashes(\%chromosomeII_clones, "CHROMOSOME_II.clone_ends.gff");
  &UpdateHashes(\%chromosomeIII_clones, "CHROMOSOME_III.clone_ends.gff");
  &UpdateHashes(\%chromosomeIV_clones, "CHROMOSOME_IV.clone_ends.gff");
  &UpdateHashes(\%chromosomeV_clones, "CHROMOSOME_V.clone_ends.gff");
  &UpdateHashes(\%chromosomeX_clones, "CHROMOSOME_X.clone_ends.gff");
}
open (LOG,">$log");
print LOG "$0 start at $runtime on $rundate\n----------------------------------\n\n";
#get allele info from database
my $db = Ace->connect(-path  => $database) || do { print  "$database Connection failure: ",Ace->error; die();};
my @alleles = $db->fetch(-query =>'Find Allele;flanking_sequences');

#my @slinks = $db->fetch(sequence =>'SUPER*');
#foreach (@slinks){
#  print "$_ ", $_->name," ",$_->Source,"\n";
#}



my %allele_data;   #  allele => [ (0)type, (1)5'flank_seq , (2)3'flank_seq, (3)DNA_matched_to, (4)end of 5'match, (5)start of 3'match , (6)chromosome]
my $count = 0;
my $scoreCutoff = 29;  #SCAN score is no of matching bases
my $sequence;
my $source;
my $OLR;
my $OLL;
my $sourceSeq;
my $OLLseq = "";
my $OLRseq = "";
my $SEQspan;
my $chromosome;

foreach my $allele (@alleles)
  {
    my $name = $allele->name;print "$name\n";
    my $type = $allele->Allelic_difference->name;
    my $left = $allele->Flanking_sequences->name;
    my $right = $allele->Flanking_sequences->right->name;

    $allele_data{$name}[0] = $type;
    $allele_data{$name}[1] = $left;
    $allele_data{$name}[2] = $right;

    #get to the sequence object
    $sequence = $allele->Sequence;

    if ( $sequence )
      {
	#examine the sequence name and extract clone name if required
	my $clone;
	if( $sequence->name =~ m/(\w+)\.\w+/ ) {
	  $clone = $1;
	}
	else {
	  $clone = $sequence->name;
	}
	$source = $db->fetch(Sequence => "$clone");
	
	#findout which CHROMOSOME this clone is on for later mappings.
	$chromosome = $source;
	while( $chromosome )
	  {
	    if( $chromosome->name =~ m/CHROMOSOME_(\w+)/ ) {
	      $chromosome = $1;
	      last;
	    }
	    else {
	      my $newchrom = $chromosome;
	      $chromosome = $newchrom->Source;
	      print "\t$chromosome\n";
	    }
	  }
	$allele_data{$name}[6] = $chromosome;
      }
    else {
      print LOG "$allele has no sequence\n";
      next;
    }

    if( $source )
      {
	$sourceSeq = $source->asDNA();
	$allele_data{$name}[3] = $source->name;  #set 5' DNA to allele source in case of no OLL

	$OLR = $source->Overlap_Right;
	$OLL = $source->Overlap_Left;
	
	if( $OLL ) {
	  $OLLseq = $OLL->asDNA();
	  $allele_data{$name}[3] = $OLL->name;   #set 5' DNA to clone left of allele clone
	}
	else {
	  print LOG "$allele $source has no OLL\n";
	}
	
	if( $OLR ) {
	  $OLRseq = $OLR->asDNA();    
	}
	else {
	  print LOG "$allele $source has no OLR\n";
	}   

	$SEQspan = $OLLseq.$sourceSeq.$OLRseq;
	$SEQspan = &FASTAformat( $SEQspan );
      }
    else {
      print LOG "$allele $sequence has no Source\n";
      next;
    }

    #make the genomic seq fa file
    open (GFA,">$genome_fa_file" ) or die "cant open genome file for $allele\n";
    print GFA ">$name\n$SEQspan";
    close GFA;

    #make the allele.fa file
    open (AFA, ">$allele_fa_file") or die "cant open $allele_fa_file\n";
    print AFA "\>$name\_5\n$allele_data{$allele}[1]\n";
    print AFA "\>$name\_3\n$allele_data{$allele}[2]\n";
    close AFA;

    print "Starting $name SCAN with score cutoff $scoreCutoff\n";

    my $scan = "/usr/local/pubseq/bin/scan -f -t";
    open (SCAN, ">$scan_file") or die "cant write $scan_file\n";
    print SCAN `$scan $scoreCutoff $genome_fa_file $allele_fa_file` or die "SCAN failed\n";
    close SCAN;

    #read in the results
    my @data;
    open (SCAN, "<$scan_file") or die "cant read $scan_file\n";
    while (<SCAN>)
      {
	# scan output-  gk2_3       1     30  56424 56453     30  100 % match, 0 copies, 0 gaps
	@data = split(/\s+/,$_);
	if( defined( $data[1]) )
	  {
	    if( $data[1] =~ /(\w+)_([3,5])/ )
	      {
		if( "$2" eq "5" ) {
		  if( $data[5] ){
		    $allele_data{$1}[4] = $data[5];     #end of 5'
		  }
		  else{ print " 5\' stall - @data\n"; sleep; }
		}
		else{
		  if( $data[4] ){
		    $allele_data{$1}[5] = $data[4];     #start of 3'
		  }
		  else{ print "3\' stall - @data\n "; sleep; }
		}
	      }
	  }
      }
    close SCAN;
    last if $count++ > 1;
  }
my $output = glob("~ar2/allele_mapping/map_results");
open (OUT,">$output");
foreach (keys %allele_data )
  {
    if( $allele_data{$_}[3]) { 
      print OUT "$_ is a $allele_data{$_}[0] from $allele_data{$_}[4] to $allele_data{$_}[5] ( ",$allele_data{$_}[5] - $allele_data{$_}[4]," bp )in $allele_data{$_}[3] chromosome $allele_data{$_}[6]\n";

      #map position on genome
      #(0)type, 
      #(1)5'flank_seq ,
      #(2)3'flank_seq
      #(3)DNA_matched_to,
      #(4)end of 5'match
      #(5)start of 3'match
      #(6)chromosome]

      my $csome_lookup = $allele_data{$_}[6];
      my $map_hash; 
      if( "$csome_lookup" eq "I" ){
	$map_hash = \%chromosomeI_clones;
      }
      elsif ( "$csome_lookup" eq "II" ){
	$map_hash = \%chromosomeII_clones;
      }
      elsif ( "$csome_lookup" eq "III" ){
	$map_hash = \%chromosomeIII_clones;
      }
      elsif ( "$csome_lookup" eq "IV" ){
	$map_hash = \%chromosomeIV_clones;
      }
      elsif ( "$csome_lookup" eq "V" ){
	$map_hash = \%chromosomeV_clones;
      }
      elsif ( "$csome_lookup" eq "X" ){
	$map_hash = \%chromosomeX_clones;
      }
      else {
	warn "allele $_ has no chromsome\n";
	next;
      }

      print "$_ map position is ",$$map_hash{$allele_data{$_}[3]} + $allele_data{$_}[4],"\n";
	


    }
  }
print "$count alleles\n";
$db->close;
$runtime = `date +%H:%M:%S`; chomp $runtime;
print LOG "$0 end at ",$runtime," \n-------------- END --------------------\n\n";
close LOG;
`less $output`;

# Always good to cleanly exit from your script
exit(0);

sub FASTAformat
  {
    print "formatting sequence\n";
    my $seq = shift;
    my @chars = split( //,$seq);
    my $fasta;
    my $position = 1;
    my $charcount = 0;
    foreach ( @chars )
      {
	if( /[a-z]/ )
	  {
	    $fasta .= $_;
	    if( ($charcount % 10000) == 0 ){ print ".";}
	    $position++;
	    $charcount++;
	    if( $position == 60 ){ 
	      $position = 1;
	      $fasta .= "\n";
	    }
	  }
      }
    print "done . .($charcount bases) \n";
    return $fasta;
  }

sub UpdateHashes #(hash, file)
  {
    my $dir = "/wormsrv2/autoace/GFF_SPLITS/WS$ver/";
    my $hash = shift;
    my $file = shift;
    $file = $dir.$file;
    my @data;
    open (CLONES,"<$file") or die "cant open $file";
    while (<CLONES>)
      {
	if( $_ =~ m/left/ )
	  {
	    @data = split(/\s+/, $_);
	    unless( $data[0] =~m/\#/ )
	      {
		$data[8] =~ s/\"//g;
		$$hash{$data[8]} = $data[2];
	      }
	  }
      }
  }


# Add perl documentation in POD format
# This should expand on your brief description above and add details of any options
# that can be used with the program.  Such documentation can be viewed using the perldoc
# command.


__END__

=pod

=head2 NAME - map_alleles.pl

=head1 USAGE

=over 4

=item map_alleles.pl  [-options]

=back

This script:


map_alleles.pl MANDATORY arguments:

=over 4

=item none

=back

map_alleles.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers ( ar2@sanger.ac.uk)

=back

=cut
