#!/usr/local/bin/perl5.6.0 -w                    # perl5.6.0 and -w flag
#
# map_alleles.pl
#
# by Anthony Rogers
#
# This maps alleles to the genome based on their flanking sequence
#
# Last updated by: $Author: ar2 $                      # These lines will get filled in by cvs and helps us
# Last updated on: $Date: 2002-10-22 13:08:19 $        # quickly see when script was last changed and by whom


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

my (%chromosomeI_clones, %chromosomeII_clones, %chromosomeIII_clones, %chromosomeIV_clones, %chromosomeV_clones, %chromosomeX_clones);

my @chromo_clones_refs = (
			  \%chromosomeI_clones, 
			  \%chromosomeII_clones, 
			  \%chromosomeIII_clones, 
			  \%chromosomeIV_clones,
			  \%chromosomeV_clones,
			  \%chromosomeX_clones
			 );

if( defined($opt_d) ) {
  $allele_fa_file = glob("~ar2/allele_mapping/alleles.fa");
  $genome_fa_file = glob("~ar2/allele_mapping/genome.fa");
  $scan_file = glob("~ar2/allele_mapping/alleles.scan");
  $log        = glob("~ar2/allele_mapping/map_alleles.$rundate");
 # $database = glob("~ar2/testace/");
  $database = "/wormsrv2/current_DB/";
  $ver--;
  our $count;
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



my %allele_data;   #  allele => [ (0)type, (1)5'flank_seq , (2)3'flank_seq, (3)CDS, (4)end of 5'match, (5)start of 3'match , (6)clone]
my $count = 0;
my $scoreCutoff = 29;  #SCAN score is no of matching bases
my $source;
my $OLR;
my $OLL;
my $sourceSeq;
my $OLLseq = "";
my $OLRseq = "";
my $SEQspan;
my $sequence;
my $clone;
my $chromosome;
my $name;
my $type;
my $left;
my $right;
my $allele;

foreach $allele (@alleles)
  {
    if( $opt_d ) {
      last if $count++ > 3;
    }
    $name = $allele->name;print "$name\n";
    $type = $allele->Allelic_difference->name;
    $left = $allele->Flanking_sequences->name;
    $right = $allele->Flanking_sequences->right->name;
    
    $allele_data{$name}[0] = $type;
    $allele_data{$name}[1] = $left;
    $allele_data{$name}[2] = $right;

    #get to the sequence object
    $sequence = $allele->Sequence;
    #if seq is not valid )ie no source check for other sequences
    unless ($sequence->Source) {
      $sequence = $allele->Sequence(1)->down;
      next unless( $sequence->Source );
    }      
      
    if ( $sequence )
      {
	$allele_data{$name}[7] = &findChromosome( $sequence );
	
	#make the allele.fa file
	open (AFA, ">$allele_fa_file") or die "cant open $allele_fa_file\n";
	print AFA "\>$name\_5\n$allele_data{$allele}[1]\n";
	print AFA "\>$name\_3\n$allele_data{$allele}[2]\n";
	close AFA;

	$allele_data{$name}[3] = "unknown";

	#try and map to clone
	#examine the sequence name and extract clone name if required
	if( $sequence->name =~ m/(\w+)\.\w+/ ) {
	  $clone = $1;
	  $allele_data{$name}[3] = "$&";
	}
	else {
	  $clone = $sequence->name;
	}
	
	$source = $db->fetch(Sequence => "$clone");	
	$SEQspan = $source->asDNA;
	
	if( &MapAllele == 1){
	  #allele mapped to clone
	  $allele_data{$name}[6] = $source->name;
	  print "Did $allele with clone\n";
	}
	
	#try and map to OverlapLeft - Clone - OverLapRight concatenated sequence
	else
	  {
	    #$source is a clone obj
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
		if( &MapAllele == 0 ) {
		  print LOG "$name cannot be mapped\n";
		  
		}
	      }
	    else {
	      print LOG "$allele $sequence has no Source\n";
	      next;
	    }
	  }
      }
  }

#parse GFF to determine where the genes are
my (%chromI_genes, %chromII_genes, %chromIII_genes, %chromIV_genes, %chromV_genes, %chromX_genes);
my @hashrefs = (

		\%chromI_genes, 
		\%chromII_genes,
		\%chromIII_genes,
		\%chromIV_genes,
		\%chromV_genes,
		\%chromX_genes
	       );

my @chromosomes = qw( I II III IV V X );
my $i = 0;
foreach (@chromosomes) {
  my $gff = "/wormsrv2/autoace/GFF_SPLITS/WS$ver/CHROMOSOME_$_.genes.gff";
  my $chrom2gene = $hashrefs[$i];
  open (GFF,"<$gff") or die "cant open $gff\n";
  while(<GFF>) {
   # CHROMOSOME_I    curated Sequence        222722  223159  .       +       .       Sequence "Y48G1BM.3" wp_acc=CE26120
    if( /curated/ ) {
      my @data = split;
      $data[9] =~ s/\"//g;
      my $gene = $data[9];
      my $end5 = $data[3];
      my $end3 = $data[4];
      $$chrom2gene{$gene} = [$end5, $end3];
    }
  }
  $i++;
}
my (@chromI_genes,  @chromII_genes, @chromIII_genes, @chromIV_genes, @chromV_genes, @chromX_genes);
my @arrayrefs = (

		\@chromI_genes, 
		\@chromII_genes,
		\@chromIII_genes,
		\@chromIV_genes,
		\@chromV_genes,
		\@chromX_genes
	       );

#need to sort these in to an ordered array - ordered on 5' coord of gene
$i = 0;
foreach (@chromosomes) {
  my $chrom_array = $arrayrefs[$i];
  my $chrom2gene = $hashrefs[$i];
  foreach (sort { $$chrom2gene{$a}[1]<=>$$chrom2gene{$b}[1]} keys  %$chrom2gene )  {
    my @temp = ($_,$$chrom2gene{$_}[0],$$chrom2gene{$_}[1]);
    push(@$chrom_array,[ @temp ]);
  } 
  $i++;
}

#this hash is to get from roman numeral chromosome to array index to get references for
# @chromo_clones_refs @arrayrefs @hashrefs

my %roman2num = ( I   => '0',
		  II  => '1',
		  III => '2',
		  IV  => '3',
		  V   => '4',
		  X   => '5' 
		);

  my %allele2gene;
#now find any genes that overlap with the alleles
foreach $allele (keys %allele_data) {
  my $chromosome = $allele_data{$allele}[7];
  my $clone =  $allele_data{$allele}[6];
  my $clone_pos =  $chromo_clones_refs[ $roman2num{"$chromosome"} ];
  my $genelist= $arrayrefs[ $roman2num{"$chromosome"} ];

  my $chromosomal_coords_5 = $$clone_pos{"$clone"} + $allele_data{$allele}[4];
  my $chromosomal_coords_3 = $$clone_pos{"$clone"} + $allele_data{$allele}[5];
  
  foreach my $gene ( @$genelist ) {
    if(
       ( $chromosomal_coords_5 > $$gene[1] and $chromosomal_coords_5 < $$gene[2] ) or  # 5' end of allele is in gene
       ( $chromosomal_coords_3 > $$gene[1] and $chromosomal_coords_3 < $$gene[2] ) or  # 3' end of allele is in gene
       ( $chromosomal_coords_5 < $$gene[1] and $chromosomal_coords_3 > $$gene[2] )     # whole gene removed
      ) {
      
      $allele2gene{$allele} .= "$$gene[0] ";
    }
  }
}
  
my $output = glob("~ar2/allele_mapping/map_results");
open (OUT,">$output");
foreach (keys %allele_data )
  {
    if( $allele_data{$_}[3] and $allele_data{$_}[4] and  $allele_data{$_}[5]) { 
      #print OUT "$_ is a $allele_data{$_}[0] from $allele_data{$_}[4] to $allele_data{$_}[5] ( ",$allele_data{$_}[5] - $allele_data{$_}[4]," bp )in $allele_data{$_}[3]\n";

      print OUT "Sequence : \"$allele_data{$_}[6]\"\nAllele $_ $allele_data{$_}[4] $allele_data{$_}[5] $allele_data{$_}[0] \n\n";
      print "$_ $allele2gene{$_}\n";
      #map position on genome
      #(0)type, 
      #(1)5'flank_seq ,
      #(2)3'flank_seq
      #(3)CDS
      #(4)end of 5'match
      #(5)start of 3'match
      #(6)clone
      #(7)chromosome

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

sub findChromosome
  {
    my $csome;
    my $seqObj = shift;
    my $sourceName;
    while(! $csome ) {
      $sourceName = $seqObj->name;
      if( $sourceName =~ /CHROMOSOME_(\w+)/ ) {
	$csome = "$1";
      }
      else {
	$seqObj = $seqObj->Source;
      }
    }
    return "$csome";
  }



sub MapAllele
  {
    # This routine converts whatever is in $SEQspan in to "genome_fa_file
    # then performs a SCAN analysis against whatever is in it
    # then checks the results to see if both flanking sequences have been matched.  If they have
    # it fills in the allele_data fields and returns 1; otherwise 0;


    &makeTargetfile;
    print "Starting $name SCAN with score cutoff $scoreCutoff\n";

    my $scan = "/usr/local/pubseq/bin/scan -f -t";
    open (SCAN, ">$scan_file") or die "cant write $scan_file\n";
    print SCAN `$scan $scoreCutoff $genome_fa_file $allele_fa_file` or die "SCAN failed\n";
    close SCAN;

    #read in the results
    my @data;
    my $end5;
    my $end3;
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
		    $end5 = $data[5]+1;     #end of 5'
		  }
		}
		else{
		  if( $data[4] ){
		    $end3 = $data[4]-1;     #start of 3'
		  }
		}
	      }
	  }
      }
    close SCAN;

    if( defined($end5) and defined($end3) ) {
      $allele_data{$name}[4] = $end5;
      $allele_data{$name}[5] = $end3;
      return 1;
    }
    else { return 0; }
  }

sub makeTargetfile
  {
    open (GFA,">$genome_fa_file" ) or die "cant open genome file for $allele\n";
    print GFA "$SEQspan";
    close GFA;
  }

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
