#!/usr/local/bin/perl5.6.1 -w                    
#
# map_alleles.pl
#
# by Anthony Rogers
#
# This maps alleles to the genome based on their flanking sequence
#
# Last updated by: $Author: ar2 $                      # These lines will get filled in by cvs and helps us
# Last updated on: $Date: 2003-01-30 11:45:15 $        # quickly see when script was last changed and by whom


use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Ace;

use Getopt::Long;
#######################################
# command-line options                #
#######################################

my ($debug, $update, $limit);
my $database;
my $ver;
my $verbose;
my $restart = "go";
my $help; { `perldoc $0`;};
my $no_geneace;

# $debug   -  all output goes to ar/allele_mapping

GetOptions( "debug"     => \$debug,
	    "limit=s"   => \$limit,
	    "database=s"=> \$database,
	    "WS=s"      => \$ver,
	    "verbose"   => \$verbose,
	    "help"      => \$help,
	    "restart=s" => \$restart,
	    "no_geneace"=> \$no_geneace
	  );

if ($help) { print `perldoc $0`;exit;}

##############
# variables  #
##############

# Most checking scripts should produce a log file that is a) emailed to us all
# and b) copied to /wormsrv2/logs

my $maintainers = "All";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
my $log;
$ver = &get_wormbase_version unless defined $ver;


my (%chromosomeI_clones, %chromosomeII_clones, %chromosomeIII_clones, %chromosomeIV_clones, %chromosomeV_clones, %chromosomeX_clones);

my @chromo_clones_refs = (
			  \%chromosomeI_clones, 
			  \%chromosomeII_clones, 
			  \%chromosomeIII_clones, 
			  \%chromosomeIV_clones,
			  \%chromosomeV_clones,
			  \%chromosomeX_clones
			 );

my $data_dump_dir;

my $allele_fa_file;
my $genome_fa_file;
my $ace_file;
my $strain_file;

if( $debug )  {
  
  $data_dump_dir = glob("~ar2/allele_mapping");
  $allele_fa_file = "$data_dump_dir/alleles.fa";
  $genome_fa_file = "$data_dump_dir/genome.fa";
  $log        = "$data_dump_dir/map_alleles.$rundate";
  $database = "/wormsrv2/current_DB/" unless $database;
  $ver--;
  $ace_file = "$data_dump_dir/mapped_alleles.ace";
}
else { 
  $log        = "/wormsrv2/logs/map_allelesWS$ver.$rundate.$$";
  $allele_fa_file = "/wormsrv2/tmp/alleles.fa";
  $genome_fa_file = "/wormsrv2/tmp/genome.fa";
  $ace_file = "/wormsrv2/autoace/MAPPINGS/allele_mapping.WS$ver.ace";
  $database = "/wormsrv2/autoace/" unless $database;
}

#update clone position hashes
&UpdateHashes(\%chromosomeI_clones, "CHROMOSOME_I.clone_path.gff");
&UpdateHashes(\%chromosomeII_clones,  "CHROMOSOME_II.clone_path.gff");
&UpdateHashes(\%chromosomeIII_clones, "CHROMOSOME_III.clone_path.gff");
&UpdateHashes(\%chromosomeIV_clones, "CHROMOSOME_IV.clone_path.gff");
&UpdateHashes(\%chromosomeV_clones, "CHROMOSOME_V.clone_path.gff");
&UpdateHashes(\%chromosomeX_clones, "CHROMOSOME_X.clone_path.gff");


##########  File handles etc #############
open (LOG,">$log") or die "cant open $log\n\n";
print LOG "$0 start at $runtime on $rundate\n----------------------------------\n\n";

my $geneace_update = "/wormsrv2/autoace/MAPPINGS/map_alleles_geneace_update$ver.ace";
my $geneace_update_delete = "/wormsrv2/autoace/MAPPINGS/map_alleles_geneace_update_delete$ver.ace";
unless ($no_geneace) {
  open (GENEACE,">$geneace_update") or die "cant open $geneace_update: $!\n";
  open (GEN_DEL,">$geneace_update_delete") or die "cant open $geneace_update_delete\n";
}
open (OUT,">$ace_file") or die "cant open $ace_file\n";
#open (STR,">$strain_file") or die "cant open $strain_file\n";


########### database accesss ####################

#get allele info from database
my $db = Ace->connect(-path  => $database) || do { print  "$database Connection failure: ",Ace->error; die();};
my @alleles = $db->fetch(-query =>'Find Allele;flanking_sequences');


########## parse GFF to determine where the genes are ##############
my (%chromI_genes, %chromII_genes, %chromIII_genes, %chromIV_genes, %chromV_genes, %chromX_genes);
my @hashrefs = (

		\%chromI_genes, 
		\%chromII_genes,
		\%chromIII_genes,
		\%chromIV_genes,
		\%chromV_genes,
		\%chromX_genes
	       );

my %allele2gene;
my @chromosomes = qw( I II III IV V X );
my $i = 0;
foreach (@chromosomes) {
  my $gff = "/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS/CHROMOSOME_$_.genes.gff";
  if ($debug ) {
    $gff = "/wormsrv2/autoace/GFF_SPLITS/WS94/CHROMOSOME_$_.genes.gff";
  }
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

my %superlink_coords;
&UpdateSuperlinkCoords;

####### allele mapping loop ######################

my %allele_data;   #  allele => [ (0)type, (1)5'flank_seq , (2)3'flank_seq, (3)CDS, (4)end of 5'match, (5)start of 3'match , (6)clone, (7)chromosome, (8)strains]
my $count = 0;
my $scoreCutoff;  #SCAN score is no of matching bases
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
my $onSuperlink;
my $go = 0;
ALLELE:
foreach $allele (@alleles)
  {
    undef $onSuperlink;
    $name = $allele->name;

    # debug facility - this bit is so that it can restart from given allele name
    unless ("$restart" eq "go"){
      if ("$restart" eq "$name") {
	$restart = "go";
      }
      else { 
	print "skipping $name\n";
	next;
      }
    }

    # debug facility - after the restart means you can specify where to start and how many to do 
    if( $limit ) {
      last if $count++ > $limit;
    }

    $left = $allele->Flanking_sequences->name;
    $right = $allele->Flanking_sequences->right->name;

    next unless ($left and $right);
    #$type = $allele->Allele_type->name;
    print "mapping $name\n";
    #$allele_data{$name}[0] = $type;
    $allele_data{$name}[1] = $left;
    $allele_data{$name}[2] = $right;
    
    #allow 2 bp mismatch in mapping (using smaller flank seq)
    if( length $right < length $left ) {
      $scoreCutoff = (length $right) - 2;
    }
    else {
      $scoreCutoff = (length $left) - 2;
    }

    #get to the sequence object

    $sequence = $allele->Sequence;
    
    unless( $sequence ) {
      # allele has no Sequence Tag - should not happen
      # for 1st run get seq from Predicted_gene
      my $p_gene = $allele->Predicted_gene;
      
      unless($p_gene ) {
	#this allele has no Sequence or Predicted_gene
	# can we get it via the locus 

	if(  $allele->Gene ) {
	  $p_gene = $allele->Gene->Genomic_sequence;
	}
      }

      unless ( $p_gene ){
	print LOG "$allele has no predicted gene $p_gene\n";
	next ALLELE;
      }

      $p_gene =~ /(\S+)\./;
      $sequence = $1;
      $sequence = $db->fetch(Sequence => "$sequence");
      unless( $sequence ) {
	print LOG "$name - cant get a Sequence from Sequence tag or Predicted_gene\n";
	next ALLELE;
      }
    }     


    #if seq is not valid )ie no source check for other sequences
    unless ($sequence->Source) {
      $sequence = $allele->Sequence(1)->down;
      next unless $sequence;
      next unless( $sequence->Source );
    }      
    
    #find strains that contain this allele
    my @strains = $allele->Strain;
    foreach (@strains) {
      $allele_data{$name}[8] .= $_->name."***";
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
	
	$onSuperlink = 1 if ($sequence =~ /SUPER/);
	if ($onSuperlink != 1) 
	  {
	    if(  &MapAllele == 1 ) { #only try and map to clone if sequence isn't a SUPERLINK
	      # allele mapped to clone
	      print "Did $allele with clone\n" if $verbose;
	    }
	    
	    #try and map to reverse strand
	    else
	      {
		# try mapping to reverse strand
		if( &MapAllele('r') == 1 ) {
		  # need to make some coordinate adjustments to account for - strand
		  # get clone details by getting csome and using chromo_clones_refs array to retrieve correct hash of clone positions
		  my $csome_index = $roman2num{ $allele_data{$name}[7] };
		  my $clone_info  = $chromo_clones_refs[ $csome_index ]{"$clone"};
		  my $clone_length = $$clone_info[1] - $$clone_info[0] + 2; # add 2 to deal with aceDB seq starting at -1 (i think thats why!) 
		  
		  my $new_5 = $clone_length - $allele_data{$name}[5];
		  my $new_3 = $clone_length - $allele_data{$name}[4];
		  
		  $allele_data{$name}[4] = $new_5;
		  $allele_data{$name}[5] = $new_3;
		}
		else {
		  print LOG "$name failed mapping\n";
		  print  "$name failed mapping\n";
		}
	      }
	  }
	else {
	  #$source is a clone obj
	  if( $source )
	    {	    
	      if( &MapAllele == 1 ) {
		#allele mapped to superlink
		$allele_data{$name}[6] = $source->name;
		print "Did $allele with superlink\n" ;
	      }
	    }	
	  else {
	    print LOG "$allele $sequence has no Source\n";
	    next;
	  }
	}

	&findOverlapGenes($name);
	&outputAllele($name);
      }
    else {
      print LOG "$allele - cant find valid sequence\n";
    }
  }

print "$count alleles\n";

print LOG "Update files for geneace allele mappings are available - \n$geneace_update\n$geneace_update_delete\n\n";

$db->close;


close OUT;
close GEN_DEL unless $no_geneace;
close GENEACE unless $no_geneace;

##############################
# read acefiles into autoace #
##############################
unless ( $debug ) {
  print LOG "\nStart parsing $ace_file in to $database\n\n";
  my $command =<<END;
pparse $geneace_update_delete
pparse $ace_file
save
quit
END
  my $tace = &tace;
  open (TACE,"| $tace -tsuser map_allele $database") || die "Couldn't open tace connection to $database\n";
  print TACE $command;
  close (TACE);
  print LOG "finished parsing\n";
}
# close LOG and send mail
close LOG;
print LOG "$0 end at ",&runtime," \n-------------- END --------------------\n\n";
&mail_maintainer("map_alleles","$maintainers","$log");

exit(0);

#######################################
#                                     #
#          SUB ROUTINES               #
#                                     #
#######################################


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

sub outputAllele
  {
    my $to_dump = shift;
    if( $allele_data{$to_dump}[3] and $allele_data{$to_dump}[4] and  $allele_data{$to_dump}[5]) { 
      
      print OUT "\nSequence : \"$allele_data{$to_dump}[6]\"\nAllele $to_dump $allele_data{$to_dump}[4] $allele_data{$to_dump}[5]\n";

      if( $allele2gene{$to_dump} ) {
	my @myStrains;
	my @affects_genes = split(/\s/,"$allele2gene{$to_dump}");
	
	# in Predicted gene object
	foreach my $ko (@affects_genes) {
	  print OUT "\nSequence : \"$ko\"\nHas_allele $to_dump\n";
	}

	# in Allele object
	print OUT "\nAllele : $to_dump\n";
	print GEN_DEL "\nAllele : $to_dump\n-D Predicted_gene\n-D Sequence\n" unless $no_geneace;# remove current sequence and predicted genes from Geneace

	print GENEACE "\nAllele : \"$to_dump\"\nSequence \"$allele_data{$to_dump}[6]\"\n" unless $no_geneace;# allele -> sequence
	if( $allele_data{$to_dump}[8] ) {
	  @myStrains = split(/\*\*\*/,"$allele_data{$to_dump}[8]");
	}
	foreach my $ko (@affects_genes) {
	  #allele - seq connection
	  print OUT "Predicted_gene $ko\n";
	  print GENEACE "Predicted_gene $ko\n" unless $no_geneace;# update geneace with allele -> Predicted_genes


#         NOT DOING THIS ANY MORE
#	  foreach my $str (@myStrains) {
#	    #strain - seq connection
#	    print STR "Strain : \"$str\"\n";
#	    print STR "Sequence $ko\n\n";
#	  }


	}
      }
      else {
	print "no overlapping gene for $to_dump\n" if $verbose;
      }

      #map position on genome
      #(0)type, 
      #(1)5'flank_seq ,
      #(2)3'flank_seq
      #(3)CDS
      #(4)end of 5'match
      #(5)start of 3'match
      #(6)clone
      #(7)chromosome
      #(8)strains containing this allele

    }
  }


sub MapAllele
  {
    # This routine converts whatever is in $SEQspan in to "genome_fa_file
    # then performs a SCAN analysis against whatever is in it
    # then checks the results to see if both flanking sequences have been matched.  If they have
    # it fills in the allele_data fields and returns 1; otherwise 0;


    &makeTargetfile;
    print "Starting $name SCAN with score cutoff $scoreCutoff\n" if $verbose;
    
    my $reverse = shift;
    if( $reverse ) {
      print "revcomping the target sequence file $genome_fa_file\n" if $verbose;
      my $rev = $genome_fa_file."_rev";
      `revcomp $genome_fa_file > $rev`;
      `mv $rev $genome_fa_file`;
    }
    # SCAN 
    my $scan = "/usr/local/pubseq/bin/scan -f -t";
    my @data;
    my $end5;
    my $end3;
    open (SCAN, "$scan $scoreCutoff $genome_fa_file $allele_fa_file |") or warn "cant pipe SCAN : $!\n";
    while (<SCAN>)
      {
	# scan output-  gk2_3       1     30  56424 56453     30  100 % match, 0 copies, 0 gaps
	s/^\s+//;  #remove any space at begining - some have it some dont f's up where things split to! 
	@data = split(/\s+/,$_);
	if( defined( $data[0]) )
	  {
	    print "@data\n" if $verbose;
	    if( $data[0] =~ /(\S+)_([3,5])/ )
	      {
		if( "$2" eq "5" ) {
		  if( $data[4] ){
		    $end5 = $data[4]+1;     #end of 5'
		  }
		}
		else{
		  if( $data[3] ){
		    $end3 = $data[3]-1;     #start of 3'
		  }
		}
	      }
	  }
      }
    close SCAN;

    if( defined($end5) and defined($end3) ) {
      $allele_data{$name}[4] = $end5;
      $allele_data{$name}[5] = $end3;
      $allele_data{$name}[6] = $source->name;
      return 1;
    }
    else { return 0; }
  }

sub makeTargetfile
  {
    open (GFA,">$genome_fa_file" ) or die "cant open genome file for $allele\n";
    print GFA uc "$SEQspan";
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


#CHROMOSOME_IV   Genomic_canonical       Sequence        16392472        16417376        .       +       .       Sequence "Y65A5A"
sub UpdateHashes #(hash, file)
  {
    my $dir = "/wormsrv2/autoace/yellow_brick_road/";
    my $hash = shift;
    my $file = shift;
    $file = $dir.$file;
    my @data;
    open (CLONES,"<$file") or die "cant open $file";
    while (<CLONES>)
      {
	@data = split(/\s+/, $_);
	unless( $data[0] =~m/\#/ )
	  {
	    $data[9] =~ s/\"//g;
	    $$hash{$data[9]}[0] = $data[3];
	    $$hash{$data[9]}[1] = $data[4];
	  }
      }
  }

sub UpdateSuperlinkCoords
  {
    #this creates a hash of SUPERLINK => [start, end] for determining overlapping genes#
    my @chromosomes = $db->fetch(Sequence => 'CHROMO*');
    foreach (@chromosomes) {
      my @superlinks = $_->Subsequence;
      foreach (@superlinks) {
	my $SLstart = $_->right->name;
	my $SLend = $_->right(2)->name;
	$superlink_coords{$_->name} = [$SLstart, $SLend];
	print "@{$superlink_coords{$_}}\n" ;
      }
    }
  }


sub findOverlapGenes
  {
    my $allele = shift;
    #now find any genes that overlap with the alleles
    my $chromosome = $allele_data{$allele}[7];
    my $clone =  $allele_data{$allele}[6];
    my $clone_pos =  $chromo_clones_refs[ $roman2num{"$chromosome"} ];
    my $genelist= $arrayrefs[ $roman2num{"$chromosome"} ];
    my $chromosomal_coords_5;
    my $chromosomal_coords_3;

    if( $onSuperlink ) {
     #get coords from the SUPERLINKS
      my $Superlink_start = $superlink_coords{$clone}[0];
      $chromosomal_coords_5 = $Superlink_start + $allele_data{$allele}[4];
      $chromosomal_coords_3 = $Superlink_start + $allele_data{$allele}[5];
    }
    else {
      # get from clone hashes
      $chromosomal_coords_5 = $$clone_pos{"$clone"}[0] + $allele_data{$allele}[4];
      $chromosomal_coords_3 = $$clone_pos{"$clone"}[0] + $allele_data{$allele}[5];
    }
    
    # if allele of > 1bp on - strand 3' coord will be bigger than 5'
    if ($chromosomal_coords_5 < $chromosomal_coords_3) {
      my $tmp = $chromosomal_coords_3;
      $chromosomal_coords_3 = $chromosomal_coords_5;
      $chromosomal_coords_5 = $tmp;
    }
      
    
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






# Add perl documentation in POD format
# This should expand on your brief description above and add details of any options
# that can be used with the program.  Such documentation can be viewed using the perldoc
# command.



__END__

=pod

=head2 NAME - map_alleles.pl

=head1 USAGE

=over 4

=item map_alleles.pl  [-debug -limit -database=s -WS=s -verbose -restart=s]

=back

This script:

Gets alleles with flanking sequence from designated database and maps them to a clone or superlink using SCAN, then checks which if any CDSs they overlap with.

Also requires that the allele has as a "seed" sequence either a Sequence or Predicted_gene (a locus with an associated
genomic_sequence will also work).

Also writes two files for updating allele->Sequence and Allele->Predicted_gene in geneace, one to remove the current connections and one to enter the new ones.

Outputs acefiles which are loaded in to the same database.

map_alleles.pl MANDATORY arguments:

=over 4

=item none

=back

map_alleles.pl  OPTIONAL arguments:

=over 4

=item -help

this stuff

=back

=over 4

=item -debug 
   
output goes to ar2 and uses current_DB

=back

=over 4

=item -limit 

limt the number of alleles mapped (debug tool)

=back

=over 4

=item -database 

specify which database to read info from and load mapping results in to

=back

=over 4

=item -restart 

choose which allele to start with. all preceding (alphabetically) alleles will be skipped

=back

=over 4

=item -verbose 

greater indication of what procedured are being used to map the allele

=back

=over 4

=item -ver 

select a version of the database other than that being built

=back




=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=item

=item Must be run AFTER gff splitting has produced CHROMOSOME_*.genes.gff

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers ( ar2@sanger.ac.uk)

=back

=cut
