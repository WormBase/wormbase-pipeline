#!/usr/local/bin/perl5.8.0 -w
#
# cds2gene.pl
# 
# by Keith Bradnam, aged 12 and half
#
# A script to create ?Gene objects for CDSs, Pseudogenes, and Transcript objects not yet
# linked to an existing ?Gene object
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2004-05-27 14:46:40 $   

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Getopt::Long;


################################################################
# set up input/output streams
################################################################

my $cds; # flag to process CDS class
my $transcript; # flag to process transcript class
my $pseudogene; # flag to process pseudogene class
my $database; # target database to query
my $id_limit; # what is the highest gene ID already used so far

GetOptions ("cds"        => \$cds,
	    "transcript" => \$transcript,
	    "pseudogene" => \$pseudogene,
	    "database=s" => \$database,
	    "id_limit=i" => \$id_limit);


# grab hash of CDS/Transcript/Pseudogene identifiers to gene IDs
my %worm_gene2cgc = &FetchData('worm_gene2cgc_name');


&process_cds_class if ($cds);
&process_transcript_class if ($transcript);
&process_pseudogene_class if ($pseudogene);



##############################################################################
# processes CDS class to create Gene objects for CDS not attached to loci
##############################################################################

sub process_cds_class{

  # start from highest ID that was used previously
  my $id = $id_limit;

  open(OUT,">processed_cds_class.ace") || die "Couldn't open output file\n";

  my $db_path = "$database";
  my $db = Ace->connect(-path  => $db_path) || do { print "Connection failure: ",Ace->error; die();};
  
  # we want all genes, all species??? But what about Ensembl, twinscan etc.
  # could start with just elegans CDS (ignore briggsae and elegans history genes for now?)
  my $query = "Find CDS WHERE Method = curated AND Species = \"Caenorhabditis elegans\" AND NOT Gene";
  push(my @CDSs, $db->find($query));

  # device for tracking multiple isoforms of same gene
  my $last_sequence_name = "";
  my $last_gene = "";

  foreach my $gene (@CDSs){

    # check that this CDS doesn't already have a gene ID in geneace
    if ($worm_gene2cgc{$gene}){
      print "ERROR: $gene already exists in worm_genes2cgc hash, skipping this gene\n"; 
      # can now ignore this gene
      next;
    }

    # set sequence name, i.e. chop off trailing a,b,c etc.
    my $sequence_name = $gene;
    $sequence_name =~ s/[a-z]$//;
    
    # does this isoform belong to the same gene as the last one
    # if so, we can append the CDS name to the last gene
    if($sequence_name eq $last_sequence_name){
      print OUT "Gene : \"$last_gene\"\n";
      print OUT "CDS $gene\n\n";
    }
    # if not, it's a new gene, write output
    else{

      $id++;
      my $id_padded = sprintf "%08d" , $id;
      my $name = "WBGene$id_padded";  
      
      #grab species info
      my $species = $gene->Species;
      
      $last_sequence_name = $sequence_name;
      $last_gene = $name;
      
      print OUT "Gene : \"$name\"\n";
      print OUT "Version 1\n";
      print OUT "Sequence_name $sequence_name\n";
      print OUT "Public_name $sequence_name\n";
      print OUT "Species \"$species\"\n";
      print OUT "Version_change 1 now \"WBPerson1971\" Imported \"Initial conversion from CDS class of WS125\"\n";
      print OUT "Live\n";
      print OUT "CDS $gene\n\n";

    }      

  }
  close(OUT);
  $db->close;
}


###########################################################################################
# processes Transcript class to create Gene objects for Transcript not attached to Genes
# can only do this *after* CDS class has been processed and loaded to test database
# this is because a transcript might be an isoform of an existing CDS in which case
# we must look up the CDS Gene ID first
###########################################################################################

sub process_transcript_class{
  
  # unused gene IDs
  my @unused_ids = (7475,7539,8104,8135,8179,9206,9527,9555,9936,10276,308,309,310,311,10376,10430,10431,10733,10777,11346,11470,12292,12350,12377,12653,12752,13010,13386,13918,13948,13993,14026,14126,14147);

  # start from highest ID that was used previously
  my $id = $id_limit;

  open(OUT,">processed_transcript_class.ace") || die "Couldn't open output file\n";

  my $db_path = "$database";

  # connect to source database, but also to geneace for extra checking
  my $db  = Ace->connect(-path  => $db_path) || do { print "Connection failure: ",Ace->error; die();};
  my $db2 = Ace->connect(-path  => "/wormsrv1/geneace") || do { print "Connection failure: ",Ace->error; die();};  

  # we want all genes, all species??? Don't want coding_transcripts
  # leave other species for later
  my $query = "Find Transcript WHERE NOT Method = coding_transcript AND Species = \"Caenorhabditis elegans\" AND NOT Gene AND NOT Method = history";
  push(my @transcripts, $db->find($query));

  # device for tracking multiple isoforms of same gene
  my $last_sequence_name = "";
  my $last_gene = "";

  foreach my $gene (@transcripts){

    # set sequence name, i.e. chop off trailing a,b,c etc.
    my $sequence_name = $gene;
    $sequence_name =~ s/[a-z]$//;
    
    # first check to see that this name isn't already in geneace, i.e. where an other isoform is a CDS
    my ($obj) = $db2->fetch(-class=>'Gene_name',-name=>"$sequence_name");
    if ($obj){
      my $gene_id = $obj->Sequence_name_for;
      print "WARNING: Transcript $gene ($sequence_name) in camace already exists in geneace under Gene $gene_id\n";
      print OUT "Gene : \"$gene_id\"\n";
      print OUT "Transcript $gene\n\n";
      # can now skip to next gene
      next;
    }
    

    # does this isoform belong to the same gene as the last one
    # if so, we can append the Transcript name to the last gene
    if($sequence_name eq $last_sequence_name){
      print OUT "Gene : \"$last_gene\"\n";
      print OUT "Transcript $gene\n\n";
    }
    # if not, it's a new gene, write output
    else{

      # use unused ids first if possible
      my $final_id;
      my $id_padded;

      if($unused_ids[0]){
	my $tmp_id = $unused_ids[0];
	shift(@unused_ids);
	$id_padded = sprintf "%08d" , $tmp_id;
      }
      else{
	$id++;
	$id_padded = sprintf "%08d" , $id;

      }
      my $name = "WBGene$id_padded";  

      
      #grab species info
      my $species = $gene->Species;
      
      $last_sequence_name = $sequence_name;
      $last_gene = $name;
      
      print OUT "Gene : \"$name\"\n";
      print OUT "Version 1\n";
      print OUT "Sequence_name $sequence_name\n";
      print OUT "Public_name $sequence_name\n";
      print OUT "Species \"$species\"\n";
      print OUT "Version_change 1 now \"WBPerson1971\" Imported \"Initial conversion from Transcript class of WS125\"\n";
      print OUT "Live\n";
      print OUT "Transcript $gene\n\n";

    }      

  }
  close(OUT);
  $db->close;
}


###########################################################################################
# processes Pseudogene class to create Gene objects for Pseudogenes not attached to Genes
# can only do this *after* CDS class has been processed and loaded to test database
# this is because a Pseudogene might be an isoform of an existing CDS in which case
# we must look up the CDS Gene ID first
###########################################################################################

sub process_pseudogene_class{
  
  # start from highest ID that was used previously
  my $id = $id_limit;

  open(OUT,">processed_pseudogene_class.ace") || die "Couldn't open output file\n";

  my $db_path = "$database";

  # connect to source database, but also to geneace for extra checking
  my $db  = Ace->connect(-path  => $db_path) || do { print "Connection failure: ",Ace->error; die();};
  my $db2 = Ace->connect(-path  => "/wormsrv1/geneace") || do { print "Connection failure: ",Ace->error; die();};  

  # use subclass to get just elegans pseudogenes
  my $query = "Find elegans_pseudogenes";
  push(my @pseudogenes, $db->find($query));

  # device for tracking multiple isoforms of same gene
  my $last_sequence_name = "";
  my $last_gene = "";

  foreach my $gene (@pseudogenes){

    # set sequence name, i.e. chop off trailing a,b,c etc.
    my $sequence_name = $gene;
    $sequence_name =~ s/[a-z]$//;
    
    # first check to see that this name isn't already in geneace, i.e. where an other isoform is a CDS
    my ($obj) = $db2->fetch(-class=>'Gene_name',-name=>"$sequence_name");
    if ($obj){
      my $gene_id = $obj->Sequence_name_for;
      print "WARNING: Pseudogene $gene ($sequence_name) in camace already exists in geneace under Gene $gene_id\n";
      print OUT "Gene : \"$gene_id\"\n";
      print OUT "Pseudogene $gene\n\n";
      # can now skip to next gene
      next;
    }
    

    # does this isoform belong to the same gene as the last one
    # if so, we can append the Transcript name to the last gene
    if($sequence_name eq $last_sequence_name){
      print OUT "Gene : \"$last_gene\"\n";
      print OUT "Pseudogene $gene\n\n";
    }
    # if not, it's a new gene, write output
    else{

      $id++;
      my $id_padded = sprintf "%08d" , $id;

      
      my $name = "WBGene$id_padded";  

      
      #grab species info
      my $species = $gene->Species;
      
      $last_sequence_name = $sequence_name;
      $last_gene = $name;
      
      print OUT "Gene : \"$name\"\n";
      print OUT "Version 1\n";
      print OUT "Sequence_name $sequence_name\n";
      print OUT "Public_name $sequence_name\n";
      print OUT "Species \"$species\"\n";
      print OUT "Version_change 1 now \"WBPerson1971\" Imported \"Initial conversion from Pseudogene class of WS125\"\n";
      print OUT "Live\n";
      print OUT "Pseudogene $gene\n\n";

    }      

  }
  close(OUT);
  $db->close;
}






exit(0);


 
