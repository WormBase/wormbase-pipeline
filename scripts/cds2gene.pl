#!/usr/local/bin/perl5.8.0 -w
#
# cds2gene.pl
# 
# by Keith Bradnam, aged 12 and half
#
# A script to create ?Gene objects for CDSs, Pseudogenes, and Transcript objects not yet
# linked to an existing ?Locus object
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2004-03-24 16:55:21 $   

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
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
  my $query = "Find CDS WHERE Method = curated AND Species = \"Caenorhabditis elegans\" AND NOT Locus";
  push(my @CDSs, $db->find($query));

  # device for tracking multiple isoforms of same gene
  my $last_sequence_name = "";
  my $last_gene = "";

  foreach my $gene (@CDSs){

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
      print OUT "Public_name $gene\n";
      print OUT "Species \"$species\"\n";
      print OUT "Version_change 1 now \"WBPerson1971\" Imported \"Initial conversion from CDS class of WS121\"\n";
      print OUT "Live\n";
      print OUT "CDS $gene\n\n";

    }      

  }
  close(OUT);
  $db->close;
}


###########################################################################################
# processes Transcript class to create Gene objects for Transcript not attached to loci
# can only do this *after* CDS class has been processed and loaded to test database
# this is because a transcript might be an isoform of an existing CDS in which case
# we must look up the CDS Gene ID first
###########################################################################################

sub process_transcript_class{

  # start from highest ID that was used previously
  my $id = $id_limit;

  open(OUT,">processed_transcript_class.ace") || die "Couldn't open output file\n";

  my $db_path = "$database";
  my $db = Ace->connect(-path  => $db_path) || do { print "Connection failure: ",Ace->error; die();};
  
  # we want all genes, all species??? Don't want coding_transcripts
  # leave other species for later
  my $query = "Find Transcript WHERE NOT Method = coding_transcript AND Species = \"Caenorhabditis elegans\"";
  push(my @transcripts, $db->find($query));

  # device for tracking multiple isoforms of same gene
  my $last_sequence_name = "";
  my $last_gene = "";

  foreach my $gene (@transcripts){

    # set sequence name, i.e. chop off trailing a,b,c etc.
    my $sequence_name = $gene;
    $sequence_name =~ s/[a-z]$//;
    
    # does this isoform belong to the same gene as the last one
    # if so, we can append the Transcript name to the last gene
    if($sequence_name eq $last_sequence_name){
      print OUT "Gene : \"$last_gene\"\n";
      print OUT "Transcript $gene\n\n";
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
      print OUT "Public_name $gene\n";
      print OUT "Species \"$species\"\n";
      print OUT "Version_change 1 now \"WBPerson1971\" Imported \"Initial conversion from Transcript class of WS121\"\n";
      print OUT "Live\n";
      print OUT "Transcript $gene\n\n";

    }      

  }
  close(OUT);
  $db->close;
}


exit(0);


 
