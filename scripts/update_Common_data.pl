#!/usr/local/bin/perl5.8.0 -w
#
# update_Common_data.pl
# 
# by Anthony Rogers et al
#
# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2006-03-28 14:23:18 $

#################################################################################
# Initialise variables                                                          #
#################################################################################

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Data::Dumper;
use Getopt::Long;
use Storable;
use Log_files;

##############################
# command-line options       #
##############################


my $build;             # for when you want to query autoace, i.e. you are building.  Otherwise defaults to current_DB
my $test;              # test mode, uses ~wormpub/TEST_BUILD
my $all;               # performs all of the below options:
my $debug;
     
my $clone2accession;   # Hash: %clone2accession     Key: Genomic_canonical                 Value: GenBank/EMBL accession
                       # Hash: %accession2clone     Key: GenBank/EMBL accession            Value: Genomic_canonical
my $clone2size;        # Hash: %clonesize           Key: Genomic_canonical                 Value: DNA length
my $cds2wormpep;       # Hash: %cds2wormpep         Key: CDS name                          Value: Wormpep ID
                       # Hash: %wormpep2cds         Key: Wormpep ID                        Value: CDS name
my $cds2protein_id;    # Hash: %cds2protein_id      Key: CDS name                          Value: Protein_ID
                       # Hash: %protein_id2cds      Key: Protein_ID                        Value: CDS name
my $cds2cgc;           # Hash: %cds2cgc             Key: CDS name                          Value: CGC name
my $cds2status;        # Hash: %cds2status          Key: CDS name                          Value: Prediction status 
my $clone2seq;         # Hash: %clone2seq           Key: Genomic_canonical                 Value: DNA sequence (lower case)
my $clone2sv;          # Hash: %clone2sv            Key: Genomic_canonical                 Value: Sequence version (integer)
my $clone2type;        # Hash: %clone2type          Key: Genomic_canonical                 Value: Type information (Cosmid, Fosmid, YAC, Plasmid)
my $genes2lab;         # Hash: %worm_gene2lab       Key: Gene (CDS|Transcript|Pseudogene)  Value: From_laboratory (HX, RW, DRW)
my $worm_gene2cgc;     # Hash: %worm_gene2cgc_name  Key: CGC name                          Value: Gene ID, plus molecular name (e.g. AH6.1), also a hash of cgc_name2gene
my $worm_gene2geneID;  # Hash: %worm_gene2geneID    Key: Gene (CDS|Transcript|Pseudogene)  Value: Gene ID
my $worm_gene2class;   # Hash: %worm_gene2class     Key: CDS/Transcript/Pseudogene name    Value: 'CDS', 'Transcript', or 'Pseudogene'
my $estdata;           # Hash: %NDBaccession2est    Key: GenBank/EMBL accession            Value: EST name (WormBase)  
                       # Hash: %estorientation      Key: EST name (WormBase)               Value: EST_5 = 5, EST_3 = 3
my $est2feature;       # Hash: %est2feature         Key: EST name (WormBase)               Value: Feature name (WormBase)
my $cds2gene_id;       # Hash: %cds2gene_id         Key: CDS name                          Value: WBGene_id


my %Table_defs = (
		  'cds2protein'      => 'CommonData:CDS_proteinID.def',
		  'clone2sv'         => 'CommonData:Clone_SequenceVersion.def',
		  'clone2accession'  => 'CommonData:Clone_Accession.def', 
		  'clone2size'       => 'CommonData:Clone_Size.def',
                  'clone2type'       => 'CommonData:Clone_Type.def',
		  'cds2status'       => 'CommonData:CDS_Status.def',
                  'cds2cgc'          => 'CommonData:CDS_CGCname.def',
		  'est2feature'      => 'CommonData:EST_Feature.def',
		  'estdata'          => 'CommonData:EST_data.def',
		  'cds2lab'          => 'CommonData:CDS_Lab.def',
		  'pseudogene2lab'   => 'CommonData:Pseudogene_Lab.def',
		  'RNAgene2lab'      => 'CommonData:RNAgene_Lab.def',
		  'wormgene2cgc'     => 'CommonData:WormGene_CGCname.def',
		  'wormgene2geneid'  => 'CommonData:WormGene_GeneID.def',
		  'cds2wormpep'      => 'CommonData:CDS2wormpep'
		  );

my $store;
GetOptions (
	    "build"              => \$build,
	    "clone2acc"          => \$clone2accession,
	    "clone2size"         => \$clone2size,
	    "cds2wormpep"        => \$cds2wormpep,
	    "cds2pid"            => \$cds2protein_id,
	    "cds2status"         => \$cds2status,
	    "clone2seq"          => \$clone2seq,
	    "clone2sv"           => \$clone2sv,
	    "genes2lab"          => \$genes2lab,
	    "worm_gene2cgc"      => \$worm_gene2cgc,
	    "worm_gene2geneID"   => \$worm_gene2geneID,
	    "worm_gene2class"    => \$worm_gene2class,
	    "est"                => \$estdata,
	    "est2feature"        => \$est2feature,
	    "gene_id"            => \$cds2gene_id,
	    "all"                => \$all,
	    "test"               => \$test,
	    "store:s"            => \$store,
	    "debug:s"            => \$debug,
	    "clone2type"         => \$clone2type,
	    "cds2cgc"            => \$cds2cgc
	   );

my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}

my $log = Log_files->make_build_log( $wormbase );
##########################################
# Set up database paths                  #
##########################################

# Set up top level base directory which is different if in test mode
# Make all other directories relative to this

my $basedir    = $wormbase->basedir;
my $data_dir   = $wormbase->common_data;
my $wquery_dir = $wormbase->autoace."/wquery";
my $ace_dir    = $wormbase->autoace;

$log->write_to("Updating COMMON_DATA in $data_dir\n");
##############################
# ACEDB executables          #
##############################

our $tace = $wormbase->tace;


# run the various options depending on command line arguments
&write_cds2protein_id   if ($cds2protein_id   || $all);
&write_clone2accession  if ($clone2accession  || $all);
&write_clonesize        if ($clone2size       || $all);
&write_cds2wormpep      if ($cds2wormpep      || $all);
&write_cds2status       if ($cds2status       || $all);
&write_cds2cgc          if ($cds2cgc          || $all);
&write_clones2seq       if ($clone2seq        || $all);
&write_clones2sv        if ($clone2sv         || $all);
&write_clone2type       if ($clone2type       || $all);
&write_genes2lab        if ($genes2lab        || $all);
&write_worm_gene2class  if ($worm_gene2class  || $all);
&write_EST              if ($estdata          || $all);
&write_Feature          if ($est2feature      || $all);
&write_Gene_id          if ($cds2gene_id      || $all);
&write_worm_gene2geneID if ($worm_gene2geneID || $all);
&write_worm_gene2cgc    if ($worm_gene2cgc    || $all);

# hasta luego

$log->mail;
exit(0);


###################################
##########  SUBROUTINES  ##########
###################################


#######################################################################
# Data writing routines - actually create and dump the data           #
#######################################################################

sub write_cds2protein_id {

  $log->write_to("Updating cds2protein\n");
  my %cds2protein_id;
  my %protein_id2cds;
    
  ####################################################################
  # connect to AceDB using TableMaker,
  # populating %accession2name (maps embl accession to contig name)
  ####################################################################

  my $command = "Table-maker -p $wquery_dir/$Table_defs{'cds2protein'}\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  while (<TACE>) {
      chomp;
      s/\"//g;
      next if ($_ eq "");
      next if (/acedb\>/);
      last if (/\/\//);

      if (/(\S+)\s+(\S+)\s+(\d+)/) {
	  my $protein_id = "$2".".$3";
	  my $gene = $1;
	  $cds2protein_id{"$gene"}       = $protein_id;
	  $protein_id2cds{"$protein_id"} = $gene;
      }
  }
  close TACE;
  
  # now dump data to file
  
  open (G2P, ">$data_dir/cds2protein_id.dat") or die "cant write $data_dir/cds2protein_id.dat :$!";
  open (P2G, ">$data_dir/protein_id2cds.dat") or die "cant write $data_dir/protein_id2cds.dat :$! ";
  
  print G2P Data::Dumper->Dump([\%cds2protein_id]);
  print P2G Data::Dumper->Dump([\%protein_id2cds]);
  
  close G2P;
  close P2G;
    
}


########################################################################################################

sub write_clone2accession  {   
  $log->write_to("Updating clone2accession\n");

  my %clone2accession;
  my %accession2clone;
    
  # connect to AceDB using TableMaker,
  my $command="Table-maker -p $wquery_dir/$Table_defs{'clone2accession'}\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  while (<TACE>) {
      chomp;
      s/\"//g;
      next if ($_ eq "");
      next if (/acedb\>/);
      last if (/\/\//);
      if (/(\S+)\s+(\S+)/) {
	  $clone2accession{$1} = $2;
	  $accession2clone{$2} = $1;
      }
  }
  close TACE;
  
  # now dump data to file
  
  open (C2A, ">$data_dir/clone2accession.dat") or die "cant write $data_dir/clone2accession.dat :$!";
  open (A2C, ">$data_dir/accession2clone.dat") or die "cant write $data_dir/accession2clone.dat :$! ";
  
  print C2A Data::Dumper->Dump([\%clone2accession]);
  print A2C Data::Dumper->Dump([\%accession2clone]);
  
  close C2A;
  close A2C;
  
}

########################################################################################################

sub write_clones2sv  {   

    my %clone2sv;
    
    # connect to AceDB using TableMaker,
    my $command="Table-maker -p $wquery_dir/$Table_defs{'clone2sv'}\nquit\n";
    
    open (TACE, "echo '$command' | $tace $ace_dir |");
    while (<TACE>) {
	chomp;
	s/\"//g;
	next if ($_ eq "");
	next if (/acedb\>/);
	if (/^(\S+)\s+(\S+)\.(\d+)/) {
	    $clone2sv{$1} = $3;
	}
    }
    close TACE;
	
    # now dump data to file

    open (C2SV, ">$data_dir/clone2sv.dat") or die "cant write $data_dir/clone2sv.dat :$!";
    
    print C2SV Data::Dumper->Dump([\%clone2sv]);
    
    close C2SV;
    
}

########################################################################################################

sub write_clone2type  {   

    my %clone2type;
    
    # connect to AceDB using TableMaker,
    my $command="Table-maker -p $wquery_dir/$Table_defs{'clone2type'}\nquit\n";
    
    open (TACE, "echo '$command' | $tace $ace_dir |");
    while (<TACE>) {
	chomp;
	s/\"//g;
	next if ($_ eq "");
	next if (/acedb\>/);
	if (/^(\S+)\s+(\S+)/) {
	    $clone2type{$1} = $2;
	}
    }
    close TACE;
	
    # now dump data to file

    open (C2TYPE, ">$data_dir/clone2type.dat") or die "cant write $data_dir/clone2type.dat :$!";
    
    print C2TYPE Data::Dumper->Dump([\%clone2type]);
    
    close C2TYPE;
    
}


########################################################################################################

sub write_clonesize  {   

  my %clonesize;
    
  # connect to AceDB using TableMaker,
  my $command="Table-maker -p $wquery_dir/$Table_defs{'clone2size'}\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  while (<TACE>) {
      chomp;
      s/\"//g;
      next if ($_ eq "");
      next if (/acedb\>/);
      if (/(\S+)\s+(\d+)/) {
	  $clonesize{$1} = $2;
      }
  }
  close TACE;
  
  # now dump data to file

  open (CLONESIZE, ">$data_dir/clonesize.dat") or die "cant write $data_dir/clonesize.dat :$!";
  
  print CLONESIZE Data::Dumper->Dump([\%clonesize]);

  close CLONESIZE;
  
}


########################################################################################################

sub write_cds2wormpep  {   

  $log->write_to("Updating cds2wormpep\n");
  my $WPver = $wormbase->get_wormbase_version;

  # connect to AceDB using TableMaker,
  my $command="Table-maker -p $wquery_dir/$Table_defs{'cds2wormpep'}\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  my %cds2wormpep;
  my %wormpep2cds;
  while (<TACE>) {
    chomp;
    s/\"//g;
    next if ($_ eq "");
    next if (/acedb\>/);
    my @data = split;
    my $gene = $data[0];
    my $pep = $data[1];
    $pep =~ s/WP://;
    $cds2wormpep{$gene} = $pep;
    $wormpep2cds{$pep} .= "$gene ";
  }

  #now dump data to file
  open (C2G, ">$data_dir/wormpep2cds.dat") or die "$data_dir/wormpep2cds.dat";
  open (G2C, ">$data_dir/cds2wormpep.dat") or die "$data_dir/cds2wormpep.dat";

  print C2G Data::Dumper->Dump([\%wormpep2cds]);
  print G2C Data::Dumper->Dump([\%cds2wormpep]);

  close C2G;
  close G2C;
}

# Hash: %cds2status             Key: CDS name                          Value: Prediction status 
#
# The prediction status {Confirmed|Partially_confirmed|Predicted) based on mapping of transcript data
# to the coding exons.


sub write_cds2status  {   
    
    my %cds2status;
    my @f;

    # connect to AceDB using TableMaker,
    my $command="Table-maker -p $wquery_dir/$Table_defs{'cds2status'}\nquit\n";
    
    open (TACE, "echo '$command' | $tace $ace_dir |");
    while (<TACE>) {
	chomp;
	s/\"//g;
	next if ($_ eq "");
	next if (/acedb\>/);
	last if (/\/\//);

	@f = split /\t/;
	next unless ( $f[1] or $f[2] or $f[3] ); # this will be the case if update is done when no data available
	
	if    ($f[1] eq "Predicted")           {$cds2status{$f[0]} = $f[1];}
	elsif ($f[2] eq "Partially_confirmed") {$cds2status{$f[0]} = $f[2];}
	elsif ($f[3] eq "Confirmed")           {$cds2status{$f[0]} = $f[3];}
    }
    close TACE;
    
    # now dump data to file

    open (CDS, ">$data_dir/cds2status.dat") or die "Can't open file: $data_dir/cds2status.dat";
    print CDS Data::Dumper->Dump([\%cds2status]);
    close CDS;
}

########################################################################################################
########################################################################################################

sub write_cds2cgc {

  # connect to AceDB using TableMaker,
  my $command="Table-maker -p $wquery_dir/$Table_defs{'cds2cgc'}\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  my %cds2cgc;

  while (<TACE>) {
      chomp;
      s/\"//g;
      next if ($_ eq "");
      next if (/acedb\>/);

      if (/^\S+\s+(\S+)\s+(\S+)/) {
	  $cds2cgc{$2} = $1;
      }
  }

  # now dump data to file
  open (cds2cgc, ">$data_dir/cds2cgc.dat") or die "$data_dir/cds2cgc.dat";

  print cds2cgc Data::Dumper->Dump([\%cds2cgc]);
  
  close cds2cgc;
  }


########################################################################################################

sub write_Feature  {   

  $log->write_to("Updating write_Feature\n");
  my %est2feature;
  my $EST;
  my $feature;

  # connect to AceDB using TableMaker,
  my $command="Table-maker -p $wquery_dir/$Table_defs{'est2feature'}\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  while (<TACE>) {
      chomp;
      s/\"//g;
      next if ($_ eq "");
      next if (/acedb\>/);
      last if (/\/\//);
      if (/^(\S+)\s+(\S+)/) {
	  $est2feature{$1} = $2;
      }
  }
  close TACE;
  
  # now dump data to file

  open (CDS, ">$data_dir/est2feature.dat") or die "Can't open file: $data_dir/est2feature.dat";
  print CDS Data::Dumper->Dump([\%est2feature]);
  close CDS;
}


########################################################################################################

sub write_EST  {   

  $log->write_to("Updating write_EST\n");
  my %NDBaccession2est;
  my %estorientation;
  my @f;

  # connect to AceDB using TableMaker,
  my $command="Table-maker -p $wquery_dir/$Table_defs{'estdata'}\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  while (<TACE>) {
    chomp;
    next if ($_ eq "");    # shortcut at empty lines
    next if (/acedb\>/);
    last if (/\/\//);      # end when you get to the end
    s/\"//g;               # remove speech marks
    @f = split (/\t/);
    
    # NDB_accession to WormBase name
    $NDBaccession2est{$f[1]} = $f[0];

    # EST orientation
    $estorientation{$f[0]} = 5 if ($f[2]);
    $estorientation{$f[0]} = 3 if ($f[3]);
  }
  close TACE;
  
  # now dump data to file
  open (EST, ">$data_dir/NDBaccession2est.dat") or die "Can't open file: $data_dir/NDBaccession2est.dat";
  print EST Data::Dumper->Dump([\%NDBaccession2est]);
  close EST;

  open (ESTorient, ">$data_dir/estorientation.dat") or die "Can't open file: $data_dir/estorientation.dat";
  print ESTorient Data::Dumper->Dump([\%estorientation]);
  close ESTorient;

}

####################################################################################

sub write_clones2seq  {
 
  $log->write_to("Updating clone2seq\n");
  my %clone2seq;
  my $seq; my $newname; my $seqname;
 
  # connect to AceDB using tace
  my $command = "query find Genome_sequence\nDNA\nquit\n";
 
  open (TACE, "echo '$command' | $tace $ace_dir | ");
  while (<TACE>) {
    chomp;
    next if ($_ eq "");
    next if (/acedb\>/);
    next if (/\/\//);
    
    if (/^>(\S+)/) {
      
      $newname = $1;
      unless (defined $seqname) {
	$seqname = $newname;
	next;
      }
      
      # assign sequence
      $clone2seq{$seqname} = $seq;
      
      # reset vars
      $seqname = $newname;
      $seq = "";
      next;
      
    }
    $seq .= $_;
  }
  close TACE;
  $clone2seq{$seqname} = $seq;
  
  #now dump data to file
  open (CDS, ">$data_dir/clone2sequence.dat") or die "Can't open file: $data_dir/clone2sequence.dat";
  print CDS Data::Dumper->Dump([\%clone2seq]);
  close CDS;
}

####################################################################################

sub write_genes2lab  {  
 
  $log->write_to("Updating genes2lab\n");
  my %genes2lab;

  # connect to AceDB using TableMaker, 3 separate queries for CDS, Transcript and Pseudogene
  my @commands;
  $commands[0] = "Table-maker -p $wquery_dir/$Table_defs{'cds2lab'}\nquit\n";
  $commands[1] = "Table-maker -p $wquery_dir/$Table_defs{'pseudogene2lab'}\nquit\n";
  $commands[2] = "Table-maker -p $wquery_dir/$Table_defs{'RNAgene2lab'}\nquit\n";

  # add to hash for each of the three queries
  foreach my $command (@commands){
    open (TACE, "echo '$command' | $tace $ace_dir |");
    while (<TACE>) {
      chomp;
      next if ($_ eq "");
      next if (/acedb\>/);
      last if (/\/\//);
      if (/^\"(\S+)\"\s+\"(\S+)\"/) {
	$genes2lab{$1} = $2;
      }
      
    }
    close TACE;
  }


  #now dump data to file
  open (GENES2LAB, ">$data_dir/worm_gene2lab.dat") or die "Can't open file: $data_dir/worm_gene2lab.dat";
  print GENES2LAB Data::Dumper->Dump([\%genes2lab]);
  close GENES2LAB;
}

####################################################################################



sub write_worm_gene2cgc  {   

  $log->write_to("Updating gene2cgc\n");
  my %worm_gene2cgc;
  my %cgc_name2gene;

  # connect to AceDB using TableMaker, but use /wormsrv2/geneace for Table-maker definition
  my $command="Table-maker -p $wquery_dir/$Table_defs{'wormgene2cgc'}\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  while (<TACE>) {
    chomp;
    s/\"//g;
    next if ($_ eq "");
    next if (/acedb\>/);
    last if (/\/\//);

    # split the line into various fields
    my ($gene,$cds,$transcript,$pseudogene,$cgc_name) = split(/\t/, $_) ;

    # populate first hash with CGC name as key
    $cgc_name2gene{$cgc_name} = $gene;

    # add to hash, molecular name is key.  Value consists of three parts:
    # 1) class of molecular name, 2) CGC name 3) ?Gene name

    if ($cds) {
	$worm_gene2cgc{$cds} = "CDS $cgc_name $gene";
    }
    elsif ($transcript) {
	$worm_gene2cgc{$transcript} = "Transcript $cgc_name $gene";
    }
    elsif ($pseudogene) {
	$worm_gene2cgc{$pseudogene} = "Pseudogene $cgc_name $gene";
    }
  }
  close TACE;
  
  # now dump data to file

  open (CGC, ">$data_dir/worm_gene2cgc_name.dat") or die "Can't open file: $data_dir/worm_gene2cgc_name.dat";
  print CGC Data::Dumper->Dump([\%worm_gene2cgc]);
  close CGC;

  #now dump data to file
  open (CGC, ">$data_dir/cgc_name2gene.dat") or die "Can't open file: $data_dir/cgc_name2gene.dat";
  print CGC Data::Dumper->Dump([\%cgc_name2gene]);
  close CGC;

}


sub write_worm_gene2geneID  {   

  $log->write_to("Updating gene2geneID\n");
  my %worm_gene2geneID;

  # connect to AceDB using TableMaker, but use /wormsrv2/geneace for Table-maker definition
  my $command="Table-maker -p $wquery_dir/$Table_defs{'wormgene2geneid'}\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  while (<TACE>) {
      chomp;
      s/\"//g;
      next if ($_ eq "");
      next if (/acedb\>/);
      last if (/\/\//);
      
      # split the line into various fields
      my ($gene,$cds,$transcript,$pseudogene) = split(/\t/, $_) ;
      
      # add to hash. CDS, Pseudogene, or Transcript name is key, gene ID is value
      
      if ($cds) {
	  $worm_gene2geneID{$cds} = "$gene";
      }
      elsif ($transcript) {
	  $worm_gene2geneID{$transcript} = "$gene";
      }
      elsif ($pseudogene) {
	  $worm_gene2geneID{$pseudogene} = "$gene";
      }
  }
  close TACE;
  
  # now dump data to file

  open (OUT, ">$data_dir/worm_gene2geneID_name.dat") or die "Can't open file: $data_dir/worm_gene2geneID_name.dat";
  print OUT Data::Dumper->Dump([\%worm_gene2geneID]);
  close OUT;

}

sub write_Gene_id{

  $log->write_to("Updating Gene_id\n");
  my %CDS2gene;
  my %gene2CDS;
  
  my $query = "select CDS, CDS->Gene from CDS in class CDS where CDS->method = \"curated\"";
  open (TACE, "echo 'select CDS, CDS->Gene from CDS in class CDS where CDS->method = \"curated\"' | $tace $ace_dir |") or die "cant open tace connection :Gene_id\t$!\n";
  while( <TACE> ) {
    next if ($_ eq "");
    next if (/acedb\>/);
    last if (/\/\//);
    
    s/\"//g;
    my ($cds,$gene) = split;
    if( defined $gene and $gene ne "NULL" ) {
      $cds =~ s/CDS://;
      $gene =~ s/Gene://;
      
      $CDS2gene{$cds} = $gene;
      push(@{$gene2CDS{$gene}},$cds);
    }
  }
  
  # now dump data to file

  open (C2G, ">$data_dir/cds2wbgene_id.dat") or die "cant write $data_dir/cds2wbgene_id.dat :$!";
  open (G2C, ">$data_dir/wbgene_id2cds.dat") or die "cant write $data_dir/wbgene_id2cds.dat :$! ";
  
  print C2G Data::Dumper->Dump([\%CDS2gene]);
  print G2C Data::Dumper->Dump([\%gene2CDS]);
  
  close C2G;
  close G2C;
}

###########################################################################################################



sub write_worm_gene2class  {   

  $log->write_to("Updating worm_gene2class\n");
  my %worm_gene2class;

  # loop through three subclasses using AQL to get list of objects in that class
  foreach my $class ("elegans_CDS","elegans_pseudogenes", "elegans_RNA_genes"){
    my $command = "select i from i in class $class";
    open (TACE, "echo '$command' | $tace $ace_dir |");
    while (<TACE>) {
      chomp;
      next if ($_ eq "");
      next if (/acedb\>/);
      last if (/\/\//);
      
      # get rid of quote marks
      s/\"//g;
      
      # populate hash with CDS/Pseudogene/Transcript name as key
      $worm_gene2class{$_} = "CDS"        if ($class eq "elegans_CDS");
      $worm_gene2class{$_} = "Pseudogene" if ($class eq "elegans_pseudogenes");
      $worm_gene2class{$_} = "Transcript" if ($class eq "elegans_RNA_genes");
    }
    close TACE;
  }

  # now dump data to file

  open (DAT, ">$data_dir/worm_gene2class.dat") or die "Can't open file: $data_dir/worm_gene2class.dat";
  print DAT Data::Dumper->Dump([\%worm_gene2class]);
  close DAT;

}

#################################################################################################################


__END__

=pod

=head1 NAME - update_Common_data.pl

=head2 DESCRIPTION

The update_Common_data.pl updates the common data sets retrieved by the Fetch_data routine in Wormbase.pm gives 
quick easy acces to a variety of data frequently used in Wormbase scripts. It comprises of one part.

Part One generates the data and writes it to a file using the Data::Dumper module.

This module updates the following data sets:

=over 4

=item *

%accession2clone

=item *

%clone2accession

=item *

%clonesize - list of DNA lengths for each genome sequence

=item *

%cds2status - list of CDSs with confirmed status as value (Confirmed, Predicted or Partially_confirmed)

=item *

%wormpep2cds  -  bear in mind that a peptide may have multiple genes. If so the genes are concatenated, separated by a space.

=item *

%cds2wormpep

=back

=item *

%NDBaccession2est - NDB accession is key, WormBase EST name is value

=back

=item *

%estorientation - WormBase EST name is key, EST orientation is value (EST_5 = 5, EST_3 = 3)

=back

=item *

%clone2seq - genomic clone name is key, DNA sequence is value

=back

=item *

%cds2protein_id - connections between CDS name and protein_ID, also outputs the reverse

=back

=item *

%worm_gene2lab - connections between C. elegans CDS, Transcript, or Pseudogene and corresponding lab designation (RW, HX, DRW)

=back

=itme *
%worm_gene2class - connections between worm_genes composite class and the actual class of each object (e.g. AH6.1 -> CDS)

=back

=item *

%worm_gene2cgc - connections between C. elegans CDS, Transcript, or Pseudogene and CGC name...the value actually stores three fields
(separated by spaces). 1) class name of the key (CDS, Transcript, or Pseudogene) 2) CGC name 3) Gene ID

=back

=item *

%worm_gene2geneID - simplification of above: CDS, Transcript, or Pseudogene is key, gene ID is value

=item *

%cgc_name2gene - generated at same time as above.  Simple list of CGC approved names as keys,
and Gene ID as value

=back




=over 4

=head2 EXAMPLE USAGE

=over 4

=item UPDATING THE DATA

=over 4

We need to be very careful about updating this data.  Depending on whether it is being updated during the build or 
otherwise we need to use autoace or current_DB. 

=item 

update_Common_data.pl -build -all

However if the build is underway and you want to write out the CDS 2 wormpep info but that is not yet in the 
building database you need to use current_DB. so don't include -build ie 

=item

update_Common_data.pl -cds2wormpep


At the end of the build ( in finish_build.pl) all data will be refreshed to be in synch with the current
release ie

=item

update_Common_data.pl -build -all.

There shouldn't really be any need to alter this apart from actually during the next build itself.  
Scripts that generate the data that is included in Common_data should call the updating routines 
themselves - so you wont have to ! 
