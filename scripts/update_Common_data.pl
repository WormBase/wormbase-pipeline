#!/usr/bin/env perl
#
# update_Common_data.pl
# 
# by Anthony Rogers et al
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2015-03-23 10:29:14 $

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
use Modules::WormSlurm;

##############################
# command-line options       #
##############################


my $test;              # test mode, uses ~wormpub/TEST_BUILD
my $all;               # performs all of the below options:
my $debug;
my $species;
my $verbose;
my $database;
     
my $clone2accession;   # Hash: %clone2accession     Key: Genomic_canonical                 Value: GenBank/EMBL accession
                       # Hash: %accession2clone     Key: GenBank/EMBL accession            Value: Genomic_canonical
my $clone2size;        # Hash: %clonesize           Key: Genomic_canonical                 Value: DNA length
my $cds2wormpep;       # Hash: %cds2wormpep         Key: CDS name                          Value: Wormpep ID
                       # Hash: %wormpep2cds         Key: Wormpep ID                        Value: CDS name
my $cds2protein_id;    # Hash: %cds2protein_id      Key: CDS name                          Value: Protein_ID
                       # Hash: %protein_id2cds      Key: Protein_ID                        Value: CDS name
my $cds2cgc;           # Hash: %cds2cgc             Key: CDS name                          Value: CGC name
my $rna2cgc;           # Hash: %rna2cgc             Key: transcript name                   Value: CGC name
my $rna2briefID;        # Hash: %rna2briefID         Key: transcript name                   Value: Brief_ID
my $pseudo2cgc;        # Hash: %pseudo2cgc          Key: Pseudogene name                   Value: CGC name
my $cds2status;        # Hash: %cds2status          Key: CDS name                          Value: Prediction status 
my $clone2seq;         # Hash: %clone2seq           Key: Genomic_canonical                 Value: DNA sequence (lower case)
my $clone2sv;          # Hash: %clone2sv            Key: Genomic_canonical                 Value: Sequence version (integer)
my $clone2type;        # Hash: %clone2type          Key: Genomic_canonical                 Value: Type information (Cosmid, Fosmid, YAC, Plasmid)
my $clone2centre;      # Hash: %clone2type          Key: Genomic_canonical                 Value: From_laboratory (HX, RW, DRW)
my $clone2dbid;        # Hash: %clone2dbid          Key: Genomic_canonical                 Value: EMBL dbID (_*)
my $genes2lab;         # Hash: %worm_gene2lab       Key: Gene (CDS|Transcript|Pseudogene)  Value: From_laboratory (HX, RW, DRW)
my $worm_gene2cgc;     # Hash: %worm_gene2cgc_name  Key: CGC name                          Value: Gene ID, plus molecular name (e.g. AH6.1), also a hash of cgc_name2gene
my $worm_gene2geneID;  # Hash: %worm_gene2geneID    Key: Gene (CDS|Transcript|Pseudogene)  Value: Gene ID
my $worm_gene2class;   # Hash: %worm_gene2class     Key: CDS/Transcript/Pseudogene name    Value: 'CDS', 'Transcript', or 'Pseudogene'
my $estdata;           # Hash: %NDBaccession2est    Key: GenBank/EMBL accession            Value: EST name (WormBase)  
                       # Hash: %estorientation      Key: EST name (WormBase)               Value: EST_5 = 5, EST_3 = 3
my $est2feature;       # Hash: %est2feature         Key: EST name (WormBase)               Value: Feature name (WormBase)
my $cds2gene_id;       # Hash: %cds2gene_id         Key: CDS name                          Value: WBGene_id


my $store;
GetOptions (
	    "clone2acc"          => \$clone2accession,
	    "clone2size"         => \$clone2size,
	    "cds2wormpep"        => \$cds2wormpep,
	    "cds2pid"            => \$cds2protein_id,
	    "cds2status"         => \$cds2status,
	    "clone2seq:s"        => \$clone2seq, # 'all' or species the sequence comes from (i.e. restrict to this Species tag)
	    "clone2sv"           => \$clone2sv,
	    "clone2centre"       => \$clone2centre,
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
            "clone2dbid"         => \$clone2dbid,
	    "cds2cgc"            => \$cds2cgc,
	    "rna2cgc"            => \$rna2cgc,
	    "rna2briefID"        => \$rna2briefID,
	    "pseudo2cgc"         => \$pseudo2cgc,
	    "species:s"		 => \$species,
	    "verbose"            => \$verbose,
            "database:s"         => \$database,
	   );

my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			     -organism => $species
			   );
}

my $log = Log_files->make_build_log( $wormbase );

$species = $wormbase->species;

my %Table_defs = (
		  'cds2protein'      => 'CommonData:CDS_proteinID.def',
		  'clone2sv'         => 'CommonData:Clone_SequenceVersion.def',
		  'clone2accession'  => 'CommonData:Clone_Accession.def', 
		  'clone2size'       => "CommonData:Clone_Size_${\$wormbase->species}.def",
                  'clone2type'       => 'CommonData:Clone_Type.def',
		  'clone2centre'     => 'CommonData:Clone_Lab.def',
                  'clone2dbid'       => 'CommonData:clone2dbid.def',
		  'cds2status'       => 'CommonData:CDS_Status.def',
                  'cds2cgc'          => 'CommonData:CDS_CGCname.def',
		  'rna2cgc'          => 'CommonData:RNA_CGCname.def',
		  'rna2BriefID'      => 'CommonData:RNA_BriefID.def',
		  'pseudo2cgc'       => 'CommonData:Pseudogene_CGCname.def',
		  'est2feature'      => 'CommonData:EST_Feature.def',
		  'estdata'          => 'CommonData:EST_data.def',
		  'cds2lab'          => 'CommonData:CDS_Lab.def',
		  'pseudogene2lab'   => 'CommonData:Pseudogene_Lab.def',
		  'RNAgene2lab'      => 'CommonData:RNAgene_Lab.def',
		  'wormgene2cgc'     => 'CommonData:WormGene_CGCname.def',
		  'wormgene2geneid'  => 'CommonData:WormGene_GeneID.def',
		  'cds2wormpep'      => 'CommonData:CDS2wormpep'
		  );



##########################################

# Set up database paths                  #
##########################################

my $data_dir   = (defined $database) ? "$database/COMMON_DATA" : $wormbase->common_data;
my $wquery_dir = (defined $database) ? "$database/wquery" : $wormbase->autoace."/wquery";
my $ace_dir    = (defined $database) ? $database : $wormbase->autoace;


$log->write_to("Updating COMMON_DATA in $data_dir\n");
##############################
# ACEDB executables          #
##############################

our $tace= $wormbase->tace;

# run '-all' Common_data dumps under Slurm
if ($all) { 

  my @all_args = qw( clone2acc clone2size cds2wormpep cds2pid
	    cds2status clone2sv clone2centre genes2lab
	    worm_gene2cgc worm_gene2geneID worm_gene2class est
	    est2feature gene_id clone2type cds2cgc rna2cgc pseudo2cgc clone2dbid rna2briefID);

  # Deal with clone2seq up front, because that routine submits its own Slurm jobs
  &write_clones2seq();

  $wormbase->check_slurm;
  my $store_file = $wormbase->build_store; # make the store file to use in all commands
  my $scratch_dir = $wormbase->logs;
  my $job_name = "worm_".$wormbase->species."_commondata";

  my %slurm_jobs;
  foreach my $arg (@all_args) {
    my $err = "$scratch_dir/update_Common_data.pl.slurm.${arg}.err";
    my $out = "$scratch_dir/update_Common_data.pl.slurm.${arg}.out";
    
    my $cmd = "update_Common_data.pl -${arg}";
    $cmd .= " -database $database" if defined $database;
    $cmd = $wormbase->build_cmd_line($cmd, $store_file);
    $log->write_to("Submitting Slurm job: $cmd\n");
    my $job_id = WormSlurm::submit_job_with_name($cmd, 'production', '6g', '3:00:00', $err, $out, $job_name);
    $slurm_jobs{$job_id} = $cmd;
  }

  WormSlurm::wait_for_jobs(keys %slurm_jobs);
  $log->write_to("All Common_data jobs have completed!\n");
  my @problem_cmds;
  for my $job_id (keys %slurm_jobs) {
    if (WormSlurm::get_exit_code($job_id) != 0) {
      $log->error("Job $job_id (" . $slurm_jobs{$job_id} . ") exited non zero\n");
      push @problem_cmds, $slurm_jobs{$job_id};
    }
  }
} else {
  # run the various options depending on command line arguments
  &write_cds2protein_id   if ($cds2protein_id);
  &write_clone2accession  if ($clone2accession);
  &write_clonesize        if ($clone2size);
  &write_cds2wormpep      if ($cds2wormpep);
  &write_cds2status       if ($cds2status);
  &write_cds2cgc          if ($cds2cgc);
  &write_rna2cgc          if ($rna2cgc);
  &write_pseudo2cgc       if ($pseudo2cgc);
  &write_clones2seq       if ($clone2seq);
  &write_clones2sv        if ($clone2sv);
  &write_clone2type       if ($clone2type);
  &write_clone2centre     if ($clone2centre);
  &write_genes2lab        if ($genes2lab);
  &write_worm_gene2class  if ($worm_gene2class);
  &write_EST              if ($estdata);
  &write_Feature          if ($est2feature);
  &write_Gene_id          if ($cds2gene_id);
  &write_worm_gene2geneID if ($worm_gene2geneID);
  &write_worm_gene2cgc    if ($worm_gene2cgc);
  &write_clone2dbid       if ($clone2dbid);
  &write_rna2BriefID      if ($rna2briefID);
}

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
      next if ($_ eq "");
      next if (/^acedb\>/);
      next if (/^\/\//);

      if (/\"(\S+)\"\s+\"(\S+)\"\s+(\d+)/) {
	  my $protein_id = "$2".".$3";
	  my $gene = $1;
	  $cds2protein_id{"$gene"}       = $protein_id;
	  $protein_id2cds{"$protein_id"} = $gene;
      }
  }
  close TACE;
  
  # now dump data to file
  
  open (G2P, ">$data_dir/cds2protein_id.dat") or $log->log_and_die("cant write $data_dir/cds2protein_id.dat :$!");
  open (P2G, ">$data_dir/protein_id2cds.dat") or $log->log_and_die("cant write $data_dir/protein_id2cds.dat :$! ");
  
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
      next if ($_ eq "");
      next if (/^acedb\>/);
      next if (/^\/\//);
      if (/\"(\S+)\"\s+\"(\S+)\"/) {
	  $clone2accession{$1} = $2;
	  $accession2clone{$2} = $1;
      }
  }
  close TACE;
  
  # now dump data to file
  
  open (C2A, ">$data_dir/clone2accession.dat") or $log->log_and_die("cant write $data_dir/clone2accession.dat :$!");
  open (A2C, ">$data_dir/accession2clone.dat") or $log->log_and_die( "cant write $data_dir/accession2clone.dat :$! ");
  
  print C2A Data::Dumper->Dump([\%clone2accession]);
  print A2C Data::Dumper->Dump([\%accession2clone]);
  
  close C2A;
  close A2C;
  
}

########################################################################################################

sub write_clone2centre{   
  $log->write_to("Updating clone2centre\n");

  my %clone2centre;
    
  # connect to AceDB using TableMaker,
  my $command="Table-maker -p $wquery_dir/$Table_defs{'clone2centre'}\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  while (<TACE>) {
      chomp;
      next if ($_ eq "");
      next if (/^acedb\>/);
      next if (/^\/\//);
      if (/\"(\S+)\"\s+\"(\S+)\"/) {
	  $clone2centre{$1} = $2;
      }
  }
  close TACE;
  
  # now dump data to file
  
  open (C2c, ">$data_dir/clone2centre.dat") or $log->log_and_die("cant write $data_dir/clone2centre.dat :$!");
  print C2c Data::Dumper->Dump([\%clone2centre]);
  close C2c;
  
}

########################################################################################################


sub write_clones2sv  {   
  $log->write_to("Updating clones2sv\n");
    my %clone2sv;
    
    # connect to AceDB using TableMaker,
    my $command="Table-maker -p $wquery_dir/$Table_defs{'clone2sv'}\nquit\n";
    
    open (TACE, "echo '$command' | $tace $ace_dir |");
    while (<TACE>) {
	chomp;
	next if ($_ eq "");
	next if (/^acedb\>/);
        next if /^\/\//;
	if (/^\"(\S+)\"\s+\"(\S+)\.(\d+)\"/) {
	    $clone2sv{$1} = $3;
	}
    }
    close TACE;
	
    # now dump data to file

    open (C2SV, ">$data_dir/clone2sv.dat") or $log->log_and_die("cant write $data_dir/clone2sv.dat :$!");
    
    print C2SV Data::Dumper->Dump([\%clone2sv]);
    
    close C2SV;
    
}

########################################################################################################

sub write_clone2type  {   
  $log->write_to("Updating clone2type\n");
    my %clone2type;
    
    # connect to AceDB using TableMaker,
    my $command="Table-maker -p $wquery_dir/$Table_defs{'clone2type'}\nquit\n";
    
    open (TACE, "echo '$command' | $tace $ace_dir |");
    while (<TACE>) {
	chomp;
	next if ($_ eq "");
	next if (/acedb\>/);
        next if /^\/\//;

	if (/^\"(\S+)\"\s+\"(\S+)\"/) {
	    $clone2type{$1} = $2;
	}
    }
    close TACE;
	
    # now dump data to file

    open (C2TYPE, ">$data_dir/clone2type.dat") or $log->log_and_die("cant write $data_dir/clone2type.dat :$!");
    
    print C2TYPE Data::Dumper->Dump([\%clone2type]);
    
    close C2TYPE;
    
}


########################################################################################################

sub write_clonesize  {   
  $log->write_to("Updating clone2size\n");

  my %clonesize;
  
  unless ($wormbase->species eq 'elegans') { # needs to be beautified
	  my $command="perl $ENV{CVS_DIR}/get_clone_sizes.pl -species ".$wormbase->species;
	  open(TACE, "$command|");
  }
  else {
	  # connect to AceDB using TableMaker,
	  my $command="Table-maker -p $wquery_dir/$Table_defs{'clone2size'}\nquit\n";
	  open (TACE, "echo '$command' | $tace $ace_dir |");
  }
  while (<TACE>) {
    chomp;
    s/\"//g;
    next if ($_ eq "");
    next if (/acedb\>/);
    next if /^\/\//;
    if (/^(\S+)\s+(\d+)$/) {
      $clonesize{$1} = $2;
    }
  }
  close TACE;

  # now dump data to file

  open (CLONESIZE, ">$data_dir/clonesize.dat") or $log->log_and_die("cant write $data_dir/clonesize.dat :$!");
  
  print CLONESIZE Data::Dumper->Dump([\%clonesize]);

  close CLONESIZE;
  
}


########################################################################################################

sub write_cds2wormpep  {   
  $log->write_to("Updating cds2wormpep\n");

  use File::Temp qw /:POSIX/;
  my $fname = File::Temp::tmpnam('/tmp/');

  $log->write_to("Updating cds2wormpep\n");
  my $WPver = $wormbase->get_wormbase_version;

  # connect to AceDB using TableMaker,
  my $species=$wormbase->full_name;
  my $tablemakerFile="$wquery_dir/$Table_defs{'cds2wormpep'}";

  $wormbase->run_command("perl -pne 's/Caenorhabditis elegans/$species/' $tablemakerFile > $fname",$log);
  my $command="Table-maker -p $fname\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  my %cds2wormpep;
  my %wormpep2cds;
  while (<TACE>) {
    chomp;
    next if ($_ eq "");
    next if (/^acedb\>/);
    next if /^\/\//;
    
    if (/\"(\S+)\"\s+\"(\S+)\"/) {
      my $gene = $1;
      my $pep = $2;
      $pep =~ s/\w+://;
      $cds2wormpep{$gene} = $pep;
      $wormpep2cds{$pep} .= "$gene ";
    }
  }

  unlink $fname;

  #now dump data to file
  open (C2G, ">$data_dir/wormpep2cds.dat") or $log->log_and_die("$data_dir/wormpep2cds.dat");
  open (G2C, ">$data_dir/cds2wormpep.dat") or $log->log_and_die("$data_dir/cds2wormpep.dat");

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
  $log->write_to("Updating cds2status\n");

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
	next if (/^\/\//);

	@f = split /\t/;
	next unless ( $f[1] or $f[2] or $f[3] ); # this will be the case if update is done when no data available
	
	if    ($f[1] eq "Predicted")           {$cds2status{$f[0]} = $f[1];}
	elsif ($f[2] eq "Partially_confirmed") {$cds2status{$f[0]} = $f[2];}
	elsif ($f[3] eq "Confirmed")           {$cds2status{$f[0]} = $f[3];}
    }
    close TACE;
    
    # now dump data to file

    open (CDS, ">$data_dir/cds2status.dat") or $log->log_and_die("Can't open file: $data_dir/cds2status.dat");
    print CDS Data::Dumper->Dump([\%cds2status]);
    close CDS;
}

########################################################################################################
########################################################################################################

sub write_cds2cgc {
  $log->write_to("Updating cds2cgc\n");

  # connect to AceDB using TableMaker,
  my $command="Table-maker -p $wquery_dir/$Table_defs{'cds2cgc'}\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  my %cds2cgc;

  while (<TACE>) {
      chomp;
      next if ($_ eq "");
      next if (/acedb\>/);
      next if /^\/\//;
      
      if (/^\S+\s+\"(\S+)\"\s+\"(\S+)\"/) {
	  $cds2cgc{$2} = $1;
      }
  }

  # now dump data to file
  open (cds2cgc, ">$data_dir/cds2cgc.dat") or $log->log_and_die("$data_dir/cds2cgc.dat");

  print cds2cgc Data::Dumper->Dump([\%cds2cgc]);
  
  close cds2cgc;
  }

########################################################################################################

sub write_rna2cgc {
  $log->write_to("Updating rna2cgc\n");

  # connect to AceDB using TableMaker,
  my $command="Table-maker -p $wquery_dir/$Table_defs{'rna2cgc'}\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  my %rna2cgc;

  while (<TACE>) {
      chomp;
      next if ($_ eq "");
      next if (/acedb\>/);
      next if /^\/\//;

      if (/^\S+\s+\"(\S+)\"\s+\"(\S+)\"/) {
	  $rna2cgc{$2} = $1;
      }
  }

  # now dump data to file
  open (rna2cgc, ">$data_dir/rna2cgc.dat") or $log->log_and_die("$data_dir/rna2cgc.dat");

  print rna2cgc Data::Dumper->Dump([\%rna2cgc]);
  
  close rna2cgc;
}

########################################################################################################

sub write_pseudo2cgc {
  $log->write_to("Updating pseudo2cgc\n");

  # connect to AceDB using TableMaker,
  my $command="Table-maker -p $wquery_dir/$Table_defs{'pseudo2cgc'}\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  my %pseudo2cgc;

  while (<TACE>) {
      chomp;
      next if ($_ eq "");
      next if (/acedb\>/);
      next if /^\/\//;

      if (/^\S+\s+\"(\S+)\"\s+\"(\S+)\"/) {
	  $pseudo2cgc{$2} = $1;
      }
  }

  # now dump data to file
  open (ps2cgc, ">$data_dir/pseudo2cgc.dat") or $log->log_and_die("$data_dir/pseudo2cgc.dat");

  print ps2cgc Data::Dumper->Dump([\%pseudo2cgc]);
  
  close ps2cgc;
}


########################################################################################################

sub write_rna2BriefID {
  $log->write_to("Updating rna2Brief_ID\n");

  # connect to AceDB using TableMaker,
  my $command="Table-maker -p $wquery_dir/$Table_defs{'rna2BriefID'}\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  my %rna2BriefID;

  while (<TACE>) {
      chomp;
      next if ($_ eq "");
      next if (/acedb\>/);
      next if /^\/\//;

      if (/^\"(\S+)\"\s+\"(.+)\"/) {
	  $rna2BriefID{$1} = $2;
      }
  }

  # now dump data to file
  open (rna2BriefID, ">$data_dir/rna2BriefID.dat") or die "$data_dir/rna2BriefID.dat";

  print rna2BriefID Data::Dumper->Dump([\%rna2BriefID]);
  
  close rna2BriefID;
  }


########################################################################################################

sub write_Feature  {   

  $log->write_to("Updating est2feature\n");
  my %est2feature;
  my $EST;
  my $feature;

  # connect to AceDB using TableMaker,
  my $command="Table-maker -p $wquery_dir/$Table_defs{'est2feature'}\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  while (<TACE>) {
      chomp;
      next if ($_ eq "");
      next if (/acedb\>/);
      next if (/^\/\//);
      if (/^\"(\S+)\"\s+\"(\S+)\"/) {
	push @{$est2feature{$1}}, $2; # make a list of features defined by this EST
      }
  }
  close TACE;
  
  # now dump data to file

  open (CDS, ">$data_dir/est2feature.dat") or $log->log_and_die("Can't open file: $data_dir/est2feature.dat");
  print CDS Data::Dumper->Dump([\%est2feature]);
  close CDS;
}


########################################################################################################

sub write_EST  {   

  $log->write_to("Updating estdata\n");
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
    next if (/^\/\//);      # end when you get to the end
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
  open (EST, ">$data_dir/NDBaccession2est.dat") or $log->log_and_die("Can't open file: $data_dir/NDBaccession2est.dat");
  print EST Data::Dumper->Dump([\%NDBaccession2est]);
  close EST;

  open (ESTorient, ">$data_dir/estorientation.dat") or $log->log_and_die("Can't open file: $data_dir/estorientation.dat");
  print ESTorient Data::Dumper->Dump([\%estorientation]);
  close ESTorient;

}

####################################################################################
# The -clone2seq parameter holds the name of the species to dump
# sequences for or it takes the value 'all' to dump out genomic
# sequences from all species. If the value is 'all' then Slurm jobs are
# fired off, each doing a different species. If the value is the short
# name of a species, then that species is dumped to a file.

sub write_clones2seq  {
 
  $log->write_to("Updating clone2seq\n");

  # populate the fullnames and sort out whether the species have a supercontig file that we can read in or not 
  my %full_names;		# key=species value=full species name
  my %supercontigs;		# key=species value=supercontigs file (if it exists)
  my %accessors = ($wormbase->species_accessors);
  $accessors{elegans} = $wormbase;
  foreach my $spDB (values %accessors) {
    my $species = $spDB->species;
    $full_names{$species} = $spDB->full_name();
    if ($spDB->assembly_type eq 'contig') {
      my $chromdir = $spDB->sequences;
      $supercontigs{$species} = $spDB->genome_seq;
    }
  }

  if ($clone2seq eq 'all' or not $clone2seq) { 
  # set up Slurm jobs each doing a different species, otherwise TACE
  # can crash or it takes a VERY long time
    $wormbase->check_slurm;
    my $store_file = $wormbase->build_store; # make the store file to use in all commands
    my $scratch_dir = $wormbase->logs;
    my $job_name = "worm_".$wormbase->species."_commondata";

    my %slurm_jobs;
    my $c2sname = $clone2seq ? $clone2seq : 'all';
    foreach my $this_species (keys %full_names) {
      my $err = "$scratch_dir/update_Common_data.pl.slurm.clone2seq.$c2sname.err";
      my $out = "$scratch_dir/update_Common_data.pl.slurm.clone2seq.$c2sname.out";
      my $cmd = "update_Common_data.pl -clone2seq $this_species";
      $cmd = $wormbase->build_cmd_line($cmd, $store_file);
      $log->write_to("Submitting Slurm job: $cmd\n");
      my $job_id = WormSlurm::submit_job_with_name($cmd, 'production', '8000m', '02:00:00', $err, $out, $job_name);
      $slurm_jobs{$job_id} = $cmd;
    }

    WormSlurm::wait_for_jobs(keys %slurm_jobs);
    $log->write_to("All Common_data jobs have completed!\n");
    for my $job_id (keys %slurm_jobs) {
      if (WormSlurm::get_exit_code($job_id) != 0) {
	$log->write_to("Job $job_id (" . $slurm_jobs{$job_id} . ") exited non zero\n");
      }
    }
    

  } else {			# dump files for a single species
    $log->write_to("Dumping clone2seq for $clone2seq\n");
    my %clone2seq;
    my $seq; 
    my $newname; 
    my $seqname;
 
    # do we have a file of sequences ready to read?
    if (exists $supercontigs{$clone2seq}) {
      my $supercontig_file = $supercontigs{$clone2seq};
      $log->write_to("Reading data from $supercontig_file\n");
      open (CONTIGS, "< $supercontig_file") || $log->log_and_die("Can't open file $supercontig_file : $!\n");
      while (<CONTIGS>) {
	next if ($_ eq "");
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
      close (CONTIGS);
      $clone2seq{$seqname} = $seq;

    } else {
      # connect to AceDB using tace
      $log->write_to("Reading data from tace\n");
      my $full_name = $full_names{$clone2seq};
      my $command = "query find Sequence where Genomic_canonical AND Species = \"$full_name\"\nDNA\nquit\n";
 
      open (TACE, "echo '$command' | $tace $ace_dir | ");
      while (<TACE>) {
	chomp;
	next if ($_ eq "");
	next if (/acedb\>/);
	next if (/^\/\//);
	
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
    }

    #now dump data to file for this species
    open (CDS, ">$data_dir/clone2sequence_$clone2seq.dat") or $log->log_and_die("Can't open file: $data_dir/clone2sequence_$clone2seq.dat");
    print CDS Data::Dumper->Dump([\%clone2seq]);
    close CDS;
  }
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
      next if (/^\/\//);
      if (/^\"(\S+)\"\s+\"(\S+)\"/) {
	$genes2lab{$1} = $2;
      }
      
    }
    close TACE;
  }


  #now dump data to file
  open (GENES2LAB, ">$data_dir/worm_gene2lab.dat") or $log->log_and_die("Can't open file: $data_dir/worm_gene2lab.dat");
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
    next if (/^\/\//);

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

  open (CGC, ">$data_dir/worm_gene2cgc_name.dat") or $log->log_and_die("Can't open file: $data_dir/worm_gene2cgc_name.dat");
  print CGC Data::Dumper->Dump([\%worm_gene2cgc]);
  close CGC;

  #now dump data to file
  open (CGC, ">$data_dir/cgc_name2gene.dat") or $log->log_and_die("Can't open file: $data_dir/cgc_name2gene.dat");
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
      next if (/^\/\//);
      
      # split the line into various fields
      my ($gene,$cds,$transcript,$pseudogene) = split(/\t/, $_) ;
      
      # add to hash. CDS, Pseudogene, or Transcript name is key, gene ID is value
      
      if ($cds) {
	  		$worm_gene2geneID{$cds} = "$gene";
      }
      if ($transcript) {
	  		$worm_gene2geneID{$transcript} = "$gene";
      }
      if ($pseudogene) {
	  		$worm_gene2geneID{$pseudogene} = "$gene";
      }
  }
  close TACE;
  
  # now dump data to file

  open (OUT, ">$data_dir/worm_gene2geneID_name.dat") or $log->log_and_die("Can't open file: $data_dir/worm_gene2geneID_name.dat");
  print OUT Data::Dumper->Dump([\%worm_gene2geneID]);
  close OUT;

}

sub write_Gene_id{

  $log->write_to("Updating Gene_id\n");
  my %CDS2gene;
  my %gene2CDS;
  
  my $query = "select CDS, CDS->Gene from CDS in class CDS where CDS->method = \"curated\"";
  open (TACE, "echo 'select CDS, CDS->Gene from CDS in class CDS where CDS->method = \"curated\"' | $tace $ace_dir |") 
      or $log->log_and_die("cant open tace connection :Gene_id\t$!\n");
  while( <TACE> ) {
    next if ($_ eq "");
    next if (/acedb\>/);
    next if (/^\/\//);
    
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

  open (C2G, ">$data_dir/cds2wbgene_id.dat") or $log->log_and_die("cant write $data_dir/cds2wbgene_id.dat :$!");
  open (G2C, ">$data_dir/wbgene_id2cds.dat") or $log->log_and_die("cant write $data_dir/wbgene_id2cds.dat :$! ");
  
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
      next if (/^\/\//);
      
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

  open (DAT, ">$data_dir/worm_gene2class.dat") or $log->log_and_die("Can't open file: $data_dir/worm_gene2class.dat");
  print DAT Data::Dumper->Dump([\%worm_gene2class]);
  close DAT;

}

sub write_clone2dbid {
  
  $log->write_to("Updating clone2dbid\n");
  my %clone2dbid;
    
  # connect to AceDB using TableMaker,
  my $command="Table-maker -p $wquery_dir/$Table_defs{'clone2dbid'}\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  while (<TACE>) {
    chomp;
    s/\"//g;
    next if ($_ eq "");
    next if (/acedb\>/);
    next if (/^\/\//);
    if (/(\S+)\s+(\S+)/) {
      $clone2dbid{$1} = $2;
    } elsif (/^(\S+)/) {
      $clone2dbid{$1} = "";
    }
  }
  close TACE;
  
  # now dump data to file
  my $datafile = "$data_dir/clone2dbid.dat";

  open (my $c2dbid, ">$datafile") or $log->log_and_die("Cant write $datafile $!\n");
  print $c2dbid Data::Dumper->Dump([\%clone2dbid]);
  close($c2dbid);
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

%clone2seq all | species_name - This writes the file clone2sequence_<species>.dat for the specified species or for elegans and all Tier II species. Genomic clone name is key, DNA sequence is value

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




