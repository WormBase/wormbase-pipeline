#!/usr/local/bin/perl5.8.0 -w
#
# update_Common_data.pl
# 
# by Anthony Rogers
#
# Last updated by: $Author: dl1 $
# Last updated on: $Date: 2004-04-28 10:08:51 $

#################################################################################
# Initialise variables                                                          #
#################################################################################

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Data::Dumper;
use Getopt::Long;

##############################
# command-line options       #
##############################


my $build;          # for when you want to query autoace, i.e. you are building.  Otherwise defaults to current_DB
my $test;           # test mode, uses ~wormpub/TEST_BUILD
my $all;            # performs all of the below options:
     
my $clone2accession;# Hash: %clone2accession     Key: Genomic_canonical                 Value: GenBank/EMBL accession
                    # Hash: %accession2clone     Key: GenBank/EMBL accession            Value: Genomic_canonical
my $cds2wormpep;    # Hash: %cds2wormpep         Key: CDS name                          Value: Wormpep ID
                    # Hash: %wormpep2cds         Key: Wormpep ID                        Value: CDS name
my $cds2protein_id; # Hash: %cds2protein_id      Key: CDS name                          Value: Protein_ID
                    # Hash: %protein_id2cds      Key: Protein_ID                        Value: CDS name
my $CDS_list;       # Hash: %CDSlist             Key: CDS name                          Value: Confirmed status (confirmed = 1, not confirmed = 0) 
my $clone2seq;      # Hash: %clone2seq           Key: Genomic_canbonical                Value: DNA sequence (lower case)
my $genes2lab;      # Hash: %worm_gene2lab       Key: Gene (CDS|Transcript|Pseudogene)  Value: From_laboratory (HX, RW, DRW)
my $worm_gene2cgc;  # Hash: %worm_gene2cgc_name  Key: CGC name                          Value: Gene ID, plus molecular name (e.g. AH6.1), also a hash of cgc_name2gene
my $estdata;        # Hash: %NDBaccession2est    Key: EST name (WormBase)               Value: GenBank/EMBL accession
                    # Hash: %estorientation      Key: EST name (WormBase)               Value: EST_5 = 5, EST_3 = 3


GetOptions("build"         => \$build,
	   "clone2acc"     => \$clone2accession,
	   "cds2wormpep"   => \$cds2wormpep,
	   "cds2pid"       => \$cds2protein_id,
	   "CDS_list"      => \$CDS_list,
	   "clone2seq"     => \$clone2seq,
	   "genes2lab"     => \$genes2lab,
	   "worm_gene2cgc" => \$worm_gene2cgc,
	   "est"           => \$estdata,
	   "all"           => \$all,
	   "test"          => \$test
	   );

##########################################
# Set up database paths                  #
##########################################

# Set up top level base directory which is different if in test mode
# Make all other directories relative to this
my $basedir   = "/wormsrv2";
$basedir      = glob("~wormpub")."/TEST_BUILD" if ($test); 

my $data_dir   = "$basedir/autoace/COMMON_DATA";
my $wquery_dir = "$basedir/autoace/wquery";
my $ace_dir = "/nfs/disk100/wormpub/DATABASES/current_DB";


##############################
# ACEDB executables          #
##############################

our $tace = &tace;

# use autoace if -build specified, else use current_DB
if($build) {
  $ace_dir = "$basedir/autoace";
  print "during build so using $ace_dir - ensure that the data you are updating is actually in the database.\n";
}
else {
  print "- NOT as part of build so using $ace_dir. If this is part of the build data MAY be stale\n";
}


# run the various options depending on command line arguments
&write_cds2protein_id  if ($cds2protein_id || $all );
&write_clone2accession if ($clone2accession || $all );
&write_cds2wormpep     if ($cds2wormpep || $all );
&write_CDSlist         if ($CDS_list || $all );
&write_clones2seq      if ($clone2seq || $all);
&write_genes2lab       if ($genes2lab || $all);
&write_worm_gene2cgc   if ($worm_gene2cgc);
&write_EST             if ($estdata);

# hasta luego

exit(0);

#######################################################################
# Data writing routines - actually create and dump the data           #
#######################################################################

sub write_cds2protein_id {

  my %cds2protein_id;
  my %protein_id2cds;
    
  ####################################################################
  # connect to AceDB using TableMaker,
  # populating %accession2name (maps embl accession to contig name)
  ####################################################################

  my $command="Table-maker -p $wquery_dir/gene2pid.def\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  while (<TACE>) {
    #gene pid version
    chomp;
    if (/\"(\S+)\"\s+\"(\S+)\"\s+(\d)/) {
      my $protein_id = "$2".".$3";
      my $gene = $1;
      $cds2protein_id{"$gene"} = $protein_id;
      $protein_id2cds{"$protein_id"} = $gene;
    }
  }
  close TACE;
  
  #now dump data to file
  open (G2P, ">$data_dir/cds2protein_id.dat") or die "cant write $data_dir/cds2protein_id.dat :$!";
  open (P2G, ">$data_dir/protein_id2cds.dat") or die "cant write $data_dir/protein_id2cds.dat :$! ";
  
  print G2P Data::Dumper->Dump([\%cds2protein_id]);
  print P2G Data::Dumper->Dump([\%protein_id2cds]);
  
  close G2P;
  close P2G;
    
}


########################################################################################################

sub write_clone2accession  {   

  my %clone2accession;
  my %accession2clone;
    
  # connect to AceDB using TableMaker,
  my $command="Table-maker -p $wquery_dir/accession2clone.def\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  while (<TACE>) {
    chomp;
    if (/\"(\S+)\"\s+\"(\S+)\"\s*\d*/) {
      $clone2accession{$1} = $2;
      $accession2clone{$2} = $1;
    }
    }
  close TACE;
  
  #now dump data to file
  open (C2A, ">$data_dir/clone2accession.dat") or die "cant write $data_dir/clone2accession.dat :$!";
  open (A2C, ">$data_dir/accession2clone.dat") or die "cant write $data_dir/accession2clone.dat :$! ";
  
  print C2A Data::Dumper->Dump([\%clone2accession]);
  print A2C Data::Dumper->Dump([\%accession2clone]);
  
  close C2A;
  close A2C;
  
}


########################################################################################################

sub write_cds2wormpep  {   

  my $WPver = &get_wormbase_version;
  open (FH,"<$basedir/WORMPEP/wormpep$WPver/wormpep$WPver") or die "cant open wormpep$WPver\n";
  my %cds2wormpep;
  my %wormpep2cds;
  while(<FH>) {
    if( />/ ) {
      chomp;
      my @data = split;
      # >2L52.1 CE32090   Zinc finger, C2H2 type status:Predicted TR:Q9XWB3
      my $pep = $data[1];
      my $gene = substr("$data[0]",1);
      $cds2wormpep{$gene} = $pep;
      $wormpep2cds{$pep} .= "$gene ";
    }
  }
  
  #now dump data to file
  open (C2G, ">$data_dir/wormpep2cds.dat") or die "$data_dir/wormpep2cds.dat";
  open (G2C, ">$data_dir/cds2wormpep.dat") or die "$data_dir/cds2wormpep.dat";
  
  print C2G Data::Dumper->Dump([\%wormpep2cds]);
  print G2C Data::Dumper->Dump([\%cds2wormpep]);
  
  close C2G;
  close G2C;
}

########################################################################################################


sub write_CDSlist  {   

  my %CDSlist;
  my $CDS;

  # connect to AceDB using TableMaker,
  my $command="Table-maker -p $wquery_dir/CDSlist.def\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  while (<TACE>) {
    chomp;
    next if ($_ eq "");
    last if (/\/\//);
    ($CDS) = (/^\"(\S+)\"/);
    # Flag confirmed genes with a 1
    if (/Confirmed_by/) {$CDSlist{$CDS} = 1;}
    else {$CDSlist{$CDS} = 0;}
    print "assigned $CDS with status '$CDSlist{$CDS}'\n" if ($test);
  }
  close TACE;
  
  #now dump data to file
  open (CDS, ">$data_dir/CDS_list.dat") or die "Can't open file: $data_dir/CDS_list.dat";
  print CDS Data::Dumper->Dump([\%CDSlist]);
  close CDS;
}
########################################################################################################


sub write_EST  {   

  my %NDBaccession2est;
  my %estorientation;
  my @f;

  # connect to AceDB using TableMaker,
  my $command="Table-maker -p $wquery_dir/ESTdata.def\nquit\n";
  
  open (TACE, "echo '$command' | $tace $ace_dir |");
  while (<TACE>) {
    chomp;
    next if ($_ eq "");    # shortcut at empty lines
    last if (/\/\//);      # end when you get to the end
    s/\"//g;               # remove speech marks
    @f = split /\t/;
    
    # NDB_accession to WormBase name
    $NDBaccession2est{$f[0]} = $f[1];

    # EST orientation
    $estorientation{$f[1]} = 5 if ($f[3]);
    $estorientation{$f[1]} = 3 if ($f[4]);

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
 
  my %clone2seq;
  my $seq; my $newname; my $seqname;
 
  # connect to AceDB using tace
  my $command = "query find Genome_sequence\nDNA\nquit\n";
 
  open (TACE, "echo '$command' | $tace $ace_dir | ");
  while (<TACE>) {
    chomp;
    next if ($_ eq "");
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

  my %genes2lab;

  # connect to AceDB using TableMaker, 3 separate queries for CDS, Transcript and Pseudogene
  my @commands;
  $commands[0] = "Table-maker -p $wquery_dir/cds2lab.def\nquit\n";
  $commands[1] = "Table-maker -p $wquery_dir/pseudogene2lab.def\nquit\n";
  $commands[2] = "Table-maker -p $wquery_dir/RNAgene2lab.def\nquit\n";

  # add to hash for each of the three queries
  foreach my $command (@commands){
    open (TACE, "echo '$command' | $tace $ace_dir |");
    while (<TACE>) {
      chomp;
      next if ($_ eq "");
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

  my %worm_gene2cgc;
  my %cgc_name2gene;

  # connect to AceDB using TableMaker, but use /wormsrv2/geneace for Table-maker definition
  my $command="Table-maker -p /wormsrv2/geneace/wquery/cgc_names_for_worm_genes.def\nquit\n";
  
  open (TACE, "echo '$command' | $tace /wormsrv2/geneace |");
  while (<TACE>) {
    chomp;
    next if ($_ eq "");
    last if (/\/\//);

    # get rid of quote marks
    s/\"//g;

    # split the line into various fields
    my ($gene,$cds,$transcript,$pseudogene,$cgc_name) = split(/\t/, $_) ;

    # populate first hash with CGC name as key
    $cgc_name2gene{$cgc_name} = $gene;

    # add to hash, molecular name is key.  Value consists of three parts:
    # 1) class of molecular name, 2) CGC name 3) ?Gene name

    if($cds){
      $worm_gene2cgc{$cds} = "CDS $cgc_name $gene";
    }
    if($transcript){
      $worm_gene2cgc{$transcript} = "Transcript $cgc_name $gene";
    }
    if($pseudogene){
      $worm_gene2cgc{$pseudogene} = "Pseudogene $cgc_name $gene";
    }
  }
  close TACE;
  
  #now dump data to file
  open (CGC, ">$data_dir/worm_gene2cgc_name.dat") or die "Can't open file: $data_dir/worm_gene2cgc_name.dat";
  print CGC Data::Dumper->Dump([\%worm_gene2cgc]);
  close CGC;

  #now dump data to file
  open (CGC, ">$data_dir/cgc_name2gene.dat") or die "Can't open file: $data_dir/cgc_name2gene.dat";
  print CGC Data::Dumper->Dump([\%cgc_name2gene]);
  close CGC;

}

###########################################################################################################



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

%CDS_list - list of CDSs with confirmed status as value (1 = confirmed, 0 = predicted or partially confirmed)

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

=item *

%worm_gene2cgc - connections between C. elegans CDS, Transcript, or Pseudogene and CGC name...the value actually stores three fields
(separated by spaces). 1) class name of the key (CDS, Transcript, or Pseudogene) 2) CGC name 3) Gene ID

=back

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
