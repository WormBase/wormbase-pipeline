#!/usr/local/bin/perl5.6.1 -w
          
# 
# by Anthony Rogers                             
#
# Last updated by: $Author: ar2 $               
# Last updated on: $Date: 2002-11-28 12:10:23 $         

use strict;                    
use lib "/wormsrv2/scripts/";
use Wormbase;
use Data::Dumper;
use Getopt::Long;

my $update;
my $build;
my $a2c;
my $c2g;
my $g2p;
my $all;
my $test;

GetOptions("update"     => \$update,
	   "in_build"   => \$build,
	   "accession"  => \$a2c,
	   "ce"         => \$c2g,
	   "pid"        => \$g2p,
	   "all"        => \$all,
	   "test"       => \$test
	  );

if( $all ) {
  $c2g = 1;
  $a2c = 1;
  $g2p = 1;
}

my $this_file = "/wormsrv2/scripts/Common_data.pm";
our $data_dir = "/wormsrv2/autoace/COMMON_DATA";
$data_dir = "/wormsrv2/tmp" if $test;
our $ace_dir;
our $wquery_dir = "/wormsrv2/autoace/wquery";
our $tace = &tace;

if( $update ) {
  print "Updating $data_dir data files ";

  # AceDB data
  if( $build ) {
    $ace_dir = "/wormsrv2/autoace";
    print "during build so using $ace_dir - ensure that the data you are updating is actually in the database.\n";
  }
  else {
    $ace_dir = "/wormsrv2/current_DB";
    print "- NOT as part of build so using $ace_dir. If this is part of the build data MAY be stale\n";
  }
  &Common_data_update;
}



our %sub2file = ( 'gene2CE' => "$data_dir/gene2CE.dat",
		  'CE2gene' => "$data_dir/CE2gene.dat",
		  'clone2acc' => "$data_dir/clone2acc.dat",
		  'acc2clone' => "$data_dir/acc2clone.dat",
		  'gene2pid' => "$data_dir/gene2pid.dat",
		  'pid2gene' => "$data_dir/pid2gene.dat"
		);


sub Common_data_update
  {
    &write_gene2pid if $g2p;
    &write_clone2acc if $a2c;
    &write_gene2CE if $c2g;
  }

# Data writing routines - actually create and dump the data
##################################
sub write_gene2pid
  {
    unless( $update ) {
      print "please update using the script ie Common_data.pm -update -pid\n";
      `perldoc $this_file`;
      exit(1);
    }

    my %gene2pid;
    my %pid2gene;
    
    ####################################################################
    # connect to AceDB using TableMaker,
    # populating %accession2name (maps embl accession to contig name)
    ####################################################################
    my $command="Table-maker -p $wquery_dir/gene2pid.def\nquit";
    
    open (TACE, "echo '$command' | $tace $ace_dir |");
    while (<TACE>) {
      #gene pid version
      print;
      chomp;
      if (/\"(\S+)\"\s+\"(\S+)\"\s+(\d)/) {
	my $pid = "$2".".$3";
	my $gene = $1;
	$gene2pid{"$gene"} = $pid;
	$pid2gene{"$pid"} = $gene;
      }
    }
    close TACE;

    #now dump data to file
    open (G2P, ">$sub2file{'gene2pid'}") or die "cant write $sub2file{'gene2pid'} :$!";
    open (P2G, ">$sub2file{'pid2gene'}") or die "cant write $sub2file{'pid2gene'} :$! ";

    print G2P Data::Dumper->Dump([\%gene2pid]);
    print P2G Data::Dumper->Dump([\%pid2gene]);
    
    close G2P;
    close P2G;
        
  }

sub write_clone2acc
  {   
    unless( $update ) {
      print "please update using the script ie Common_data.pm -update -accession\n";
      `perldoc $this_file`;
      exit(1);
    }

    my %clone2acc;
    my %acc2clone;
    
    ####################################################################
    # connect to AceDB using TableMaker,
    # populating %gene2pid
    ####################################################################
    my $command="Table-maker -p $wquery_dir/accession2clone.def\nquit";
    
    open (TACE, "echo '$command' | $tace /wormsrv2/autoace |");
    while (<TACE>) {
      chomp;
      if (/\"(\S+)\"\s+\"(\S+)\"\s*\d*/) {
	$clone2acc{$1} = $2;
	$acc2clone{$2} = $1;
      }
    }
    close TACE;

    #now dump data to file
    open (C2A, ">$sub2file{'clone2acc'}") or die "cant write $sub2file{'clone2acc'} :$!";
    open (A2C, ">$sub2file{'acc2clone'}") or die "cant write $sub2file{'acc2clone'} :$! ";

    print C2A Data::Dumper->Dump([\%clone2acc]);
    print A2C Data::Dumper->Dump([\%acc2clone]);
    
    close C2A;
    close A2C;
        
  }

sub write_gene2CE
  {   
    unless( $update ) {
      print "please update using the script ie Common_data.pm -update -ce\n";
      exec ('perldoc', $this_file);
      exit(1);
    }
    my $WPver = &get_wormbase_version;
    open (FH,"</wormsrv2/WORMPEP/wormpep$WPver/wormpep$WPver") or die "cant open wormpep$WPver\n";
    my %gene2CE;
    my %CE2gene;
    while(<FH>)
      {
	if( />/ ) {
	  chomp;
	  my @data = split;
	  # >2L52.1 CE32090   Zinc finger, C2H2 type status:Predicted TR:Q9XWB3
	  my $pep = $data[1];
	  my $gene = substr("$data[0]",1);
	  $gene2CE{$gene} = $pep;
	  $CE2gene{$pep} .= "$gene ";
	}
      }
    
    #now dump data to file
    open (C2G, ">$sub2file{'CE2gene'}") or die "$sub2file{'CE2gene'}";
    open (G2C, ">$sub2file{'gene2CE'}") or die "$sub2file{'gene2CE'}";

    print C2G Data::Dumper->Dump([\%CE2gene]);
    print G2C Data::Dumper->Dump([\%gene2CE]);
    
    close C2G;
    close G2C;
  }


# Data retrieval routines - all work thru same sub but pass different files

sub gene2pid 
  {
    my $ref = shift;
    my $file = $sub2file{'gene2pid'};
    &getData($ref, $file);
  }

sub clone2acc 
  {
    my $ref = shift;
    my $file = $sub2file{'clone2acc'};
    &getData($ref, $file);

  }

#data for this dumped from write_clone2acc sub
sub acc2clone
  {
    my $ref = shift;
    my $file = $sub2file{'acc2clone'};
    &getData($ref, $file);

  }

sub gene2CE
  {
    my $ref = shift;
    my $file = $sub2file{'gene2CE'};
    &getData($ref, $file);
  }

sub CE2gene
  {
    my $ref = shift;
    my $file = $sub2file{'CE2gene'};
    &getData($ref, $file);
  }

# Actually gets the data from the .dat files and populated the hash whose reference is passed.
sub getData
  {
    my $ref = shift;
    my $file = shift;
    open (FH, "<$file") or die "cant open $file";
    undef $/;
    my $VAR1;
    my $data = <FH>;
    eval $data;
    die if $@;
    $/ = "\n";
    close FH;
    %$ref = (%$VAR1);    
  }

####################
#Return a true value
####################

1;


__END__

=pod

=head1 NAME - Common_data.pm

=head2 DESCRIPTION

The Common_data.pm gives quick easy acces to a variety of data frequently used in Wormbase scripts.
It comprises of two parts

Part One generates the data and writes it to a file using the Data::Dumper module.
Part Two is the rapid retrieval of this data ( rather than recreating it multiple times in different scripts)


This module provides access to the following data sets:

=over 4

=item *

%acc2clone

=item *

%clone2acc

=item *

%CE2gene  -  bear in mind that a peptide may have multiple genes. If so the genes are concatenated, separated by a space.

=item *

%gene2CE

=back

=over 4

=head2 EXAMPLE USAGE

=over4

=item include this file

  use Common_data.pm;

=item

#declare the hash you want to use

 my %hash_of_genes_2_CEs;

=item

#call the relavent sub-routine passing a reference to the hash you just declared

&gene2CE(\%hash_of_genes_2_CEs);

=item

#use your hash as if you'd filled it yourself

my $histone = $hash_of_genes_2_CEs{'F08G2.2'}

=item

returns 'CE04501'

=item UPDATING THE DATA

=over4

We need to be very careful about updating this data.  Depending in wether it is being updated during the build or otherwise we need to use autoace or current_DB. 

If you have just generated new accession to clone information in the build, autoace should be used by adding the -build option ie

=item 

Common_data.pm -update -build -a2c

However if the build is underway and you want to write out the gene 2 CE info but that is not yet in the building database you need to use current_DB. so dont include -build ie 

=item

Common_data.pm -update -a2c


At the end of the build ( in finish_build.pl) all data will be refreshed to be in synch with the current release ie

=item

Common_data.pm -update -build -all.

There shouldn't really be any need to alter this apart from actually during the next build itself.  Scripts that generate the data that is included in Common_data will call the updating routines themselves - so you wont have to ! 

BTW If you try and call the updating routines directly they will complain.
