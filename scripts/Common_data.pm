#!/usr/local/bin/perl5.6.1 -w
          
# 
# by Anthony Rogers                             
#
# Last updated by: $Author: ar2 $               
# Last updated on: $Date: 2002-11-19 14:33:26 $         

use strict;                    
use lib "/wormsrv2/scripts/";
use Wormbase;
use Data::Dumper;
my $data_dir = "/wormsrv2/autoace/COMMON_DATA";

our %sub2file = ( 'gene2CE' => "$data_dir/gene2CE.dat",
		  'CE2gene' => "$data_dir/CE2gene.dat",
		  'clone2acc' => "$data_dir/clone2acc.dat",
		  'acc2clone' => "$data_dir/acc2clone.dat"
		 # '' => "$data_dir/.dat",
		 # '' => "$data_dir/.dat",
		);


# Data writing routines - actually create and dump the data
##################################
sub write_clone2acc
  {
    # AceDB database
    my $ace_dir = "/wormsrv2/autoace";
    my $wquery_dir = "/wormsrv2/autoace/wquery";
    my $tace = &tace;

    my %clone2acc;
    my %acc2clone;
    
    ####################################################################
    # connect to AceDB using TableMaker,
    # populating %accession2name (maps embl accession to contig name)
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
