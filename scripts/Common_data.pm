#!/usr/local/bin/perl5.8.0 -w
#
# Common_data.pm          
# 
# by Anthony Rogers                             
#
# Last updated by: $Author: dl1 $               
# Last updated on: $Date: 2004-03-16 12:46:12 $         

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

use vars qw($update $build $a2c $c2g $g2p $cloneseq $clonesize $all $predicted_CDS $test);

GetOptions("update"        => \$update,
	   "in_build"      => \$build,
	   "accession"     => \$a2c,
	   "ce"            => \$c2g,
	   "pid"           => \$g2p,
	   "cloneseq"      => \$cloneseq,
	   "clonesize"     => \$clonesize,
	   "all"           => \$all,
	   "predicted_CDS" => \$predicted_CDS,
	   "test"          => \$test
	   );


##############################
# database paths             #
##############################

# Set up top level base directory which is different if in test mode
# Make all other directories relative to this
my $basedir   = "/wormsrv2";
$basedir      = glob("~wormpub")."/TEST_BUILD" if ($test); 
my $db_path   = "$basedir/autoace";

my $this_file   = "$basedir/scripts/Common_data.pm";
our $data_dir   = "$basedir/autoace/COMMON_DATA";
our $wquery_dir = "$basedir/autoace/wquery";
our $ace_dir;

our %sub2file = ( 'gene2CE'   => "$data_dir/gene2CE.dat",
		  'CE2gene'   => "$data_dir/CE2gene.dat",
		  'clone2acc' => "$data_dir/clone2acc.dat",
		  'acc2clone' => "$data_dir/acc2clone.dat",
		  'gene2pid'  => "$data_dir/gene2pid.dat",
		  'cloneseq'  => "$data_dir/clone2seq.dat",
		  'clonesize' => "$data_dir/clonesize.dat",
		  'CDSlist'   => "$data_dir/CDSlist.dat",
		  'pid2gene'  => "$data_dir/pid2gene.dat"
		);

##############################
# ACEDB executables          #
##############################

our $tace = &tace;

##############################
# debug mode                 #
##############################

my $debug = 0;

# update mode 
if( $update ) {
  print "Updating $data_dir data files ";

  # AceDB data
  if( $build ) {
    $ace_dir = "$basedir/autoace";
    print "during build so using $ace_dir - ensure that the data you are updating is actually in the database.\n";
  }
  else {
    $ace_dir = "/nfs/disk100/wormpub/DATABASES/current_DB";
    print "- NOT as part of build so using $ace_dir. If this is part of the build data MAY be stale\n";
  }
  &Common_data_update;
}
# checks to stop you running the writeable outside update mode
else {
  if ( ($g2p) || ($a2c) || ($c2g) ) {
    print "please update using the script ie Common_data.pm -update -pid\n";
    `perldoc $this_file`;
    exit(1);
  }
}

sub Common_data_update {
    &write_gene2pid  if $g2p;
    &write_clone2acc if $a2c;
    &write_gene2CE   if $c2g;
    &write_CDSlist   if $predicted_CDS;
    &write_cloneseq  if $cloneseq;
    &write_clonesize if $clonesize;

}


#######################################################################
# Data writing routines - actually create and dump the data           #
#######################################################################

sub write_gene2pid {

    my %gene2pid;
    my %pid2gene;
    
    ####################################################################
    # connect to AceDB using TableMaker,
    # populating %accession2name (maps embl accession to contig name)
    ####################################################################
    my $command="Table-maker -p $wquery_dir/gene2pid.def\nquit\n";
    
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

sub write_clone2acc  {   

    my %clone2acc;
    my %acc2clone;
    
    ####################################################################
    # connect to AceDB using TableMaker,
    # populating %gene2pid
    ####################################################################
    my $command="Table-maker -p $wquery_dir/accession2clone.def\nquit\n";
    
    open (TACE, "echo '$command' | $tace $ace_dir |");
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

sub write_gene2CE  {   

    my $WPver = &get_wormbase_version;
    $WPver = "666" if ($test);

    open (FH,"<$basedir/WORMPEP/wormpep$WPver/wormpep$WPver") or die "cant open wormpep$WPver\n";
    my %gene2CE;
    my %CE2gene;
    while(<FH>) {
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

sub write_CDSlist  {   

    my %CDSlist;
    my $CDS;

    ####################################################################
    # connect to AceDB using TableMaker,
    # populating %gene2pid
    ####################################################################
    my $command="Table-maker -p $wquery_dir/CDSlist.def\nquit\n";
    
    open (TACE, "echo '$command' | $tace $ace_dir |");
    while (<TACE>) {
	chomp;
	next if ($_ eq "");
	last if (/\/\//);
	($CDS) = (/^\"(\S+)\"/);
	if (/Confirmed_by/) {$CDSlist{$CDS} = 1;}
	else {$CDSlist{$CDS} = 0;}
	print "assigned $CDS with status '$CDSlist{$CDS}'\n" if ($debug);
    }
    close TACE;
    
    #now dump data to file
    open (CDS, ">$sub2file{'CDSlist'}") or die "Can't open file: $sub2file{'CDSlist'}";
    print CDS Data::Dumper->Dump([\%CDSlist]);
    close CDS;
}

sub write_cloneseq  {   

    my %clone2seq;
    my $seq; my $newname; my $seqname;

    ####################################################################
    # connect to AceDB using tace
    ####################################################################

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
    open (CDS, ">$sub2file{'cloneseq'}") or die "Can't open file: $sub2file{'cloneseq'}";
    print CDS Data::Dumper->Dump([\%clone2seq]);
    close CDS;
}



####################################################################################

####################
#Return a true value
####################

1;


__END__

=pod

=head1 NAME - Common_data.pm

=head2 DESCRIPTION

The Common_data.pm gives quick easy acces to a variety of data frequently used in Wormbase scripts.
It comprises of one part.

Part One generates the data and writes it to a file using the Data::Dumper module.

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

=over 4

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

=over 4

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
