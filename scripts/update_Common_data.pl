#!/usr/local/bin/perl5.8.0 -w
#
# update_Common_data.pl
# 
# by Anthony Rogers
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2003-12-05 16:05:12 $

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

use vars qw($update $build $a2c $c2g $g2p $all $predicted_CDS $test);

GetOptions("update"        => \$update,
	   "in_build"      => \$build,
	   "accession"     => \$a2c,
	   "ce"            => \$c2g,
	   "pid"           => \$g2p,
	   "all"           => \$all,
	   "predicted_CDS" => \$predicted_CDS,
	   "test"          => \$test
	   );

# do all of the data sets if -all
#if ($all) {
#    $c2g = 1; $a2c = 1; $g2p = 1;
#}


##########################################
# Set up database paths                  #
##########################################

# Set up top level base directory which is different if in test mode
# Make all other directories relative to this
my $basedir   = "/wormsrv2";
$basedir      = glob("~wormpub")."/TEST_BUILD" if ($test); 

our $data_dir   = "$basedir/autoace/COMMON_DATA";


our $wquery_dir = "$basedir/autoace/wquery";
our $ace_dir;

our %sub2file = ( 'gene2CE'   => "$data_dir/gene2CE.dat",
		  'CE2gene'   => "$data_dir/CE2gene.dat",
		  'clone2acc' => "$data_dir/clone2acc.dat",
		  'acc2clone' => "$data_dir/acc2clone.dat",
		  'gene2pid'  => "$data_dir/gene2pid.dat",
		  'CDSlist'   => "$data_dir/CDSlist.dat",
		  'pid2gene'  => "$data_dir/pid2gene.dat"
		);

##############################
# ACEDB executables          #
##############################

our $tace = &tace;

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
	`perldoc $0`;
	exit(1);
    }
}

sub Common_data_update {
    &write_gene2pid  if ( $g2p || $all );
    &write_clone2acc if ( $a2c || $all );
    &write_gene2CE   if ( $c2g || $all );
    &write_CDSlist   if ( $predicted_CDS || $all );
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
	print "assigned $CDS with status '$CDSlist{$CDS}'\n" if ($test);
    }
    close TACE;
    
    #now dump data to file
    open (CDS, ">$sub2file{'CDSlist'}") or die "Can't open file: $sub2file{'CDSlist'}";
    print CDS Data::Dumper->Dump([\%CDSlist]);
    close CDS;
}

####################################################################################

####################
#Return a true value
####################

1;


__END__

=pod

=head1 NAME - update_Common_data.pl

=head2 DESCRIPTION

The update_Common_data.pl updates the common data sets retrieved by the Fetch_data routine in Wormbase.pmgives quick easy acces to a variety of data frequently used in Wormbase scripts.
It comprises of one part.

Part One generates the data and writes it to a file using the Data::Dumper module.

This module updates the following data sets:

=over 4

=item *

%acc2clone

=item *

%clone2acc

=item *

%predictedCDSs

=item *

%CE2gene  -  bear in mind that a peptide may have multiple genes. If so the genes are concatenated, separated by a space.

=item *

%gene2CE

=back

=over 4

=head2 EXAMPLE USAGE

=over 4

=item UPDATING THE DATA

=over 4

We need to be very careful about updating this data.  Depending in wether it is being updated during the build or otherwise we need to use autoace or current_DB. 

=item 

update_Common_data.pm -update -build -all

However if the build is underway and you want to write out the gene 2 CE info but that is not yet in the building database you need to use current_DB. so dont include -build ie 

=item

Common_data.pm -update -ce


At the end of the build ( in finish_build.pl) all data will be refreshed to be in synch with the current
 release ie

=item

update_Common_data.pl -update -build -all.

There shouldn't really be any need to alter this apart from actually during the next build itself.  Scripts that generate the data that is included in Common_data should call the updating routines themselves - so you wont have to ! 
