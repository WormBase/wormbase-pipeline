#!/usr/local/bin/perl5.8.0 -w 
#
# reorder_exons.pl
#
# by Darin Blasiar
#
# This script checks the exon order and corrects them if needed
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2005-12-16 11:18:55 $



#############################################
# USAGE: reorder_exons.pl <database>        #
#############################################

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Ace;
use Cwd;
use Log_files;
use Storable;

#################################
# Command-line options          #
#################################

my $database; # which database to test
my $debug;
my $test;     # use test environment?
my $quicktest;# use test environment but gets less clones?
my $store;


GetOptions ("database=s"   => \$database,
	    "debug:s"      => \$debug,
            "test"         => \$test,
	    "quicktest"    => \$quicktest,
	    "store:s"      => \$store
	   );

$test = 1 if $quicktest;
my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}
#################################
# Set paths                     #
#################################

# set root level
my $basedir     = $wormbase->basedir;

# check database name and path, default to autoace if -database not specified
my $dbpath;
my $CWD = cwd;

if(!$database){
  print "No database specified, using ".$wormbase->autoace." by default\n\n";
  $dbpath = $wormbase->autoace;
}
elsif($database =~ /^(\~\w+)\//){ # do we need to expand path if ~path specified?
  $dbpath = glob("$1");
  $dbpath =~ s/\/tmp_mnt//;
  my $filename = "$'";
  $dbpath = "$dbpath"."/"."$filename";
} 
elsif($database =~ /^(\w+)/) { # for incomplete paths, expand using CWD
  $dbpath="$CWD"."/"."$database";
}  
elsif ($database =~ /\/\w+/) { # else assume path is ok
  $dbpath=$database;
} 
else {
  die "-database does not specify a valid path to a database\n";
}


# open output file
open (ACE, ">$dbpath/acefiles/sorted_exons.ace") || die $!;


#########################
# connect to database   #
#########################

my $tace = $wormbase->tace;

my $db = Ace->connect (-path => "$dbpath",
		       -program => $tace) || 
    die "cannot connect to $dbpath\n";
$db->auto_save(0);



####################################
# get genomic_canonical sequences  #
####################################
my @sequences;
if( $quicktest ) {
  @sequences = $db->fetch(-query => 'Find Sequence AH6');
}
else {
  @sequences = $db->fetch(-query => 'Find Sequence Properties == Genomic_canonical & Species="Caenorhabditis elegans"');
}

unless($sequences[0]) {
    die print "database contained no genomic_canonical sequences\n";
}


#################################################################
# set up hashes of subsequence objects and transcript objects   #
#################################################################


my %subsequencehash;
my %transcripthash;
my @seqoutoforder;
my @tranoutoforder;

foreach (@sequences) {
  my $clone = $_;
  my $obj = $db->fetch(Sequence => $clone);
  my @CDSs = $obj->CDS_child;
  my @transcripts = $obj->Transcript(1);
  
#######################################
# check objects of %subsequencehash   #
#######################################

  foreach(@CDSs) {
    my $cds = $_;
    my $obj = $db->fetch(CDS => $cds);
    my @exons = $obj->Source_exons(1);
    
    my $oldleft = 0;
    foreach(@exons) {
      my ($a, $b) = $_->row;
      my @ends = ($a, $b);
      if($a < $oldleft) {
#		print "exons out of order in $cds\n";
	push(@seqoutoforder, $cds);
      }
      $oldleft = $a;
      push(@{$subsequencehash{$cds}}, @ends);
    }
  }
  
  
######################################
# check objects of %transcripthash   #
######################################
  
  foreach(@transcripts) {
    
    my $transcript = $_;
    my $obj = $db->fetch(Transcript => $transcript);
    my @exons = $obj->Source_exons(1);
    
    my $oldleft = 0;
    foreach(@exons) {
      my ($a, $b) = $_->row;
      my @ends = ($a, $b);
      if($a < $oldleft) {
	#	print "exons out of order for $transcript\n";
	push(@tranoutoforder, $transcript);
      }
      $oldleft = $a;
      push(@{$transcripthash{$transcript}}, @ends);
    }
  }
}



################################
# fix exon order if necessary  #
################################


if($seqoutoforder[0]) {

  my $seqoutoforder_ref = \@seqoutoforder;
  my $subsequence_ref = \%subsequencehash;
  
  &seq_fix_order($seqoutoforder_ref, $subsequence_ref);
}

if($tranoutoforder[0]) {

  my $tranoutoforder_ref = \@tranoutoforder;
  my $transcript_ref = \%transcripthash;
  
  &tran_fix_order($tranoutoforder_ref, $transcript_ref);
}

##############################################
# Upload to autoace, tidy up and exit cleanly
##############################################
$db->close;
close(ACE);

# Only upload if working with autoace
if($dbpath =~ m/autoace$/){
  my $command = "pparse $dbpath/acefiles/sorted_exons.ace";
  $command .= "\nsave\nquit\n";
  print "\nUploading sorted exon info to autoace\n";
  open (WRITEDB, "| $tace -tsuser reorder_exons $dbpath") || die "Couldn't open pipe to autoace\n";
  print WRITEDB $command;
  close WRITEDB;

}

print "\nFinished\n";
exit(0);

#################################
# SUBROUTINES    ################
#################################

sub seq_fix_order {
    my ($order_ref, $sub_ref) = @_;
    foreach(@{$order_ref}) {
	print ACE "\nCDS : \"$_\"\n";
	print ACE "-D Source_exons\n\n";
	print ACE "CDS : \"$_\"\n";
#	print ACE "$$sub_ref{$_}->[0] $$sub_ref{$_}->[1] $$sub_ref{$_}->[2] $$sub_ref{$_}->[3] \n";
	my @thisarray = @{$$sub_ref{$_}}; 
	my @sorted = sort {$a <=> $b} (@thisarray);
	foreach(my $ii = 0; $ii <= $#sorted; $ii+=2) {
	    print ACE "Source_exons $sorted[$ii] $sorted[$ii+1]\n"
	}

    }

}

sub tran_fix_order {
    my ($order_ref, $tran_ref) = @_;
    foreach(@{$order_ref}) {
	print ACE "\nTranscript : \"$_\"\n";
	print ACE "-D Source_exons\n\n";
	print ACE "Transcript : \"$_\"\n";
	my @thisarray = @{$$tran_ref{$_}}; 
	my @sorted = sort {$a <=> $b} (@thisarray);
	foreach(my $ii = 0; $ii <= $#sorted; $ii+=2) {
	    print ACE "Source_exons $sorted[$ii] $sorted[$ii+1]\n"
	}

    }

}






__END__

=pod

=head1 NAME - reorder_exons.pl

=head2 USAGE

=over 4

=item reorder_exons.pl -[options]

=back

=head1 DESCRIPTION

Sorts all exons in an acedb database to be in a linear order.  This is useful when new
gene curation adds a new 5' exon to a gene but by default acedb adds this to the end
of the list of exons.

Usually run as part of the build but can be run against any valid acedb database that has
'Source_exons'

=back

=head1 MANDATORY arguments:

=over 4

=item none

=back

=head1 OPTIONAL arguments: -database, -test

=over 4

=item -database <path to database>

Specifies location of target database, will default to using autoace if not specified

=item -test

Uses test environment in ~wormpub/TEST_BUILD/

=back


=head1 AUTHOR Darin Blasiar (dblasiar@watson.wustl.edu) 

=back

=cut

