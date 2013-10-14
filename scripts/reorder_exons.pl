#!/usr/local/bin/perl5.8.0 -w 
#
# reorder_exons.pl
#
# by Darin Blasiar
#
# This script checks the exon order and corrects them if needed
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-10-14 10:16:24 $



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
my $debug;    # Log messages will only go to the specified user.
my $test;     # use test environment?
my $quicktest;# use test environment but gets the defined clone.
my $store;    # Provide the Wormbase store object to the script
my $out;      # allows an output file to be specified to avoid writing to BUILD
my $verbose;  # gives additional output if interested
my $noload;

GetOptions ("database=s"     => \$database,
	    "debug:s"        => \$debug,
            "test"           => \$test,
	    "quicktest:s"    => \$quicktest,
	    "store:s"        => \$store,
	    "out:s"          => \$out,
	    "noload"         => \$noload,
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
my $acefile;
if (defined $out) {
  $acefile = $out;
} else {
  $acefile = $dbpath."/acefiles/sorted_exons.ace";
}

print "output will be written to $acefile\n";
open (ACE, ">$acefile") || die $!;


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
my (@sequences,$query);
if( $quicktest ) {
  $query = "Find Sequence $quicktest";
  print "only fetching data from $quicktest\n\n" if ($verbose);
}
else {
  $query = "Find Sequence Properties == Genomic_canonical & Species=\"Caenorhabditis elegans\"";
}

@sequences = $db->fetch(-query => $query);

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
		print "exons out of order in $cds\n" if ($verbose);
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
		print "exons out of order for $transcript\n" if ($verbose);
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

unless ($noload) {
  my $command = "pparse $acefile";
  $command .= "\nsave\nquit\n";
  print "\nUploading sorted exon info to $database\n";
  open (WRITEDB, "| $tace -tsuser reorder_exons $dbpath") || die "Couldn't open pipe to $database\n";
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
	print "$$sub_ref{$_}->[0] $$sub_ref{$_}->[1] $$sub_ref{$_}->[2] $$sub_ref{$_}->[3] \n" if ($verbose);
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


=head1 MANDATORY arguments:

=over 4

=item none

=back

=head1 OPTIONAL arguments: -database, -test, -out, -debug, -verbose, -quicktest

=over 4

=item -database <Database path>

Specifies location of target database, will default to using autoace if not specified

=item -test

Uses test environment in ~wormpub/TEST_BUILD/

=item -quicktest <Clone>

This option allows you to specify a clone that you know has an inconsistent CDS/Transcript

=item -out <File destination>

If you specify a file, the data will be output into this file instead of touching the build env.

=item -debug <User id>

Log messages will only go to the specified user.

=item -verbose

gives additional output if interested

=back

=head1 USAGE Example:

=over 4

=item reorder_exons.pl -database ~/DATABASES/camace -out ~/DATABASES/camace/reordered_exons.ace -debug pad -test 

This will check the database camace for inconsistent exon order and output a patch file to "~/DATABASES/camace/reordered_exons.ace". 

=back

=head1 AUTHOR Darin Blasiar (dblasiar@watson.wustl.edu) 

=cut

