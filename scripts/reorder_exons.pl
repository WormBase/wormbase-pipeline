#!/usr/local/bin/perl5.8.0 -w 
#
# reorder_exons.pl
#
# by Darin Blasiar
#
# This script checks the exon order and corrects them if needed
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2003-07-24 13:45:20 $



#############################################
# USAGE: reorder_exons.pl <database>        #
#############################################

use strict;
use Ace;
use lib "/wormsrv2/scripts";
use Wormbase;

# database paths at the Sanger
my %db_list = (
	    'autoace',  '/wormsrv2/autoace',
            'curr_db',  '/nfs/disk100/wormpub/DATABASES/current_DB'  
	);

open (ACE, ">/wormsrv2/autoace/acefiles/sorted_exons.ace") || die $!;

#########################
# connect to database   #
#########################

my $tace = &tace;
my $db_name = $ARGV[0];
my $db_path = $db_list{$db_name};
my $db = Ace->connect (-path => "$db_path",
		       -program => $tace) || 
    die "cannot connect to acedb\n";
$db->auto_save(0);



####################################
# get genomic_canonical sequences  #
####################################

my @sequences;

if ($db_name eq 'autoace') {
    @sequences = $db->fetch(-query => 'Find Sequence Properties == Genomic_canonical & Species="Caenorhabditis elegans"');
}
if ($db_name eq 'curr_db') {
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
  my @subsequences = $obj->Subsequence(1);
  my @transcripts = $obj->Transcript_child(1);
  
#######################################
# check objects of %subsequencehash   #
#######################################

  foreach(@subsequences) {
    my $clone = $_;
    my $obj = $db->fetch(Sequence => $clone);
    my @exons = $obj->Source_Exons(1);
    
    my $oldleft = 0;
    foreach(@exons) {
      my ($a, $b) = $_->row;
      my @ends = ($a, $b);
      if($a < $oldleft) {
#		print "exons out of order in $clone\n";
	push(@seqoutoforder, $clone);
      }
      $oldleft = $a;
      push(@{$subsequencehash{$clone}}, @ends);
    }
  }
  
  
######################################
# check objects of %transcripthash   #
######################################
  
  foreach(@transcripts) {
    
    my $clone = $_;
    my $obj = $db->fetch(Transcript => $clone);
    my @exons = $obj->Source_Exons(1);
    
    my $oldleft = 0;
    foreach(@exons) {
      my ($a, $b) = $_->row;
      my @ends = ($a, $b);
      if($a < $oldleft) {
	#	print "exons out of order for $clone\n";
	push(@tranoutoforder, $clone);
      }
      $oldleft = $a;
      push(@{$transcripthash{$clone}}, @ends);
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
if($db_name eq "autoace"){
  my $command = "pparse /wormsrv2/autoace/acefiles/sorted_exons.ace";
  $command .= "save\nquit\n";
  print "\nUploading sorted exon info to autoace\n";
  open (WRITEDB, "| $tace -tsuser reorder_exons /wormsrv2/autoace |") || die "Couldn't open pipe to autoace\n";
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
	print ACE "\nSequence : \"$_\"\n";
	print ACE "-D Source_Exons\n\n";
	print ACE "Sequence : \"$_\"\n";
#	print ACE "$$sub_ref{$_}->[0] $$sub_ref{$_}->[1] $$sub_ref{$_}->[2] $$sub_ref{$_}->[3] \n";
	my @thisarray = @{$$sub_ref{$_}}; 
	my @sorted = sort {$a <=> $b} (@thisarray);
	foreach(my $ii = 0; $ii <= $#sorted; $ii+=2) {
	    print ACE "Source_Exons $sorted[$ii] $sorted[$ii+1]\n"
	}

    }

}

sub tran_fix_order {
    my ($order_ref, $tran_ref) = @_;
    foreach(@{$order_ref}) {
	print ACE "\nTranscript : \"$_\"\n";
	print ACE "-D Source_Exons\n\n";
	print ACE "Transcript : \"$_\"\n";
	my @thisarray = @{$$tran_ref{$_}}; 
	my @sorted = sort {$a <=> $b} (@thisarray);
	foreach(my $ii = 0; $ii <= $#sorted; $ii+=2) {
	    print ACE "Source_Exons $sorted[$ii] $sorted[$ii+1]\n"
	}

    }

}





