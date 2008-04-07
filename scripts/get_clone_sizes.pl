#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  get_clone_sizes.pl
#
#        USAGE:  ./get_clone_sizes.pl 
#
#  DESCRIPTION:  prints the clone sizes of briggsae based on the chromosome
#                subsequences. Used by update_Common_data.pl for C.briggsae
#
#      OPTIONS:  -store BRIGGSAE -debug USERNAME -test -species Briggsae
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:   (Michael Han), <mh6@sanger.ac.uk>
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  02/18/08 10:48:17 GMT
#     REVISION:  $revision:$
#===============================================================================
use lib $ENV{'CVS_DIR'};
use Ace;
use Wormbase;
use Getopt::Long;
use strict;

my ($store,$debug,$test,$species);

GetOptions(
	'debug=s'   => \$debug,
	'store=s'   => \$store,
	'test'      => \$test,
	'species=s' => \$species,
)||die(@!);

my $wormbase;

if ($store) {
    $wormbase = Storable::retrieve($store)
      or croak("cannot restore wormbase from $store");
} else { 
    $wormbase = Wormbase->new(-debug => $debug, -test => $test,-organism => $species);
}

my $db=Ace->connect(-path => $wormbase->autoace);
my @chromosomes=$wormbase->get_chromosome_names(-prefix => 1);

# builds a table:
# cb25.NA_364     2750
foreach my $chromosome(@chromosomes){
	 my $chr=$db->fetch(Sequence => $chromosome);
	 map {printf("$_\t%s\n",$_->right->right - $_->right+1)}
	     $chr->at("Structure.Subsequence");
}
$db->close;
