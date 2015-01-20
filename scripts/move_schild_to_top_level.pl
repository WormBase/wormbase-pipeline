#!/usr/bin/env perl
#
# move_schild_to_topleve.pl
#
# In a typical Acedb database containing a multi-level assembly
# (e.g. chromsome, superlink and contig for C. elegans), features
# will be stored at a mixture of the multiple levels. This script
# will move features from non-toplevel layers (superlink, clone)
# to the top level (chromosome)
#
# In principle, the script can be used to project any SChild feature,
# but by default only CDS, Transcript and Pseudogene are projected
# (since this is needed for EMBL dumping)
#

use strict;
use Ace;
use Getopt::Long;
use lib $ENV{CVS_DIR};

use Coords_converter;

my (
  $wormbase,
  $species,
  $store,
  $debug,
  $test,
  $database, 
  $embl_feat_only,
  $acefile,
  $acefh,
  $noload);

&GetOptions(
  "species=s"    => \$species,
  "store=s"      => \$store,
  "debug=s"      => \$debug,
  "test"         => \$test,
  "database=s"   => \$database,
  "acefile=s"    => \$acefile,
  "emblfeatonly" => \$embl_feat_only,
  "noload"       => \$noload,);

my (@schild_list) = ("CDS_child", "Transcript", "Pseudogene");

############################
# recreate configuration   #
############################
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("cant restore wormbase from $store\n"); 
}
else { 
  $wormbase = Wormbase->new(-debug => $debug, 
                            -test => $test,
                            -organism => $species,
      );
}

$database = $wormbase->autoace if not defined $database;

if (defined $acefile) {
  open($acefh, ">$acefile") or die "Could not open $acefile for writing\n";
} else {
  $acefh = \*STDOUT;
}

my $full_name = $wormbase->full_name;

my $cc = Coords_converter->invoke( $database, undef, $wormbase );
my $db = Ace->connect(-path => $database) or die "Could not connect to $database\n";

my (%tl_seqs, %good_method_cache, %embl_method_cache);

my $iter = $db->fetch_many(-query => "Find Method EMBL_feature");
while(my $meth = $iter->next) {
  $embl_method_cache{$meth->name} = 1;
}

my $method_clause = "(" . join(" OR ", map { "Method = \"$_\"" } sort keys %embl_method_cache) . ")"; 

foreach my $class ("CDS", "Transcript", "Pseudogene") {
  my $query =  "Find $class Species = \"$full_name\" AND $method_clause";

  my $iter = $db->fetch_many(-query => "$query");
  while(my $obj = $iter->next()) {    
    $good_method_cache{$obj->name} = 1;
  }
}

foreach my $schild (@schild_list) {
  my (%add_seqs, %del_seqs);
  
  my $seqs_it = $db->fetch_many(-query => "Find Sequence $schild AND Species = \"$full_name\"");
  while(my $seq = $seqs_it->next()) {
    foreach my $row ($seq->at("SMap.S_child.${schild}")) {
      if ($embl_feat_only) {
        next if not $good_method_cache{$row->name};
      }
      if ($row->right and $row->right->right) {
        
        my ($sname, $start, $end) = ($seq->name, $row->right->name, $row->right->right->name);
        my ($mapped_seq, $mapped_start, $mapped_end) = $cc->LocateSpanUp($sname, $start, $end);
        
        $add_seqs{$mapped_seq}->{$row->name} = [$mapped_start, $mapped_end];
        $del_seqs{$sname}->{$row->name} = [$start, $end];
        $tl_seqs{$mapped_seq} = 1;
      }
    }
  }

  foreach my $seq (keys %add_seqs) {
    foreach my $k (sort keys %{$add_seqs{$seq}}) {
      my ($st, $en) = @{$add_seqs{$seq}->{$k}};
      print $acefh "\nSequence : \"$seq\"\n";
      print $acefh "$schild \"$k\" $st $en\n";
    }

  }
}

if ($acefile) {
  close($acefh) or die "Could not close $acefile after writing\n";
  $wormbase->load_to_database($database, $acefile, "move_to_toplevel") unless $noload;
}


exit(0);
