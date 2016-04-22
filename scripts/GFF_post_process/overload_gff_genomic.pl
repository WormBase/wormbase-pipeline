#!/usr/bin/env perl
#
# overload_gff_genomic.pl
# 
# overloads Genomic_canonical lines with GenBank accesssions
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2014-08-27 21:50:10 $      

use strict;                                      
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use Ace;
use strict;

######################################
# variables and command-line options # 
######################################

my ( $debug, $test, $verbose, $store, $wormbase, $species, $database );
my ( $infile, $outfile, $gff3, %TF, $changed_lines, $added_lines);

GetOptions (
  "debug=s"    => \$debug,
  "test"       => \$test,
  "store:s"    => \$store,
  "species:s"  => \$species,
  "gff3"       => \$gff3,
  "infile:s"   => \$infile,
  "outfile:s"  => \$outfile,
  "database:s" => \$database,
    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( '-debug'   => $debug,
                             '-test'    => $test,
                             '-organism' => $species,
			     );
}

$database = $wormbase->autoace if not defined $database;
$species = $wormbase->species;

my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}

my %acc = $wormbase->FetchData('clone2accession');

open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  

while( <$gff_in_fh> ) {
  unless (/Genomic_canonical\tregion\t/ or
          /Genomic_canonical\tassembly_component\t/) {
    print $gff_out_fh $_;
    next;
  };
  chomp;
  my @f = split(/\t/, $_);
  
  my ($fid);
  if ($f[8] =~ /Sequence\s+\"(\S+)\"/) {
    $fid = $1;
  } elsif ($f[8] =~ /Sequence:([^;]+)/) {
    $fid = $1;
  } else {
    $log->log_and_die("Could not get sequence name from line: @f\n");
  }
  
  if (exists $acc{$fid}) {
    if ($gff3) {
      $f[8] .= ";Note=Clone:$fid,GenBank:$acc{$fid}";
    } else {
      $f[8] .= " ; Note \"Clone $fid\" ; Genbank \"$acc{$fid}\"";
    }
    $changed_lines++;
    # add in a duplicate line with Genbank source
    my $new_line_attr = ($gff3) ? "genbank=$acc{$fid}" : "genbank \"$acc{$fid}\"";
    print $gff_out_fh join("\t", $f[0], "Genbank", @f[2..7], $new_line_attr), "\n";
    $added_lines++;
  }
  print $gff_out_fh join("\t", @f), "\n";
}
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n");
$log->write_to("Finished processing : $changed_lines lines modified and $added_lines added\n");
$log->mail();
exit(0);


1;
