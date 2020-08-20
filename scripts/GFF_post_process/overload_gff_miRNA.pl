#!/usr/bin/env perl
#
# overload_gff_miRNA.pl
#
# Overloads miRNA lines with type info
#
# Last updated by: $Author: pad $
# Last updated on: $Date: 2020-07-15 18:50:10 $

use lib $ENV{CVS_DIR};
use Wormbase;
use Storable;
use IO::File;
use Getopt::Long;
use strict;
use Ace;

my ($debug,$test,$species,$store,$wormbase,$database);
my ($gff3,$infile,$outfile, $changed_lines);

GetOptions(
  'debug=s'       => \$debug,
  'test'          => \$test,
  'species:s'     => \$species,
  'store:s'       => \$store,
  'infile:s'      => \$infile,
  'outfile:s'     => \$outfile,
  'gff3'          => \$gff3,
  'database:s'    => \$database,
    );

if ($store) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species
			     );
}

$database = $wormbase->autoace if not defined $database;
my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}

open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  

my $db = Ace->connect(-path => $database);



# fill some reference hashes
#my ($miRNA_type) = &get_miRNA_data();

while (<$gff_in_fh>){
#    print $_;
    next unless (/^\S+/);
    unless ((/WormBase\s+miRNA\s+/)||(/miRNA_mature\s+miRNA\s+/)){
        print $gff_out_fh $_;
        next;
    }
    chomp;
    my $Transcript_name;
    #Name=Y71F9AL.23;
    #  if (m/Name\=(\S+)?\;/){
    if ($gff3) {
        if (m/Name\=(\w+\.\w+)\;?/){
            $Transcript_name = $1
        }
    }
    else {
        if (m/Transcript\s+\"(\S+)\"/){
            $Transcript_name = $1
        } 
    }
    my $transobj = $db->fetch(Transcript => $Transcript_name);
    print "Transcript : $Transcript_name"; 
    my $trans_type =  $transobj->Transcript->right->name;
    print " type $trans_type\n";
    if ($gff3) {
        print $gff_out_fh $_.";type=$trans_type\n";
    }
    else {
        print $gff_out_fh $_." ; type \"$trans_type\"\n";
    }
    $changed_lines++;
}
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n");
$log->write_to("Finished processing : $changed_lines lines modified\n");
$log->mail();
exit(0);

