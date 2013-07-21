#!/usr/bin/env perl
#
# overload_gff_sage.pl
#
# Overloads the SAGE_tag GFF lines with additional info
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2013-07-21 11:07:59 $      
#

use lib $ENV{CVS_DIR};
use strict;
use Ace;
use Storable;
use Wormbase;
use Log_files;
use Getopt::Long;

my ($help, $debug, $test, $store, $wormbase );
my ($gff3, $infile, $outfile, $changed_lines);

GetOptions (
  "debug=s"      => \$debug,
  "test"         => \$test,
  "store:s"      => \$store,
  "gff3"         => \$gff3,
  "infile:s"     => \$infile,
  "outfile:s"    => \$outfile,
    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

my $log = Log_files->make_build_log($wormbase);
if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}

open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  


$log->write_to("Looking up SAGE ids for $infile\n");
my (%sage_info);
open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
while(<$gff_in_fh>) {
  /(SAGE:[a-z]+)/ and $sage_info{$1} = {};
}

$log->write_to("Reading SAGE info from DB from objects found in $infile...\n");
&markup(\%sage_info);

$log->write_to("Processing $infile...\n");
open($gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
while(<$gff_in_fh>) {
  if( /^\#/ or 
      not /SAGE_tag/) {
    print $gff_out_fh $_;
    next;
  }
  
  my ($ref,$src,$met,$start,$stop,$scr,$strand,$phase,$att) = split(/\t+/, $_);
  
  my ($sequence) = $att =~ /(SAGE:[a-z]+)/;
  
  # infer strandedness from targets
  if ($gff3) {
    ($strand) = ($att =~ /Target=\S+ \d+ \d+ (\S)/);
  } else {
    if ($att =~ /Target\s+\S+\s+(\d+)\s+(\d+)/) {
      if ($2 < $1) {
        $strand = '-';
      }
      else {
        $strand = '+';
      }
    }
  }
  $strand =~ tr/\./+/;
  
  my @atts = (["Sequence", $sequence]);
  
  if (ref($sage_info{$sequence})) {
    push @atts, @{$sage_info{$sequence}};
  }
  
  my $old_att = $att;
  if ($gff3) {
    $att = join(";", map { join("=", @$_) } @atts);
  } else {
    $att = join(" ; ", map { "$_->[0] \"$_->[1]\"" } @atts);
  }    
  print $gff_out_fh join("\t",$ref,$src,$met,$start,$stop,$scr,$strand,$phase,$att), "\n";
  
  $changed_lines++ if $att ne $old_att;
}
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n");
$log->write_to("Finished processing : $changed_lines lines modified\n");
$log->mail;
exit(0);


####################################################
sub markup {
  my ($hashr) = @_;

  my $db = Ace->connect( -path => $wormbase->autoace );

  foreach my $tag (keys %$hashr) {
    my $r = $db->fetch(SAGE_tag => $tag);
    next if not defined $r;

    my @genes    = ($r->Gene,$r->Transcript,$r->Pseudogene,$r->Predicted_CDS);
    my @unambig  = grep { scalar $_->get('Unambiguously_mapped')} @genes;
    my @three_p  = grep { scalar $_->get('Most_three_prime') } @genes;
    if (!@unambig && !@three_p) {
      my @unique = grep {$_->class eq 'Gene' || $_->class eq 'Pseudogene'} @genes;
      @unambig = @unique if @unique == 1;
    }
    
    my @count = $r->Results;
    my $count;
    for my $exp (@count) {
      my $cnt = $exp->get('Frequency')->right;
      $count += $cnt;
    }
    for (@unambig, @three_p) {
      my $name = eval{ $_->CGC_name} || eval{$_->Sequence_name} || $_->name;
      $_ = [$_->class,$name];
    }

    $hashr->{$tag} = [ ["count", $count], @unambig, @three_p ];
  }                       

  $db->close();
}
