#!/usr/bin/env perl
#
# restructure_miRNAs
#
# Restructures miRNA-type transcript features into something resembling the correct hierarchy, 
# As required by AGR, the structure is:
#   Gene => miRNA_primary_transcript => exon
#   Gene => pre_miRNA => miRNA (no exon)
#
# Note that this is only applicable to GFF3; in GFF2 mode, the file is left as is

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use Ace;
use strict;

######################################
# variables and command-line options # 
######################################

my ( $debug, $test, $store, $wormbase,$database);
my ($datavase, $species, $gff3, $infile, $outfile, %var, $changed_lines, $removed_lines);

GetOptions (
  "debug=s"    => \$debug,
  "test"       => \$test,
  "store:s"    => \$store,
  "species:s"  => \$species,
  "database:s" => \$database,
  "gff3"       => \$gff3,
  "infile:s"   => \$infile,
  "outfile:s"  => \$outfile,
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
my $sp_full_name = $wormbase->full_name;

my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}

my ($remove_exons_h, $mirna_parent_h) = &get_mirna_data();


open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  

while (<$gff_in_fh>) {

  if (/^\#/ or not $gff3) {
    print $gff_out_fh $_;
    next;
  }


  # - Remove exon children of miRNA and pre_miRNA features
  # - Change parent of miRNA feature to point to pre_miRNA of same gene (if exists)

  chomp;
  my @l = split(/\t/, $_);

  if ($l[2] ne 'exon' and $l[2] ne 'miRNA') {
    print $gff_out_fh join("\t", @l), "\n"; next;
  }

  if ($l[2] eq 'exon') {
    my ($parent) = $l[8] =~ /Parent=Transcript:([^;]+)/;
    if (not exists $remove_exons_h->{$parent}) {
      print $gff_out_fh join("\t", @l), "\n";
    } else {
      $removed_lines++;
    }
    next;
  } elsif ($l[2] eq 'miRNA') {
    my ($id) = $l[8] =~ /ID=Transcript:([^;]+)/;
    my ($old_parent) = $l[8] =~ /Parent=([^;]+)/;
    if (exists $mirna_parent_h->{$id}) {
      my $new_parent = "Transcript:" . $mirna_parent_h->{$id};

      $l[8] =~ s/Parent=$old_parent/Parent=$new_parent/;
      $changed_lines++;
    }
    print $gff_out_fh join("\t", @l), "\n";
  } else {
    $log->log_and_die("Unexpected line: " . join("\t", @l));
  }
}
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n");

$log->write_to("Finished processing : $changed_lines lines modified, $removed_lines removed\n");
$log->mail();
exit(0);


##############################################################
#
# Subroutines
#
##############################################################

sub get_mirna_data {
  
  my (%genes);

  my $table = $wormbase->table_maker_query($database, &write_tm_def_file);
  
  while(<$table>) {
    chomp;
    s/\"//g;
    next if (not defined $_);
    next if (/acedb/ or /\/\//);
    my ($gene_id, $trans_id, $method) = split(/\s+/,$_);
    next if not defined $method;

    $genes{$gene_id}->{$method}->{$trans_id} = 1;
  } 
  close($table);

  my (%remove_exons, %mirna_parents);

  foreach my $g (keys %genes) {
    my (@mirna, @pre_mirna, @prim_mirna);

    @mirna = keys %{$genes{$g}->{miRNA}} if exists $genes{$g}->{miRNA};
    @pre_mirna = keys %{$genes{$g}->{pre_miRNA}} if exists $genes{$g}->{pre_miRNA};
    @prim_mirna = keys %{$genes{$g}->{miRNA_primary_transcript}} if exists $genes{$g}->{miRNA_primary_transcript};
    
    map { $remove_exons{$_} = 1 } @mirna;
    map { $remove_exons{$_} = 1 } @pre_mirna;
    
    if (scalar(@pre_mirna) == 1) {
      my ($pre_mirna) =  @pre_mirna;
      
      foreach my $mirna (@mirna) {
        $mirna_parents{$mirna} = $pre_mirna;
      }
    }
  }

  return (\%remove_exons, \%mirna_parents);

}



sub write_tm_def_file {
  my $def = '/tmp/mirna_restructure.def';
  open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
  my $txt = <<END;
Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Gene 
From 1 
Condition Species = "$sp_full_name" AND Corresponding_transcript AND NOT Corresponding_CDS

Colonne 2 
Width 12 
Optional 
Visible 
Class 
Class Transcript 
From 1 
Tag Corresponding_transcript 

Colonne 3 
Width 12 
Mandatory 
Visible 
Class 
Class Method 
From 2 
Tag Method 
Condition "miRNA" OR "pre_miRNA" OR "miRNA_primary_transcript"


END

  print TMP $txt;
  close TMP;
  return $def;
}
