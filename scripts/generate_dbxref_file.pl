#!/usr/bin/env perl
#
# generate_dbxrefs_file.pl
#
# Generates a table of xrefs that should be generally useful for
# Uniprot, ENA and Ensembl. Columns are:
#
# 1  Gene sequence name
# 2  WBGene id
# 3  CGC name
# 4  Transcript name (not CDS name for coding transcripts, but the Coding_transcript name)
# 5  Pep name
# 6  ENA clone accession
# 7  ENA protein_id
# 8  Uniprot accession
#
# 
#  Last updated on: $Date: 2014-11-13 13:48:43 $
#  Last updated by: $Author: klh $

use strict;
use Getopt::Long;
use Storable;

use lib $ENV{'CVS_DIR'};
use Wormbase;

my ($test,
    $debug,
    $store,
    $species,
    $database,
    $wormbase,
    $outfile,
    $no_header,
    );

GetOptions (
  "test"            => \$test,
  "debug=s"         => \$debug,
  "store:s"         => \$store,
  "species:s"       => \$species,
  "database:s"      => \$database,
  "outfile:s"       => \$outfile,
  "noheader"        => \$no_header,
    );


if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species,
                             -autoace  => $database,
      );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

$species = $wormbase->species;
my $tace = $wormbase->tace;
my $full_species_name = $wormbase->full_name;
my $wormbase_version = $wormbase->get_wormbase_version_name;
my $dbdir = ($database) ? $database : $wormbase->autoace;
$outfile = $wormbase->acefiles . "/DBXREFs.txt" if not defined $outfile;

my (%wbgene, %gene, %cds, %transcds,%clone2acc, $out_fh);

$wormbase->FetchData('clone2accession', \%clone2acc, "$dbdir/COMMON_DATA");

$log->write_to("Generating protein-coding table\n");

my $query = &generate_coding_query($full_species_name);
my $command = "Table-maker -p $query\nquit\n";

open (TACE, "echo '$command' | $tace $dbdir |");
while (<TACE>) {
  chomp; s/\"//g;

  my ($cds, $gene, $trans, $prot, $clone, $pid, $pid_version, $uniprot ) = split(/\t/, $_);

  next if $gene !~ /^WBGene/;

  $wbgene{$gene}->{cds}->{$cds} = 1;
  $wbgene{$gene}->{transcript}->{$trans} = 1;
  $transcds{$trans} = $cds;
  
  $prot =~ s/\S+://;
  $cds{$cds}->{protein} = $prot;
  if ($clone and $pid and $pid_version) {
    $cds{$cds}->{pid}->{"$clone:$pid:$pid_version"} = 1;
  }
  $cds{$cds}->{uniprot} = $uniprot if $uniprot;

}
close TACE;
unlink $query;

foreach my $class ('Transcript', 'Pseudogene') {
  $log->write_to("Generating non-coding $class table\n");

  $query = &generate_noncoding_query($full_species_name, $class);
  $command = "Table-maker -p $query\nquit\n";
  open (TACE, "echo '$command' | $tace $dbdir |");
  while (<TACE>) {
    chomp; s/\"//g;
    
    my ($trans, $gene, $parent ) = split(/\t/, $_);
    next if $gene !~ /^WBGene/;
            
    $wbgene{$gene}->{transcript}->{$trans} = 1;
    $wbgene{$gene}->{sequence} = $parent;
  }
  close TACE;
  unlink $query;
}  


$log->write_to("Generating gene table\n");

$query = &generate_gene_query($full_species_name);
$command = "Table-maker -p $query\nquit\n";

open (TACE, "echo '$command' | $tace $dbdir |");
while (<TACE>) {
  chomp; s/\"//g;
  
  my ($wbgene, $sequence_name, $cgc_name, $locus_tag) = split(/\t/, $_);
  next if $wbgene !~ /^WBGene/;

  $gene{$sequence_name}->{$wbgene} = 1;
  if ($cgc_name) {
    $wbgene{$wbgene}->{cgc} = $cgc_name;
  }
  if ($locus_tag) {
    $wbgene{$wbgene}->{locus_tag} = $locus_tag;
  }
}
close TACE;
unlink $query;

open($out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");

&write_header($out_fh) unless $no_header;

foreach my $g (sort keys %gene) {
  foreach my $wbgeneid (keys %{$gene{$g}}) {
    my $cgc_name = (exists $wbgene{$wbgeneid}->{cgc}) ? $wbgene{$wbgeneid}->{cgc} : ".";
    my $locus_tag = (exists $wbgene{$wbgeneid}->{locus_tag}) ? $wbgene{$wbgeneid}->{locus_tag} : ".";
    foreach my $trans (keys %{$wbgene{$wbgeneid}->{transcript}}) {
      my @pid_list;
      my ($cds, $pepid, $uniprot) = (".", ".", ".");
      
      if (exists $transcds{$trans}) {
        # coding
        $cds = $transcds{$trans};        
        $pepid = $cds{$cds}->{protein};
        $uniprot = $cds{$cds}->{uniprot}  if exists $cds{$cds}->{uniprot};
        
        if (exists $cds{$cds}->{pid}) {
          foreach my $str (keys %{$cds{$cds}->{pid}}) {
            my ($clone, $pid) = split(/:/, $str);
            push @pid_list, [$clone2acc{$clone}, $pid];
          }
        } else {
          @pid_list = (['.', '.']);
        }
      } else {
        my ($pid, $clone) = (".", ".");

        if (exists $wbgene{$wbgeneid}->{sequence}) {
          $clone = $wbgene{$wbgeneid}->{sequence};
          if (exists $clone2acc{$clone}) {
            $clone = $clone2acc{$clone};
          } else {
            $clone = "$clone:NOACC";
          }
        }
        @pid_list =  ([$clone, $pid]);
      }

      foreach my $pidpair (@pid_list) {
        my ($clone, $pid) = @$pidpair;
        
        print $out_fh join("\t", 
                           $g, 
                           $wbgeneid, 
                           $cgc_name, 
                           $trans,
                           $pepid,
                           $clone,
                           $locus_tag,
                           $pid, 
                           $uniprot), "\n";
      }
    }
  }
}

close($out_fh) or $log->log_and_die("Could not cleanly close output file\n");

$log->mail();
exit(0);

##########################################
sub write_header {
  my ($out_fh) = @_;

  my $header = <<"HERE";
//
// WormBase $full_species_name XREFs for $wormbase_version
//
// Columns (tab separated) are:
//    1. WormBase Gene sequence name
//    2. WormBase Gene accession
//    3. WormBase Gene CGC name
//    4. WormBase Transcript sequence name
//    5. WormPep protein accession
//    6. INSDC parent sequence accession
//    7. INSDC locus_tag id
//    8. INSDC protein_id
//    9. UniProt accession
//
// Missing or not applicable data (e.g. protein identifiers for non-coding RNAs) is denoted by a "."
//
HERE

  print $out_fh $header;

}

##########################################
sub generate_coding_query {
  my ($full_species) = @_;

  my $tmdef = "/tmp/cod_tmquery.$$.def";
  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  

  my $tablemaker_template = <<"EOF";


Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class CDS 
From 1 
Condition Method = "curated" AND Species = "$full_species"
 
Colonne 2 
Width 12 
Optional 
Visible 
Class 
Class Gene 
From 1 
Tag Gene 
 
Colonne 3 
Width 12 
Optional 
Visible 
Class 
Class Transcript 
From 1 
Tag Corresponding_transcript 
 
Colonne 4
Width 12 
Optional 
Visible 
Class 
Class Protein
From 1 
Tag Corresponding_protein

Colonne 5
Width 12 
Optional 
Visible 
Class 
Class Sequence 
From 1 
Tag Protein_id 
 
Colonne 6 
Width 12 
Optional 
Visible 
Text 
Right_of 5 
Tag  HERE  
 
Colonne 7 
Width 12 
Optional 
Visible 
Integer 
Right_of 6 
Tag  HERE  
 
Colonne 8 
Width 12 
Optional 
Hidden 
Class 
Class Database 
From 1 
Tag Database 
Condition UniProt
 
Colonne 9 
Width 12 
Optional 
Hidden 
Class 
Class Database_field 
Right_of 8 
Tag  HERE  
Condition UniProtAcc
 
Colonne 10 
Width 12 
Optional 
Visible 
Class 
Class Accession_number 
Right_of 9 
Tag  HERE  

EOF

  print $qfh $tablemaker_template;
  return $tmdef;


}


sub generate_noncoding_query {
  my ($full_species, $class) = @_;

  my $tmdef = "/tmp/nc_tmquery.$$.def";
  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  

  my $condition = "";
  if ($class eq 'Transcript') {
    $condition = "NOT Method = \"Coding_transcript\" AND NOT Method = \"history_transcript\" AND Species = \"$full_species\"";
  } elsif ($class eq 'Pseudogene') {
    $condition = "Method = \"Pseudogene\" AND Species = \"$full_species\"";
  } else {
    $log->log_and_die("Unrecognised non-coding class: $class\n");
  }

  my $tablemaker_template = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class $class
From 1 
Condition $condition 
 
Colonne 2 
Width 12 
Optional 
Visible 
Class 
Class Gene 
From 1 
Tag Gene 
 
Colonne 3 
Width 12 
Optional 
Visible 
Class 
Class Sequence 
From 1 
Tag Sequence 
 
EOF

  print $qfh $tablemaker_template;
  return $tmdef;

}



sub generate_gene_query {
  my ($full_species) = @_;

  my $tmdef = "/tmp/gene_tmquery.$$.def";
  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  

  my $condition = "";

  my $tablemaker_template = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Gene 
From 1 
Condition Species = "$full_species"
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class Gene_name 
From 1 
Tag Sequence_name 
 
Colonne 3 
Width 12 
Optional 
Visible 
Class 
Class Gene_name 
From 1 
Tag CGC_name 

Colonne 4 
Width 12 
Optional 
Hidden 
Class 
Class Database 
From 1 
Tag Database 
Condition NDB
 
Colonne 5 
Width 12 
Optional 
Hidden 
Class 
Class Database_field 
Right_of 4 
Tag  HERE  
Condition locus_tag
 
Colonne 6 
Width 12 
Optional 
Visible 
Class 
Class Accession_number 
Right_of 5 
Tag  HERE  
 
EOF

  print $qfh $tablemaker_template;
  return $tmdef;


}
__END__
