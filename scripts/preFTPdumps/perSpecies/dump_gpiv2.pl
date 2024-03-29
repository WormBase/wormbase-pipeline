#!/usr/bin/env perl
#
# exporter to dump gene / transcript / protein info as GPI file
#   specs: http://www.geneontology.org/page/gene-product-information-gpi-format
#   version 2: https://github.com/geneontology/go-annotation/blob/master/specs/gpad-gpi-2-0.md
#
# usage:
#   perl dump_gpi.pl -species elegans

use strict;
use lib $ENV{CVS_DIR};
use Wormbase;
use Log_files;
use Getopt::Long;
use Storable;
use Ace;
use IO::File;

my ($debug,$test,$species,$store,$output,$database);

GetOptions(
   'debug=s'   => \$debug,   # send log mails only to one person
   'species=s' => \$species, # specify the species to run on
   'test'      => \$test,    # use the test database instead of the live one
   'store=s'   => \$store,   # pass a storable (for the build)
   'outfile=s' => \$output,  # write somewhere else and not to REPORTS/
   'database=s'=> \$database,# specify a different database than BUILD/$species
)||die(@!);

my $wormbase;
if ($store) {
  $wormbase = retrieve($store) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

$output||=$wormbase->reports . '/' . $wormbase->species . '.gene_product_info.gpi2';
my $outfile = IO::File->new($output,'w')||$log->log_and_die(@!);

$log->write_to("creating a GPI file at $output for ${\$wormbase->long_name}\n");

my $db = Ace->connect(-path => ($database||$wormbase->autoace));

my $genesh = $db->fetch_many(-query => 'Find Gene;Species="'.$wormbase->long_name.'";Live')
            or $log->log_and_die(Ace->error);

print $outfile "!gpi-version: 2.0\n";
print $outfile "!namespace: WB\n";

my (%genes, %coding_trans, %nc_trans, %proteins);

while (my $g = $genesh->next){
    next unless $g->name =~ /^WBGene\d+$/;
  # Gene block
  my ($desc) = $g->Gene_class ? $g->Gene_class->Description : '';
  
  print STDERR "processing $g\n" if $debug;
  push @{$genes{$g}}, [
    "WB:$g", # 1 id
    $g->Public_name, # 2 symbol
    $desc, # 3 Name
    join('|',$g->Other_name), # 4 synonyms
    'gene', # 5 SO
    'NCBITaxon:' . $g->Species->NCBITaxonomyID, # 6 taxon
    '', # 7 gene
    '', # 8 parent _protein
    '', # 9 complex members
    join("|", &get_xrefs($g, 'gene')), # 10 db_xrefs
    '', # 11 properties
    ]; 
  
  foreach my $c($g->Corresponding_CDS){
    # Transcript/CDS block
    push @{$coding_trans{$c}}, [
      "WB:$c",
      $g->Public_name,
      $desc,
      join('|',$g->Other_name),
      'transcript',
      'NCBITaxon:' . $c->Species->NCBITaxonomyID, # 6
      "WB:$g",
      '',
      '',
      '',
      ''];
    
    foreach my $p ($c->Corresponding_protein) {
      push @{$proteins{$p}}, [
        "WB:$p",
        uc($g->Public_name),
        $desc,
        join('|',map { uc($_) } $g->Other_name),
        'protein',
        'NCBITaxon:' .$p->Species->NCBITaxonomyID,
        "WB:$g",
	'',
	'',
        join('|', &get_xrefs($p, 'protein')),
        '',
        ];
    }
  }
  
  foreach my $t($g->Corresponding_transcript){
    next if "${\$t->Method}" eq 'Coding_transcript';
    # ncRNA transcript block
    push @{$nc_trans{$t}}, [
      "WB:$t",
      $g->Public_name,
      $desc,
      join('|',$g->Other_name),
      $t->Method->GFF_SO->SO_name,
      'NCBITaxon:' . $t->Species->NCBITaxonomyID,
      "WB:$g",
      '',
      '',
      join("|". &get_xrefs($t, 'transcript')), 
      '',
      ];
  }
}

foreach my $g (sort keys %genes) {
  my @g = @{$genes{$g}};
  $log->log_and_die("ERROR: multiple gene entries with id $g\n") if scalar(@g) != 1;

  check_for_missing_data($log, $g[0]);
  print $outfile join("\t", @{$g[0]}), "\n";
}

foreach my $ct (sort keys %coding_trans) {
  my @ct = @{$coding_trans{$ct}};
  $log->log_and_die("ERROR: multiple gene entries with id $ct\n") if scalar(@ct) != 1;

  check_for_missing_data($log, $ct[0]);
  print $outfile join("\t", @{$ct[0]}), "\n";
}

foreach my $nct (sort keys %nc_trans) {
  my @nct = @{$nc_trans{$nct}};
  $log->log_and_die("ERROR: multiple gene entries with id $nct\n") if scalar(@nct) != 1;

  check_for_missing_data($log, $nct[0]);
  print $outfile join("\t", @{$nct[0]}), "\n";
}

foreach my $p (sort keys %proteins) {
  my @p = @{$proteins{$p}};

  if (scalar(@p) > 1) {
    my @all_t;
    foreach my $op (@p) {
      push @all_t, $op->[7];
    }
    $p[0]->[7] = join("|", @all_t);
  }
  check_for_missing_data($log, $p[0]);
  print $outfile join("\t", @{$p[0]}), "\n";
}

$log->mail;

# returns status 13 when finishing, so try to clean up
$outfile->close;
$db->close;

#############################################################
sub check_for_missing_data {
    my ($log, $columns) = @_;

    for my $ix (0, 1, 4, 5) {
	my $nr = $ix + 1;
	$log->log_and_die("Required field missing in column $nr:\n\t" .
			  join("\t", @$columns) . "\n")
	    unless length $columns->[$ix] > 0;
    }

    return;
}


sub get_xrefs {
  my ($obj, $type) = @_;

  my (%accs_by_source, @ret);
  
  foreach my $db ($obj->Database) {
    my $db_name = lc($db->name);
    if ($db_name eq "swissprot" or $db_name eq 'trembl' or $db_name eq 'uniprot_gcrp' or $db_name eq 'rnacentral') {
      foreach my $subtype ($db->col) {
        foreach my $acc ($subtype->col) {
          push @{$accs_by_source{$db_name}}, $acc->name;
        }
      }
    }
  }

  if ($type eq 'gene') {
    # for protein-coding genes, get GCRP xref(s) if one exists; if not, get SwissProt/Trembl/RNAcentral if there is only one

    if (exists $accs_by_source{uniprot_gcrp}) {
      @ret = map { "UniProtKB:$_" } @{$accs_by_source{uniprot_gcrp}};
    } else {
      my @list;
      if (exists $accs_by_source{swissprot}) {
        push @list, map { "UniProtKB:$_" } @{$accs_by_source{swissprot}};
      }
      if (exists $accs_by_source{trembl}) {
        push @list, map { "UniProtKB:$_" } @{$accs_by_source{trembl}};
      }
      if (exists $accs_by_source{rnacentral}) {
        push @list, map { "RNAcentral:$_" } @{$accs_by_source{rnacentral}};
      }
      
      if (scalar(@list) == 1) {
        @ret = @list;
      }
    }
  } elsif ($type eq 'protein') {
    # transcript or protein
    if (exists $accs_by_source{swissprot}) {
        push @ret, map { "UniProtKB:$_" } @{$accs_by_source{swissprot}};
    }
    if (exists $accs_by_source{trembl}) {
      push @ret, map { "UniProtKB:$_" } @{$accs_by_source{trembl}};
    }
  } else {
    # transcripts and ncRNAs; only going to add RNAcentral xrefs for these
    if (exists $accs_by_source{rnacentral}) {
      push @ret, map { "RNAcentral:$_" } @{$accs_by_source{rnacentral}};
    }
  }


  return @ret;
}

exit(0);
