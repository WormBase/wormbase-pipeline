#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

use lib $ENV{CVS_DIR};

use Wormbase;

my $reg_conf;
&GetOptions(
  'reg_conf=s'      => \$reg_conf,
    ) or die "Bad options?";

my $WB_JBROWSE_BASE = 'https://wormbase.org/tools/genome/jbrowse-simple/?data=data%2F[WB_GENOME]';
my $PARASITE_JBROWSE_BASE = 'http://parasite.wormbase.org/jbrowse/index.html?data=%2Fjbrowse-data%2F[PARASITE_GENOME]%2Fdata';

die "You must supply a valid Registry file\n" if not defined $reg_conf or not -e $reg_conf;

my $wb = Wormbase->new();
my %accessors = ($wb->all_species_accessors);
my %wb_genomes;
foreach my $acc ($wb, values %accessors) {
  my $spec = $acc->long_name;
  $spec =~ s/\s+/_/; 
  my $prod_name = lc(join("_", $spec, $acc->ncbi_bioproject)); 
  $wb_genomes{$prod_name} = sprintf("%s_%s", $acc->gspecies_name, $acc->ncbi_bioproject);
}

Bio::EnsEMBL::Registry->load_all($reg_conf);

my $dbs = Bio::EnsEMBL::Registry->get_all_DBAdaptors( -group => "core" );

my %list;
for my $mode (qw/nematode flatworm/){
 foreach my $db (@$dbs) {

  my $div = $db->get_MetaContainer->get_division();
  next if $div !~ /Parasite/;

  my $classification = $db->get_MetaContainer->get_classification();
  my $class;
  
  if (grep { $_ eq 'Nematoda' } @$classification) {
    # Nematode
    next if $mode eq 'flatworm';

    $class = $db->get_MetaContainer->single_value_by_key("species.nematode_clade");
    die "species.nematode_clade missing for $db\n" if not defined $class;

  } elsif (grep { $_ eq 'Platyhelminthes' } @$classification) {
    # Platyhelminth
    next if $mode eq 'nematode';

    for(my $i=0; $i< @$classification; $i++) {
      if ($classification->[$i] eq 'Platyhelminthes') {
        $class = $classification->[$i-1];
        last;
      }
    }
    die "species.nematode_clade missing for $db\n" if not defined $class;
  } else {
    die "$db is neither nematode not Platyhelminth\n";
  }

  my $genus_species = $db->get_MetaContainer->get_scientific_name();
  my $prod_name = $db->get_MetaContainer->get_production_name();
  my $tax_id = $db->get_MetaContainer->get_taxonomy_id();

  my $provider = $db->get_MetaContainer->single_value_by_key("provider.name");
  my $provider_url = $db->get_MetaContainer->single_value_by_key("provider.url");
  my $bioproject = $db->get_MetaContainer->single_value_by_key("species.bioproject_id");
  my $url_name = $db->get_MetaContainer->single_value_by_key("species.url");
  my $strain = $db->get_MetaContainer->single_value_by_key("species.strain");
  $strain = "Not specified" if not $strain;

  my $gspecies_field = sprintf("=HYPERLINK(\"https://parasite.wormbase.org/%s\",\"%s\")", $url_name, $genus_species);
  my $tax_field = sprintf("=HYPERLINK(\"http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s\",\"%s\")", $tax_id, $tax_id);
  my $bp_field  = sprintf("=HYPERLINK(\"http://www.ncbi.nlm.nih.gov/bioproject/%s\",\"%s\")", $bioproject, $bioproject);
  my $provider_field = sprintf("=HYPERLINK(\"%s\",\"%s\")", $provider_url, $provider);

  my $parasite_jbrowse = $PARASITE_JBROWSE_BASE;
  $parasite_jbrowse =~ s/\[PARASITE_GENOME\]/${prod_name}/;
  $parasite_jbrowse = "=HYPERLINK(\"$parasite_jbrowse\", \"Link\")";

  my $clade =$mode eq 'nematode' && exists $wb_genomes{$prod_name} ?  do {
      my $gname = $wb_genomes{$prod_name};
      my $wb_jbrowse = $WB_JBROWSE_BASE;
      $wb_jbrowse =~ s/\[WB_GENOME\]/$gname/; 
      $wb_jbrowse = "=HYPERLINK(\"$wb_jbrowse\", \"Link\")";
      $wb_jbrowse
    } : ""; 
  push @{$list{$mode}}, [$gspecies_field, $class, $tax_field, $bp_field, $strain, $provider_field, $parasite_jbrowse, $clade];
 }
}
for my $mode (qw/nematode flatworm/){
  foreach my $l (sort { $a->[1] cmp $b->[1] or $a->[0] cmp $b->[0] } @{$list{$mode}}) {
    print join("\t", @$l), "\n";
  }
}
