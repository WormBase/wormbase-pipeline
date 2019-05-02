use strict;
use warnings;
use ProductionMysql;
use Test::More;
use LWP::Simple;

subtest "Species regex ok" => sub {
  ok($ProductionMysql::GOLDEN_SPECIES_REGEX_MATCH && $ProductionMysql::GOLDEN_SPECIES_REGEX_REPLACEMENT, "Export the regex strings");

  my %core_db_to_biomart_names = map {$_ => ProductionMysql::core_db_to_biomart_name($_) } ProductionMysql->staging->core_databases;

  my %seen;
  my @dups = grep {$seen{$_} > 1} map {$seen{$_}++; $_} values %core_db_to_biomart_names;

  ok(not(@dups), "No duplicate biomart names") or diag explain join("\n\t", 
    "",
    map {
    my $biomart_name = $_;
    join("\t", $biomart_name, grep {$core_db_to_biomart_names{$_} eq $biomart_name} keys %core_db_to_biomart_names)
  } @dups);

  my @too_long_names = grep {length $_ > 17} values %core_db_to_biomart_names;
  ok(not(@too_long_names), "Biomart names under 18 chars") or diag explain join("\n\t", 
    "",
    map {
    my $biomart_name = $_;
    join("\t", $biomart_name, length $biomart_name,  grep {$core_db_to_biomart_names{$_} eq $biomart_name} keys %core_db_to_biomart_names)
  } @too_long_names);
};

subtest "Taxon tree JS" => sub {
  my @biomart_names_for_core_dbs = sort map {ProductionMysql::core_db_to_biomart_name($_) } ProductionMysql->staging->core_databases;

  my $taxon_tree_js = get "http://test.parasite.wormbase.org/taxon_tree_data.js";
  ok ($taxon_tree_js, "Site serving taxon tree");
  my @biomart_names_tree = sort $taxon_tree_js =~ m{"biomart"\s*:\s*"(\S+)"}g;
  is_deeply(\@biomart_names_for_core_dbs, \@biomart_names_tree);
};
done_testing;
