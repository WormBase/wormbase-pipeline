use strict;
use warnings;
use ProductionMysql;
use Test::More;

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

#If you want to see all matches then uncomment:
# my $key;
# foreach $key (keys %core_db_to_biomart_names)
# {
#   print "$key -> $core_db_to_biomart_names{$key}\n";
# };

done_testing;
