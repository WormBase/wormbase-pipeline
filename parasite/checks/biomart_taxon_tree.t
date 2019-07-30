use strict;
use warnings;
use ProductionMysql;
use Test::More;
use LWP::Simple;
use File::Slurp qw/read_file/;

my ($doc_path) = @ARGV;
$doc_path //= "http://test.parasite.wormbase.org/taxon_tree_data.js";
my $doc = -f $doc_path ? read_file $doc_path : get $doc_path;

like($doc, qr/^taxonTreeData =/, "Doc with taxonTreeData");

subtest "All species mentioned in the doc" => sub {
  ok ($doc =~ /$_/, "Doc has: $_") for ProductionMysql->staging->species;
};

subtest "All biomart names mentioned in the doc" => sub {
  ok ($doc =~ /$_/, "Doc has: $_") for sort map {ProductionMysql::core_db_to_biomart_name($_) } ProductionMysql->staging->core_databases; 
};

done_testing;
