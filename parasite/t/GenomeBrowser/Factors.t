use Test::More;

use File::Temp qw/tempdir/;
use GenomeBrowser::ArrayExpressMetadata;
use GenomeBrowser::RnaseqerMetadata;
use GenomeBrowser::Factors;

sub studies_grouped_by_factor_as_expected {
    my (
        $rnaseqer_characteristics_per_study,
        $array_express_factors_per_study,
        $expected, $description
    ) = @_;
    my $root_dir = tempdir( CLEANUP => 1 );
    my $species  = "schistosoma_mansoni";
    my $assembly = "schisto_v7.2";
    my $mock_rm  = bless { metadata => { $assembly => $rnaseqer_characteristics_per_study }},
      'GenomeBrowser::RnaseqerMetadata';
    my $mock_aem = bless {
        primary_accession_to_factor_type => $array_express_factors_per_study,
        secondary_to_primary_accession   => {},
      },
      'GenomeBrowser::ArrayExpressMetadata';

    my $subject =
      GenomeBrowser::Factors->new( $root_dir, $species, $assembly, $mock_rm,
        $mock_aem );
    is_deeply( $subject, bless( $expected, 'GenomeBrowser::Factors' ),
        $description ) or diag explain $subject, $expected;
}

studies_grouped_by_factor_as_expected( {}, {}, [], "Null case" );
studies_grouped_by_factor_as_expected(
    { study => { run => { "type" => "value" } } },
    {}, ["type"], "One type one run" );
studies_grouped_by_factor_as_expected(
    { study => { run => { "type" => "value", "type_2" => "value_2" } } },
    {}, ["type","type_2"], "Two types one run" );
studies_grouped_by_factor_as_expected(
    { study => { run => { "type" => "value" , "synonym"=> "value of blacklisted type"} } },
    {}, ["type"], "Blacklisted type one run" );
studies_grouped_by_factor_as_expected(
    {
        study => {
            run   => { "type" => "value" },
            run_2 => { "type" => "value_2" }
        }
    },
    {},
    ["type"],
    "One type two runs"
);
studies_grouped_by_factor_as_expected(
    {
        study => {
            run   => { "type" => "value", "type_2"=>"value_2" },
            run_2 => { "type" => "value_3" , "type_2"=>"value_4" }
        }
    },
    {},
    ["type", "type_2"],
    "Two types two runs"
);
studies_grouped_by_factor_as_expected(
    {
        study => {
            run   => { "type" => "value", "type_2"=>"value_2" },
            run_2 => { "type" => "value_3" , "type_2"=>"value_4" }
        }
    },
    {"study"=>["type"]},
    ["type"],
    "Two types two runs AE override"
);
studies_grouped_by_factor_as_expected(
    {
        study => {
            run   => { "type" => "value", "type_2"=>"value_2" },
            run_2 => { "type" => "value_3" , "type_2"=>"value_4" }
        }
    },
    {"study"=>["type", "type 2 that got renamed by AE"]},
    ["type", "type_2"],
    "Two types two runs do not use AE if they renamed values"
);
studies_grouped_by_factor_as_expected(
    {
        study => {
            run   => { "type" => "value", "common_type"=>"common_value" },
            run_2 => { "type" => "value_2" , "common_type"=>"common_value" }
        }
    },
    {},
    ["type"],
    "One type two runs extra common type"
);
studies_grouped_by_factor_as_expected(
    {
        study => {
            run   => { "type" => "value" },
        },
        study_2 => {
            run_2   => { "type" => "value" },
        }
    },
    {},
    ["type"],
    "One type two studies"
);

studies_grouped_by_factor_as_expected(
    {
        study => {
            run   => { "type" => "value" },
            run_2 => { "type_2" => "value_2" }
        }
    },
    {},
    ["type", "type_2"],
    "Two types"
);
studies_grouped_by_factor_as_expected(
    {
        study => {
            run   => { "type" => "value" },
        },
        study_2 => {
            run_2   => { "type_2" => "value" },
        }
    },
    {},
    ["type", "type_2"],
    "Two types two studies"
);

studies_grouped_by_factor_as_expected(
    {
        study => {
            run   => { "type" => "value" , "common_type"=>"common_value"},
            run_2 => { "type_2" => "value_2" ,"common_type"=>"common_value" }
        }
    },
    {},
    ["type", "type_2"],
    "Two types common type"
);
studies_grouped_by_factor_as_expected(
    {
        study => {
            run   => { "type" => "value" , "common_type"=>"common_value"},
        },
        study_2 => {
            run_2   => { "type_2" => "value" , "common_type"=>"common_value"},
        }
    },
    {},
    ["common_type", "type", "type_2"],
    "Two types two studies common type across studies (Debatable)"
);

done_testing();
