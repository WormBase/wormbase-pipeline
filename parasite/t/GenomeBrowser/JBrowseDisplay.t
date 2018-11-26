use Test::More;
use List::Util qw/pairmap/;
use File::Temp qw/tempdir/;
use GenomeBrowser::JBrowseDisplay;

sub studies_grouped_by_factor_as_expected {
    my ( $studies, $expected, $description ) = @_;
    my $actual =
      GenomeBrowser::JBrowseDisplay::column_headers_for_studies( @{$studies} );
    is_deeply( $actual, $expected, $description )
      or diag explain $actual, $expected;
}

studies_grouped_by_factor_as_expected(
    [
        {
            'runs' => [
                {
                    'characteristics' =>
                      { 'type_2' => 'value_2', 'type' => 'value' },
                    'run_id' => 'run'
                }
            ],
            'study_id' => 'study'
        }
    ],
    [ 'type', 'type_2' ],
    'Two types one run'
);
studies_grouped_by_factor_as_expected(
    [
        {
            'study_id' => 'study',
            'runs'     => [
                {
                    'run_id'          => 'run',
                    'characteristics' => {
                        'synonym' => 'value of blacklisted type',
                        'type'    => 'value'
                    }
                }
            ]
        }
    ],
    ['type'],
    'Blacklisted type one run'
);
studies_grouped_by_factor_as_expected(
    [
        {
            'runs' => [
                {
                    'run_id'          => 'run_2',
                    'characteristics' => { 'type' => 'value_2' }
                },
                {
                    'run_id'          => 'run',
                    'characteristics' => { 'type' => 'value' }
                }
            ],
            'study_id' => 'study'
        }
    ],
    ['type'],
    'One type two runs'
);
studies_grouped_by_factor_as_expected(
    [
        {
            'study_id' => 'study',
            'runs'     => [
                {
                    'characteristics' =>
                      { 'type' => 'value', 'type_2' => 'value_2' },
                    'run_id' => 'run'
                },
                {
                    'characteristics' =>
                      { 'type_2' => 'value_4', 'type' => 'value_3' },
                    'run_id' => 'run_2'
                }
            ]
        }
    ],
    [ 'type', 'type_2' ],
    'Two types two runs'
);
studies_grouped_by_factor_as_expected(
    [
        {
            'study_id' => 'study',
            'runs'     => [
                {
                    'characteristics' =>
                      { 'common_type' => 'common_value', 'type' => 'value_2' },
                    'run_id' => 'run_2'
                },
                {
                    'run_id' => 'run',
                    'characteristics' =>
                      { 'common_type' => 'common_value', 'type' => 'value' }
                }
            ]
        }
    ],
    ['type'],
    'One type two runs extra common type'
);
studies_grouped_by_factor_as_expected(
    [
        {
            'study_id' => 'study',
            'runs'     => [
                {
                    'characteristics' => { 'type' => 'value' },
                    'run_id'          => 'run'
                }
            ]
        },
        {
            'runs' => [
                {
                    'run_id'          => 'run_2',
                    'characteristics' => { 'type' => 'value' }
                }
            ],
            'study_id' => 'study_2'
        }
    ],
    ['type'],
    'One type two studies'
);
studies_grouped_by_factor_as_expected(
    [
        {
            'study_id' => 'study',
            'runs'     => [
                {
                    'run_id'          => 'run_2',
                    'characteristics' => { 'type_2' => 'value_2' }
                },
                {
                    'characteristics' => { 'type' => 'value' },
                    'run_id'          => 'run'
                }
            ]
        }
    ],
    [ 'type', 'type_2' ],
    'Two types'
);
studies_grouped_by_factor_as_expected(
    [
        {
            'study_id' => 'study_2',
            'runs'     => [
                {
                    'characteristics' => { 'type_2' => 'value' },
                    'run_id'          => 'run_2'
                }
            ]
        },
        {
            'study_id' => 'study',
            'runs'     => [
                {
                    'characteristics' => { 'type' => 'value' },
                    'run_id'          => 'run'
                }
            ]
        }
    ],
    [ 'type', 'type_2' ],
    'Two types two studies'
);
studies_grouped_by_factor_as_expected(
    [
        {
            'runs' => [
                {
                    'run_id'          => 'run_2',
                    'characteristics' => {
                        'type_2'      => 'value_2',
                        'common_type' => 'common_value'
                    }
                },
                {
                    'characteristics' =>
                      { 'type' => 'value', 'common_type' => 'common_value' },
                    'run_id' => 'run'
                }
            ],
            'study_id' => 'study'
        }
    ],
    [ 'type', 'type_2' ],
    'Two types common type'
);
studies_grouped_by_factor_as_expected(
    [
        {
            'runs' => [
                {
                    'characteristics' =>
                      { 'type' => 'value', 'common_type' => 'common_value' },
                    'run_id' => 'run'
                }
            ],
            'study_id' => 'study'
        },
        {
            'runs' => [
                {
                    'run_id' => 'run_2',
                    'characteristics' =>
                      { 'type_2' => 'value', 'common_type' => 'common_value' }
                }
            ],
            'study_id' => 'study_2'
        }
    ],
    [ 'common_type', 'type', 'type_2' ],
    'Two types two studies common type across studies (Debatable)'
);
done_testing();
