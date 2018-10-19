use Test::More;
use Test::MockModule;

use File::Temp qw/tempdir/;
use JSON;
use GenomeBrowser::Resources::ArrayExpressMetadata;

my $module = new Test::MockModule('GenomeBrowser::Resources::LocallyCachedResource');
$module->mock('get_text', do {local $/; <DATA>});

my $root_dir = tempdir( CLEANUP => 1 );
my $species = "schistosoma_mansoni";

my $subject = GenomeBrowser::Resources::ArrayExpressMetadata->new($root_dir, $species);

my $expected_factor_types = [
  "organism_part",
  "sample_description",
  "sex"
];

is_deeply($subject->factor_types("E-ERAD-516"),
  $expected_factor_types,
  "Can get factor types from AE (primary acc.)"
) or diag explain $subject;

is_deeply($subject->factor_types("ERP016356"),
  $expected_factor_types,
  "Can get factor types from AE (secondary acc.)"
);

done_testing();

__DATA__
{
    "experiments": {
        "api-version": 3,
        "api-revision": "091015",
        "version": 3.0,
        "revision": "091015",
        "total": 45,
        "total-samples": 1145,
        "total-assays": 1578,
        "experiment": [{
                "id": 569232,
                "accession": "E-ERAD-516",
                "secondaryaccession": ["ERP016356"],
                "name": "S  mansoni male and female adult worm reproductive organs pre  and post pairing",
                "releasedate": "2016-07-05",
                "lastupdatedate": "2018-03-14",
                "organism": ["Schistosoma mansoni"],
                "experimenttype": ["RNA-seq of coding RNA"],
                "description": [{
                    "id": null,
                    "text": "Study description:This data is part of a pre-publication release. For information on the proper use of pre-publication data shared by the Wellcome Trust Sanger Institute (including details of any publication moratoria), please see <a href=\"http:\/\/www.sanger.ac.uk\/datasharing\/\" target=\"_blank\">http:\/\/www.sanger.ac.uk\/datasharing\/<\/a> Schistosomes are blood dwelling digenean trematodes that mature as adults in the intestinal or urinary veins, predominantly of mammals. Adult female Schistosomes produce eggs of which approximately 40% fail to pass into faeces or urine (species-dependent) but are dispersed by the blood stream into different organs where they provoke severe inflammation. This parasitic infection is known as schistosomiasis and considered by the WHO as the second most socioeconomically devastating parasitic disease, next only to malaria, with hundreds of millions infected worldwide. As the eggs represent the causative pathogenic agents, the understanding of egg forming processes, and therefore of the schistosomal reproductive biology in general, is of fundamental interest. Additionally the differences between immature (pre-mating) and mature (post-mating) adult worms are of significant interest; adults only reach maturity upon mating. Trasncriptomic sequencing of gonad-specific cellular material will help to unravel signal transduction cascades involved in e.g. gametogenesis and\/or vitellogenesis. This project aims to characterize the transcriptome profiles of immature and mature adult ovaries and adult testes."
                }],
                "provider": [{
                    "contact": " Datahose",
                    "role": "submitter",
                    "email": "datahose@sanger.ac.uk"
                }],
                "samplecharacteristic": [{
                    "category": "organism",
                    "value": ["Schistosoma mansoni"]
                }, {
                    "category": "organism part",
                    "value": ["ovaries", "testes", "whole female", "whole male"]
                }, {
                    "category": "sample description",
                    "value": ["Schistosoma mansoni ovaries, mixed sex infection", "Schistosoma mansoni ovaries, single sex infection", "Schistosoma mansoni, testes from mixed sex infection", "Schistosoma mansoni, testes from single sex infection", "Schistosoma mansoni, whole female from mixed infection", "Schistosoma mansoni, whole female from mixed sex infection", "Schistosoma mansoni, whole female from single sex infection", "Schistosoma mansoni, whole male from mixed infection", "Schistosoma mansoni, whole male from mixed sex infection", "Schistosoma mansoni, whole male from single sex infection"]
                }, {
                    "category": "sex",
                    "value": ["female", "male"]
                }, {
                    "category": "strain",
                    "value": ["Liberia"]
                }],
                "experimentalvariable": [{
                    "name": "organism part",
                    "value": ["ovaries", "testes", "whole female", "whole male"]
                }, {
                    "name": "sample description",
                    "value": ["Schistosoma mansoni ovaries, mixed sex infection", "Schistosoma mansoni ovaries, single sex infection", "Schistosoma mansoni, testes from mixed sex infection", "Schistosoma mansoni, testes from single sex infection", "Schistosoma mansoni, whole female from mixed infection", "Schistosoma mansoni, whole female from mixed sex infection", "Schistosoma mansoni, whole female from single sex infection", "Schistosoma mansoni, whole male from mixed infection", "Schistosoma mansoni, whole male from mixed sex infection", "Schistosoma mansoni, whole male from single sex infection"]
                }, {
                    "name": "sex",
                    "value": ["female", "male"]
                }],
                "protocol": [{
                    "id": 1245480,
                    "accession": "P-ERP016356-1"
                }, {
                    "id": 1245481,
                    "accession": "P-ERP016356-2"
                }],
                "bioassaydatagroup": [{
                    "id": null,
                    "name": "scan",
                    "bioassaydatacubes": 276,
                    "arraydesignprovider": null,
                    "dataformat": "scan",
                    "bioassays": 276,
                    "isderived": 0
                }]
            }

        ]
    }
}
