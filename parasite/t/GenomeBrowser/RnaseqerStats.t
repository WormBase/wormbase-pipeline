use Test::More;
use Test::MockModule;

use File::Temp qw/tempdir/;
use GenomeBrowser::RnaseqerMetadata;
use GenomeBrowser::RnaseqerStats;
our $run_id = "SRR1124914";
my $module = new Test::MockModule('GenomeBrowser::LocallyCachedResource');
our $text = do { local $/; <DATA> };
$module->mock(
    'get_text',
    sub {
        my ( $class, $url ) = @_;
        if ( $url =~ /$run_id/ ) {
            return $text;
        }
        else {
            die @_;
            return undef;
        }
    }
);

my $root_dir          = tempdir( CLEANUP => 1 );
my $species           = "ancylostoma_ceylanicum";
my $assembly          = "acey_assembly";
my $study_id          = "SRP035476";
my $rnaseqer_metadata = bless {
    metadata => {
        $assembly => {
            $study_id => {
                $run_id => {}
            },
        },
    },
    location_per_run_id => {
        SRR1124914 =>
          "ftp://invalid/SRR112/004/SRR1124914/SRR1124914.nospliced.bw"
    },
  },
  'GenomeBrowser::RnaseqerMetadata';

my $subject = GenomeBrowser::RnaseqerStats->new( $root_dir, $species,
    $rnaseqer_metadata );
is_deeply(
    $subject,
    bless(
        {
            'SRR1124914' => {
                'fraction_reads_uniquely_mapped' => '0.856',
                'library_size'                   => '59712662'
            }
        },
        'GenomeBrowser::RnaseqerStats'
    ),
    "Get properties from FTP"
) or diag explain $subject;
is_deeply(
    $subject->get_formatted_stats('SRR1124914'),
 {
    'library_size_approx' => '50-75mln',
    'library_total_amount_of_reads' => '59712662',
    'mapping_fraction_of_uniquely_mapped_reads' => '0.856',
    'mapping_quality_approx' => '80-90%'
  },
   "Can format stats" 
) or diag explain $subject->get_formatted_stats('SRR1124914');
done_testing();
__DATA__
All_entries,59712662
Valid_entries,59712662
Duplicate,0
Alignments,55347142
Spliced,26499976
ReadsSpliced,24425911
Paired,59712662
Paired_mapped,51761148
Paired_mapped_mate1,27610474
Paired_mapped_mate2,27736668
ReadsUnmapped,4365520
ReadsMapped,52431398
  NM0,46116968
  NM1,7384545
  NM2,1845629
NMGE3,0
plus_strand,21888348
minus_strand,21541988
UniquelyMappedReads,51104386
MultimapReads,1327012
  NH1,51104386
NHGE2,1327012
NHGE10,112259
NHGE20,0
Spliced1,24675787
Spliced2,1782894
Spliced3,41124
SplicedGE4,171
 MAPQ,46.2785306240384
