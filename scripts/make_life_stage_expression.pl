#!/nfs/panda/ensemblgenomes/wormbase/software/bin/perl -w

# 2017-06-28 Added output of CSV table of values


use Ace;
use strict;                                      
use Getopt::Long;
use List::Util qw(max sum);
use File::Path qw(make_path);
use lib $ENV{CVS_DIR};
use Wormbase;

my ($wormbase, $debug, $test, $store);
my ($database, $gene_name, $outfile);
GetOptions(
           "database=s"     => \$database,
           "gene=s"         => \$gene_name, # do a single gene, for testing
	   "outfile=s"      => \$outfile, # output filename, for testing
	   "test"           => \$test,
	   "debug=s"        => \$debug,
	   "store=s"        => \$store,
          );

# The embryo times
# mapping the EE_50-* life stages to WBls terms
# Remark "Library 'EE_50-0' Illumina sequencing of C. elegans N2 early
# embryo EE_50-0 polyA+ RNAseq random fragment library Make eggs hatch
# in the absence of food to get them all in L1 arrest add food and wait
# 50 hours. Eggs are then harvested from the adults -- just a few have
# an egg or two. Then the eggs are incubated for 0 to 720
# minutes. Synchronization was only approximate with a distribution of
# embryo ages. Sample may include multiple embryonic stages since
# mothers maintain eggs until the 30-cell stage leaving a possible 150
# minute range of ages (in this case 0-150m post-fertilization or -40 to
# 110 mins post-cleavage)."
#  knock off 40 mins and add 110 mins to get the range of ages post cleavage (used by WormBase)
#  add 150 to get the range post fertilization
#  My approximate mapping table from EE_50 time to WB life stage
# EE_50 time   range    approximate life-stage
# 0     -40-110   WBls:0000004
# 30    -10-140   WBls:0000004
# 60    20-170    WBls:0000004
# 90    50-200   WBls:0000004
# 120   80-270   WBls:0000010
# 150   110-300  WBls:0000010
# 180   140-290  WBls:0000010
# 210   170-320  WBls:0000013
# 240   200-350  WBls:0000013
# 270   230-380  WBls:0000014
# 300   260-410  WBls:0000014
# 330   290-440  WBls:0000014
# 360   320-470  WBls:0000014
# 390   350-500  WBls:0000015
# 420   380-530  WBls:0000015
# 450   410-560  WBls:0000015
# 480   440-590  WBls:0000015
# 510   470-620  WBls:0000015
# 540   500-650  WBls:0000020
# 570   530-680  WBls:0000020
# 600   560-710  WBls:0000020
# 630   590-740  WBls:0000021
# 660   620-770  WBls:0000021
# 690   650-800  WBls:0000021
# 720   680-830  WBls:0000021
#  life stage    time post cleavage
# WBls:0000004 = 0-350min        proliferating embryo Ce
# WBls:0000010 = 100-290min      gastrulating embryo Ce
# WBls:0000013 = 290-350min      enclosing embryo Ce
# WBls:0000014 = 210-350min      late cleavage stage embryo Ce
# WBls:0000015 = 350-620min      elongating embryo Ce
# WBls:0000019 = 460-520min      2-fold embryo Ce
# WBls:0000020 = 520-620min      3-fold embryo Ce
# WBls:0000021 = 620-800min      fully-elongated embryo Ce
# Embryo
# WBls:0000003
# embryo Ce
# 4 cell embryo
# WBls:0000008

my @etimes = (
	      "0",
	      "30",
	      "60",
	      "90",
	      "120",
	      "150",
	      "180",
	      "210",
	      "240",
	      "270",
	      "300",
	      "330",
	      "360",
	      "390",
	      "420",
	      "450", # e.g. SRX1022654 "embryo sampled at 450 minutes"
	      "480",
	      "510",
	      "540",
	      "570",
	      "600",
	      "660",
	      "690",
	      "720",
	    );

# The classical time points;
my @ctimes = (
	      "EE",
	      "LE",
	      "L1",
	      "L2",
	      "L3",
	      "L4",
	      "YA",
	      "Dauer entry",
	      "Dauer",
	      "Dauer exit",
	      "Male L4",
	      "Soma L4"
	     );


# the RNASeq SRP000401 experiments - excluding tissue-specific ones, and ribomins ones and mutants and anything else Julie doesn't like
my %experiments = (
		 SRX092477 => ['0','polyA', 'N2_EE_50-0'],
		 SRX092478 => ['0','polyA', 'N2_EE_50-0'],
		 SRX099902 => ['0','polyA', 'N2_EE_50-0'],
		 SRX099901 => ['0','polyA', 'N2_EE_50-0'],
		 SRX103649 => ['0','polyA', 'N2_EE_50-0'],
		 SRX1022600 => ['0','ribozero', '20120411_EMB-0'],
		 SRX1020637 => ['0','ribozero', '20120223_EMB-0'],
		 SRX1020636 => ['0','ribozero', '20120223_EMB-0'],

		 SRX092371 => ['30','polyA', 'N2_EE_50-30'],
		 SRX092372 => ['30','polyA', 'N2_EE_50-30'],
		 SRX099908 => ['30','polyA', 'N2_EE_50-30'],
		 SRX099907 => ['30','polyA', 'N2_EE_50-30'],
		 SRX103650 => ['30','polyA', 'N2_EE_50-30'],
		 SRX1020634 => ['30','ribozero', '20120223_EMB-30'],
		 SRX1022610 => ['30','ribozero', '20120419_EMB-30'],
		 SRX1020635 => ['30','ribozero', '20120223_EMB-30'],

		 SRX085112 => ['60','polyA', 'N2_EE_50-60'],
		 SRX085111 => ['60','polyA', 'N2_EE_50-60'],
		 SRX1022599 => ['60','ribozero', '20120411_EMB-60'],
		 SRX1020638 => ['60','ribozero', '20120223_EMB-60'],
		 SRX1020639 => ['60','ribozero', '20120223_EMB-60'],

		 SRX092480 => ['90','polyA', 'N2_EE_50-90'],
		 SRX092479 => ['90','polyA', 'N2_EE_50-90'],
		 SRX099915 => ['90','polyA', 'N2_EE_50-90'],
		 SRX103651 => ['90','polyA', 'N2_EE_50-90'],
		 SRX1022605 => ['90','ribozero', '20120411_EMB-90'],
		 SRX1020640 => ['90','ribozero', '20120223_EMB-90'],
		 SRX1020641 => ['90','ribozero', '20120223_EMB-90'],
		 SRX1022611 => ['90','ribozero', '20120419_EMB-90'],

		 SRX085217 => ['120','polyA', 'N2_EE_50-120'],
		 SRX085218 => ['120','polyA', 'N2_EE_50-120'],
		 SRX1022602 => ['120','ribozero', '20120411_EMB-120'],
		 SRX1022645 => ['120','ribozero', '20120419_EMB-120'],
		 SRX1020630 => ['120','ribozero', '20120223_EMB-120'],
		 SRX1020631 => ['120','ribozero', '20120223_EMB-120'],

		 SRX099995 => ['150','polyA', 'N2_EE_50-150'],
		 SRX1022601 => ['150','ribozero', '20120411_EMB-150'],
		 SRX1020632 => ['150','ribozero', '20120223_EMB-150'],
		 SRX1020633 => ['150','ribozero', '20120223_EMB-150'],
		 SRX1022646 => ['150','ribozero', '20120419_EMB-150'],

		 SRX099985 => ['180','polyA', 'N2_EE_50-180'],
		 SRX1022603 => ['180','ribozero', '20120411_EMB-180'],
		 SRX1022584 => ['180','ribozero', '20120223_EMB-180'],
		 SRX1022585 => ['180','ribozero', '20120223_EMB-180'],
		 SRX1022647 => ['180','ribozero', '20120419_EMB-180'],

		 SRX099996 => ['210','polyA', 'N2_EE_50-210'],
		 SRX099997 => ['210','polyA', 'N2_EE_50-210'],
		 SRX099998 => ['210','polyA', 'N2_EE_50-210'],
		 SRX103652 => ['210','polyA', 'N2_EE_50-210'],
		 SRX1022570 => ['210','ribozero', '20120223_EMB-210'],
		 SRX1022571 => ['210','ribozero', '20120223_EMB-210'],

		 SRX099986 => ['240','polyA', 'N2_EE_50-240'],
		 SRX099987 => ['240','polyA', 'N2_EE_50-240'],
		 SRX103653 => ['240','polyA', 'N2_EE_50-240'],
		 SRX1022604 => ['240','ribozero', '20120411_EMB-240'],
		 SRX1022566 => ['240','ribozero', '20120223_EMB-240'],
		 SRX1022567 => ['240','ribozero', '20120223_EMB-240'],
		 SRX1022648 => ['240','ribozero', '20120419_EMB-240'],

		 SRX099999 => ['270','polyA', 'N2_EE_50-270'],
		 SRX100000 => ['270','polyA', 'N2_EE_50-270'],
		 SRX100001 => ['270','polyA', 'N2_EE_50-270'],
		 SRX103677 => ['270','polyA', 'N2_EE_50-270'],
		 SRX1022568 => ['270','ribozero', '20120223_EMB-270'],
		 SRX1022569 => ['270','ribozero', '20120223_EMB-270'],
		 SRX1022649 => ['270','ribozero', '20120419_EMB-270'],

		 SRX100819 => ['300','polyA', 'N2_EE_50-300'],
		 SRX1022580 => ['300','ribozero', '20120223_EMB-300'],
		 SRX1022581 => ['300','ribozero', '20120223_EMB-300'],
		 SRX1022608 => ['300','ribozero', '20120411_EMB-300'],
		 SRX1022650 => ['300','ribozero', '20120419_EMB-300'],

		 SRX099980 => ['330','polyA', 'N2_EE_50-330'],
		 SRX1022572 => ['330','ribozero', '20120223_EMB-330'],
		 SRX1022573 => ['330','ribozero', '20120223_EMB-330'],
		 SRX1022651 => ['330','ribozero', '20120419_EMB-330'],

		 SRX099981 => ['360','polyA', 'N2_EE_50-360'],
		 SRX1022574 => ['360','ribozero', '20120223_EMB-360'],
		 SRX1022575 => ['360','ribozero', '20120223_EMB-360'],
		 SRX1022607 => ['360','ribozero', '20120411_EMB-360'],
		 SRX1022652 => ['360','ribozero', '20120419_EMB-360'],

		 SRX099982 => ['390','polyA', 'N2_EE_50-390'],
		 SRX099983 => ['390','polyA', 'N2_EE_50-390'],
		 SRX1022576 => ['390','ribozero', '20120223_EMB-390'],
		 SRX1022577 => ['390','ribozero', '20120223_EMB-390'],

		 SRX099984 => ['420','polyA', 'N2_EE_50-420'],
		 SRX1022578 => ['420','ribozero', '20120223_EMB-420'],
		 SRX1022579 => ['420','ribozero', '20120223_EMB-420'],
		 SRX1022653 => ['420','ribozero', '20120419_EMB-420'],

		 SRX100002 => ['450','polyA', 'N2_EE_50-450'],
		 SRX1022582 => ['450','ribozero', '20120223_EMB-450'],
		 SRX1022583 => ['450','ribozero', '20120223_EMB-450'],
		 SRX1022654 => ['450','ribozero', '20120419_EMB-450'],

		 SRX099988 => ['480','polyA', 'N2_EE_50-480'],
		 SRX099989 => ['480','polyA', 'N2_EE_50-480'],
		 SRX099990 => ['480','polyA', 'N2_EE_50-480'],
		 SRX103672 => ['480','polyA', 'N2_EE_50-480'],
		 SRX1022586 => ['480','ribozero', '20120223_EMB-480'],
		 SRX1022587 => ['480','ribozero', '20120223_EMB-480'],

		 SRX100003 => ['510','polyA', 'N2_EE_50-510'],
		 SRX100004 => ['510','polyA', 'N2_EE_50-510'],
		 SRX100005 => ['510','polyA', 'N2_EE_50-510'],
		 SRX103673 => ['510','polyA', 'N2_EE_50-510'],
		 SRX1022588 => ['510','ribozero', '20120223_EMB-510'],
		 SRX1022589 => ['510','ribozero', '20120223_EMB-510'],

		 SRX099991 => ['540','polyA', 'N2_EE_50-540'],
		 SRX099992 => ['540','polyA', 'N2_EE_50-540'],
		 SRX099993 => ['540','polyA', 'N2_EE_50-540'],
		 SRX103669 => ['540','polyA', 'N2_EE_50-540'],
		 SRX1022592 => ['540','ribozero', '20120223_EMB-540'],
		 SRX1022593 => ['540','ribozero', '20120223_EMB-540'],

		 SRX099973 => ['570','polyA', 'N2_EE_50-570'],
		 SRX099974 => ['570','polyA', 'N2_EE_50-570'],
		 SRX103671 => ['570','polyA', 'N2_EE_50-570'],
		 SRX1022597 => ['570','ribozero', '20120223_EMB-570'],
		 SRX1022598 => ['570','ribozero', '20120223_EMB-570'],

		 SRX099975 => ['600','polyA', 'N2_EE_50-600'],
		 SRX099976 => ['600','polyA', 'N2_EE_50-600'],
		 SRX099977 => ['600','polyA', 'N2_EE_50-600'],
		 SRX103670 => ['600','polyA', 'N2_EE_50-600'],
		 SRX1022596 => ['600','ribozero', '20120223_EMB-600'],
		 SRX1022595 => ['600','ribozero', '20120223_EMB-600'],
		 SRX1022609 => ['600','ribozero', '20120411_EMB-600'],

		 SRX099978 => ['630','polyA', 'N2_EE_50-630'],

		 SRX099979 => ['660','polyA', 'N2_EE_50-660'],

		 SRX099994 => ['690','polyA', 'N2_EE_50-690'],

		 SRX100006 => ['720','polyA', 'N2_EE_50-720'],

		 SRX004863 => ['EE','polyA', 'EE_ce0128_rw005'],
		 SRX004864 => ['EE','polyA', 'EE_ce1003_rw005'],
		 SRX037186 => ['EE','polyA', 'N2_EE-2'],
		 SRX004866 => ['EE','polyA', 'EE_ce0129_rw006'], # checked with LaDeana - she says this is an early embryo
		 SRX145660 => ['EE','ribozero', 'N2_EE_RZ-54'],
		 SRX190369 => ['EE','ribozero', 'N2_EE_RZ-54'],

		 SRX004865 => ['LE','polyA', 'LE_ce0129_rw006'],
		 SRX047446 => ['LE','polyA', 'N2_LE-1'],

		 SRX004867 => ['L1','polyA', 'L1_ce0132_rw007'], # fastq files downloaded again because they were in an odd format - all ok now
		 SRX037288 => ['L1','polyA', 'N2_L1-1'],

		 SRX001872 => ['L2','polyA', 'L2_ce0109_rw001'],
		 SRX047653 => ['L2','polyA', 'N2_L2-4'],
		 SRX190370 => ['L2','ribozero', 'N2_L2_RZ-53'],
		 SRX145661 => ['L2','ribozero', 'N2_L2_RZ-53'],

		 SRX001875 => ['L3','polyA', 'L3_ce0120_rw002'],
		 SRX036881 => ['L3','polyA', 'N2_L3-1'],

		 SRX008144 => ['L4','polyA', 'L4_ce1009_rw1001'],
		 SRX001874 => ['L4','polyA', 'L4_ce0121_rw003'],

		 SRX001873 => ['YA','polyA', 'YA_ce0122_rw004'],
		 SRX047787 => ['YA','polyA', 'N2_Yad-1'],
		 SRX103986 => ['YA','ribozero', 'N2_YA_RZ-1'],
		 SRX103987 => ['YA','ribozero', 'N2_YA_RZ-1'],
		 SRX103988 => ['YA','ribozero', 'N2_YA_RZ-1'],
		 SRX103989 => ['YA','ribozero', 'N2_YA_RZ-1'],

		 SRX011569 => ['Male EM','polyA', 'EmMalesHIM8_ce1005_rw1001'],
		 SRX037198 => ['Male EM','polyA', 'EmMalesHIM8-2'],

		 SRX004868 => ['Male L4','polyA', 'L4_ce1001_rw1001'],
		 SRX047469 => ['Male L4','polyA', 'L4MALE5'],

		 SRX014010 => ['Soma L4','polyA', 'L4JK1107soma_ce1014_rw1001'],
		 SRX037200 => ['Soma L4','polyA', 'L4JK1107soma-2'],

		 SRX008139 => ['Dauer entry','polyA', 'DauerEntryDAF2_ce1007_rw1001'],
		 SRX047470 => ['Dauer entry','polyA', 'DauerEntryDAF2-2'],
		 SRX103273 => ['Dauer entry','polyA', 'DauerEntryDAF2-1-1'],
		 SRX103274 => ['Dauer entry','polyA', 'DauerEntryDAF2-1-1'],
		 SRX103275 => ['Dauer entry','polyA', 'DauerEntryDAF2-1-1'],
		 SRX103276 => ['Dauer entry','polyA', 'DauerEntryDAF2-1-1'],
		 SRX103277 => ['Dauer entry','polyA', 'DauerEntryDAF2-4-1'],

		 SRX008138 => ['Dauer','polyA', 'DauerDAF2_ce1006_rw1001'],
		 SRX103983 => ['Dauer','polyA', 'DauerDAF2-2-1'],
		 SRX103984 => ['Dauer','polyA', 'DauerDAF2-2'],
		 SRX103985 => ['Dauer','polyA', 'DauerDAF2-5-1'],

		 SRX008140 => ['Dauer exit','polyA', 'DauerExitDAF2_ce1008_rw1001'],
		 SRX037199 => ['Dauer exit','polyA', 'DauerExitDAF2-2'],
		 SRX103269 => ['Dauer exit','polyA', 'DauerExitDAF2-3-1'],
		 SRX103270 => ['Dauer exit','polyA', 'DauerExitDAF2-3-1'],
		 SRX103271 => ['Dauer exit','polyA', 'DauerExitDAF2-3-1'],
		 SRX103272 => ['Dauer exit','polyA', 'DauerExitDAF2-3-1'],
		 SRX103278 => ['Dauer exit','polyA', 'DauerExitDAF2-6-1'],
		 SRX103281 => ['Dauer exit','polyA', 'DauerExitDAF2-6-1'],
		 SRX103280 => ['Dauer exit','polyA', 'DauerExitDAF2-6-1'],
		 SRX103279 => ['Dauer exit','polyA', 'DauerExitDAF2-6-1'],

		);

my %libraries;
foreach my $experiment (keys %experiments) {
  my ($stage, $type, $library) = @{$experiments{$experiment}};
  if (! exists $libraries{$library}) {
    push @{$libraries{$library}}, ($stage, $type);
  }
}


if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                           );
}


$database = $wormbase->autoace unless $database;
print "Using database $database\n";


print "Connecting to $database ...\n";
my $db = Ace->connect (-path => "$database") || die("Cannot connect to database at $database\n");

if (!defined $outfile) {$outfile = $wormbase->spell ."/expr_graph.csv"}

open (OUT, ">$outfile") || die "Can't open $outfile";
print OUT '"Gene","Gene_name","Life-stage","Library","Protocol","FPKM value"'."\n";

my $gene;
if (defined $gene_name) { # for testing
  $gene = $db->fetch(Gene => $gene_name);
  # get data
  my ($gene_data, $cds_data, $trans_data, $pseud_data) = get_gene_data($gene);
  if (!defined $gene_name) {$outfile = "$gene.png"}

  # draw graphs  
    my ($name) = keys %{$gene_data};
    my $data = $gene_data->{$name};
    if (scalar keys %{$data} == 0) {print "No data found for $gene\n"; next} # no data - usually because it is a transposon
    draw_graph('Gene', $data, $name, $gene_name);
  
    # draw CDS graphs
#    foreach my $name (keys %{$cds_data} ) {
#      if (scalar keys %{$cds_data->{$name}} == 0) {print "No data found for $name\n"; next} # no data
#      draw_graph('CDS', $cds_data->{$name}, $name, $gene->name);
#    }

} else {
  my $iterator = $db->fetch_many(-query => 'Find Gene where Status = "Live" AND Species = "Caenorhabditis elegans"');
  while ($gene = $iterator->next) {
	  unless ($gene=~/WBGene/) {next}; # Make sure not to read VC2010 genes
    print "Doing $gene ...\n";
#    if ($gene->name eq 'WBGene00001495') {
#      next; # lacks a Public_name in WS273
#    }
    # get data
    my ($gene_data, $cds_data, $trans_data, $pseud_data) = get_gene_data($gene);
    #       $half_length = $full_length / 2;

    # draw gene graph
    my ($name) = keys %{$gene_data};
    my $data = $gene_data->{$name};
    if (scalar keys %{$data} == 0) {print "No data found for $gene\n"; next} # no data - usually because it is a transposon
    draw_graph('Gene', $data, $name, $gene->name);
    
    # draw CDS graphs
#    foreach my $name (keys %{$cds_data} ) {
#      if (scalar keys %{$cds_data->{$name}} == 0) {print "No data found for $name\n"; next} # no data
#      draw_graph('CDS', $cds_data->{$name}, $name, $gene->name);
#    }

    # draw Transcript graphs
#    foreach my $name (keys %{$trans_data} ) {
#      if (scalar keys %{$trans_data->{$name}} == 0) {print "No data found for $name\n"; next} # no data
#      draw_graph('Transcript', $trans_data->{$name}, $name, $gene->name);
#    }
    # draw Pseudogene graphs
#    foreach my $name (keys %{$pseud_data} ) {
#      if (scalar keys %{$pseud_data->{$name}} == 0) {print "No data found for $name\n"; next} # no data
#      draw_graph('Pseudogene', $pseud_data->{$name}, $name, $gene->name);
#    }
  }
}

close(OUT);
$db->close;

##################################

sub draw_graph {
  my ($type, $data, $name, $gene_id) = @_;
  # type of ace object - 'gene', 'cds', 'trans', 'pseud'
  # ref to hash of data points to plot, keyed by library name
  # title = name to display at top or ''
  # gene_id


  # get width of image - the classical stages are displayed with 5 times the width of the embryonic time series
  # we have several potential parts of the x-axis:
  # embryonic time series 0-720 minutes
  # classical stages (EE, LE, L1-4, YA)
  # Male EM, L4 stage
  # Soma
  # Dauer entry/exit = optional dauer stages


  # get values by stage
  my %stages;
  foreach my $library (keys %{$data}) {
    my $fpkm_value = $data->{$library};
    my ($stage, $type) = @{$libraries{$library}};
    #print "$library $stage, $type\n";
    push @{$stages{$stage}}, [$fpkm_value, $type];
  }


  # make histograms
  foreach my $stage (keys %stages) {
    my @flist;
    my $median;
    foreach my $value (@{$stages{$stage}}) {
      my ($fpkm, $type) = @{$value};
      push @flist, $fpkm;
    }
    $median = median(@flist);

    print OUT "\"$gene_id\", \"$name\", \"$stage\", \"Median\", \"Median\", $median\n";
  }
  


  # plot points
  # $y_axis_y is the zero position of the y axis
  # $scale is the number of y-axis pixels per 1 data point
  foreach my $library (keys %{$data}) {
    my $fpkm_value = $data->{$library};
    my ($stage, $type) = @{$libraries{$library}};
    print OUT "\"$gene_id\", \"$name\", \"$stage\", \"$library\", \"$type\", $fpkm_value\n";
  }



}

############################################################################
# my ($name, %data) = get_gene_data($gene);
# Args: $gene - ace object of gene to process

sub get_gene_data {
  my ($gene) = @_;

  my @CDS =  $gene->Corresponding_CDS;
  my @Transcript = $gene->Corresponding_transcript;
  my @Pseudogene = $gene->Corresponding_pseudogene;
  
  my ($gene_data) = get_data('gene',  ($gene));
  my ($cds_data) = get_data('cds', @CDS);
  my ($trans_data);# = get_data('trans', @Transcript);
  my ($pseud_data);# = get_data('pseud', @Pseudogene);

  return ($gene_data, $cds_data, $trans_data, $pseud_data);
}

############################################################################
# my ($name, $data) = get_data('gene', $gene);
# Args: type of ace object - 'gene', 'cds', 'trans', 'pseud'
#       @objs - array of ace objects of gene to process
# Return: ref to hash (keyed by object name) of objects' expression values (hash keyed by library name containing FPKM values)
#         

sub get_data {
  my ($type, @objs) = @_;
  my %data;

  foreach my $obj (@objs) {
    my %objdata;
    my %data_by_library;
    my $dummy;
    my $name;

    if ($type eq 'gene') {
      $name = $obj->Public_name->name;
    } else {
      $name = $obj->name;
    }

    my @FPKM = $obj->RNASeq_FPKM;
    my $prev_value = '';
#    print ">>>>>START $type $name\n";
    foreach my $F (@FPKM) {
      my $ace_text = $F->asAce;
      my @lines = split /\n/, $ace_text;
      foreach my $line (@lines) {
	chomp $line;
	$line =~ s/\t$//; # some lines end with a spurious TAB
	if ($line eq '') {next}
#	print "$line (prev=$prev_value)\n";
	# '11.8211	From_analysis	"RNASeq.elegans.N2.WBls:0000038.Hermaphrodite.WBbt:0007833.SRP035479.SRX435700"'
	my ($value, $analysis) = ($line =~ /(\S+)\s+\S+\s+(\S+)/); 

	# '	From_analysis	"RNASeq.elegans.PS4730.WBls:0000038.Male.WBbt:0005062.SRP015688.SRX185680"'
	if (!defined $value || $value eq '') {
	  $value = 0; ($dummy, $analysis) = ($line =~ /(\S+)\s+(\S+)/);
	} 
	# '""	""	"RNASeq.elegans.PS4730.WBls:0000038.Male.WBbt:0005062.SRP015688.SRX185663"'
	if ($value eq '""') {
	  $value = $prev_value;
	} 

	$prev_value = $value;
	my ($srx) = ($analysis =~ /\.(\w+)"$/);
#	print "val=$value\tanal=$analysis\tprev=$prev_value\n";
	if (exists $experiments{$srx}) {
	  my ($stage, $type, $library) = @{$experiments{$srx}};
	  push @{$data_by_library{$library}}, $value;
	}
      }
    }
    
    # now get median value of technical replicates in each library
    foreach my $library (keys %data_by_library) {
      $objdata{$library} = median(@{$data_by_library{$library}});
    }

    $data{$name} = {%objdata};
  }

  return (\%data);
}
############################################################################
# return the median value of a list of values
sub median {

    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}
############################################################################
# return the mean value of a list of values
# this expects there to be some values in the input array!
sub mean {
  return sum(@_)/@_;
}
############################################################################
