use strict;                                      
use warnings;
use Path::Class;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Ace;
use Modules::Curate;
use FileHandle;

my ($debug, $outfile, $idfile, $species, $database, $test, $class);

GetOptions (
    "debug=s"    => \$debug,
    "outfile=s"  => \$outfile,
    "species=s"  => \$species,
    "test"       => \$test,
    "ids=s"      => \$idfile,
    "class=s"    => \$class,
    "database=s" => \$database
    );

my $wormbase = Wormbase->new(
    -debug    => $debug,
    -test     => $test,
    -organism => $species,
    -autoace  => $database
    );

my $out = FileHandle->new("> $outfile");
my $log = Log_files->make_build_log($wormbase);

if (!defined $class) {$class = 'CDS';}
if (!defined $species) {$species = $wormbase->species;}
if (!defined $database) {
    if ($species eq 'elegans') {
	$database = $wormbase->wormpub."/CURATION_DATABASES/camace_$ENV{USER}";
    } else {
	$database = $wormbase->wormpub."/CURATION_DATABASES/${species}_curation";
    }
}
if (!-e $database || !-d $database) {die "Can't find database $database\n"}

my $ace = Ace->connect (-path => $database) || die "cannot connect to database at $database\n";
my $curate = Curate->new($out, $ace, $wormbase, $log);


my $ids = parse_id_file($idfile);

for my $id (@$ids) {
    my $gene = $curate->SeqName2Gene($id);
    if (!defined $gene) {
	$log->error("ERROR Can't make a history object from $id - it is not attached to a Gene\n");
    } else {
	$curate->Make_history($class, $id);
    }
    $curate->check_force_flag();
}

$ace->close;
$log->mail();
exit(0);


sub parse_id_file {
    my $file = shift;

    my @ids;
    my $fh = file($file)->openr;
    while (my $line = $fh->getline()) {
	chomp $line;
	push @ids, $line;
    }

    return \@ids;
}
