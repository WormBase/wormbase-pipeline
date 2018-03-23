use CoreCreation::Conf;
use File::Spec;
use File::Basename;
 
my $conf = CoreCreation::Conf->new(shift) or die "usage: $0 <species conf>";

my $fasta = $conf -> {fasta};
my $gff3 = $conf -> {gff3};
my $db_name = $conf -> {core_database}->{dbname};
my $host = $conf -> {core_database}->{host};
my $port = $conf -> {core_database}->{port};

my $healthcheck_script = File::Spec->catfile(dirname($0) , "healthcheck.sh"); 
system("$healthcheck_script -d $db_name -g $fasta -a $gff3 -u ensro -p $port -s $host");
 
