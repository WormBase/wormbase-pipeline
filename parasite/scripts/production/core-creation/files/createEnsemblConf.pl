#Experimental
#Inflate the little config we need to provide into a massive one
#Uses env vars for database passwords based on the right module load
#Writes config to stdout
#TODO: let config values override what this script would rather do
use YAML;
use File::Basename;
use Getopt::Long qw(GetOptionsFromArray);

my ($generics) = YAML::Load(<<'...');
cvsdir: /nfs/panda/ensemblgenomes/wormbase/software/packages/ensembl/branches/branch-89
taxonomy_database:
  host: mysql-eg-pan-prod.ebi.ac.uk
  port: 4276
  dbname: ncbi_taxonomy
production_database:
  host: mysql-eg-pan-prod.ebi.ac.uk
  user: ensro
  port: 4276
  dbname: ensembl_production_parasite
analysisconf:  /nfs/panda/ensemblgenomes/wormbase/parasite/config/ensembl-config/ParaSite/pipe_conf/analysis.conf
...
$generics->{'cvsdir'} = "/nfs/panda/ensemblgenomes/wormbase/software/packages/ensembl/branches/branch-$ENV{ENSEMBL_VERSION}";

my $path_to_conf = shift or die "Usage: $0 <path to species conf>";

my $staging_parse = `\${PARASITE_STAGING_MYSQL}-ensrw details script`;
my @a = split /\s+/, $staging_parse;
my ($host, $password, $port, $user); 
GetOptionsFromArray(\@a ,
  "host=s"=>\$host,
  "pass=s"=>\$password,
  "port=s"=>\$port,
  "user=s"=>\$user
);
$user eq "ensrw" or die "Parsing db details didn't work:( Got $staging_parse from system call";

my ($filename, $dirs, $suffix) = fileparse($path_to_conf, qw/conf/);
$filename =~ s/\.$//;
$dirs=~s/\/$//;
die "Expected: <species alias>.conf, got $filename.$suffix" unless $suffix eq "conf";

my ($spe, $cies, $bioproject) = split "_", $filename, 3;
my $species = "${spe}_${cies}";
 
my $species_conf = YAML::LoadFile($path_to_conf);

$species_conf-> { "core_database" } = {
    host => $host,
    password => $password,
    port => $port,
    user => $user,
    dbname => "${species}_${bioproject}_core_$ENV{PARASITE_VERSION}_$ENV{ENSEMBL_VERSION}_1"
};
my $fasta_location="$dirs/${filename}.fa";
my $gff3_location="$dirs/${filename}.gff3";
warn "We expect a good FASTA at $fasta_location" unless -f $fasta_location;
warn "We expect a good GFF3 at $gff3_location" unless -f $gff3_location;

$species_conf->{fasta} = $fasta_location;
$species_conf->{gff3} = $gff3_location;

my $taxon_id = $species_conf->{"taxon_id"};
$taxon_id =~ /\d+/ or die "taxon_id should be an NCBI string, got: $taxon_id";

my $assembly_version = $species_conf->{"assembly_version"} or die "Required: assembly_version";

$species_conf->{meta}->{"species.db_name"}="${species}_${bioproject}";
$species_conf->{meta}->{"species.display_name"}=sprintf("%s %s (%s)" , ucfirst($spe) , $cies, uc($bioproject));
$species_conf->{meta}->{"species.division"}="EnsemblParasite";
$species_conf->{meta}->{"species.alias"} = sprintf("%s %s", ucfirst($spe) , $cies);
$species_conf->{meta}->{"species.bioproject_id"} = uc($bioproject);
$species_conf->{meta}->{"species.ftp_genome_id"} = uc($bioproject);
$species_conf->{meta}->{"species.production_name"} = "${species}_${bioproject}";
$species_conf->{meta}->{"species.scientific_name"} = sprintf("%s %s", ucfirst($spe) , $cies);
$species_conf->{meta}->{"species.species_taxonomy_id"} = $taxon_id;
$species_conf->{meta}->{"species.taxonomy_id"} = $taxon_id;
$species_conf->{meta}->{"species.url"} = "${species}_${bioproject}";
$species_conf->{meta}->{"assembly.default"} = $assembly_version;
$species_conf->{meta}->{"assembly.name"} = $assembly_version;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $gene_build =sprintf("%s-%02d-WormBase", $year + 1900, $mon+1); 
$species_conf->{meta}->{"genebuild.start_date"} = $gene_build;
$species_conf->{meta}->{"genebuild.version"} = $gene_build;

print Dump({generics => $generics,"${species}_${bioproject}" => $species_conf});
