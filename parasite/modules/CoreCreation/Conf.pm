#Experimental
#Inflate the little config we need to provide into a massive one
#Uses env vars for database passwords based on the right module load
#Writes config to stdout
#TODO: let config values override what this script would rather do
package CoreCreation::Conf;
use YAML;
use File::Basename;
use File::Spec;
use Getopt::Long qw(GetOptionsFromArray);

sub new {
    my ( $class, $path_to_conf ) = @_;
    die "Required: <path to species conf>" unless -f $path_to_conf;
    my $config = &read_config($path_to_conf);
    bless $config, $class;
    return $config;
}
sub db_name {
   return shift->{'meta'}->{"species.db_name"};
}

sub generics {
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
    return $generics;
}

sub dump {
    my $self = shift;
    my $db_name = $self->db_name();
    my $c={};
    while ( ( $key, $value ) = each %{$self} ) {
        $c->{$key} = $value;
    }
    return Dump({'generics' => $self->generics(),$db_name => $c});
}
sub read_config {
    my $path_to_conf = shift;

    my $staging_parse = `\${PARASITE_STAGING_MYSQL}-ensrw details script`;
    my @a = split /\s+/, $staging_parse;
    my ( $host, $password, $port, $user );
    GetOptionsFromArray(
        \@a,
        "host=s" => \$host,
        "pass=s" => \$password,
        "port=s" => \$port,
        "user=s" => \$user
    );
    $user eq "ensrw"
      or die
      "Parsing db details didn't work:( Got $staging_parse from system call";

    my ( $filename, $dirs, $suffix ) = fileparse( $path_to_conf, qw/conf/ );
    $filename =~ s/\.$//;
    $dirs =~ s/\/$//;
    die "Expected: <species alias>.conf, got $filename.$suffix"
      unless $suffix eq "conf";

    my ( $spe, $cies, $bioproject ) = split "_", $filename, 3;
    my $species = "${spe}_${cies}";

    my $species_conf = YAML::LoadFile($path_to_conf);

    $species_conf->{"core_database"} = {
        host     => $host,
        password => $password,
        port     => $port,
        user     => $user,
        dbname =>
"${species}_${bioproject}_core_$ENV{PARASITE_VERSION}_$ENV{ENSEMBL_VERSION}_1"
    };
    my $fasta_location = File::Spec->rel2abs("$dirs/${filename}.fa");
    my $gff3_location  = File::Spec->rel2abs("$dirs/${filename}.gff3");
    warn "We expect a good FASTA at $fasta_location" unless -f $fasta_location;
    warn "We expect a good GFF3 at $gff3_location"   unless -f $gff3_location;

    $species_conf->{fasta} = $fasta_location;
    $species_conf->{gff3}  = $gff3_location;

    my $taxon_id = $species_conf->{"taxon_id"};
    $taxon_id =~ /\d+/
      or die "taxon_id should be an NCBI string, got: $taxon_id";

    my $assembly_version = $species_conf->{"assembly_version"}
      or die "Required: assembly_version";

    $species_conf->{meta}->{"species.db_name"} = "${species}_${bioproject}";
    $species_conf->{meta}->{"species.display_name"} =
      sprintf( "%s %s (%s)", ucfirst($spe), $cies, uc($bioproject) );
    $species_conf->{meta}->{"species.division"} = "EnsemblParasite";
    $species_conf->{meta}->{"species.alias"} =
      sprintf( "%s %s", ucfirst($spe), $cies );
    $species_conf->{meta}->{"species.bioproject_id"} = uc($bioproject);
    $species_conf->{meta}->{"species.ftp_genome_id"} = uc($bioproject);
    $species_conf->{meta}->{"species.production_name"} =
      "${species}_${bioproject}";
    $species_conf->{meta}->{"species.scientific_name"} =
      sprintf( "%s %s", ucfirst($spe), $cies );
    $species_conf->{meta}->{"species.species_taxonomy_id"} = $taxon_id;
    $species_conf->{meta}->{"species.taxonomy_id"}         = $taxon_id;
    $species_conf->{meta}->{"species.url"}      = "${species}_${bioproject}";
    $species_conf->{meta}->{"assembly.default"} = $assembly_version;
    $species_conf->{meta}->{"assembly.name"}    = $assembly_version;

    my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) =
      localtime(time);
    my $gene_build = sprintf( "%s-%02d-WormBase", $year + 1900, $mon + 1 );
    $species_conf->{meta}->{"genebuild.start_date"} = $gene_build;
    $species_conf->{meta}->{"genebuild.version"}    = $gene_build;

    # Let the values from original file win
    my %c = %{ YAML::LoadFile($path_to_conf) };
    while ( ( $key, $value ) = each %c ) {
        next if $key eq "meta";
	$species_conf->{$key} = $value;
    }
    while ( ( $key, $value ) = each %{$c{"meta"}} ) {
       $species_conf->{'meta'}->{$key} = $value;
    }
    return $species_conf;
}

1;
