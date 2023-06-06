use strict;
use Carp;
package ProductionMysql;
use Bio::EnsEMBL::Registry;
# Ensembl has a nice production environment.
# The databases we need are all in PATH so we can alias them.
# To get parameters in command line you can e.g. `$PARASITE_STAGING_MYSQL details host`
# Sometimes we need to prepare configs and stuff, and we want the same from Perl.
# This module interfaces the production env to Perl. The format is right for passing into Ensembl's pipelines as config, hopefully :)
# Usage: 
#  ProductionMysql->staging->url
#  ProductionMysql->new("mysql-ps-prod")->conn->{host}
#  ProductionMysql->new("mysql-pan-1")->conn("ncbi_taxonomy")
sub new { 
 my ($class, $db_cmd) = @_;
 return bless {db_cmd => $db_cmd}, $class;
}
sub staging {
  my $v = $ENV{PARASITE_STAGING_MYSQL} or die "PARASITE_STAGING_MYSQL not in env. You need to module load parasite_prod_relx";
  return new(shift, $v);
}
sub previous_staging {
  my $v = $ENV{PREVIOUS_PARASITE_STAGING_MYSQL} or die "PREVIOUS_PARASITE_STAGING_MYSQL not in env. You need to module load parasite_prod_relx";
  return new(shift, $v);
}
sub staging_writable {
  my $v = $ENV{PARASITE_STAGING_MYSQL} or die "PARASITE_STAGING_MYSQL not in env. You need to module load parasite_prod_relx";
  return new(shift,"$v-w"); 
}
sub core_databases {
  my $db_cmd= shift -> {db_cmd};
  my @patterns = @_;
  my @all_core_dbs;
  open(my $fh, " { $db_cmd  -Ne 'show databases like \"%core%\" ' 2>&1 1>&3 | grep -v \"can be insecure\" 1>&2; } 3>&1 |") or Carp::croak "ProductionMysql: $db_cmd not in your PATH\n";
  while(<$fh>) {
   chomp;
   push @all_core_dbs, $_ if $_;
  }
  return @all_core_dbs unless @patterns;

  my @result;
  for my $core_db (@all_core_dbs){
    my $include;
    for my $pat (@patterns){
      chomp $pat;
      $include = 1 if $core_db =~ /$pat/;
    }
    push @result, $core_db if $include;
  }
  return @result;
}

sub variation_databases {
  my $db_cmd= shift -> {db_cmd};
  my @patterns = @_;
  my @all_variation_dbs;
  open(my $fh, " { $db_cmd  -Ne 'show databases like \"%variation%\" ' 2>&1 1>&3 | grep -v \"can be insecure\" 1>&2; } 3>&1 |") or Carp::croak "ProductionMysql: $db_cmd not in your PATH\n";
  while(<$fh>) {
   chomp;
   push @all_variation_dbs, $_ if $_;
  }
  return @all_variation_dbs unless @patterns;

  my @result;
  for my $var_db (@all_variation_dbs){
    my $include;
    for my $pat (@patterns){
      chomp $pat;
      $include = 1 if $var_db =~ /$pat/;
    }
    push @result, $var_db if $include;
  }
  return @result;
}

sub core_db {
  my ($self, @patterns) = @_;
  my ($core_db, @others) = $self->core_databases(@patterns);
  Carp::croak "ProductionMysql $self->{db_cmd}: no db for: @patterns" unless $core_db;
  Carp::croak "ProductionMysql $self->{db_cmd}: multiple dbs for: @patterns" if @others;
  return $core_db;
}
sub species_for_core_db {
  my ($spe, $cies, $bp ) = split "_", shift;
  
  return join "_", $spe, $cies, ($bp eq 'core' ? () : $bp);
}
sub one_species {
  my ($self, @patterns) = @_;
  return species_for_core_db($self->core_db(@patterns));
}
sub species {
  my ($self, @patterns) = @_;
  my %h;
  for ($self->core_databases(@patterns)) {
    $h{&species_for_core_db($_)}++;
  }
  return sort keys %h;
}
sub meta_value {
  my ($self, $core_db_pattern, $pattern) = @_;
  my $db_cmd = $self->{db_cmd};
  my $core_db = $self->core_db($core_db_pattern);
  my @result;
  open(my $fh, "$db_cmd $core_db -Ne 'select meta_value from  meta where meta_key like \"$pattern\" ' |") or Carp::croak "ProductionMysql: $db_cmd not in your PATH\n";
  while(<$fh>) {
   chomp;
   push @result, $_ if $_;
  }
  return @result if wantarray;
  return $result[0] if @result;
  return undef;
}

sub set_meta_value {
  my ($self, $core_db_pattern, $meta_key, $new_value) = @_;
  my $db_cmd = $self->{db_cmd};
  my $core_db = $self->core_db($core_db_pattern);
  open(my $fh_1, " { $db_cmd  $core_db -Ne 'delete from meta where meta_key=\"$meta_key\" ' 2>&1 1>&3 | grep -v \"can be insecure\" 1>&2; } 3>&1 |") or Carp::croak "ProductionMysql: $db_cmd not in your PATH\n";
  while(<$fh_1>) {
   chomp;
   print STDERR "$_\n" if $_;
  }
  close $fh_1;
  open(my $fh_2, " { $db_cmd  $core_db -Ne 'insert into meta(meta_key, meta_value) values (\"$meta_key\", \"$new_value\") ' 2>&1 1>&3 | grep -v \"can be insecure\" 1>&2; } 3>&1 |") or Carp::croak "ProductionMysql: $db_cmd not in your PATH\n"; 
  while(<$fh_2>) {
   chomp;
   print STDERR "$_\n" if $_;
  }
  close $fh_2;
}

sub conn {
  my $db_cmd= shift -> {db_cmd};
  my $db_name = shift;
  open(my $fh, "$db_cmd details script |") or Carp::croak "ProductionMysql: $db_cmd not in your PATH\n";
  my %conn;
  while(<$fh>) {
    /host\s+(\S+)/ and $conn{host}     = $1;
    /port\s+(\d+)/ and $conn{port}     = $1;
    /user\s+(\S+)/ and $conn{user}     = $1;
    /pass\s+(\S+)/ and $conn{password} = $1;
  }
  $conn{dbname}=$db_name if $db_name;
  return \%conn;
}
sub url {
  my $db_cmd= shift -> {db_cmd};

  my $url;

  open(my $fh, "$db_cmd details url |") or Carp::croak "ProductionMysql: $db_cmd not in your PATH\n";
  while(<$fh>) {
    /^(mysql:\S+)/ and $url = $1;
  }

  $url =~ s/\/$//;

  return $url;
}
sub adaptor {
   my ($self, $pattern , @adaptor_args) = @_;
   my $registry = 'Bio::EnsEMBL::Registry';
   $registry->load_registry_from_url( $self->url );
   return $registry->get_adaptor($self->species($pattern), 'core', @adaptor_args);
}
sub dbc {
   my ($self, $pattern , @adaptor_args) = @_;
   my $registry = 'Bio::EnsEMBL::Registry';
   $registry->load_registry_from_url( $self->url );
   return $registry->get_DBAdaptor($self->species($pattern), 'core', @adaptor_args)->dbc;
}
# Needs to map each core db to something unique and within a character limit
# Not sure exactly what the true limit is - it's a limitation in mysql schema - ceiling is 17 chars
# capturing groups across alternatives, that are separated by |
# Stuff in alternative groups should be empty
# Special case for t. pseudospiralis so that the names don't get too long
# Joint treatment of cores and comparators
#  "Up to n letters" with .{1,n}
# It needs to be executed from Java, and the two implementations are slightly different: https://docs.oracle.com/javase/7/docs/api/java/util/regex/Pattern.html#jcc
our $GOLDEN_SPECIES_REGEX_MATCH = join ("|",
  "^trichinella_pseudospiralis_(iss[0-9]+prjna257433)_core_$ENV{PARASITE_VERSION}_$ENV{ENSEMBL_VERSION}_[0-9]*\$", # $1
  "^globodera_rostochiensis_(l[0-9]+prjna695196)_core_$ENV{PARASITE_VERSION}_$ENV{ENSEMBL_VERSION}_[0-9]*\$", # $2
  "^heterobilharzia_americana_(td[0-9]+prjeb44434)_core_$ENV{PARASITE_VERSION}_$ENV{ENSEMBL_VERSION}_[0-9]*\$", # $3
  "^aphelenchoides_besseyi_(a[a-z]{3}prjna834627)_core_$ENV{PARASITE_VERSION}_$ENV{ENSEMBL_VERSION}_[0-9]*\$", # $4
  "^caenorhabditis_briggsae_([a-z]{2}[0-9]{2,4}+prjna784955)_core_$ENV{PARASITE_VERSION}_$ENV{ENSEMBL_VERSION}_[0-9]*\$", # $5
  "^(sch)istosoma_([a-z]{3})[a-z]{2,9}_(td[1-2])prjeb44434_core_$ENV{PARASITE_VERSION}_$ENV{ENSEMBL_VERSION}_[0-9]*\$", # $6,7,8
  "^(s)chmidtea_(med)iterranea_s2f19(h[1-2]prjna885486)_core_$ENV{PARASITE_VERSION}_$ENV{ENSEMBL_VERSION}_[0-9]*\$", # $9,10,11
  "^meloidogyne_chitwoodi_(r[a-z]{3}[1-2]?prjna666745)_core_$ENV{PARASITE_VERSION}_$ENV{ENSEMBL_VERSION}_[0-9]*\$", # $12
  "^(.{1,2}).*?_(.{1,2}).*?_(v[1-9].*?)_core_$ENV{PARASITE_VERSION}_$ENV{ENSEMBL_VERSION}_[0-9]*\$", # $13,14,15
  "^(.{1,2}).*?_(.{1,4}).*?_(.*?)_core_$ENV{PARASITE_VERSION}_$ENV{ENSEMBL_VERSION}_[0-9]*\$", # $16,17,18
  "^([a-z])[^_]+_([^_]+)_core_.*\$" # $19,20
);


# $GOLDEN_SPECIES_REGEX_REPLACEMENT should have as many $ as the capturing groups above! So if you
# add any regexes make sure to add some $ here:
our $GOLDEN_SPECIES_REGEX_REPLACEMENT = '$1$2$3$4$5$6$7$8$9$10$11$12$13$14$15$16$17$18$19$20';

sub core_db_to_biomart_name {
  my ($core_db) = @_;
  eval "\$core_db =~ s/$GOLDEN_SPECIES_REGEX_MATCH/$GOLDEN_SPECIES_REGEX_REPLACEMENT/"; 
  die $@ if $@;
  return $core_db;
}

sub core_db_to_local_ftp_filename {
  my ($db_cmd, $core_db) = @_;
  my $ftp_id = meta_value($db_cmd, $core_db, 'species.ftp_genome_id');
  my ($spe, $cies, $bp) = split "_", $core_db;
  my $species_name = $spe . "_" . $cies;
  return $species_name . "." . $ftp_id . "." . 'WBPS' . $ENV{PARASITE_VERSION}
}

sub core_db_to_local_ftp_path_n_filename {
  my ($db_cmd, $core_db) = @_;
  my $ftp_id = meta_value($db_cmd, $core_db, 'species.ftp_genome_id');
  my ($spe, $cies, $bp) = split "_", $core_db;
  my $species_name = $spe . "_" . $cies;
  return $species_name . "/" . $ftp_id . "/" . $species_name . "." . $ftp_id . "." . 'WBPS' . $ENV{PARASITE_VERSION}
}



1;
