use strict;
use Carp;
package ProductionMysql;
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
sub staging_writable {
  my $v = $ENV{PARASITE_STAGING_MYSQL} or die "PARASITE_STAGING_MYSQL not in env. You need to module load parasite_prod_relx";
  return new(shift,"$v-ensrw"); 
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

1;
