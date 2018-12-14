package IdMapping::TemporaryIdGenerator;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::IdMapping::StableIdGenerator::EnsemblGeneric);

sub increment_stable_id {
  my ( $self, $lastId ) = @_;

  $lastId =~ /^Tmp([GTEP])(\d{11})$/ or  die "invalid stable ID: $lastId.";
  my $num = $2+1;
  return sprintf("Tmp%s%011d",$1,$num);

}

sub is_valid {
  my ( $self, $stableId ) = @_;

  if ( defined($stableId) && $stableId =~ /^Tmp([GTEP])(\d{11})$/)
     {return 1} else {return undef}
}

1;
