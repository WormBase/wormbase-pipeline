package Bio::DB::IdI;

use strict;
use Carp 'croak';

# connect to database
# pass host, username, password
sub connect ($;$$) {
  my $class            = shift;
  my ($name,$password) = @_;
  $class->_abstract_death;
  return bless {},$class;
}

# namespaces are selected by "domain" and "species"
sub getAllDomains () {
  my $self = shift;
  my @domains;
  $self->_abstract_death;
  return \@domains;
}

# namespaces are selected by "domain" and "species"
sub getAllSpecies () {
  my $self = shift;
  my @species;
  $self->_abstract_death;
  return \@species;
}

# get/select current domain
sub getDomain {
  my $self = shift;
  my $domain;
  $self->_abstract_death;
  return $domain;
}

# get/select current domain
sub setDomain {
  my $self = shift;
  my $new_domain = shift;
  my $domain;
  $self->_abstract_death;
  return $domain;
}

# get/select current species
sub getSpecies {
  my $self = shift;
  my $species;
  $self->_abstract_death;
  return $species;
}

sub setSpecies {
  my $self = shift;
  my $new_species = shift;
  my $species;
  $self->_abstract_death;
  return $species;
}

# create new domain
sub newDomain {
  my $self = shift;
  my $new_domain = shift;
  $self->_abstract_death;
  return 1;
}

# create new species
sub newSpecies {
  my $self = shift;
  my $new_species = shift;
  $self->_abstract_death;
  return 1;
}

# assign a new id creation template to current combination of species and domain
sub assignTemplate {
  my $self     = shift;
  my ($template,$domain,$species) = @_;
  $self->_abstract_death;
  return 1;
}

# alternative name spaces are assigned one at a time
sub addNameType {
  my $self = shift;
  my ($nametype,$arity) = @_;
  $self->_abstract_death;
  return 1;
}

# get namespaces
# this is an array of array in which first element is namespace name and second is arity
sub getNametypes {
  my $self = shift;
  $self->_abstract_death;
  return [];
}

# create a new object
sub idCreate {
  my $self = shift;
  my ($type,$name) = @_;
  my $idObj;  # create new unique object here
  $self->_abstract_death;
  return $idObj;
}

# get an existing object
sub idGetByID {
  my $self = shift;
  my $id = shift;
  my $idObj;  # create new unique object here
  $self->_abstract_death;
  return $idObj;
}

# get an existing object by its name
sub idGetObjsByTypedName {
  my $self           = shift;
  my ($type,$name)   = shift;
  my @objs;  # fetch objects with this name here
  $self->_abstract_death;
  return \@objs;
}

sub idGetObjsByAnyName {
  my $self = shift;
  my $anyName = shift;
  my @objs;  # create new object here
  $self->_abstract_death;
  return \@objs;
}

sub _abstract_death {
  croak "unimplemented method";
}

1;


