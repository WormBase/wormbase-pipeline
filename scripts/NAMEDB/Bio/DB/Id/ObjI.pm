package Bio::DB::Id::ObjI;
use strict;
use Carp 'croak';


# add a name to the object
sub addName {
  my $self = shift;
  my ($type,$name) = @_;
  my $success;
  $self->_abstract_death;
  return $success;
}

# delete a name
sub delName {
  my $self = shift;
  my ($type,$name) = @_;
  my $success;
  $self->_abstract_death;
  return $success;
}

# merge one object into another
sub merge {
  my $self = shift;
  my $idObj = shift;
  my $success;
  $self->_abstract_death;
  return $success;
}

# split out a new object
sub split {
  my $self = shift;
  my ($nametype,$name) = @_;
  my $idObj;
  $self->_abstract_death;
  return $idObj;
}

# kill an object
sub kill {
  my $self = shift;
  my $success;
  $self->_abstract_death;
  return $success;
}

# resurrect an object
sub resurrect {
  my $self = shift;
  my $success;
  $self->_abstract_death;
  return $success;
}

# get public name
sub publicName {
  my $self = shift;
  my $name;
  $self->_abstract_death;
  return $name;
}

# get a typed name
sub typedNames {
  my $self = shift;
  my $type = shift;
  my @names;
  $self->_abstract_death;
  return \@names;
}

# get an iterator over names
sub getNamesIterator {
  my $self = shift;
  my $iterator;  # return a Bio::DB::Id::Obj::NameIteratorI;
  $self->_abstract_death;
  return $iterator;
}

# get the ultimate merged descendent(s) of this object, but only if live
sub getDescendents {
  my $self = shift;
  my $obj;  # look up descendent here, following merges, etc
  $self->_abstract_death;
  return $obj;
}

# get the version
sub version {
  my $self = shift;
  my $version;
  $self->_abstract_death;
  return $version;
}

sub live {
  my $self = shift;
  my $live;
  $self->_abstract_death;
  return $live;
}

sub getHistory {
  my $self = shift;
  my @history_items;
  $self->_abstract_death;
  return \@history_items;
}

sub _abstract_death {
  croak "unimplemented method";
}

1;
