package NameDBI;

use strict;
use DBI;
use Carp 'croak';

#--------------------- connect, disconnect -----------------
# connect to database
# pass host, username, password
sub connect ($$;$$) {
  my $class            = shift;
  my $dsn              = @_;
  my ($name,$password) = @_;
  $class->_abstract_death;
  return bless {},$class;
}

sub disconnect {
  my $self = shift;
  $self->_abstract_death;
}

#-------------------- basic accessors --------------------------

sub dsn {
  shift->_abstract_death;
}
sub user {
  shift->_abstract_death;
}

sub password {
  shift->_abstract_death;
}

#-------------------- domains ---------------------------------

sub getDomains {
  my $self = shift;
  my @result;
  $self->_abstract_death;
  return \@result;
}

sub getDomain {
  shift->_abstract_death;
}

sub setDomain {
  my ($self,$new_domain) = @_;
  $self->_abstract_death;
}

sub addDomain {
  my $self = shift;
  my ($new_domain,$template) = @_;
  $self->_abstract_death;
}

#-------------------- public_id templates--------------------------

sub assignTemplate {
  my $self     = shift;
  my ($template,$domain) = @_;
  $self->_abstract_death;
}

sub getTemplate {
  my $self   = shift;
  my $domain = shift;
  $self->_abstract_death;
}

#-------------------- domain-specific name types --------------------------
sub addNameType {
  my $self = shift;
  my ($nametype,$multivalued,$unique) = @_;
  $self->_abstract_death;
}

sub getNameTypes {
  my $self = shift;
  my $domain = shift || $self->getDomain or croak "no domain defined";
  $self->_abstract_death;
}

# return array ($multivalued,$unique)
sub getNameTypeAttributes {
  my $self     = shift;
  my $nametype = shift;
  my $domain   = shift || $self->getDomain or croak "no domain defined";
  $self->_abstract_death;
}

#-------------------- public ID creation-----------------------
sub idCreate {
  my $self = shift;
  my ($type,$name,$domain) = @_;
  $domain  ||= $self->getDomain  or croak "no domain defined";
  $self->_abstract_death;
}

#-------------------- public_id retrieval-----------------------
sub idGetByID {
  my $self = shift;
  my $public_id = shift;
  my $domain    = shift;
  $domain  ||= $self->getDomain  or croak "no domain defined";
  $self->_abstract_death;
}

sub idGetByTypedName {
  my $self = shift;
  my ($type,$name,$domain) = @_;
  $domain  ||= $self->getDomain  or croak "no domain defined";
}

sub idGetByAnyName {
  my $self = shift;
  my ($name,$domain) = @_;
  $domain  ||= $self->getDomain  or croak "no domain defined";
  $self->_abstract_death;
}

# get the live merged/split descendent(s) of this object
sub idGetDescendents {
  my $self = shift;
  my ($public_name,$domain) = @_;
  $domain  ||= $self->getDomain  or croak "no domain defined";
  $self->_abstract_death;
}

sub idPublicName {
  my $self = shift;
  my ($public_id,$domain) = @_;
  $domain  ||= $self->getDomain  or croak "no domain defined";
  $self->_abstract_death;
}

# returns array of array where each event is:
# version,what,who,date,nametype,namevalue,related_public_name
sub idGetHistory {
  my $self = shift;
  my ($public_id,$after_version,$domain) = @_;
  $domain  ||= $self->getDomain  or croak "no domain defined";
  $self->_abstract_death;
}

#-------------------- add a typed name -----------------------
sub addName {
  my $self = shift;
  my ($public_id,$nametype,$name,$domain) = @_;
  $domain  ||= $self->getDomain  or croak "no domain defined";
  $self->_abstract_death;
}

#-------------------- delete a typed name -----------------------
sub delName {
  my $self = shift;
  my ($public_id,$nametype,$name,$domain) = @_;
  $domain  ||= $self->getDomain  or croak "no domain defined";
  $self->_abstract_death;
}

sub idKill {
  my $self        = shift;
  my ($public_name,$domain) = @_;
  $domain  ||= $self->getDomain  or croak "no domain defined";
  $self->_abstract_death;
}

sub idResurrect {
  my $self        = shift;
  my ($public_name,$domain) = @_;
  $domain  ||= $self->getDomain  or croak "no domain defined";
  $self->_abstract_death;
}

# merge eaten_name into eater_name
# eaten_name identifier is killed, and eater_name inherits
# all of eaten_name's names
sub idMerge {
  my $self = shift;
  my ($eaten_name,$eater_name,$domain) = @_;
  my $dbh = $self->dbh;
  $domain  ||= $self->getDomain  or croak "no domain defined";
  $self->_abstract_death;
}

# split out a new object from current one
# and assign indicated name (like objCreate)
sub idSplit {
  my $self = shift;
  my $dbh  = $self->dbh;
  my ($public_name,$nametype,$name,$domain) = @_;
  $domain  ||= $self->getDomain  or croak "no domain defined";
  $self->_abstract_death;
}

sub idLive {
  my $self        = shift;
  my ($public_name,$domain) = @_;
  $domain  ||= $self->getDomain  or croak "no domain defined";
  $self->_abstract_death;
}

# return array of names of given type
sub idTypedNames {
  my $self = shift;
  my ($public_id,$nametype,$domain) = @_;
  $domain  ||= $self->getDomain  or croak "no domain defined";
  $self->_abstract_death;
}

# return hash of all names, where key is type name
sub idAllNames {
  my $self = shift;
  my ($public_id,$domain) = @_;
  $domain  ||= $self->getDomain  or croak "no domain defined";
  $self->_abstract_death;
}

sub _abstract_death {
  croak "unimplemented method";
}

1;
