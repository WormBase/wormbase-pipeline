package Bio::DB::Id;

use strict;
use Bio::DB::IdI;
use Bio::DB::Id::Obj;
use Bio::DB::Id::Domain;

use vars '@ISA';
@ISA = qw(Bio::DB::IdI);

sub new {
  my $class            = shift;
  my ($name,$password) = @_;
  return bless {
		user      => $name,
		password  => $password,
		domains   => {},   # indexed by domain
		species   => {},   # indexed by species
		current_domain  => undef,
		current_species => undef,
	       },$class;
}

sub getAllDomains {
  my $self = shift;
  my @domains = keys %{$self->{domains}};
  \@domains;
}

sub getAllSpecies {
  my $self = shift;
  my @species = keys %{$self->{species}};
  \@species;
}

sub getDomain {
  shift->{current_domain};
}

sub setDomain {
  my ($self,$new_domain) = @_;
  die "Invalid domain: $new_domain" unless $self->{domains}{$new_domain};
  $self->{current_domain} = $new_domain;
}

sub getSpecies {
  shift->{current_species};
}

sub setSpecies {
  my $self = shift;
  my $new_species = shift;
  my $species_hash = $self->{species}{$new_species} or die "Invalid species: $new_species";
  die "invalid species: $new_species" unless $species_hash;
  $self->{current_species} = $new_species;
}

sub newDomain {
  my $self = shift;
  my $new_domain = shift;
  $self->{domains}{$new_domain} and die "domain already exists: $new_domain";
  $self->{domains}{$new_domain} = Bio::DB::Id::Domain->new($new_domain);
  1;
}

sub newSpecies {
  my $self = shift;
  my $new_species = shift;
  $self->{species}{$new_species} and die "species already exists: $new_species";
  $self->{species}{$new_species} = {};
  1;
}

sub assignTemplate {
  my $self = shift;
  my ($template,$domain,$species) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $species = $self->getSpecies unless defined $species;
  $self->{domains}{$domain}->addTemplate($template,$species);
}

sub addNameType {
  my $self = shift;
  my ($nametype,$arity) = @_;
  $self->_get_domain_info->addNameType($nametype,$arity);
}

sub getNameTypes {
  my $self = shift;
  my @types = $self->_get_domain_info->nameTypes;
  \@types;
}

sub idCreate {
  my $self = shift;
  my ($type,$name) = @_;

  my $domain = $self->_get_domain_info;
  my $user   = $self->{user};
  my $species = $self->getSpecies;
  $domain->createObj($type,$name,$user,$species);
}

sub idGetByID {
  my $self = shift;
  my $id   = shift;
  my $domain = $self->_get_domain_info or die "No current domain";
  $domain->getObjById($id);
}

sub idGetObjsByTypedName {
  my $self           = shift;
  my ($type,$name)   = @_;
  my $domain = $self->_get_domain_info or die "No current domain";
  my @objs = $domain->getObjsByName($type,$name);
  \@objs;
}

sub idGetObjsByAnyName {
  my $self           = shift;
  my $name           = shift;
  my $domain = $self->_get_domain_info or die "No current domain";
  my @objs = $domain->getObjsByName(undef,$name);
  \@objs;
}

sub _get_domain_info {
  my $self           = shift;
  my $domain = $self->getDomain or die "No current domain";
  return $self->{domains}{$domain};
}

1;

