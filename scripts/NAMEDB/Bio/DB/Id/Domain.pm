package Bio::DB::Id::Domain;

use strict;

sub new {
  my $class = shift;
  my $name  = shift;
  return bless {
		name        => $name,
		nametypes   => { },
		prec        => 0,
		template    => {default => '%010d'},
		id_counter  => 0,
		obj_by_id   => {},
		oby_by_name => {},
		},$class;
}

sub name { shift->{name} }

sub nameTypes {
  my $types = shift->{nametypes};
  sort { $types->{$a}{prec} <=> $types->{$b}{prec} } keys %$types
}

sub addNameType {
  my $self = shift;
  my ($nametype,$arity) = @_;
  my $types = $self->{nametypes};
  return if exists $types->{$nametype};
  $arity    ||= 1;
  $types->{$nametype}{arity} = $arity;
  $types->{$nametype}{prec}  = ++$self->{prec};
}

sub addTemplate {
  my $self = shift;
  my ($template,$species) = @_;
  $species ||= 'default';
  $self->{template}{$species} = $template;
}

sub newId {
  my $self = shift;
  my $species  = shift;
  my $template = $self->getTemplate($species);
  sprintf($template,++$self->{id_counter});
}

sub species {
  my $self    = shift;
  my @species = keys %{$self->{template}};
  grep {!/^default$/} @species;
}

sub getTemplate {
  my $self = shift;
  my $species = shift || 'default';
  $self->{template}{$species} || $self->{template}{default};
}

sub createObj {
  my $self = shift;
  my ($nametype,$name,$user,$species) = @_;
  $user = (getpwuid($<))[0] unless defined $user;

  my $arity = $self->arity($nametype) or die "no arity for $nametype";;

  if ($arity == 1) {  # if this is an arity=1 type, then we check uniqueness
    my @objs  = $self->getObjsByName($nametype,$name);
    @objs && 
      die "There is already an object of type $nametype, name $name, but $nametype should be unique";
  }

  my $template = $self->getTemplate;
  my $newId    = sprintf($template,++$self->{id_counter});
  my $obj      = Bio::DB::Id::Obj->new($newId,$self,$user);
  $obj->addName($nametype,$name);
  $self->{obj_by_id}{$newId} = $obj;
}

sub arity {
  my $self = shift;
  my $nametype = shift;
  my $types = $self->{nametypes}{$nametype} or die "unknown name type: $nametype";
  return $types->{arity};
}

sub addObjName {
  my $self = shift;
  my ($nametype,$name,$obj) = @_;
  $self->{obj_by_name}{$nametype}{$name}{$obj} = $obj;
}

sub delObjName {
  my $self = shift;
  my ($nametype,$name,$obj) = @_;
  delete $self->{obj_by_name}{$nametype}{$name}{$obj};
  delete $self->{obj_by_name}{$nametype}{$name} unless
    keys %{$self->{obj_by_name}{$nametype}{$name}} > 0;
}

sub getObjsByName {
  my $self = shift;
  my ($nametype,$name) = @_;

  my @nametypes = defined $nametype ? $nametype : $self->nameTypes;
  my @objs;
  for my $type (@nametypes) {
    my $objs = $self->{obj_by_name}{$type}{$name};
    push @objs,values %$objs;
  }
  # uniquefy
  my %seen;
  return grep {$_->live && !$seen{$_}++} @objs;
}

sub getObjById {
  my $self = shift;
  my $id   = shift;
  return $self->{obj_by_id}{$id};
}

1;
