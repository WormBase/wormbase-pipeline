package Bio::DB::Id::Obj;

use strict;
use Bio::DB::Id::ObjI;
use Bio::DB::Id::Obj::History;

use vars '@ISA';
@ISA = 'Bio::DB::Id::ObjI';

sub new {
  my $class = shift;
  my ($id,$domain,$user) = @_;
  my $self = bless {
		    id       => $id,
		    domain   => $domain,
		    user     => $user,
		    version  => -1,  # creation is version 0
		    names    => {},
		    history  => [],
		    live     => 1,
		   },ref($class) || $class;
  $self->_add_history('created');
  $self;
}

sub addName {
  my $self = shift;
  my ($type,$name) = @_;
  my $arity = $self->_domain->arity($type) or die "Invalid type $type";
  if ($arity == 1 && @{$self->typedNames($type)} ) {  # replace
    my $oldname = $self->typedNames($type)->[0];
    $self->{names}{$type} = {$name=>1};
    $self->_domain->delObjName($type,$oldname,$self);
    $self->_domain->addObjName($type,$name,$self);
    $self->_add_history('changeName',$type,$name);
  } else {
    $self->{names}{$type}{$name}++;
    $self->_domain->addObjName($type,$name,$self);
    $self->_add_history('addName',$type,$name);
  }
  return 1;
}

sub delName {
  my $self = shift;
  my ($type,$name) = @_;
  $self->{names}{$type}{$name} or die "Name $name of type $type does not exist";
  delete $self->{names}{$type}{$name};
  $self->_domain->delObjName($type,$name,$self);
  $self->_add_history('delName',$type,$name);
}

sub kill {
  my $self = shift;
  $self->{live} = 0;
  $self->_add_history('kill');
}

sub resurrect {
  my $self = shift;
  $self->{live} = 1;
  $self->_add_history('resurrect');
}

sub live { shift->{live} }

sub version { shift->{version} }

sub id   { shift->{id}   }

sub typedNames {
  my $self = shift;
  my $type = shift;
  my @names;
  if (my $names = $self->{names}{$type}) {
    @names = keys %{$names};
  }
  \@names;
}

sub publicName {
  my $self = shift;
  my $domain = $self->{domain};
  my @types  = $domain->nameTypes;
  foreach (@types) {
    my $names = $self->typedNames($_);
    return $names->[0] if @$names;
  }
  return;
}


sub merge {
  my $self = shift;
  my $obj = shift;
  $self->{live} = 0;
  my $domain = $self->_domain;

  # this copies the names over
  $obj->_merge_from($self);

  # this unmaps our names
  my $nameIterator = $self->getNamesIterator;
  while (my($type,$name) = $nameIterator->nextName) {
    $domain->delObjName($type,$name,$self);
  }
  $self->_add_history('mergeTo',undef,undef,$obj);
}

sub split {
  my $self = shift;
  my ($nametype,$name) = @_;

  my $domain = $self->_domain;
  my $newID  = $domain->newId;
  my $newObj = $self->new($newID,$domain,$self->{user});
  $newObj->_add_history('splitFrom',undef,undef,$self);

  # no names are inherited... I think?
  $newObj->addName($nametype,$name);

  $self->_add_history('splitTo',undef,undef,$newObj);
  $newObj;
}

sub getNamesIterator {
  my $self = shift;
  return Bio::DB::Id::Obj::NameIterator->new($self->{names});
}

sub getHistory {
  my $self = shift;
  return $self->{history};
}

sub _merge_from {
  my $self = shift;
  my $obj  = shift;
  my $nameIterator = $obj->getNamesIterator or die "name iterator failure";
  while (my($type,$name) = $nameIterator->nextName) {
    $self->addName($type,$name) unless $self->_domain->arity($type) == 1;
  }
  $self->_add_history('mergeFrom',undef,undef,$obj);
}

sub _version_bump {
  shift->{version}++;
}

sub _domain {
  my $self = shift;
  return $self->{domain};
}

# get the ultimate merged descendent(s) of this object, but only if live
sub getDescendents {
  my $self = shift;
  return [$self] if $self->live;
  # follow all merges/splits
  my @descendents;
  my $history = $self->getHistory;
  for my $ev (@$history) {
    if ($ev->event eq 'mergeTo' or $ev->event eq 'splitTo') {
      my $id  = $ev->relatedId;
      my $obj = $self->_domain->getObjById($id);
      die ("data inconsistency, ",
	   $self->id,
	   " was merged into $id, but object is gone") unless $obj;
      push @descendents,@{$obj->getDescendents};
    }
  }
  \@descendents;
}

sub _add_history {
  my $self = shift;
  my ($what,$nametype,$namevalue,$related_object) = @_;
  $self->_version_bump;
  push @{$self->{history}},
    Bio::DB::Id::Obj::History->new(
				   $self->version,
				   $what,
				   $nametype,
				   $namevalue,
				   $related_object ? $related_object->id : undef,
				   $self->{user},
				   time,
				   );
}

