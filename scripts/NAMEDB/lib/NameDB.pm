package NameDB;

=head1 NAME

NameDB - interface to biological name tracking database

=head1 SYNOPSIS

  use NameDB;

  #connecting to database
  my $db = NameDB->connect('name_test','lstein',undef);

  # define domains
  $db->addDomain('Protein' => 'P%10d');
  $db->addDomain('Gene');
  $db->addDomain('Sequence');

  # get domains
  my @domains = $db->getDomains;

  # set current domain
  $db->setDomain('Protein');

  #name creation templates
  $db->setTemplate('P%010d');
  my $t =  $db->getTemplate('Protein');

  $db->setDomain('Protein');
  $db->addNameType(CGC,$isUnique,$isPrimary);
  $db->addNameType(GenBank,$isUnique,$isPrimary);
  $db->addNameType(EMBL,$isUnique,$isPrimary);
  $db->addNameType(WormPep,$isUnique,$isPrimary);
  $db->addNameType(WTP,$isUnique,$isPrimary);
  my @types = $db->getNameTypes;

  ($isUnique,$isPrimary) = $db->getNameTypeAttributes('EMBL','Protein');

  my $id = $db->idCreate;
  $db->addName($id,WormPep=>'wp123');
  $db->addName($id,WTP => 'I23Z31');

  $db->delName($id,WTP => 'I23Z31');

  # test history
  $db->setDomain('Gene');
  $result = $db->idGetByTypedName(CGC=>'unc-1');
  $unc1   = $result->[0];
  my $history = $db->idGetHistory($unc1,[$after,[,$version]]);
  for my $event (@$history) {
    print join "\t",'version',@$event,"\n";
  }

  # search 'em all
  my @result = $db->idGetByAnyName('ace');

  # wildcard search
  @names = $db->idSearch(CGC => '*');

  # living, dead
  print $db->idExists($id),"\n";
  print $db->idLive($unc1),"\n";
  $db->idKill($unc1);
  $db->idResurrect($unc1);

  # merges
  $db->idMerge($unc1=>$unc2);

   # splits
  ($ace) = $db->idGetByTypedName(EMBL=>'ace');
  my $ace2 = $db->idSplit($ace);

  # union of split and merge
  my @children = $db->idGetChildren($ace);
  print join "\n",@descendents,"\n";

  # ultimate descendent following merges
  my $child = $db->idGetUltimate($ace);

  # get public name
  my $name = $db->idPublicName($id[,$domain]);

  # everything that has changed since timestamp
  my @ids  = $db->idChanged($timeStampChanged);
    # NOTE: timestamp in acedb format


=head1 DESCRIPTION

This is an interface to a MySQL-based name tracking system for
biological identifiers.  Each identifier has the following attributes:

  1. biological domain, e.g. "Protein"
  2. public ID.  The combination of public ID and domain is unique
  3. public name.  The public name is the "best" name for the object.
  4. typed name(s). One or more qualified names.  Qualified names may be
          unique, or may be sured among multiple identifiers.
  5. a live flag.  An identifier may be live or dead.
  6. a version number
  7. a modification history

Operations on identifiers include creating them, merging and splitting
them, and assigning them names.

=head1 METHODS

This section describes the methods available on this module.

=head2 Constructors

=cut

use constant SUCCESS => 0;
use constant NOPERM  => 1;
use constant BADCONN => 2;
use constant BADDOM  => 3;
use constant BADID   => 4;
use constant BADTARG => 5;
use constant BADTYPE => 6;
use constant BADNAME => 7;
use constant DUPNAME => 8;
use constant BADPASS => 9;
use constant UNIMPL  => 10;
use constant DBERR   => 11;
use constant INTERR  => 12;

use Crypt::CBC;
use constant SECRET_PASSWORD => 'ViVaAce!';
use constant DB_TYPE         => 'InnoDB';

my $CIPHER = Crypt::CBC->new({key => SECRET_PASSWORD,
			      cipher => 'DES'}) or die;

my %exceptions = (
   SUCCESS() => 'success',
   NOPERM()  => 'no permission',
   BADCONN() => 'bad connection',
   BADDOM()  => 'bad domain',
   BADID()   => 'bad id',
   BADTARG() => 'bad target',
   BADTYPE() => 'bad name type',
   DUPNAME() => 'duplicate name',
   BADPASS() => 'bad password',
   UNIMPL()  => 'unimplemented',
   DBERR()   => 'database error',
   INTERR()  => 'integrity error',
);

# ERROR CODES FROM RD:
# 0 success
# 1 no permission
# 2 bad connection
# 3 bad domain
# 4 bad ID
# 5 bad target (merge/split)
# 6 bad name type
# 7 bad name
# 8 duplicate_name
# 9 bad password
# 10 unimplemented
# 11 database error
# 12 integrity error

use strict;
use NameDBI;
use DBI;
use Carp 'croak';

use vars '@ISA','$VERSION';
@ISA = qw(NameDBI);
$VERSION = '1.10';

# package globals
my %HANDLES;
my @SCHEMA;

=over 4

=item $db = NameDB->connect($dsn [,$name,$password])

Connect to the database indicated by $dsn (using the DBI DSN
nomenclature), and optionally the indicated username and password.
Returns a NameDB object if successful.

=cut

sub ping  { return "Bonjour!" }

#--------------------- connect, disconnect -----------------
sub connect {
  my $class             = shift;
  my $dsn               = shift;
  my ($name,$password)  = @_;
  $dsn = "dbi:mysql:$dsn" unless $dsn =~ /^dbi/;
  my $self = bless { domain   => undef,
		     session  => undef,
		   },$class;
  my $db = $self->dbh($dsn,$name,$password) or $class->throw(BADCONN);
  $self;
}

=item $db->disconnect

Sever the connection.

=cut

sub disconnect {
  my $self = shift;
  $self->disconnect_dbh;
}

=item $db->initialize($overwrite_flag);

Load the schema.  If $overwrite_flag is set, then the entire database
will be zeroed!  (access permissions permitting)

=back

=cut

sub initialize {
  my $self = shift;
  my $overwrite = shift;
  my $dbh = $self->dbh;
  if ($overwrite) {
    foreach (qw(domain name_template name_type
		last_identifier primary_identifier secondary_identifier
		identifier_log)) {
      $dbh->do("drop table $_");
    }
  }
  local $/ = ";\n";
  my $table_type = DB_TYPE;
  chomp(@SCHEMA = <DATA>) unless @SCHEMA;
  foreach (@SCHEMA) {
    next unless /\S/;
    s/\$DBTYPE/$table_type/g;
    $dbh->do($_) or $self->throw(DBERR);
  }
  _scalar_result(1);
}

#-------------------- basic accessors --------------------------

=head2 Basic Accessors

=over 4

=item $dsn = $db->dsn

Return the database DSN.

=cut

sub dsn {(shift->_session)[0]}

=item $user = $db->user

Return the username.

=cut

sub user {(shift->_session)[1]}

=item $pass = $db->password

Return the password.

=back

=cut

sub password {(shift->_session)[3]}

sub _session {
  my $self = shift;
  if (@_ == 3) {
    my($dsn,$name,$password) = @_;
    $self->{session} = pack('u',$CIPHER->encrypt(join $;,$dsn,$name,$password));
    return $self->{session};
  }
  return split $;,$CIPHER->decrypt(unpack('u',$self->{session}));
}

#-------------------- domains ---------------------------------

=head2 Manipulating Domains

=over 4

=item $db->addDomain($new_domain [,$template])

Add a new domain to the database.  If unsuccessful, throws a "Can't
insert domain" exception.  Example:

   eval { $db->addDomain('Protein') } or die "can't add domain: $@";

Templates are sprintf()-style strings used by the idCreate() method to
generate unique IDs.  If not specified, a default template of
"%010d" will be used.  The template can be changed with setTemplate().

=cut

sub addDomain {
  my $self = shift;
  my ($new_domain,$template) = @_;
  my $dbh = $self->dbh;
  my $query = 'insert ignore into domain (domain_name) values (?)';
  $dbh->do($query,undef,$new_domain) or $self->throw(DBERR);
  $self->setTemplate($template,$new_domain) if defined $template;
  _scalar_result(1);
}

=item $domains = $db->getDomains

In a scalar context, returns an arrayref containing names of all
defined domains.  In a list context, returns a list of all domains.

=cut

sub getDomains {
  my $self = shift;
  my $dbh  = $self->_getDomains;
  _list_result(@$dbh);
}

sub _getDomains {
  my $self = shift;
  my $dbh  = $self->dbh;
  my $query = 'select domain_name from domain';
  my $result = $dbh->selectcol_arrayref($query) or $self->throw(DBERR);
  $result;
}

=item $db->setDomain($new_domain)

Set the current domain.  This domain will be the default in subsequent
operations if not explicitly overridden.

=cut

sub setDomain {
  my ($self,$new_domain) = @_;
  my %domains = map {$_=>1} @{$self->_getDomains};
  $self->throw(BADDOM) unless $domains{$new_domain};
  $self->{domain} = $new_domain;
  _scalar_result(1);
}

=item $domain = $db->getDomain

Get the current domain.

=back

=cut

sub getDomain {
  _scalar_result(shift->{domain});
}

#-------------------- public_id templates--------------------------

=head2 Manipulating Templates

Templates are sprintf()-style strings used by the idCreate() method to
generate unique IDs.

=over 4

=item $db->setTemplate($template [,$domain])

Assign an ID creation template to a domain.  The current domain is
used if not explicitly given.  The ID should contain the sprintf code
"%010d" which will be replaced with a unique integer.  For example:

  $db->setTemplate('P%010d','Protein');

=cut

sub setTemplate {
  my $self     = shift;
  my ($template,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $dbh = $self->dbh;
  my $query =<<END;
replace into name_template (name_template.domain_id,template_template)
  select domain.domain_id,? from domain
  where domain.domain_name=?
END
  ;
  my $result = $dbh->do($query,undef,$template,$domain) or $self->throw(DBERR);
  _scalar_result($result > 0);
}

=item $template = $db->getTemplate([$domain])

Get the id creation template assigned to $domain or the current domain
if not specified.

=back

=cut

sub getTemplate {
  my $self   = shift;
  my $domain = shift;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $query =<<END;
select template_template from name_template,domain
  where domain.domain_id = name_template.domain_id
      and domain_name=?
END
  ;
  my $dbh = $self->dbh;
  my $arrayref = $dbh->selectcol_arrayref($query,undef,$domain)
    or $self->throw(DBERR);
  _scalar_result($arrayref->[0] || "%011d");
}

#-------------------- domain-specific name types --------------------------

=head2 Domain-specific name types

Each domain can have one or more name types associated with it.  Each
name type can be restricted to occur only once within its identifier
(in which case it is "isUnique"), or can be set so that within the
entire domain only one object can have a particular value for that
type (in which case it is "globally unique").

=over 4

=item $db->addNameType($nametype,$isUnique, $isPrimary  [,$domain])

Add a name type to the domain with the indicated unique flags.  If
$domain is undefined, then the current domain is used.  The unique
flags default to false.

If the type already exists in the database, this will raise a DBERR
exception.

=cut

sub addNameType {
  my $self = shift;
  my ($nametype,$isUnique,$isPrimary,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(DBERR);

  my $dbh    = $self->dbh;
  my $domainq = 'select domain_id from domain where domain_name=?';
  my $result  = $dbh->selectcol_arrayref($domainq,undef,$domain) or $self->throw(DBERR);
  my $domain_id = $result->[0];

  my $query =<<END;
insert into name_type (domain_id,name_type_name,
                       name_type_unique,name_type_primary)
  values (?,?,?,?)
END
;
  $dbh->do($query,undef,$domain_id,$nametype,$isUnique||0,$isPrimary||0,$domain)
    or $self->throw(DBERR);
  _scalar_result(1);
}

=item $types = $db->getNameTypes([$domain])

In a scalar context return an arrayref containing all the name types
for the indicated domain.  In a list context, returns an array with
the same information.  If $domain is not provided, uses the current
domain.

=cut

sub getNameTypes {
  my $self = shift;
  my $domain = shift;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $dbh    = $self->dbh;
  my $query = <<END;
select name_type_name from name_type,domain
    where domain.domain_id = name_type.domain_id 
    and domain_name=?
END
;
  my $arrayref = $dbh->selectcol_arrayref($query,undef,$domain)
    or $self->throw(DBERR);
  _list_result(@$arrayref);
}

=item ($isUnique,$isPrimary) = $db->getNameTypeAttributes($nametype,[$domain])

Returns the attributes of the indicated name type.  If $domain is not
provided, uses the current domain.

Returns undef if the name type is invalid.

In SOAP context, returns a hash with the keys "isUnique" and "isPrimary".

=back

=cut

sub getNameTypeAttributes {
  my $self     = shift;
  my $nametype = shift;
  my $domain   = shift;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $dbh      = $self->dbh;
  my $query =<<END; 
select name_type_unique,name_type_primary from name_type,domain
   where domain.domain_id = name_type.domain_id 
   and domain_name=?
   and name_type_name=?
END
  my $arrayref = $dbh->selectall_arrayref($query,undef,$domain,$nametype)
    or $self->throw(DBERR);
  return unless @$arrayref;
  return (caller(2) && caller(2) =~ /SOAP/) ? (SOAP::Data->name(isUnique  => $arrayref->[0][0]),
					       SOAP::Data->name(isPrimary => $arrayref->[0][1]))
                                            : @{$arrayref->[0]};
}

#ar2----------------- import pre-existing data -------------------

sub import {
  my $self   = shift;
  my $domain = shift;
  my $object = shift; #ref to hash containing data
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);



  #get name template and confirm passed id
  my $name_template = $self->getTemplate($domain);
  
  
}

#-------------------- public ID creation-----------------------
=head2 Identifier Creation

=over 4

=item $new_id = $db->idCreate([$domain])

Create a new unique identifier in the indicated domain and returns the
ID.  The nametype and name are assigned to the identifier.

Example:

  $new_id = $db->idCreate;

=back

=cut

sub idCreate {
  my $self   = shift;
  my $domain = shift;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $dbh = $self->dbh;
  my $new_id;

  local $dbh->{RaiseError} = 1;
  $dbh->begin_work;
  eval {
    $new_id = $self->_idCreate(undef,$domain);
    $dbh->commit;
  };
  if ($@) {
    $dbh->rollback;
    die $@;
  }

  _scalar_result($new_id);
}

#-------------------- public_id retrieval-----------------------
=head2 Identifier Lookup

=over 4

=item $id = $db->idGetByID($id [,$domain])

Look up the ID in the database and return it if present.

=cut

sub idExists {
  my $self = shift;
  my $public_id = shift;
  my $domain    = shift;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $id        = $self->_internal_id($public_id,$domain);
  _scalar_result($id);
}

=item $ids = $db->idGetByTypedName($nametype,$name [,$domain])

Look up identifier(s) using a typed name.  There may be more than one
identifier with the indicated name. Returns an arrayref in a scalar
context, a list in a list context.

Examples:

  $ids = $db->idGetByTypedName(Genbank=>'A12345');
  @ids = $db->idGetByTypedName(Genbank=>'A12345');

=cut

sub idGetByTypedName {
  my $self = shift;
  my ($type,$name,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $query =<<END;
select object_public_id
  from primary_identifier,secondary_identifier,domain,name_type
    where domain.domain_id                    = primary_identifier.domain_id
      and primary_identifier.object_id        = secondary_identifier.object_id
	and secondary_identifier.name_type_id = name_type.name_type_id
        and object_live=1
        and domain_name=?
        and name_type_name=?
        and object_name=?
END
  ;
  my $arrayref = $self->dbh->selectcol_arrayref($query,undef,$domain,$type,$name)
    or $self->throw(DBERR);
  _list_result(@$arrayref);
}

=item $ids = $db->idGetByAnyName($name [,$domain])

As for the previous method, except that the search will return a match
on any of the nametypes defined for the domain:

  @ids = $db->idGetByAnyName('A12345');

=cut

sub idGetByAnyName {
  my $self = shift;
  my ($name,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $query =<<END;
select object_public_id
  from primary_identifier,secondary_identifier,domain
    where domain.domain_id             = primary_identifier.domain_id
      and primary_identifier.object_id = secondary_identifier.object_id
      and object_live=1
      and domain_name=?
      and object_name=?
END
  ;
  my $arrayref = $self->dbh->selectcol_arrayref($query,undef,$domain,$name)
    or $self->throw(DBERR);
  _list_result(@$arrayref);
}

=item @ids = $db->idSearch($type => $wildcardpattern [,$domain])

=item @ids = $db->idSearch($wildcardpattern)

Perform a wildcard search, optionally constrained by name type.

=cut

sub idSearch {
  my $self = shift;
  my ($type,$pattern,$domain);
  if (@_ == 1) {
    $pattern = shift;
  } else {
    ($type,$pattern,$domain) = @_;
  }
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  $pattern = '*' unless defined $pattern;
  $pattern =~ tr/*?/%_/;

  my $arrayref;

  if ($type) {
    my $query = <<END;
select distinct object_public_id
  from primary_identifier,secondary_identifier,domain,name_type
    where domain.domain_id                    = primary_identifier.domain_id
      and primary_identifier.object_id        = secondary_identifier.object_id
	and secondary_identifier.name_type_id = name_type.name_type_id
        and object_live=1
        and domain_name=?
        and name_type_name=?
        and object_name LIKE ?
END
  ;
  $arrayref = $self->dbh->selectcol_arrayref($query,undef,$domain,$type,$pattern)
    or $self->throw(DBERR);
}

else {
  my $query = <<END;
select distinct object_public_id
  from primary_identifier,secondary_identifier,domain
    where domain.domain_id             = primary_identifier.domain_id
      and primary_identifier.object_id = secondary_identifier.object_id
      and object_live=1
      and domain_name=?
      and object_name LIKE ?
END
;
  $arrayref = $self->dbh->selectcol_arrayref($query,undef,$domain,$pattern)
    or $self->throw(DBERR);
  }
  _list_result(@$arrayref);
}

=item $ids = $db->idGetChildren($id [,$domain])

Return all live descendents of the identifier, following splits and
merges.

=cut

# get the live merged/split descendent(s) of this object
sub idGetChildren {
  my $self = shift;
  my ($public_name,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $id = $self->_internal_id($public_name,$domain)
    or $self->throw(BADID);

  my @descendents;

  my $history = $self->_idGetHistory($public_name,-1,$domain);
  for my $event (@$history) {
    next unless $event->{event} eq 'mergedTo' or $event->{event} eq 'splitTo';
    my $related_obj_name = $event->{related_id};
    my $related_obj_id   = $self->_internal_id($related_obj_name,$domain);
    $self->throw(INTERR) unless $related_obj_id;
    push @descendents,$related_obj_name if $self->idLive($related_obj_name,$domain);
    push @descendents,$self->idGetChildren($related_obj_name);
  }

  _list_result(@descendents);
}

=item $ids = $db->idGetUltimate($id [,$domain])

Return the ultimate descendent of the identifier, following splits and
merges.

=back

=cut

# get the ultimate living child of this ID
sub idGetUltimate {
  my $self = shift;
  my ($public_name,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);
  return _scalar_result($public_name) if $self->idLive($public_name,$domain);

  my $id = $self->_internal_id($public_name,$domain)
    or $self->throw(BADID);

  my @descendents;

  my $history = $self->_idGetHistory($public_name,-1,$domain);
  for my $event (@$history) {
    next unless $event->{event} eq 'mergedTo' or $event->{event} eq 'splitTo';
    my $related_obj_name = $event->{related_id};
    my $related_obj_id   = $self->_internal_id($related_obj_name,$domain);
    $self->throw(INTERR) unless $related_obj_id;
    return _scalar_result($related_obj_name) if $self->idLive($related_obj_name,$domain);
    return _scalar_result($self->idGetUltimate($related_obj_name,$domain));
  }

  return;
}

=head2 Methods on Individual Identifiers

These methods return a true value on success, and raise a variety of
exceptions on failure.  Exceptions include:

   no domain defined
   no ID $id in domain $domain
   method methodName() failed: details
   can't execute query: details

The most common exception is the use of an identifier that does not
exist in the database.  To catch these errors, wrap the method call
in an eval{} block:

  my $id = eval {$db->idPublicName($id)};
  warn $@ if $@;

=over 4

=item $name = $db->idPublicName($id [,$domain])

Given an identifier's ID, return its "best" name.  The best name is
the first defined name type, where name types that are added to the
database earlier have priority.

=cut

sub idPublicName {
  my $self = shift;
  my ($public_id,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $query = <<END;
select b.object_name,c.name_type_name,a.object_id
  from primary_identifier as a,secondary_identifier as b,name_type as c,domain as d
  where a.object_id      = b.object_id
      and b.name_type_id = c.name_type_id
      and a.domain_id    = d.domain_id
      and a.object_public_id = ?
      and d.domain_name      = ?
  order by c.name_type_id
END
;
  my $arrayref = $self->dbh->selectall_arrayref($query,undef,$public_id,$domain)
    or $self->throw(DBERR);
  # ignore everything but the first column
  return _scalar_result($public_id) unless @$arrayref;
  return _scalar_result($arrayref->[0][0]);
}

=item $flag = $db->idLive($id [,$domain])

Return true if the indicated ID is alive, false otherwise.

=cut

sub idLive {
  my $self        = shift;
  my ($public_name,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $query = <<END;
select object_live
  from primary_identifier,domain
  where primary_identifier.domain_id = domain.domain_id
    and object_public_id   = ?
    and domain.domain_name = ?
END
;
  my $arrayref = $self->dbh->selectcol_arrayref($query,undef,$public_name,$domain)
    or $self->throw(DBERR);
  return $arrayref->[0] unless wantarray;
  return _scalar_result($arrayref->[0]);
}

=item $version = $db->idVersion($id [,$domain])

Return the current version of the indicated ID.  

When an ID is first created its version is zero.  It is immediately
assigned a name, bumping its version to 1.  So newly-created objects
have a version number of 1.  Each subsequent operation increments the
version by 1.

=cut

sub idVersion {
  my $self      = shift;
  my $public_id = shift;
  my $domain    = shift;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my ($id,$version) = $self->_internal_id($public_id,$domain);
  return _scalar_result($version);
}


=item @names = $db->idTypedNames($id,$nametype [,$domain])

Returns the names of the indicated type that are assigned to this
identifier.  In a scalar context will return an array reference
containing zero or more assigned names.  In a list context will return
a list of zero or more names.

For unique names, the easiest way to work with this is to call in a
list context:

  ($genbank_name) = $db->idTypedNames($id,'GenBank');

=cut

# return array of names of given type
sub idTypedNames {
  my $self = shift;
  my ($public_id,$nametype,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $dbh    = $self->dbh;

  my $query = <<END;
select object_name
  from primary_identifier,secondary_identifier,domain,name_type
    where name_type.name_type_id       = secondary_identifier.name_type_id
      and primary_identifier.object_id = secondary_identifier.object_id
      and domain.domain_id             = primary_identifier.domain_id
      and primary_identifier.object_public_id = ?
      and name_type.name_type_name            = ?
      and domain.domain_name                  = ?
END
  ;
  my $names = $dbh->selectcol_arrayref($query,undef,$public_id,$nametype,$domain)
    or $self->throw(DBERR);
  return unless @$names;
  return _list_result(@$names);
}

=item %names = $db->idAllNames($id [,$domain])

Returns a hash of the names assigned to this identifier in which the
keys are the name types and the values are the names of that type.  If
a name type has more than one value, the hash value will be an array
reference.

In a list context returns the hash directly.

Here is how to iterate through all the names assigned to an
identifier:

 sub print_names {
   my ($db,$id) = @_;
   my %typed_names = $db->allNames($id);
   print "$id:\n";
   for my $type (keys %typed_names) {
     my @values = ref $typed_names{$type} ?
                   @{$typed_names{$type}}
                 : $typed_names{$type};
     print "\t",$type,"=>",join(',',@values),"\n";
   }
 }


=cut

# return hash of all names, where key is type name
sub idAllNames {
  my $self = shift;
  my ($public_id,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $dbh    = $self->dbh;

  my $query = <<END;
select name_type_name,object_name
  from primary_identifier,secondary_identifier,domain,name_type
    where name_type.name_type_id       = secondary_identifier.name_type_id
      and primary_identifier.object_id = secondary_identifier.object_id
      and domain.domain_id             = primary_identifier.domain_id
      and primary_identifier.object_public_id = ?
      and domain.domain_name                  = ?
END
  ;
  my $names = $dbh->selectall_arrayref($query,undef,$public_id,$domain)
    or $self->throw(DBERR);
  return unless @$names;
  my %hash;
  foreach (@$names) {
    my ($type,$name) = @$_;
    if (exists $hash{$type}) {
      $hash{$type} = [$hash{$type}] unless ref $hash{$type};
      push @{$hash{$type}},$name;
    } else {
      $hash{$type} = $name;
    }
  }

  # arggh.  we have to do something different for soap
  if (defined caller(2) && caller(2) =~ /^SOAP/) {
    my @result;
    foreach my $type (keys %hash) {
      my @value = ref($hash{$type}) ? @{$hash{$type}} : $hash{$type};
      push @result,SOAP::Data->type('nameType' => {type=>$type,value=>$_}) foreach @value;
    }
    return _list_result(@result);
  }

  else {
    return wantarray ? %hash : \%hash;
  }
}

=item @history = $db->idGetHistory($id,[$after, [,$domain]])

Returns the change log for the identifier.  In a scalar context,
returns an array ref, otherwise a list of log entries.  Each item of
the log is a hash reference containing the following information:

  Element   Description
  -------   -----------
     version      version number
     event        event type
     user         username
     date         date
     related_id   related identifier id
     name_type    name type
     name_value   name value

The version number is incremented by one each time the identifier is
modified.

The event type is one of the following:
 created      delName       mergedTo    splitFrom
 addName      killed        mergedFrom
 changeName   resurrected   splitTo

The username contains the name of the user who was connected to the
database at the time the change was made.

The date is the date and time of the modification, in the format
YYYY-MM-DD HH:MM:SS

For addName, changeName, and delName events, the name type is the type
of the name affected.  For other events, this field is undef.

For addName, changeName, and delName events, the name value is the
value of the affected name. For other events, this field is undef.

For mergedTo, mergedFrom, splitTo and splitFrom events, the related
identifier ID is the ID of the identifier that was involved in the
merge or split.  For other events, this field is undef.

=cut

sub idGetHistory {
  my $self = shift;
  my $history = $self->_idGetHistory(@_);
  return wantarray ? @$history : $history unless defined caller(2) && caller(2) =~ /^SOAP/;
  my @result;
  for my $entry (@$history) {
    push @result,(SOAP::Data->type(idHistory =>$entry) )
  }
  _list_result(@result);
}

# returns array of array where each event is:
# version,what,who,date,nametype,namevalue,related_public_name
sub _idGetHistory {
  my $self = shift;
  my ($public_id,$after_version,$domain) = @_;
  $after_version = -1 unless defined $after_version;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $internal_id = $self->_internal_id($public_id,$domain) or return;

  my $query =<<END;
select log_version,log_what,log_who,date_format(log_when,"%Y-%m-%d_%H:%i:%s"),
       object_public_id,name_type_name,log_name_value
  from primary_identifier,name_type,identifier_log
    where name_type.name_type_id            = identifier_log.log_name_type
      and identifier_log.log_related_object = primary_identifier.object_id
      and identifier_log.object_id = ?
      and log_version > ?
END
  ;
  my $arrayref = $self->dbh->selectall_arrayref($query,undef,$internal_id,$after_version);
  my @result;
  foreach (@$arrayref) {
    my %hash;
    @hash{qw(version event user date related_id name_type name_value)} = @$_;
    push @result,\%hash
  }
  return \@result;
}



=item @changed_ids = $db->idChanged($timestamp [,$domain])

Returns all public IDs of objects in the current or indicated domain
that have changed since the indicated timestamp (using ACEDB timestamp
format).

=cut

sub idChanged {
  my $self = shift;
  my ($acetimestamp,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);
  my ($year,$mo,$day,$hr,$min,$sec) = $acetimestamp =~ /(\d+)-(\d+)-(\d+)(?: (\d+):(\d+):(\d+))?/;
  my $when = sprintf("%04d%02d%02d%02d%02d%02d",$year,$mo,$day,$hr||0,$min||0,$sec||0);
  my $query = <<END;
select distinct object_public_id from primary_identifier,identifier_log,domain
  where primary_identifier.domain_id=domain.domain_id
    and identifier_log.object_id=primary_identifier.object_id
    and identifier_log.log_when>=?
      and domain.domain_name=?
END
  my $arrayref = $self->dbh->selectcol_arrayref($query,undef,$when,$domain)
    or $self->throw(DBERR);
  _list_result(@$arrayref);
}

=item $db->addName($id,$nametype,$name [,$domain])

Add a name to the identifier.  $nametype must be one of the valid name
types for the domain.  If the identifier already has a name of this
type, the name will be replaced if it is unique.

=cut

#-------------------- add a typed name -----------------------
sub addName {
  my $self = shift;
  my ($public_id,$nametype,$name,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $dbh = $self->dbh;
  local $dbh->{RaiseError} = 1;
  $dbh->begin_work;
  eval {
    $self->_addName($public_id,$nametype,$name,$domain);
    $dbh->commit;
  };
  if ($@) {
    $dbh->rollback;
    die $@;
  }
  _scalar_result(1);
}

=item $db->delName($id,$nametype,$name [,$domain])

Removes the indicated name from the identifier.

=cut

#-------------------- delete a typed name -----------------------
sub delName {
  my $self = shift;
  my ($public_id,$nametype,$name,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $dbh = $self->dbh;
  local $dbh->{RaiseError} = 1;
  $dbh->begin_work;
  eval {
    $self->_delName($public_id,$nametype,$name,$domain);
    $dbh->commit;
  };
  if ($@) {
    $dbh->rollback;
    die $@;
  }
  _scalar_result(1);
}

=item $db->idKill($id [,$domain])

Kill the identifier, making its "live" flag false.

=cut

sub idKill {
  my $self        = shift;
  my ($public_name,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $dbh = $self->dbh;
  local $dbh->{RaiseError} = 1;
  $dbh->begin_work;
  eval {
    $self->_idKillResurrect(0,$public_name,$domain);
    $dbh->commit;
  };
  if ($@) {
    $dbh->rollback;
    die $@;
  }
  _scalar_result(1);
}

=item $db->idResurrect($id [,$domain])

Resurrect the identifier, making its "live" flag true.

=cut

sub idResurrect {
  my $self        = shift;
  my ($public_name,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $dbh = $self->dbh;

  local $dbh->{RaiseError} = 1;
  $dbh->begin_work;
  eval {
    $self->_idKillResurrect(1,$public_name,$domain);
    $dbh->commit;
  };
  if ($@) {
    $dbh->rollback;
    die $@;
  }
  _scalar_result(1);
}

=item $db->idMerge($source_id, $target_id [,$domain])

Merge the identifier referred to by $source_id into the identifier
referred to by $destination_id.  $source_id is killed, becoming
inactive.  All names that are not globally unique are copied from
$source_id to $destination_id.

=cut

# merge eaten_name into eater_name
# eaten_name identifier is killed, and eater_name inherits
# all of eaten_name's names
sub idMerge {
  my $self = shift;
  my ($eaten_name,$eater_name,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $dbh = $self->dbh;

  local $dbh->{RaiseError} = 1;
  $dbh->begin_work;
  eval {
    $self->_idMerge($eaten_name,$eater_name,$domain);
    $dbh->commit;
  };
  if ($@) {
    $dbh->rollback;
    die $@;
  }
  _scalar_result(1);
}

=item $new_id = $db->idSplit($source_id,$nametype,$name [,$domain])

Split the identifier given by $source_id into two, returning the new
identifier.  The new identifier is given the nametype and name
indicated, but does not inherit any other names automatically.

=cut

# split out a new object from current one
# and assign indicated name (like objCreate)
sub idSplit {
  my $self = shift;
  my ($public_name,$domain) = @_;
  $domain  = $self->getDomain  unless defined $domain;
  $domain  or $self->throw(BADDOM);

  my $dbh  = $self->dbh;

  my $new_id;

  local $dbh->{RaiseError} = 1;
  $dbh->begin_work;
  eval {
    $new_id = $self->_idSplit($public_name,$domain);
    $dbh->commit;
  };
  if ($@) {
    $dbh->rollback;
    die $@;
  }
  _scalar_result($new_id);
}

=back

=head2 Utility Methods

=over 4

=item $handle = $db->dbh;

Return the DBI database handle.

=back

=cut

#---------------- fetch database handle (semi-public) -----------------------
sub dbh {
  my $self = shift;
  if (my $dbh = $self->_dbh) {
    return $dbh if $dbh->ping;
  }
  my ($dsn,$name,$password) = $self->{session} ? $self->_session() : @_;
  $name     ||= '';
  $password ||= '';
  my $db = eval { DBI->connect($dsn,$name,$password,{PrintError=>0,AutoCommit=>1})};
  my $connection_string        = $self->_session($dsn,$name,$password);
  $HANDLES{$connection_string} = $db;
  $db or DBI->errstr =~ /denied/i ? $self->throw(BADPASS) : $self->throw(BADCONN);
}

sub _dbh {
  my $self = shift;
  return unless ref $self && $self->{session};
  $HANDLES{$self->{session}};
}

sub DESTROY { shift->disconnect }

#--------------------------- internal ----------------------

sub _typeid {
  my $self     = shift;
  my $nametype = shift;
  my $domain   = shift || $self->getDomain;

  my $dbh      = $self->dbh;
  my $query =<<END;
select name_type_id 
  from name_type,domain
  where domain.domain_id = name_type.domain_id 
      and domain_name=? 
      and name_type_name=?
END
;
  my $arrayref = $dbh->selectcol_arrayref($query,undef,$domain,$nametype)
    or $self->throw(DBERR);
  $arrayref->[0];
}

sub _idCreate {
  my $self = shift;
  my ($split_from,$domain) = @_;
  my $dbh     = $self->dbh;

  # perhaps create new entry
  my $query =<<END;
insert ignore into last_identifier (last_identifier.domain_id)
  select domain_id from domain
  where domain.domain_name=?
END
  $dbh->do($query,undef,$domain);

  # get current value
  $query =<<END;
select last_identifier.object_id,domain.domain_id
  from last_identifier,domain
  where last_identifier.domain_id = domain.domain_id
      and domain_name=?
END

  my $arrayref   = $dbh->selectall_arrayref($query,undef,$domain);
  @$arrayref or $self->throw(BADID);
  my ($id_counter,$domain_id) = @{$arrayref->[0]};

  my $template = $self->getTemplate($domain);
  my $new_id   = sprintf($template,$id_counter);

  my $internal_id = $self->_createObj($new_id,$split_from,$domain);

  # bump counter -- best to do this last
  $query =<<END;
update last_identifier set last_identifier.object_id=last_identifier.object_id+1
  where last_identifier.domain_id=?
END
  ;
  $dbh->do($query,undef,$domain_id);

  return $new_id;
}

sub _idKillResurrect {
  my $self = shift;
  my $live = shift;
  my ($public_id,$domain) = @_;

  $domain = $self->getDomain unless defined $domain;
  my $event = $live ? 'resurrected' : 'killed';
  my $value = $live ? 1 : 0;

  my $internal_id = $self->_internal_id($public_id,$domain)
    or $self->throw(BADID);

  my $query = <<END;
update primary_identifier
  set object_live = $value
    where object_id = ?
END
;
  $self->dbh->do($query,undef,$internal_id)
    or $self->throw(DBERR);
  $self->_add_history($event,$public_id,undef,undef,undef,$domain);
}

sub _idMerge {
  my $self = shift;
  my ($eaten_name,$eater_name,$domain) = @_;
  my $dbh = $self->dbh;

  my $eater_id = $self->_internal_id($eater_name,$domain)
    or $self->throw(BADTARG);

  my $eaten_id = $self->_internal_id($eaten_name,$domain)
    or $self->throw(BADID);

  # kill $eaten and mark it merged
  my $query = <<END;
update primary_identifier
  set object_live = 0
    where object_id = ?
END
;
  $dbh->do($query,undef,$eaten_id)
    or $self->throw(DBERR);
  $self->_add_history('mergedTo',$eaten_name,undef,undef,$eater_name,$domain);

  # we will inherit all names from the eaten object
  # this is ugly because names may be multiple, hence
  # the hashes.
  my %eaten_names = $self->idAllNames($eaten_name,$domain);
  my %eater_names = $self->idAllNames($eater_name,$domain);
  for my $type (keys %eaten_names) {
    my ($isUnique,$isPrimary) = $self->getNameTypeAttributes($type);

    # don't replace isUnique names if already present in eater
    next if $isUnique and exists $eater_names{$type};

    my %eater = map {$_=>1} ref $eater_names{$type} ? @{$eater_names{$type}}
                                                    : $eater_names{$type};
    my %eaten = map {$_=>1} ref $eaten_names{$type} ? @{$eaten_names{$type}}
                                                    : $eaten_names{$type};
    foreach (keys %eaten) {
      $self->_delName($eaten_name,$type,$_,$domain,1);
      next if $eater{$_};  # don't bother replacing
      $self->_addName($eater_name,$type,$_,$domain,1);
    }
  }

  # Make sure that the eater is alive!!!!
  $query = <<END;
update primary_identifier
  set object_live = 1
    where object_id = ?
END
;
  $dbh->do($query,undef,$eater_id) or $self->throw(DBERR);

  $self->_add_history('mergedFrom',$eater_name,undef,undef,$eaten_name,$domain);
}

# return new object
sub _idSplit {
  my $self = shift;
  my ($public_name,$domain) = @_;
  my $dbh  = $self->dbh;

  # create new object -- this handles split entry
  my $new_name = $self->_idCreate($public_name,$domain)
    or $self->throw(DBERR);

  $self->_add_history('splitTo',$public_name,undef,undef,$new_name,$domain);

  $new_name;
}

sub _createObj {
  my $self      = shift;
  my ($public_id,$split_from,$domain) = @_;
  $domain    = $self->getDomain unless $domain;

  my $dbh = $self->dbh;
  my $query =<<END;
insert into primary_identifier (object_public_id,primary_identifier.domain_id)
  select ?,domain_id from domain
  where domain.domain_name=?
END

  $dbh->do($query,undef,$public_id,$domain)
    or $self->throw(DBERR);

  my $action = $split_from ? 'splitFrom' : 'created';
  $self->_add_history($action,$public_id,undef,undef,$split_from,$domain);
}

sub _addName {
  my $self = shift;
  my ($public_id,$nametype,$name,$domain,$no_history) = @_;

  my $id    = $self->_internal_id($public_id,$domain)
    or $self->throw(BADID);
  my $dbh   = $self->dbh;
  my ($isUnique,$isPrimary) = $self->getNameTypeAttributes($nametype,$domain)
    or $self->throw(BADTYPE);
  my $typeid = $self->_typeid($nametype,$domain);

  if ($isPrimary) {
    # check whether any other objects already have this name
    my $oldObjs = $self->idGetByTypedName($nametype,$name,$domain);
    $self->throw(DUPNAME) if @$oldObjs;
  }

  if ($isUnique && $self->idTypedNames($public_id,$nametype,$domain)) {
    return $self->_changeName($public_id,$nametype,$name,$domain);
  }

  my $query =<<END;
  insert into secondary_identifier (object_id,secondary_identifier.name_type_id,object_name)
    select primary_identifier.object_id,?,?
      from primary_identifier
	where object_public_id=?
END
  $dbh->do($query,undef,$typeid,$name,$public_id) or $self->throw(DBERR);
  $self->_add_history('addName',$public_id,$typeid,$name,undef,$domain) unless $no_history;
}

# this should only get called for unique name types!!!
sub _changeName {
  my $self = shift;
  my ($public_id,$nametype,$name,$domain) = @_;
  my $typeid = $self->_typeid($nametype,$domain)
    or $self->throw(BADTYPE);

  my $objectid = $self->_internal_id($public_id,$domain) or $self->throw(BADID);

  my $query =<<END;
update secondary_identifier
  set object_name=?
  where name_type_id = $typeid
    and object_id    = $objectid
END
;
  $self->dbh->do($query,undef,$name);
  $self->_add_history('changeName',$public_id,$typeid,$name,undef,$domain);
}

sub _delName {
  my $self = shift;
  my ($public_id,$nametype,$name,$no_history,$domain) = @_;

  my $dbh = $self->dbh;
  my $id = $self->_internal_id($public_id,$domain);
  my $typeid = $self->_typeid($nametype,$domain);

  my $query = <<END;
delete from secondary_identifier 
  where object_id=?
    and secondary_identifier.name_type_id=? 
    and object_name=?
END
;

  $dbh->do($query,undef,$id,$typeid,$name)
    or $self->throw(DBERR);

  $dbh->_add_history('delName',$public_id,$nametype,$name,undef,$domain) unless $no_history;
}

sub _internal_id {
  my $self      = shift;
  my $public_id = shift;
  my $domain    = shift || $self->getDomain;

  my $query =<<END;
select object_id,object_version 
  from primary_identifier,domain
  where domain.domain_id = primary_identifier.domain_id
      and object_public_id=?
      and domain_name=?
END
  my $arrayref = $self->dbh->selectall_arrayref($query,undef,$public_id,defined($domain) ? $domain : $self->getDomain);
  return unless @$arrayref;
  my ($id,$version) = @{$arrayref->[0]};
  return wantarray ? ($id,$version) : $id;
}

sub _add_history {
  my $self = shift;
  my ($operation,$public_name,$type,$name,$related_name,$domain) = @_;
  my $user                            = $self->user;
  my ($internal_id,$internal_version) = $self->_internal_id($public_name,$domain);
  my $related_id = $self->_internal_id($related_name,$domain);

  my $query =<<END;
insert into identifier_log (identifier_log.object_id,log_version,
			    log_what,log_who,
                            log_related_object,
			    log_name_type,log_name_value)
       values(?,?,?,?,?,?,?)
END
  ;
  my $dbh = $self->dbh;
  eval {$dbh->do($query,undef,$internal_id,$internal_version,
		 $operation,$user,$related_id||1,$type||1,$name)}
	  or $self->throw(DBERR);
  $self->_version_bump($internal_id);
}

sub _version_bump {
  my $self         = shift;
  my $internal_id = shift;
  $self->dbh->do('update primary_identifier set object_version=object_version+1 where object_id=?',
		 undef,
		 $internal_id)
    or $self->throw(DBERR);
}

sub disconnect_dbh {
  my $self = shift;
  local $^W = 0;
  my $db = $self->dbh or return;
  $db->disconnect;
  delete $HANDLES{$self->{session}};
  undef $self->{session};
  _scalar_result(1);
}

sub throw {
  my $self    = shift;
  my $errcode = shift;
  my $message = $exceptions{$errcode} || "programmer's error";
  my $syserr  = $! || '';
  my $dberr   = $DBI::errstr;
  my $dberrno = $DBI::err;
  $dberr   ||= ''; # prevent uninitialized var warnings
  $dberrno ||= '';
  my $chain = '';
  my $frame = 0;
  while (my($package,$filename,$line) = caller($frame++)) {
    $chain .= "From ($frame): $package [$filename line $line]\n";
  }
  chomp $chain;
  die <<END;
EXCEPTION
Errcode:    $errcode
Message:    $message
Errno:      $syserr
DB errcode: $dberrno
DBI err:    $dberr
$chain
END
}

sub _scalar_result {
  my $value = shift;
  my ($pack) = caller(2);
  if (defined($pack) && $pack =~ /^SOAP/ && SOAP::Data->can('name')) {
    return SOAP::Data->name('Result'=>$value);
  } else {
    return $value;
  }
}

# this detects whether we are in a SOAP context and
# returns the appropriate value
sub _list_result {
  my @values = @_;
  my ($pack) = caller(2);
  if (defined($pack) && $pack =~ /^SOAP/ && SOAP::Data->can('name')) {
    return SOAP::Data->name('Result'=>\@values);
  } else {
    return wantarray ? @values : \@values;
  }
}

1;

=head1 SCHEMA

Here's how to create an empty database:

 % mysqladmin create name_test;
 # set username and password here

 % perl -MNameDB -e "NameDB->connect('name_test','username','passwd')->initialize(1)"

=head1 AUTHOR

Lincoln Stein <lstein@cshl.org>.

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

__DATA__

CREATE TABLE IF NOT EXISTS domain (
       domain_id    int(11)  auto_increment primary key,
       domain_name  char(80) not null,
       unique(domain_name)
) type=$DBTYPE;

CREATE TABLE IF NOT EXISTS name_template (
       template_id       int(11) auto_increment primary key,
       domain_id         int(11) not null,
       template_template char(80) not null,
       unique(domain_id)
) type=$DBTYPE;

CREATE TABLE IF NOT EXISTS name_type (
       name_type_id	         int(11) auto_increment primary key,
       domain_id	         int(11) not null,
       name_type_name	         char(80) not null,
       name_type_unique          tinyint default 0,
       name_type_primary         tinyint default 0,
       unique(domain_id,name_type_name)
) type=$DBTYPE;

CREATE TABLE IF NOT EXISTS last_identifier (
       domain_id	int(11) not null,
       object_id	int(11) default 1,
       unique(domain_id)
) type=$DBTYPE;

CREATE TABLE IF NOT EXISTS primary_identifier (
       object_id	   int(11) auto_increment primary key,
       object_public_id	   char(80) not null,
       domain_id	   int(11) not null,
       object_version	   int(5)  default 0,
       object_live         tinyint default 1,
       unique(object_public_id,domain_id)
) type=$DBTYPE;

CREATE TABLE IF NOT EXISTS secondary_identifier (
       object_id	   int(11) not null,
       name_type_id	   int(11) not null,
       object_name	   char(80) not null,
       unique(object_name,name_type_id,object_id)
) type=$DBTYPE;

CREATE TABLE IF NOT EXISTS identifier_log (
       object_id    int(11) not null,
       log_version  int(11) not null,
       log_what	    enum('created','addName','changeName','delName',
			'killed','resurrected','mergedTo','splitTo',
			'mergedFrom','splitFrom'),
       log_who	    char(80) not null,
       log_when	    timestamp not null,
       log_related_object int(11) default 1,
       log_name_type      int(11) default 1,
       log_name_value	  char(80),
       index(object_id)
) type=$DBTYPE;

INSERT INTO primary_identifier values(1,'',0,0,0);
