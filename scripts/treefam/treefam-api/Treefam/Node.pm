# Author: jkh1
# 2006-08-22
#

=head1 NAME

 Treefam::Node.pm

=head1 SYNOPSIS




=head1 DESCRIPTION

 Representation of a node in a Treefam tree. Though Treefam now uses
 TFF trees, the NHX format is still supported.
 NHX tags valid in Treefam are:

    B     bootstrap value
    S     taxon name
    D     'Y' for duplication, 'N' for speciation, or not defined
    O     sequence ID
    G     gene ID
    E     gene loss stored in a string like "$-Eutheria-DROME"

=head1 CONTACT

 jkh1@sanger.ac.uk

=cut


package Treefam::Node;


use strict;
use Scalar::Util qw(weaken);
use Treefam::Tree;
use Carp;

{
  # map tags to methods
  my %_valid_tags = (O=>'sequence_id', Sd=>'', B=>'bootstrap', S=>'taxon', G=>'gene', E=>'gene_loss', D=>'is_duplication', Com=>'');
  sub _is_valid_tag {
    my $tag = shift;
    return $_valid_tags{$tag};
  }

  my $_internalID = 0;
  sub _get_new_internalID {
    return ++$_internalID;
  }
}

=head2 new

 Arg1: Treefam::DBConnection object
 Arg2: (optional) node as an unblessed hash ref,
       used by Treefam::Tree->new
 Description: Creates a new Treefam node object.
 Returntype: Treefam::Node

=cut

sub new {

  my $class = shift;
  my $dbc = shift if @_;
  my $self;
  if (@_) {
    $self = shift;
  }
  else {
    $self = {};
  }
  $self->{'DBConnection'} = $dbc;
  weaken($self->{'DBConnection'});

  $self->{'internalID'} = $class->_get_new_internalID;

  bless ($self, $class);

  return $self;
}

=head2 id

 Arg: (optional) string
 Description: Gets/sets the node id. Synonym for name.
 Returntype: string

=cut

sub id {

  my $self = shift;
  $self->{N} = shift if @_;

  return $self->{N};
}

=head2 internalID

 Description: Gets the id used internally to uniquely identify
              the node.
 Returntype: integer

=cut

sub internalID {

  my $self = shift;
  return $self->{'internalID'};
}

=head2 name

 Arg: (optional) string
 Description: Gets/sets the node name. Synonym for id.
 Returntype: string

=cut

sub name {

  my $self = shift;
  $self->{N} = shift if @_;

  return $self->{N};

}

=head2 parent

 Arg: (optional) Treefam::Node object
 Description: Gets/sets the node's parent.
 Returntype: Treefam::Node object

=cut

sub parent {

  my $self = shift;
  $self->{P} = shift if @_;

  return $self->{P};
}

=head2 taxon

 Arg: (optional) string
 Description: Gets/sets the node's taxon.
 Returntype: string

=cut

sub taxon {

  my $self = shift;
  $self->{S} = shift if @_;
  if (!$self->{S} && $self->{G}) {
    my $geneID = $self->{G};
    my $query = qq( SELECT DISTINCT s.tax_id,s.taxname
                    FROM species s, genes g
                    WHERE g.GID = ?
                    AND g.tax_id=s.tax_id);
    my $dbh = $self->{'DBConnection'}->get_DatabaseHandle;
    my $sth = $dbh->prepare($query);
    $sth->execute($geneID);
    ($self->{taxon_id},$self->{S}) = $sth->fetchrow_array;
    $sth->finish;
  }

  return $self->{S};
}

=head2 get_taxon_id

 Description: Gets the node's taxon id.
 Returntype: integer

=cut

sub get_taxon_id {

  my $self = shift;
  if (!$self->{taxon_id}) {
    my $dbh = $self->{'DBConnection'}->get_DatabaseHandle;
    my $taxon = $self->taxon;
    my $query = qq(SELECT TAX_ID FROM species WHERE SWCODE='$taxon' OR TAXNAME='$taxon');
    my $sth = $dbh->prepare($query);
    $sth->execute;
    ($self->{taxon_id}) = $sth->fetchrow_array;
    $sth->finish;
  }
  return $self->{taxon_id};
}


=head2 is_duplication

 Arg: (optional) character: Y or N
 Description: Gets/sets the node's status as a duplication event
 Returntype: 1 or 0

=cut

sub is_duplication {

  my $self = shift;
  $self->{D} = shift if @_;

  my $is_dup = uc($self->{D}) eq 'Y' ? 1 : 0;

  return $is_dup;

}

=head2 is_speciation

 Arg: (optional) character: Y or N
 Description: Gets/sets the node's status as a speciation event
 Returntype: 1 or 0

=cut

sub is_speciation {

  my $self = shift;
  $self->{D} = shift if @_;

  my $is_spec = uc($self->{D}) eq 'N' ? 1 : 0;

  return $is_spec;

}

=head2 is_leaf

 Description: Gets/sets the node's status as a leaf. A leaf is
              a node with no children
 Returntype: 1 or 0

=cut

sub is_leaf {

  my $self = shift;
  $self->{'is_leaf'} = shift if @_;
  if (!$self->{'is_leaf'}) {
    # a leaf has no children
     $self->{'is_leaf'} = $self->children ? 0:1;
  }
  return $self->{'is_leaf'};

}

=head2 gene

 Arg: (optional) Treefam::Gene
 Description: Gets/sets the gene attached to the node.
              Works only for leaves.
 Returntype: Treefam::Gene

=cut

sub gene {

  my $self = shift;
  if (@_) {
    $self->{'gene'} = shift;
    $self->{G} = ref($self->{'gene'})?$self->{'gene'}->ID:$self->{'gene'};
  }
  if (!$self->{'gene'}) {
    my ($geneID) = $self->{G};
    my $dbc = $self->{'DBConnection'};
    my $gh = $dbc->get_GeneHandle;
    if ($geneID) {
      $self->{'gene'} = $gh->get_by_id($geneID);
    }
    else {
      # try getting gene using sequence id
      my $seqID = $self->sequence_id;
      if ($seqID) {
	my $query = qq(SELECT GID FROM genes WHERE ID='$seqID');
	my $dbh = $dbc->get_DatabaseHandle;
	my $sth = $dbh->prepare($query);
	$sth->execute;
	($geneID) = $sth->fetchrow_array;
	$self->{'gene'} = $gh->get_by_id($geneID) if $geneID;
      }
    }
  }

  return $self->{'gene'};
}

=head2 children

 Arg: (optional) list of Treefam::Node objects
 Description: Gets/sets the children of the node.
 Returntype: list of Treefam::Node objects

=cut

sub children {

  my $self = shift;
  @{$self->{C}} = @_ if @_;

  return @{$self->{C}} if $self->{C};

}

=head2 sequence_id

 Arg: (optional) string
 Description: Gets/sets the id of the sequence attached to
              the node.
 Returntype: string

=cut

sub sequence_id {

  my $self = shift;
  $self->{O} = shift if @_;

  return $self->{O};

}

=head2 bootstrap

 Arg: (optional) integer
 Description: Gets/sets the bootstrap value attached to the node.
 Returntype: integer

=cut

sub bootstrap {

  my $self = shift;
  $self->{B} = shift if @_;

  return $self->{B};

}

=head2 branch_length

 Arg: (optional) double
 Description: Gets/sets the branch length value attached to
              the node.
 Returntype: double

=cut

sub branch_length {

  my $self = shift;
  $self->{dist} = shift if @_;

  return $self->{dist};

}

=head2 gene_loss

 Arg: (optional) string
 Description: Gets/sets the gene/taxon loss event attached to
              the node
 Returntype: string

=cut

sub gene_loss {

  my $self = shift;
  if (@_) {
    $self->{E} = shift;
    $self->{E} = "$-".$self->{E} unless ($self->{E}=~/^\$-/);
  }
  if ($self->{E}) {
    (my $loss = $self->{E})=~s/^\$-//;
    return $loss;
  }
  return undef;
}

=head2 get_all_tags

 Description: Gets all the tags attached to the node.
 Returntype: list of string

=cut

sub get_all_tags {

  my $self = shift;

  return grep { _is_valid_tag($_) } keys %{$self};

}

=head2 get_tag_value

 Arg: (optional) string
 Description: Gets the value associated with the given tag
              attached to the node.
 Returntype: string

=cut

sub get_tag_value {

 my $self = shift;
 my $tag = shift;
 croak "ERROR: Need a tag" unless $tag;
 croak "ERROR: Not a valid tag" unless _is_valid_tag($tag);
 return $self->{$tag};

}

=head2 get_all_ancestors

 Description: Gets list of all ancestors of this node. The list
              is ordered from parent to grand-parent, etc up to
              the root of the tree.
 Returntype: list of Treefam::Node objects

=cut

sub get_all_ancestors {

  my $self = shift;
  my $parent = $self->parent;
  my @ancestors;
  while ($parent){
    push @ancestors,$parent;
    $parent = $parent->parent;
  }
  return @ancestors;

}

=head2 height

 Arg: (optional) double
 Description: Gets/sets distance of node to root of the tree as
              cumulative sum of branch lengths
 Returntype: double

=cut

sub height {

  my $self = shift;
  $self->{'height'} = shift if @_;
  if (!$self->{'height'}) {
    my $node = $self;
    my $height = $node->branch_length;
    while( my $parent = $node->parent) {
      $height += $parent->branch_length if ($parent->branch_length);
      $node = $parent;
    }
    $self->{'height'} = $height;
  }
  return $self->{'height'};

}

=head2 species_name

 Arg1: (optional) 'latin' or 'swcode'
 Arg2: (optional) species name to attach to this node
 Description: Gets/sets species name of the node in specified
              format, default is latin name
 Returntype: string

=cut

sub species_name {

  my $self = shift;
  my $format = lc(shift) if @_;
  $format ||= 'latin';
  if ($format ne 'latin' && $format ne 'swcode') {
    croak "Argument must be 'latin' or 'swcode'";
  }
  $self->{'species'}{$format} = shift if @_;

  if (!$self->{'species'}{$format}) {
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->get_DatabaseHandle;
    my $sp = $self->get_tag_value('S');
    my $name;
    my $column = $format eq 'swcode' ? 'SWCODE' : 'TAXNAME';
    if ($self->is_leaf) {
      my $query = qq(SELECT $column FROM species WHERE TAXNAME = ? OR SWCODE = ? );
      my $sth = $dbh->prepare($query);
      $sth->execute($sp,$sp);
      ($name) = $sth->fetchrow_array;
      $column = $format eq 'swcode' ? 'SWCODE' : 'NAME';
      if (!$name) {
	$query = qq(SELECT $column FROM spec_names WHERE SWCODE = ? );
	my $sth = $dbh->prepare($query);
	$sth->execute($sp);
	($name) = $sth->fetchrow_array;
      }
      if (!$name) {
	$sp .='%';
	$query = qq(SELECT $column FROM spec_names WHERE NAME LIKE ?);
	$sth = $dbh->prepare($query);
	$sth->execute($sp);
	($name) = $sth->fetchrow_array;
      }
    }
    else {
      $name = $sp;
    }
    $self->{'species'}{$format} = $name;
  }

  return $self->{'species'}{$format};

}

=head2 sequence

 Arg1: optional, string, the type of sequence (nt: nucleotide,
       aa: amino-acid), defaults to amino-acid.
 Arg2: optional, string, the node's representative sequence
 Description: Gets/sets the node's representative sequence of given
              type used in the trees. Works only for leaves.
 Returntype: string

=cut

sub sequence {

  my $self = shift;
  my $type = shift if @_;
  $type ||='aa';
  $self->{'sequence'}{$type} = shift if @_;
  if (!defined $self->{'sequence'}{$type}) {
    my $gene = $self->gene;
    if (!$gene) {
      croak "No gene for node";
    }
    $self->{'sequence'}{$type} = $gene->sequence($type);
  }

  return $self->{'sequence'}{$type};

}


1;
