# Author: jkh1
# 2006-08-18
#

=head1 NAME

 Treefam::TreeHandle

=head1 SYNOPSIS

 use Treefam::DBConnection;

 my $dbc = new Treefam::DBConnection ();

 my $trh = $dbc->get_TreeHandle();

 my $tree = $trh->get_by_id('TF300011','FULL');

=head1 DESCRIPTION

 Enables retrieval of Treefam tree objects.

=head1 CONTACT

 jkh1@sanger.ac.uk


=cut


package Treefam::TreeHandle;


use strict;
use Carp;
use Treefam::Tree;
use Treefam::DBConnection;
use Scalar::Util qw(weaken);


=head2 new

 Arg: Treefam::DBConnection
 Description: Creates a new tree object handle
 Returntype: Treefam::TreeHandle

=cut

sub new {

  my $class = shift;
  my $self = {};
  $self->{'DBConnection'} = shift;
  weaken($self->{'DBConnection'});
  bless ($self, $class);

  return $self;

}

=head2 get_by_id

 Arg1: Treefam family database ID
 Arg2: type of tree: FULL,SEED or CLEAN
 Description: Gets tree with given Treefam database ID and type
 Returntype: Treefam::Tree object

=cut

sub get_by_id {

  my ($self,$familyID,$type) = @_;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  my $query = qq( SELECT tree FROM trees WHERE AC= ? AND type= ?);
  my $sth= $dbh->prepare ($query);
  $sth->execute($familyID,$type);
  my ($nhx) = $sth->fetchrow_array();

  if ($nhx) {
    return new Treefam::Tree($dbc,$familyID,$type,$nhx);
  }
  return undef;

}

=head2 get_by_ac

 Arg: Treefam family database ID
 Description: A synonym for get_by_id.
 Returntype: Treefam::Tree object

=cut

sub get_by_ac {

  my ($self,$familyID,$type) = @_;
  return $self->get_by_id($familyID,$type);

}

=head2 get_by_family

 Arg1: Treefam::Family object or Treefam family database ID
 Arg2: type of tree: FULL,SEED or CLEAN
 Description: Gets tree for given family and type
 Returntype: Treefam::Tree object

=cut

sub get_by_family {

  my ($self,$family,$type) = @_;
  my $familyID = ref($family) ? $family->ID() : $family;
  return $self->get_by_id($familyID,$type);

}

=head2 get_by_gene

 Arg1: Treefam::Gene object
 Arg2: type of tree: FULL, SEED or CLEAN
 Description: Gets tree of selected type containing given gene.
 Returntype: Treefam::Tree object or undef if gene doesn't
             exist in Treefam or is not in the type of tree
             requested

=cut

sub get_by_gene {

  my ($self,$gene,$type) = @_;
  my $dbc = $self->{'DBConnection'};
  my $famh = $dbc->get_FamilyHandle();
  my $geneID = ref $gene ? $gene->ID(): $gene;
  my $family = $famh->get_by_gene($geneID);
  unless ($family) {
    return undef;
  }
  my $familyID = $family->ID();
  my $tree = $self->get_by_id($familyID,$type);
  my ($node) = $tree->get_nodes_by_tag_value(-G=>$geneID) if $tree;
  if ($node) {
    return $tree;
  }
  return undef;
}

=head2 new_tree

 Arg (required): -tree => string, tree in nhx or tff format
 Arg: -acc => string, tree family ID
 Arg: -type => string, type of tree (i.e. clean, full or seed)
 Description: creates a new tree object with given attributes.
 Returntype: Treefam::Tree object

=cut

sub new_tree {

  my $self = shift;
  my %tree = @_;

  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  my $tree;
  unless ($tree{'-tree'}) {
    croak "Tree required";
  }
  $tree{'-acc'} ||= '';
  $tree{'-type'} ||= '';

  return Treefam::Tree->new($dbc,$tree{'-acc'},$tree{'-type'},$tree{'-tree'});

}

1;
