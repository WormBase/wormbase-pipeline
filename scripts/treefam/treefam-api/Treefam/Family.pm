# Author: jkh1
# 2006-08-18
#

=head1 NAME

 Treefam::Family

=head1 SYNOPSIS

 use Treefam::DBConnection;

 my $dbc = new Treefam::DBConnection ();

 my $famh = $dbc->get_FamilyHandle();

 my $family = $famh->get_by_id('TF101001');

 my $symbol= $family->symbol();
 my $date = $family->creation_date();
 my $desc = $family->description();

=head1 DESCRIPTION

 Representation of a Treefam family.

=head1 CONTACT

 jkh1@sanger.ac.uk

=cut


package Treefam::Family;


use strict;
use Treefam::DBConnection;
use Treefam::TreeHandle;
use Scalar::Util qw(weaken);

=head2 new

 Arg1: Treefam::DBConnection
 Arg2: family ID
 Arg3: type of family: familyA or familyB
 Description: Creates a new family object.
 Returntype: Treefam::Family

=cut

sub new {

  my ($class,$dbc,$familyID,$type) = @_;
  my $self = {};
  $self->{'DBConnection'} = $dbc;
  weaken($self->{'DBConnection'});
  my $dbh = $dbc->{'database_handle'};

  $self->{'ID'} = $familyID;
  $self->{'type'} = $type;

  bless ($self, $class);

  return $self;

}

=head2 ID

 Arg: optional, Treefam family ID
 Description: Gets/sets Treefam family ID
 Returntype: string

=cut

sub ID {

  my $self = shift;
  $self->{'ID'} = shift if @_;
  return $self->{'ID'};

}

=head2 AC

 Arg: optional, Treefam family ID
 Description: Synonym for ID
 Returntype: string

=cut

sub AC {

  my $self = shift;
  $self->{'ID'} = shift if @_;
  return $self->{'ID'};

}

=head2 type

 Arg: optional, family type
 Description: Gets/sets family type: familyA or familyB
 Returntype: string

=cut

sub type {

  my $self = shift;
  $self->{'type'} = shift if @_;
  return $self->{'type'};

}

=head2 symbol

 Arg: optional, family symbol
 Description: Gets/sets family's symbol.
              Only valid for type familyA.
 Returntype: string

=cut

sub symbol {

  my $self = shift;
  $self->{'symbol'} = shift if @_;
  if (!defined $self->{'symbol'}) {
    my $familyID = $self->{'ID'};
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->{'database_handle'};
    my $sth = $dbh->prepare("SELECT symbol FROM familyA WHERE AC = ?");
    $sth->execute($familyID);
    ($self->{'symbol'}) = $sth->fetchrow_array();
    $sth->finish();
  }

  return $self->{'symbol'};

}

=head2 description

 Arg: optional, family description
 Description: Gets/sets family's description.
              Only valid for type familyA.
 Returntype: string

=cut

sub description {

  my $self = shift;
  $self->{'description'} = shift if @_;
  if (!defined $self->{'description'}) {
    my $familyID = $self->{'ID'};
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->{'database_handle'};
    my $sth = $dbh->prepare("SELECT f.desc FROM familyA f WHERE f.AC = ?");
    $sth->execute($familyID);
    ($self->{'description'}) = $sth->fetchrow_array();
    $sth->finish();
  }

  return $self->{'description'};

}


=head2 creation_date

 Arg: optional, date
 Description: Gets/sets family's creation date
 Returntype: string

=cut

sub creation_date {

  my $self = shift;
  $self->{'creation_date'} = shift if @_;
  if (!defined $self->{'creation_date'}) {
    my $familyID = $self->{'ID'};
    my $type = $self->{'type'};
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->{'database_handle'};
    my $sth = $dbh->prepare("SELECT created FROM $type WHERE AC = ?");
    $sth->execute($familyID);
    ($self->{'creation_date'}) = $sth->fetchrow_array();
    $sth->finish();
  }

  return $self->{'creation_date'};

}

=head2 get_tree

 Arg: type of tree: SEED, FULL or CLEAN, defaults to FULL
 Description: Gets a tree for this family.
 Returntype: Treefam::Tree object

=cut

sub get_tree {

  my $self = shift;
  my $tree_type = shift || 'FULL';
  my $familyID = $self->{'ID'};
  my $dbc = $self->{'DBConnection'};
  my $th = $dbc->get_TreeHandle;
  my $tree = $th->get_by_id($familyID,$tree_type);

  return $tree;

}

=head2 get_genes

 Arg: optional, type of tree: FULL,SEED or CLEAN
 Description: Gets list of genes for this family.
              Limited to genes in the given tree type if any
 Returntype: list of Treefam::Gene objects

=cut

sub get_genes {

  my $self = shift;
  my $tree_type = shift if @_;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  my $familyID = $self->ID;
  my @genes;
  if (!$tree_type) {
    my $family_type = $self->type eq 'familyA' ? 'famA_gene' : 'famB_gene';
    my $sth = $dbh->prepare("SELECT DISTINCT g.GID
                             FROM genes g, $family_type t1
                             WHERE t1.AC= ? AND t1.ID=g.ID");
    $sth->execute($familyID);
    while (my ($geneID) = $sth->fetchrow_array()) {
      my $gh = $dbc->get_GeneHandle;
      my $gene = $gh->get_by_id($geneID);
      $gene->family($self);
      push @genes,$gene;
    }
  }
  else {
    my $tree = $self->get_tree(uc($tree_type));
    @genes = $tree->get_genes;
  }
  return @genes;
}

=head2 get_species

 Arg1: optional, type of species name: latin or swcode
       (5 letters species name used in Swissprot)
 Arg2: optional, type of tree: FULL,SEED or CLEAN
 Description: Gets list of species that have at least one gene
              in this family. Limited to genes in the given tree
              type if any
 Returntype: list of strings

=cut

sub get_species {

  my $self = shift;
  my $type = shift if @_;
  my $tree_type = shift if @_;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  my $familyID = $self->ID;
  my @species;
  my $taxname = lc($type) eq 'swcode' ? 'SWCODE' : 'TAXNAME';
  if (!$tree_type) {
    my $family_type = $self->type eq 'familyA' ? 'famA_gene' : 'famB_gene';
    my $sth = $dbh->prepare("SELECT DISTINCT s.$taxname
                             FROM genes g, species s, $family_type t1
                             WHERE t1.AC= ? AND t1.ID=g.ID
                             AND g.TAX_ID = s.TAX_ID");
    $sth->execute($familyID);
    while (my ($species) = $sth->fetchrow_array()) {
      push @species,$species;
    }
  }
  else {
    my $tree = $self->tree(uc($tree_type));
    @species = $tree->get_species;
  }

  return @species;
}

=head2 get_species_tree

 Arg: (optional) complete or short tree
 Description: Gets species tree for the species in this family.
              If arg is set to 'complete', the full tree up to
              the species root and with all intermediate nodes
              is returned. If arg is set to 'short', returns a
              tree with only the leaves and the internal nodes
              that have more than one child. This is the default
              behaviour.
 Returntype: Treefam::Tree object

=cut

sub get_species_tree {

  my $self = shift;
  my $type = shift if @_;
  my $dbc = $self->{'DBConnection'};
  my $species_tree = Treefam::Tree->new($dbc,undef,'species',undef);
  @{$species_tree->{_nodes}} = ();
  my %seen;
  my @node_names = $self->get_species('latin');
  my $dbh = $dbc->{'database_handle'};
  my $query = qq(SELECT t2.name FROM spec_nodes t1 LEFT JOIN spec_nodes t2 ON t1.parent_id=t2.taxon_id WHERE t1.name=?);
  my $sth = $dbh->prepare($query);
  while (my $name = shift(@node_names)) {
    my $node;
    if ($seen{$name}) {
      $node = $seen{$name};
    }
    else {
      $node = Treefam::Node->new($dbc);
      $node->name($name);
      $node->branch_length(0.1);
      $seen{$name} = $node;
    }
    push @{$species_tree->{_nodes}},$node;
    $sth->execute($name);
    my ($n) = $sth->fetchrow_array();
    if ($n) {
      my $parent;
      if ($seen{$n}) {
	$parent = $seen{$n};
      }
      else {
	$parent = Treefam::Node->new($dbc);
	$parent->name($n);
	$parent->branch_length(0.1);
	push @node_names,$n;
	$seen{$n} = $parent;
      }
      $node->parent($parent);
      push @{$parent->{C}},$node;
    }
    else { # no parent, this is the root
      $species_tree->{_root} = $node;
    }
  }

  unless (defined($type) && (lc($type) eq 'complete' || lc($type) eq 'full')) { 
    foreach my $node (@{$species_tree->{_nodes}}) {
      # delete node if it has only 1 child,
      # child becomes child of grand-parent keeping its original distance from it
      if ($node->children && scalar($node->children)==1) {
	my ($child) = $node->children;
	my $parent = $node->parent;
	my $nodeID = $node->internalID;
	my @chldrn;
	if ($parent) {
	  foreach my $n ($parent->children) {
	    if ($n->internalID eq $nodeID) {
	      push @chldrn,$child;
	      $child->{P} = $parent;
	      weaken($child->{P});
	      undef %{$n};
	      undef $n;
	    }
	    else {
	      push @chldrn,$n;
	    }
	  }
	  $parent->children(@chldrn);
	}
      }
    }

    # remove root if it has one child
    while (scalar($species_tree->{_root}->children) ==1) {
      my ($child) = $species_tree->{_root}->children;
      $species_tree->{_root} = $child;
    }
  }
  @{$species_tree->{_nodes}} = grep { $_->internalID } @{$species_tree->{_nodes}};
  return $species_tree;

}

=head2 get_domains

 Arg: (optional) e-value cut-off
 Description: Gets protein domains found in the family with
              e-value below given cut-off, default is 1e-2.
              Returned list is non-redundant.
 Returntype: list of strings (PFAM domain IDs)

=cut

sub get_domains {

  my $self = shift;
  my $cutoff = shift if @_;
  if(!$cutoff) {
    $cutoff = 1e-2;
  }
  if (!$self->{'domains'}) {
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->{'database_handle'};
    my @genes = $self->get_genes;
    my %seen;
    foreach my $gene(@genes) {
      my $geneID = $gene->sequence_id;
      my $query = qq(SELECT DISTINCT PFAM_ID FROM pfam WHERE ID= ? AND EVALUE<$cutoff);
      my $sth = $dbh->prepare($query);
      $sth->execute($geneID);
      while (my ($pfamid) = $sth->fetchrow_array()) {
	push @{$self->{'domains'}},$pfamid unless $seen{$pfamid}++;
      }
    }
  }

  return @{$self->{'domains'}} if (defined($self->{'domains'}));

}

=head2 get_hclusters

 Arg: none
 Description: Gets hclusters (a kind of superfamily) that contain
              genes from this family. Returned list is ordered by
              decreasing number of genes shared with family.
 Returntype: list of Treefam::Family objects

=cut

sub get_hclusters {

  my $self = shift;
  my @geneIDs = map { $_->sequence_id} $self->get_genes;
  my %seen;
  @seen{@geneIDs} = (1) x @geneIDs;
  my $list = "'".join("','",@geneIDs)."'";
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  my $query = qq( SELECT DISTINCT AC
                  FROM famB_gene
                  WHERE ID IN ($list)
                  AND AC>='TF500000' );
  my $sth = $dbh->prepare($query);
  $sth->execute();
  my @hclusters;
  my $famh = $dbc->get_FamilyHandle;
  while ( my ($familyID) = $sth->fetchrow_array()) {
    push @hclusters,$famh->get_by_id($familyID);
  }
  $sth->finish();
  # order hclusters by number of genes in common with the family
  my %count;
  foreach my $hcl(@hclusters) {
    my %intersect;
    my @IDs = map { $_->sequence_id} $hcl->get_genes;
    foreach my $id(@IDs) {
      $intersect{$id}++ if $seen{$id};
    }
    $count{$hcl->ID} = scalar(keys %intersect);
  }
  @hclusters = sort { $count{$b->ID}<=>$count{$a->ID} } @hclusters;

  return @hclusters;

}

1;
