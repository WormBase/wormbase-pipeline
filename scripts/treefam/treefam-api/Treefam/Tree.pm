# Author: jkh1
# 2006-08-18
#

=head1 NAME

 Treefam::Tree

=head1 SYNOPSIS

 use Treefam::DBConnection;

 my $dbc = new Treefam::DBConnection ();

 my $trh = $dbc->get_TreeHandle();

 my $tree = $trh->get_by_id('TF101001','FULL');

 # finds duplications
 my @nodes = $tree->get_nodes_by_tag_value(-D=>'y');

 # finds nodes with bootstrap>=50
 my @nodes = $tree->get_nodes_by_tag_value(-B=>'>=50');

=head1 DESCRIPTION

 Representation of a Treefam tree.
 Treefam trees are now (from Treefam 4) stored in TFF format.
 They can be exported in the previously used NHX format.
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


package Treefam::Tree;

use strict;
use Carp;
use Treefam::DBConnection;
use Scalar::Util qw(weaken);
use Treefam::Node;

=head2 new

 Arg1: Treefam::DBConnection
 Arg2: family ID or 'species' to get a species tree
 Arg3: type of tree (only if Arg2 is a family ID)
 Arg4: optional, tree in NHX or TFF format
 Description: Creates a new tree object.
 Returntype: Treefam::Tree

=cut

sub new {

  my ($class,$dbc,$familyID,$type,$tree) = @_;

  my $self = {};
  if (!$tree && $familyID) {
    # get tree from database
    my $dbh = $dbc->{'database_handle'};
    if ($dbc->get_Database eq 'treefam_3') { # using Treefam_3
      my $query = qq( SELECT tree FROM trees WHERE AC= ? AND type= ?);
      my $sth= $dbh->prepare ($query);
      $sth->execute($familyID,$type);
      ($tree) = $sth->fetchrow_array();
    }
    else { # assume Treefam_4 or newer
      my $sth = $dbh->prepare(qq(SELECT f.tree_id FROM tff_fam f WHERE f.ac=? AND f.tree_type=?));
      $sth->execute($familyID, $type);
      my ($tree_id) = $sth->fetchrow_array;
      $sth = $dbh->prepare(qq(SELECT t.node_id, t.parent_id, t.branch_len FROM tff_topology t WHERE t.tree_id=?));
      my $sth2 = $dbh->prepare(qq(SELECT n.node_key, n.node_val FROM tff_node_info n WHERE n.tree_id=? AND n.node_id=?));
      my $str = '';
      my ($node_id, $parent_id, $branch_len);
      $sth->execute($tree_id);
      while ((($node_id, $parent_id, $branch_len) = $sth->fetchrow_array)) {
	my ($siid, $key, $val);
	my $s = '';
	$sth2->execute($tree_id, $node_id);
	while ((($key, $val) = $sth2->fetchrow_array)) {
	  $s .= "$key=$val\t";
	  if ($key eq 'O') { # add gene name
	    my $gh = $dbc->get_GeneHandle;
	    my $gene = $gh->get_by_sequence_id($val);
	    if ($gene) {
	      my $t = $gene->ID;
	      $s .= "G=$t" if ($t);
	    }
	    else {
	      warn "ERROR: Can't find gene\n";
	    }
	  }
	}
	chop($s); # chop trailing "\t"
	$str .= "$node_id\t$parent_id\t$branch_len\t$s\n";
      }
    }
  }

  $self->{'DBConnection'} = $dbc;
  weaken($self->{'DBConnection'});

  $self->{'ID'} = $familyID if $familyID;
  $self->{'type'} = $type if $type;
  %{$self->{_valid_tags}} = (O=>1, Sd=>1, B=>1, S=>1, G=>1, E=>1, D=>1, Com=>1);
  # parse tree (using Li Heng's code)
  if ($tree) {
    $self->{_error} = 0;
    @{$self->{_nodes}} = ();
    my ($array, @stack);
    $array = \@{$self->{_nodes}};
    my $str = $tree;

    if ($tree=~/^\(/) { # assume nhx format
      $self->{'nhx'} = $tree;
      # parse nhx format
      my $pos;
      $_ = (($pos = index($str, ';')) >= 0)? substr($str, 0, $pos) : $str;
      s/\s//g;
      while ($_=~/(\(|((\)?[^,;:\[\]\(\)]+|\))(:[\d.eE\-]+)?(\[&&NHX[^\[\]]*\])?))/g) {
	my $st = &_parse_nhx($self,$array,\@stack,$1,$3,$4,$5);
      }
      if (@stack != 1) {
	my $count = @stack;
	my $msg = qq(ERROR: unmatched "\(" \($count\) );
	warn($msg);
	$self->{_error} = 1;
	@stack = ();
      }
      if ($self->{_error} == 0) {
	$self->{_root} = shift(@stack);
	weaken($self->{_root});
      }
      else {
	@{$self->{_nodes}} = ();
	delete($self->{_root});
      }
      if ($self->{_root}) {
	my $j = 0;
	foreach my $p (@{$self->{_nodes}}) {
	  ++$j unless ($p->{C});
	}
	$self->{_n_leaf} = $j;
      }
    }
    elsif ($tree=~/\t/) { # assume TFF format
      $self->{'tff'} = $tree;
      # parse TFF format
      my (@data, %conv, @node, $root_index);
      my $i = 0;
      $self->DESTROY if ($self->{_root});
      @{$data[$i++]} = split("\t") foreach (split("\n", $str));
      $i = 0;
      foreach my $p (sort {$a->[0]<=>$b->[0]} @data) { # initialization
	my %hash;
	push(@node, \%hash);
	$conv{$p->[0]} = $i++;
	$root_index = $p->[0] if (!$root_index && $p->[0] == $p->[1]);
      }
      if (!$root_index) {
	warn('ERROR: tree has no root');
	$self->{_error} = 1;
	return;
      }
      $root_index = $conv{$root_index};
      foreach my $p (@data) {
	my ($chi, $par) = ($conv{$p->[0]}, $conv{$p->[1]});
	my $q = \%{$node[$chi]};
	$q->{P} = $node[$par];
	push(@{$node[$par]->{C}}, $q) if ($chi != $par);
	$q->{dist} = $p->[2];
	foreach (my $i = 3; $i < @$p; ++$i) {
	  $q->{$1} = $2 if ($p->[$i] =~ /^\s*([A-Za-z][A-Za-z0-9]*)\s*=\s*([^\n\t]+)\s*$/);
	}
      }
      $self->{_root} = $node[$root_index];
      weaken($self->{_root});
      my @stack;
      my $k = 0;
      my $array = \@{$self->{_nodes}};
      my $p = \%{$stack[0]};
      @$array = ();
      $p->{P} = $self->{_root};
      $p->{i} = 0;
      for (;;) {
	while ($p->{P}{C} && $p->{i} != @{$p->{P}{C}}) {
	  $stack[++$k]{i} = 0;
	  $stack[$k]{P} = $p->{P}{C}[$p->{i}];
	  $p = \%{$stack[$k]};
	}
	push(@$array, $p->{P});
	$p = \%{$stack[--$k]};
	if ($k >= 0) { ++$p->{i}; }
	else { last; }
      }
      $i = 0;
      for my $p (@{$self->{_nodes}}) {
	if (!$p->{N} && !$p->{C}) {
	  $p->{N} = $i;
	  warn("ERROR: leaf has no name\n");
	}
	++$i;
      }
    }
    else {
      croak "ERROR: unrecognized format\n";
    }
    # make all nodes Treefam::Node objects
    foreach my $n (@{$self->{_nodes}}) {
      $n = Treefam::Node->new($dbc,$n);
    }
  }

  bless ($self, $class);

  return $self;

}

sub _parse_nhx {
  # code from li Heng
  my ($self, $array, $stack, $str, $name, $dist, $nhx) = @_;
  if ($str eq '(') {
    push(@$stack, $str);
  }
  elsif ($name) {
    my %hash;
    if ($name =~ /^\)/) {
      my (@s, $t);
      while (($t = pop(@$stack))) {
	last if (ref($t) ne 'HASH');
	push(@s, $t);
      }
      unless (defined($t)) {
	warn('ERROR: unmatched ")"');
	$self->{_error} = 1;
	return;
      }
      foreach (@s) {
	push(@{$hash{C}}, $_);
	$_->{P} = \%hash;
	weaken($_->{P});
      }
      $hash{N} = substr($name, 1) if (length($name) > 1);
    }
    else {
      $hash{N} = $name;
    }
    $hash{dist} = substr($dist, 1) if ($dist);
    while ($nhx && $nhx=~/(:([^\s=:]+)=([^:=\[\]]+))/g) {
      $hash{$2}=$3;
      $nhx =~s/$1//;
    }
 #   $nhx =~ s/:([^\s=:]+)=([^:=\[\]]+)/$hash{$1}=$2,''/eg if ($nhx);
    push(@$stack, \%hash);
    push(@$array, \%hash);
  }
  return $str;
}

=head2 ID

 Arg: optional, Treefam family ID
 Description: Gets/sets Treefam family ID for the tree
 Returntype: string

=cut

sub ID {

  my $self = shift;
  $self->{'ID'} = shift if @_;
  return $self->{'ID'};

}

=head2 type

 Arg: optional, family type
 Description: Gets/sets tree type: FULL, SEED or CLEAN
 Returntype: string

=cut

sub type {

  my $self = shift;
  $self->{'type'} = shift if @_;
  return $self->{'type'};

}

=head2 tff

 Arg: optional, string
 Description: Tree in TFF format
 Returntype: string

=cut

sub tff {

  my $self = shift;
  if (@_) {
    $self->{'tff'} = shift;
  }
  else {
    # do this each time because tree may have been modified
    # between 2 calls to this function.
    # code from Li Heng
    my ($k, $str, $s) = (0, '', '');
    # calculate _index and _lindex
    foreach my $p (@{$self->{_nodes}}) {
      $p->{_lindex} = ($p->{C})? $p->{C}[0]->{_lindex} : $k;
      $p->{_index} = $k++;
    }
    foreach my $p (@{$self->{_nodes}}) {
      $s = '';
      foreach my $q (sort keys %$p) {
	next if (ref($p->{$q}) || $q eq 'internalID' || $q eq 'dist');
	$s .= "\t$q=".$p->{$q} if ($q !~ /^_/);
      }
      $str .= "$p->{_index}\t" . ($p == $self->root ? $p->{_index} : $p->{P}->{_index});
      $str .= "\t" . (defined($p->{dist}) ? $p->{dist} : 0) . "$s\n";
    }
    $self->{'tff'} = $str;
  }

  return $self->{'tff'};
}



=head2 nhx

 Arg: optional, string
 Description: Tree in nhx format
 Returntype: string

=cut

sub nhx {

  my $self = shift;
  if (@_) {
    $self->{'nhx'} = shift;
  }
  else {
    # do this each time because tree may have been modified
    # between 2 calls to this function.
    $self->{'nhx'} = $self->_string_nhx($self->root). ";\n";
  }

  return $self->{'nhx'};

}

sub _string_nhx {
  # code from Li Heng
  my ($self, $root) = @_;
  my $str;
  if ($root->{C}) {
    $str = '(';
    for my $p (reverse @{$root->{C}}) {
      $str .= $self->_string_nhx($p) . ",\n";
    }
    chop($str); chop($str);	# chop the trailing ",\n"
    $str .= "\n)";
    $str .= $root->{N} if ($root->{N}); # node name
    $str .= ":" . $root->{dist} if ($root->{dist}); # length
    # nhx output
    $str .= '[&&NHX';
    foreach my $p (keys %{$self->{_valid_tags}}) {
      $str .= ":$p=" . $root->{$p} if ($root->{$p});
    }
    $str .= ']';
  }
  else {			# leaf
    $str = $root->{N};
    $str .= ":" . $root->{dist} if ($root->{dist});
    $str .= '[&&NHX';
    foreach my $p (keys %{$self->{_valid_tags}}) {
      $str .= ":$p=" . $root->{$p} if ($root->{$p});
    }
    $str .= ']';
  }
  return $str;

}

=head2 get_all_nodes

 Description: Gets all nodes in the tree
 Returntype: list of Treefam::Node objects

=cut

sub get_all_nodes {

  my $self = shift;

  return @{$self->{_nodes}};

}

=head2 get_last_common_ancestor

 Arg: pair of Treefam::Node or Treefam::Gene objects
 Description: Gets last common ancestor of the given nodes/genes
 Returntype: a Treefam::Node object

=cut

sub get_last_common_ancestor {

  my $self = shift;
  my @input = @_ if @_;
  my $dbc = $self->{'DBConnection'};
  my @nodes;
  foreach my $i(@input) {
    if (ref($i)=~/Gene/) {
      my $geneID = $i->ID;
      my ($node) = $self->get_nodes_by_tag_value(-G=>$geneID);
      croak "Gene not found" unless $node;
      push @nodes,$node if ($node);
    }
    else {
      push @nodes,$i;
    }
  }
  croak "ERROR: Treefam::Node objects required" unless (@nodes);

  # after code in Bio::Tree::TreeFunctionsI
  my %seen_parent;
  my $parent = $nodes[0];

  while ($parent){
    $seen_parent{$parent->internalID} = $parent;
    $parent = $parent->parent;
  }
  $parent = $nodes[1];
  while ( $parent ){
    if ( $seen_parent{$parent->internalID} ){
      return $seen_parent{$parent->internalID};
    }
    $parent = $parent->parent;
  }
  carp("Can't find last common ancestor");
  return undef;
}

=head2 root

 Arg: optional, a Treefam::Node object
 Description: Gets/sets root of the tree
 Returntype: a Treefam::Node object

=cut

sub root {

  my $self = shift;
  if (@_) {
    $self->{_root}  = shift;
    weaken($self->{_root});
  }

  return $self->{_root};
}

=head2 get_leaves

 Description: Gets leaves of the tree
 Returntype: list of Treefam::Node objects

=cut

sub get_leaves {

  my $self = shift;

  # Leaves are nodes with no children
  return grep { !$_->{C} } @{$self->{_nodes}};
}

=head2 get_node_by_id

 Arg: string, id/name of the node in the tree
 Description: Gets a node from the tree
 Returntype: Treefam::Node objects

=cut

sub get_node_by_id {

  my $self = shift;
  my $id = shift;

  my ($node) = grep { defined($_->id) && $_->id eq $id } @{$self->{_nodes}};

  return $node;
}

=head2 get_nodes_by_tag_value

 Arg: -key => value
 Description: Gets all nodes for which the 'key' tag matches
              'value'
 Returntype: list of Treefam::Node objects

=cut

sub get_nodes_by_tag_value {

  my ($self,$tag,$value) = @_;
  $tag=~s/^-//;
  my @all_nodes = $self->get_all_nodes;
  my @found_nodes;
  foreach my $node(@all_nodes) {
    my @tags = $node->get_all_tags();
    foreach my $t(@tags) {
      next unless ($t eq $tag);
      my $v = $node->get_tag_value($t);
      # look for exact matches but also allow for things like 'B>=70'
      next unless ( lc($v) eq lc($value) || ($tag eq 'B' && $t eq 'B' && $value=~/\D/ && eval("$v$value")));
      push @found_nodes,$node;
      last;
    }
  }

  return @found_nodes;

}

=head2 get_genes

 Description: Gets genes that are in the tree
 Returntype: list of Treefam::Gene objects

=cut

sub get_genes {

  my $self = shift;
  my @genes;
  foreach my $leaf($self->get_leaves) {
    my $gene = $leaf->gene;
    $gene->family($self->family);
    push @genes, $gene if $gene;
  }
  return @genes;
}

=head2 get_leaf_by_gene

 Arg: Treefam::Gene object
 Description: Gets leaf node for the given gene
 Returntype: Treefam::Node object

=cut

sub get_leaf_by_gene {

  my $self = shift;
  my $gene = shift;
  if (!$gene || ref($gene) ne 'Treefam::Gene') {
    croak "Treefam::Gene required";
  }
  my @leaves = $self->get_leaves;
  foreach my $leaf(@leaves) {
    my @tags = $leaf->get_all_tags();
    foreach my $t(@tags) {
      next unless ($t eq 'G' || $t eq 'O');
      my $v = $leaf->get_tag_value($t);
      if ($v eq $gene->ID || $v eq $gene->sequence_id) {
	return $leaf;
      }
    }
  }

  return undef;
}


=head2 score

 Arg: (optional) double
 Description: Gets/sets score for the tree, for example maximum
              likelihood
 Returntype: double

=cut

sub score {

  my $self = shift;
  $self->{'score'} = shift if @_;
  return $self->{'score'};

}

=head2 get_distance

 Arg: list of 2 Treefam::Node objects or Treefam::Gene objects
 Description: Gets distance between 2 nodes as cumulative branch
              length
 Returntype: double

=cut

sub get_distance {

  my ($self,$node1,$node2) = @_;
  if (ref($node1) eq 'Treefam::Gene') {
    $node1 = $self->get_leaf_by_gene($node1);
  }
  if (ref($node2) eq 'Treefam::Gene') {
    $node2 = $self->get_leaf_by_gene($node2);
  }
  unless (ref($node1) eq 'Treefam::Node' && ref($node2) eq 'Treefam::Node') {
    croak "ERROR: Two nodes are required";
  }
  my $distance = $node1->branch_length;
  my $lca = $self->get_last_common_ancestor($node1,$node2);
  foreach my $node ($node1->get_all_ancestors) {
    last if ($node->internalID eq $lca->internalID);
    $distance += $node->branch_length;
  }
  $distance += $node2->branch_length;
  foreach my $node ($node2->get_all_ancestors) {
    last if ($node->internalID eq $lca->internalID);
    $distance += $node->branch_length;
  }
  return $distance;

}

=head2 delete_node

 Arg: a Treefam::Node object
 Description: Removes given node (and its children) from the
              tree.
 Returntype: the modified Treefam::Tree object

=cut

sub delete_node {

  my ($self,$node) = @_;
  unless ($node) {
    croak "ERROR: Node required";
  }
  my $status = 0;
  my $parent = $node->parent;
  my $internalID = $node->internalID;
  $self->_remove_node($node);

  # delete parent node if it is left with 1 child,
  # child becomes child of grand-parent keeping its original distance from it
  if ($parent->children && scalar($parent->children)==1) {
    my ($child) = $parent->children;
    my $branch_length;
    if ($child->branch_length && $parent->branch_length) {
      $branch_length = $child->branch_length + $parent->branch_length;
      $child->branch_length($branch_length);
    }
    # add child to grand-parent's children
    my $gp = $parent->parent;
    my $parentID = $parent->internalID;
    my @chldrn;
    if ($gp) {
      foreach my $n($gp->children) {
	if  ($n->internalID eq $parentID) {
	  push @chldrn,$child;
	  $child->{P} = $gp;
	  weaken($child->{P});
	  undef %{$n};
	  undef $n;
	  $status++;
	}
	else {
	  push @chldrn,$n;
	}
      }
      $gp->children(@chldrn);
    }
  }
  # remove root if it has one child
  while (scalar($self->{_root}->children) ==1) {
    my ($child) = $self->{_root}->children;
    $self->{_root} = $child;
  }

  # remove undefined nodes from list
  @{$self->{_nodes}} = grep {$_->internalID} @{$self->{_nodes}};

  return $self;

}

sub _remove_node {

  my ($self,$node) = @_;
  unless ($node) {
    croak "ERROR: Node required";
  }
  my $parent = $node->parent;
  my $internalID = $node->internalID;
  while (my $child=shift(@{$node->{C}})) {
    $self->_remove_node($child);
  }
  # remove node from list of its parent's children
  @{$parent->{C}} = grep { $_->internalID ne $internalID } @{$parent->{C}};
  undef %{$node};
  undef $node;
}


=head2 get_subtree

 Arg1: Treefam::Node or Treefam::Gene object
 Arg2: (optional) Treefam::Node or Treefam::Gene object
 Description: Gets a subtree rooted at the given node if
              only one node is given. If 2 nodes are given
              as arguments then the subtree is rooted at the
              last common ancestor of the given nodes.
 Returntype: Treefam::Tree

=cut

sub get_subtree {

  my $self = shift;
  my @nodes = @_;
  my $dbc = $self->{'DBConnection'};
  my $nhx;
  if (scalar(@nodes)==2) {
    my $node = $self->get_last_common_ancestor(@nodes);
    $nhx = $self->_string_nhx($node);
  }
  elsif (scalar(@nodes)==1) {
    $nhx = $self->_string_nhx($nodes[0]);
  }
  else {
    croak "ERROR: One or two nodes required";
  }

  return new Treefam::Tree($dbc,$self->ID,$self->type,$nhx);

}

=head2 get_species

 Arg: (optional) type of species name: latin or swcode
      (5 letters Swissprot code)
 Description: Gets a list of species that have at least one gene
              in the tree, defaults to latin name of species
 Returntype: list of strings

=cut

sub get_species {

  my $self = shift;
  my $type = shift if @_;
  my @species;
  # Go gene by gene as the tree might be different from the database one
  my %seen;
  foreach my $gene($self->get_genes) {
    push @species,$gene->species($type) unless $seen{$gene->species($type)}++;
  }
  return @species;

}

=head2 family

 Arg: (optional) Treefam::Family object or family AC
 Description: Gets/sets family the tree belongs to
 Returntype: Treefam::Family object

=cut

sub family {

  my $self = shift;
  my $family = shift if @_;
  my $dbc = $self->{'DBConnection'};
  my $famh = $dbc->get_FamilyHandle();
  my $familyID;
  if (!$family) {
    $familyID = $self->ID;
  }
  else {
    if (ref($family)) {
      return $family;
    }
    else {
      $familyID = $family;
    }
  }
  return $famh->get_by_id($familyID) if $familyID;

}

=head2 get_species_tree

 Arg: (optional) complete or short
 Description: Gets species tree for the species in this Treefam
              tree. If arg is set to 'complete', the full tree
              up to the species root and with all intermediate
              nodes is returned. If arg is set to 'short',
              returns a tree with only the leaves and the
              internal nodes that have more than one child.
              This is the default behaviour.
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
  my $query = qq(SELECT taxon_id FROM spec_names WHERE name= ?);
  my $sth1 = $dbh->prepare($query);
  $query = qq(SELECT t2.name FROM spec_nodes t1 LEFT JOIN spec_nodes t2 ON t1.parent_id=t2.taxon_id WHERE t1.taxon_id =?);
  my $sth2 = $dbh->prepare($query);
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
    $sth1->execute($name);
    my ($taxid) = $sth1->fetchrow_array();
    $sth2->execute($taxid);
    my ($n) = $sth2->fetchrow_array();
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

  if (!defined($type)|| (lc($type) ne 'complete' && lc($type) ne 'full')) {
    foreach my $node (@{$species_tree->{_nodes}}) {
      # delete node if it has only 1 child,
      # child becomes child of grand-parent
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

=head2 get_alignment

 Arg: (optional) type of sequence: aa (for amino-acid, default)
      or nt (for nucleotide)
 Description: Gets multiple alignment used to build the tree.
              Defaults to protein alignment. Nucleotide alignment
              is guided by protein alignment.
 Returntype: string, multiple alignment in FASTA format

=cut

sub get_alignment {

  my $self = shift;
  my $type = shift if @_;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->get_DatabaseHandle;
  my $tree_aln_aa = '';
  my $tree_aln_nt = '';
  my $family = $self->family;
  my $AC = $family->AC if $family;
  if (!$AC) {
    warn "WARNING: Tree not in a family. No alignment available.\n";
    return;
  }
  if (!defined($self->{'aa_alignment'})) {
    my $query = qq( SELECT al.ID, al.CIGAR, s.SEQ
                    FROM genes g, aa_full_align al, aa_seq s
                    WHERE g.GID = ?
                    AND g.ID = al.ID
                    AND al.ID = s.ID
                    AND al.AC = ?);
    my $sth = $dbh->prepare($query);
    foreach my $gene($self->get_genes) {
      my $geneID = $gene->ID;
      my $seqID = $gene->sequence_id;
      $sth->execute($geneID,$AC);
      my ($ID,$cigar,$seq) = $sth->fetchrow_array;
      my $aln = _cigar2align($cigar,$seq,0);
      $tree_aln_aa .= ">$ID GENEID=$geneID\n";
      $tree_aln_aa .=_pretty_seq($aln)."\n";
      # convert to protein-guided nucleotide alignment
      $query = qq(SELECT s.SEQ FROM nt_seq s WHERE s.ID= ?);
      my $sth = $dbh->prepare($query);
      $sth->execute($seqID);
      my ($ntseq) = $sth->fetchrow_array;
      next unless $ntseq;
      $aln = _cigar2align($cigar,$ntseq,1);
      $tree_aln_nt .= ">$ID GENEID=$geneID\n";
      $tree_aln_nt .=_pretty_seq($aln)."\n";
    }
    $self->{'aa_alignment'} = $tree_aln_aa;
    $self->{'nt_alignment'} = $tree_aln_nt;
  }
  if (lc($type) eq 'nt' || lc($type) eq 'nucleotide') {
    return $self->{'nt_alignment'};
  }
  else {
    return $self->{'aa_alignment'};
  }
}

=head2 _cigar2align

  Arg1: string, CIGAR format
  Arg2: string, sequence
  Arg3: (optional) 0 or 1, use 1 for nucleotide sequence
  Description : Converts CIGAR string to alignment string.
  ReturnType  : string

=cut

sub _cigar2align {
  # code from Li Heng
  my $cigar = shift;
  my $seq = shift;
  my $is_nucl = (@_)? shift : 0;
  my $tmp = $cigar;
  my $start = 0;
  my $len = length($seq);
  if ($is_nucl) {
    $tmp =~ s/(\d+)D/'-'x($1*3)/eg;
    $tmp =~ s/(\d+)M/$start+=$1*3,($start<=$len)?substr($seq,$start-$1*3,$1*3):'-'x($1*3)/eg;
  }
  else {
    $tmp =~ s/(\d+)D/'-'x$1/eg;
    $tmp =~ s/(\d+)M/$start+=$1,($start<=$len)?substr($seq,$start-$1,$1):'-'x$1/eg;
  }
  return $tmp;

}

=head2 _pretty_seq

  Arg1: string, sequence
  Arg2: (optional) int, desired line length, defaults to 60
  Description : reformats string into lines of given number
                of characters
  ReturnType  : string

=cut

sub _pretty_seq {

  my $seq = shift;
  my $line_size = (@_) ? shift : 60;
  my $seqlength = length($seq);
  my $i = $line_size;
  while ($i < $seqlength) {
    substr($seq, $i, 0) = "\n";
    $seqlength++;
    $i += $line_size+1;
  }
  return $seq;
}

1;
