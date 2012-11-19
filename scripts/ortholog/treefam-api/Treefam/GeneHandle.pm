# Author: jkh1
# 2006-08-21
#

=head1 NAME

 Treefam::GeneHandle

=head1 SYNOPSIS

 use Treefam::DBConnection;

 my $dbc = new Treefam::DBConnection ();

 my $gene_handle = $dbc->get_GeneHandle();

 my $gene = $gene_handle->get_by_id('ENSG00000087586');


=head1 DESCRIPTION

 Enables retrieval of Treefam::Gene objects.

=head1 CONTACT

 jkh1@sanger.ac.uk


=cut


package Treefam::GeneHandle;


use strict;
use Carp;
use Treefam::Gene;
use Treefam::DBConnection;
use Scalar::Util qw(weaken);


=head2 new

 Arg: Treefam::DBConnection
 Description: Creates a new gene object handle
 Returntype: Treefam::GeneHandle

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

 Arg: gene ID
 Description: Gets gene with given ID
 Returntype: Treefam::Gene object

=cut

sub get_by_id {

  my ($self,$geneID) = @_;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  # check if gene exists in database
  my $query = qq( SELECT COUNT(*) FROM genes WHERE GID= ? );
  my $sth= $dbh->prepare ($query);
  $sth->execute($geneID);
  my ($count) = $sth->fetchrow_array();
  unless ($count) {
    return undef;
  }

  return new Treefam::Gene($dbc,$geneID);

}

=head2 get_by_sequence_id

 Arg: sequence/transcript ID
 Description: Gets gene with given ID
 Returntype: Treefam::Gene object

=cut

sub get_by_sequence_id {

  my ($self,$ID) = @_;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  my $query = qq( SELECT GID FROM genes WHERE ID= ? OR TID= ?);
  my $sth= $dbh->prepare ($query);
  $sth->execute($ID,$ID);
  my ($geneID) = $sth->fetchrow_array();

  return undef if (!defined($geneID));

  my $gene = Treefam::Gene->new($dbc,$geneID);
  $gene->sequence_id($ID);

  return $gene;

}

=head2 get_by_symbol

 Arg: gene symbol
 Description: Gets gene with given symbol
 Returntype: Treefam::Gene object

=cut

sub get_by_symbol {

  my ($self,$symbol) = @_;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  my $sth = $dbh->prepare("SELECT DISTINCT GID FROM genes
                           WHERE symbol=?");
  $sth->execute($symbol);
  my ($geneID) = $sth->fetchrow_array();
  $sth->finish();

  return undef if ( !defined $geneID );

  return $self->get_by_id($geneID);

}

=head2 get_all_by_name

 Arg: string
 Description: Gets genes with names matching given string
 Returntype: list of Treefam::Gene objects

=cut

sub get_all_by_name {

  my ($self,$name) = @_;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  my $sth = $dbh->prepare("SELECT DISTINCT g.GID FROM genes g WHERE g.desc LIKE '%$name%' ");
  $sth->execute();
  my @genes=();
  while (my ($g)=$sth->fetchrow_array()) {
    push (@genes,$self->get_by_id($g));
  }
  $sth->finish();

  return @genes;

}

=head2 get_all_by_species

 Arg1: string, species name
 Arg2: optional, Treefam::Family object
 Description: Gets all genes of given species present in
              the family if one is given, otherwise returns
              all genes of the requested species found in
              Treefam.
 Returntype: list of Treefam::Gene objects

=cut

sub get_all_by_species {

  my ($self,$species,$family) = @_;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  my $sth;
  if ($family) {
    my $familyID = ref($family)? $family->ID : $family;
    my $type = $family->type eq 'familyA'? 'famA_gene' : 'famB_gene';
    $sth = $dbh->prepare("SELECT DISTINCT g.GID
                          FROM genes g, species s, $type t
                          WHERE g.tax_id = s.tax_id
                          AND (s.taxname = ? OR s.swcode = ?)
                          AND g.id = t.id
                          AND t.ac = ?");
    $sth->execute($species,$species,$familyID);
  }
  else {
    $sth = $dbh->prepare("SELECT DISTINCT g.GID
                          FROM genes g, species s
                          WHERE g.tax_id = s.tax_id
                          AND (s.taxname = ? OR s.swcode = ?)");
    $sth->execute($species,$species);
  }
  my @genes=();
  while (my ($g)=$sth->fetchrow_array()) {
    push (@genes,$self->get_by_id($g));
  }
  $sth->finish();

  return @genes;

}

=head2 get_all_by_family

 Arg: Treefam::Family
 Description: Gets genes belonging to the given family
 Returntype: list of Treefam::Gene objects

=cut

sub get_all_by_family {

  my ($self,$family) = @_;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  my $familyID = $family->ID;
  my $type = $family->type eq 'familyA'? 'famA_gene' : 'famB_gene';
  my $sth = $dbh->prepare("SELECT DISTINCT g.GID FROM genes g, $type t
                           WHERE t.AC= ? AND t.ID=g.ID");
  $sth->execute($familyID);
  my @genes=();
  while (my ($g)=$sth->fetchrow_array()) {
    push (@genes,$self->get_by_id($g));
  }
  $sth->finish();

  return @genes;

}

=head2 get_all_by_domain

 Arg1: string, pfam id of domain
 Arg2: (optional) e-value cut-off
 Description: Gets all genes with the given domain with
              e-value below given cut-off (default is 1e-2).
 Returntype: list of Treefam::Gene objects

=cut

sub get_all_by_domain {

  my $self = shift;
  my $pfamid = shift;
  unless ($pfamid) {
    croak "Pfam domain id required";
  }
  my $cutoff = shift if @_;
  if (!$cutoff) {
    $cutoff = 1e-2;
  }
  my @genes;
  my %seen;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  my $query = qq(SELECT DISTINCT g.GID
                 FROM genes g, pfam p
                 WHERE p.PFAM_ID = ?
                 AND p.EVALUE <= $cutoff
                 AND p.ID = g.ID);
  my $sth = $dbh->prepare($query);
  $sth->execute($pfamid);
  while (my ($ID)=$sth->fetchrow_array()) {
    push (@genes,$self->get_by_id($ID));
  }
  $sth->finish();

  return @genes;

}

=head2 get_all_by_domain_list

 Arg: list of strings, pfam ids of domains
 Description: Gets all genes with all the given domains.
              No e-value cut-off.
 Returntype: list of Treefam::Gene objects

=cut

sub get_all_by_domain_list {

  my $self = shift;
  my @pfamids = @_;
  unless ($pfamids[0]) {
    croak "Pfam domain ids required";
  }
  my @genes;
  my %seen;
  # remove duplicates
  my @uniq_pfamids = grep { !$seen{$_}++ } @pfamids;
  %seen = ();
  foreach my $pfamid(@uniq_pfamids) {
    map { $seen{$_->ID}++ if $_ } $self->get_all_by_domain($pfamid,10000);
  }
  my $n = scalar(@uniq_pfamids);
  my @IDs = grep { $seen{$_}==$n } keys %seen;
  while (my $ID=shift @IDs) {
    push (@genes,$self->get_by_id($ID));
  }

  return @genes;

}

1;
