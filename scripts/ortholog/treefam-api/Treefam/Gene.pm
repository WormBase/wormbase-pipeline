# Author: jkh1
# 2006-08-21
#

=head1 NAME

 Treefam::Gene

=head1 SYNOPSIS

 use Treefam::DBConnection;

 my $dbc = new Treefam::DBConnection ();

 my $gene_handle = $dbc->get_GeneHandle();

 my $gene = $gene_handle->get_by_id('ENSG00000087586');

 my $symbol= $gene->symbol();
 my @names = $gene->names();
 my $family = $gene->family();
 my $hmmer = $gene->hmmer_score();
 my @orthologs = $gene->get_orthologs('Saccharomyces cerevisiae',
                                      'Schizosaccharomyces pombe',
                                      'Drosophila melanogaster',
                                      'Caenorhabditis elegans');


=head1 DESCRIPTION

 Representation of a Treefam gene.

=head1 CONTACT

 jkh1@sanger.ac.uk

=cut


package Treefam::Gene;


use strict;
use Carp;
use Treefam::DBConnection;
use Scalar::Util qw(weaken);

=head2 new

 Arg1: Treefam::DBConnection
 Arg2: optional, gene ID
 Description: Creates a new gene object.
 Returntype: Treefam::Gene

=cut

sub new {

  my ($class,$dbc,$geneID) = @_;
  my $self = {};
  $self->{'DBConnection'} = $dbc;
  weaken($self->{'DBConnection'});
  my $dbh = $dbc->{'database_handle'};

  $self->{'ID'} = $geneID;

  bless ($self, $class);

  return $self;

}

=head2 ID

 Arg: optional, Treefam gene ID
 Description: Gets/sets Treefam ID
 Returntype: string

=cut

sub ID {

  my $self = shift;
  $self->{'ID'} = shift if @_;
  return $self->{'ID'};

}

=head2 symbol

 Arg: optional, gene symbol
 Description: Gets/sets gene's symbol
 Returntype: string

=cut

sub symbol {

  my $self = shift;
  $self->{'symbol'} = shift if @_;
  if (!defined $self->{'symbol'}) {
    my $geneID = $self->{'ID'};
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->{'database_handle'};
    my $sth = $dbh->prepare("SELECT DISTINCT g.symbol FROM genes g WHERE g.GID = ?");
    $sth->execute($geneID);
    ($self->{'symbol'}) = $sth->fetchrow_array();
    $sth->finish();
  }

  return $self->{'symbol'};

}

=head2 sequence_id

 Arg: optional, id of the gene's representative sequence
 Description: Gets/sets ID of the gene's representative
              sequence used in the trees or ID for longest
              protein if gene not in a family.
 Returntype: string

=cut

sub sequence_id {

  my $self = shift;
  $self->{'sequence_id'} = shift if @_;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  my $family = $self->family;
  if ($family) {
    if (!defined $self->{'sequence_id'}) {
      my $AC = $family->AC;
      my $geneID = $self->{'ID'};
      my $sth = $dbh->prepare("SELECT DISTINCT g.ID
                               FROM genes g, aa_full_align al
                               WHERE g.GID = ?
                               AND g.ID = al.ID
                               AND al.AC = ?");
      $sth->execute($geneID,$AC);
      ($self->{'sequence_id'}) = $sth->fetchrow_array();
      # in some cases, sequence is only in seed (because of curation in older release)
      if (!$self->{'sequence_id'}) {
	$sth = $dbh->prepare("SELECT DISTINCT g.ID
                              FROM genes g, aa_seed_align al
                              WHERE g.GID = ?
                              AND g.ID = al.ID
                              AND al.AC = ?");
	$sth->execute($geneID,$AC);
	($self->{'sequence_id'}) = $sth->fetchrow_array();
      }
      $sth->finish();
    }
  }
  else {
    if (!$self->{'sequence_id'}) {
      # for orphan genes, use longest protein
      my $query = qq(SELECT s.ID
                     FROM aa_seq s LEFT JOIN genes g
                     ON s.ID=g.ID
                     WHERE g.GID= ?
                     ORDER BY LENGTH DESC
                     LIMIT 1);
      my $sth = $dbh->prepare($query);
      $sth->execute($self->ID);
      ($self->{'sequence_id'}) = $sth->fetchrow_array();
    }
  }

  return $self->{'sequence_id'};

}

=head2 sequence

 Arg1: optional, string, the type of sequence (nt: nucleotide,
       aa: amino-acid), defaults to amino-acid.
 Arg2: optional, string, the gene's representative sequence
 Description: Gets/sets the gene's representative sequence of
              given type used in the trees or sequence for
              longest protein if gene not in a family.
 Returntype: string

=cut

sub sequence {

  my $self = shift;
  my $type = shift if @_;
  $type ||='aa';
  $self->{'sequence'}{$type} = shift if @_;
  if (!defined $self->{'sequence'}{$type}) {
    my $seqID = $self->sequence_id;
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->{'database_handle'};
    my $query = "";
    if (lc($type) eq 'nt' || lc($type) eq 'nucleotide') {
      $query = qq(SELECT DISTINCT s.SEQ
                  FROM nt_seq s
                  WHERE s.ID = ?);
    }
    else {
      $query = qq(SELECT DISTINCT s.SEQ
                  FROM aa_seq s
                  WHERE s.ID = ?);
    }
    my $sth = $dbh->prepare($query);
    $sth->execute($seqID);
    ($self->{'sequence'}{$type}) = $sth->fetchrow_array();
    $sth->finish();
  }

  return $self->{'sequence'}{$type};

}

=head2 names

 Arg: optional, names as one string
 Description: Gets/sets the gene's names
 Returntype: string

=cut

sub names {

  my ($self,$names) = @_;
  $self->{'names'} = $names if $names;
  if (!defined $self->{'names'}) {
    my $geneID = $self->{'ID'};
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->{'database_handle'};
    my $sth = $dbh->prepare("SELECT DISTINCT g.desc
                             FROM genes g
                             WHERE g.gid= ?");
    $sth->execute($geneID);
    ($self->{'names'})=$sth->fetchrow_array();
    $sth->finish();
  }
  return $self->{'names'};

}

=head2 transcripts

 Arg: optional, list of transcript IDs
 Description: Gets/sets transcripts associated with gene
 Returntype: list of transcript IDs

=cut

sub transcripts {

  my ($self,@transcripts) = @_;
  push (@{$self->{'transcripts'}},@transcripts) if @transcripts;
  if (!defined @{$self->{'transcripts'}}) {
    my $geneID = $self->{'ID'};
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->{'database_handle'};
    my $sth = $dbh->prepare("SELECT DISTINCT TID
                             FROM genes
                             WHERE gid= ?");
    $sth->execute($geneID);
    while (my ($transcript)=$sth->fetchrow_array()) {
      push (@{$self->{'transcripts'}},$transcript);
    }
  }
  return @{$self->{'transcripts'}};

}

=head2 species

 Arg1: optional, type of name: latin or swcode (5 letters species
       name used in Swissprot)
 Arg2: optional, species name
 Description: Gets/sets gene's species. Returns latin name by
              default.
 Returntype: string

=cut

sub species {

  my $self = shift;
  my $type = shift if @_;
  $type ||='latin';
  $self->{'species'} = shift if @_;
  if (!defined $self->{'species'} || (defined $self->{'species_name_type'} && defined $type && lc($type) ne lc($self->{'species_name_type'}))) {  # first request or request for a different name type e.g. latin name when we previously had swcode
    my $geneID = $self->{'ID'};
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->{'database_handle'};
    my $sth;
    if (defined($type) && lc($type) eq 'swcode') {
      $sth = $dbh->prepare("SELECT DISTINCT s.swcode
                            FROM species s, genes g
                            WHERE g.GID = ?
                            AND g.tax_id=s.tax_id");
    }
    else {
      $sth = $dbh->prepare("SELECT DISTINCT s.taxname
                            FROM species s, genes g
                            WHERE g.GID = ?
                            AND g.tax_id=s.tax_id");
    }
    $sth->execute($geneID);
    ($self->{'species'}) = $sth->fetchrow_array();
    $sth->finish();
  }
  $self->{'species_name_type'} = defined($type) && lc($type) eq 'swcode' ? 'SWCODE' : 'LATIN';
  return $self->{'species'};

}

=head2 family

 Arg: optional, Treefam::Family object
 Description: Gets/sets family the gene belongs to
              (family with the highest HMMER score)
 Returntype: Treefam::Family object

=cut

sub family {

  my $self = shift;
  $self->{'family'} = shift if @_;
  if (!defined $self->{'family'}) {
    my $dbc = $self->{'DBConnection'};
    my $famh = $dbc->get_FamilyHandle();
    $self->{'family'} = $famh->get_by_gene($self);
  }

  return $self->{'family'};

}

=head2 hmmer_score

 Arg1: optional, Treefam::Family object or AC
 Arg2: optional, hmmer score
 Description: Gets/sets hmmer score for the gene in the given family,
              defaults to the family the gene belongs to
 Returntype: double

=cut

sub hmmer_score {

  my $self = shift;
  my $family = shift if @_;
  my $score = shift if @_;
  if (!$family) {
    $family = $self->family();
    return undef unless $family;
  }
  my $familyID = ref($family)? $family->ID : $family;
  if ($familyID && $score) {
    $self->{'hmmer_score'}{$familyID} = $score;
  }
  if ($familyID && !defined($self->{'hmmer_score'}{$familyID})) {
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->{'database_handle'};
    my $geneID = $self->ID;
    my $sth = $dbh->prepare("SELECT h.score
                             FROM hmmer h, genes g
                             WHERE g.gid = ?
                             AND g.id = h.id
                             AND h.ac = ?");
    $sth->execute($geneID,$familyID);
    ($self->{'hmmer_score'}{$familyID}) = $sth->fetchrow_array();
    $sth->finish();
  }

  return $self->{'hmmer_score'}{$familyID} if $familyID;

}

=head2 hmmer_evalue

 Arg1: optional, Treefam::Family object or AC
 Arg2: optional, hmmer evalue
 Description: Gets/sets hmmer evalue for the gene in the given
              family, defaults to the family the gene belongs to
 Returntype: double

=cut

sub hmmer_evalue {

  my $self = shift;
  my $family = shift if @_;
  my $evalue = shift if @_;
  if (!$family) {
    $family = $self->family();
    return undef unless $family;
  }
  my $familyID = ref($family)? $family->ID : $family;
  if ($familyID && $evalue) {
    $self->{'hmmer_evalue'}{$familyID} = $evalue;
  }
  if ($familyID && !defined $self->{'hmmer_evalue'}) {
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->{'database_handle'};
    my $geneID = $self->ID;
    my $sth = $dbh->prepare("SELECT h.evalue
                             FROM hmmer h, genes g
                             WHERE g.gid = ?
                             AND g.id = h.id
                             AND h.ac = ?");
    $sth->execute($geneID,$familyID);
    ($self->{'hmmer_evalue'}{$familyID}) = $sth->fetchrow_array();
    $sth->finish();
  }

  return $self->{'hmmer_evalue'}{$familyID};

}

=head2 get_orthologs

 Arg: optional, list of species (by latin name) from which to
      get orthologs
 Description: Gets this gene's orthologs
 Returntype: list of Treefam::Gene objects with boostrap attribute set

=cut

sub get_orthologs {

  my $self = shift;
  my @species = @_ if @_;
  if (!defined $self->{'orthologs'}) {
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->{'database_handle'};
    my $geneID = $self->ID;
    my $spec = $self->species();
    my $genequery = qq(SELECT IDX FROM genes WHERE GID= ?);
    my $sth0 = $dbh->prepare($genequery);
    # get all orthologs/paralogs
    my $orthoquery1 = qq (SELECT DISTINCT g.GID, sp.taxname, o.bootstrap
                          FROM genes g, ortholog o, species sp
                          WHERE o.idx1 = ?
                          AND o.idx2 = g.idx
                          AND g.TAX_ID = sp.TAX_ID );
    my $orthoquery2 = qq (SELECT DISTINCT g.GID, sp.taxname, o.bootstrap
                          FROM genes g, ortholog o, species sp
                          WHERE o.idx2 = ? AND o.idx1 = g.idx
                          AND g.TAX_ID = sp.TAX_ID );
    my $sth1 = $dbh->prepare($orthoquery1);
    my $sth2 = $dbh->prepare($orthoquery2);

    $sth0->execute($geneID);
    while (my ($idx) = $sth0->fetchrow_array()) {
      foreach my $tmpsth($sth1,$sth2) {
	$tmpsth->execute($idx);
	while (my($orthologID,$species,$bootstrap) = $tmpsth->fetchrow_array()) {
	  next if ($species eq $spec); # discard paralogs
	  my $gh = $dbc->get_GeneHandle();
	  my $ortholog = $gh->get_by_id($orthologID);
	  if ($ortholog) {
	    push @{$self->{'orthologs'}},$ortholog;
	    $ortholog->species('latin',$species);
	    $ortholog->bootstrap($bootstrap);
	  }
	}
      }
    }
  }
  if (defined($self->{'orthologs'})) {
    my %wanted;
    # reduce list to the requested species
    if (@species) {
      @wanted{@species} = (1) x @species;
      @{$self->{'orthologs'}} = grep { $wanted{$_->species} } @{$self->{'orthologs'}};
    }

    return @{$self->{'orthologs'}};
  }
  else {
    return ();
  }
}

=head2 get_paralogs

 Description: Gets this gene's paralogs. We currently only get
              within-species paralogs
 Returntype: list of Treefam::Gene objects with boostrap attribute set

=cut

sub get_paralogs {

  my $self = shift;
  if (!defined $self->{'paralogs'}) {
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->{'database_handle'};
    my $geneID = $self->ID;
    my $spec = $self->species();
    my $genequery = qq(SELECT IDX FROM genes WHERE GID= ?);
    my $sth0 = $dbh->prepare($genequery);
    # get all orthologs/paralogs
    my $orthoquery1 = qq (SELECT DISTINCT g.GID, sp.taxname, o.bootstrap
                          FROM genes g, ortholog o, species sp
                          WHERE o.idx1 = ?
                          AND o.idx2 = g.idx
                          AND g.TAX_ID = sp.TAX_ID );
    my $orthoquery2 = qq (SELECT DISTINCT g.GID, sp.taxname, o.bootstrap
                          FROM genes g, ortholog o, species sp
                          WHERE o.idx2 = ? AND o.idx1 = g.idx
                          AND g.TAX_ID = sp.TAX_ID );
    my $sth1 = $dbh->prepare($orthoquery1);
    my $sth2 = $dbh->prepare($orthoquery2);

    $sth0->execute($geneID);
    while(my ($idx) = $sth0->fetchrow_array()) {
      foreach my $tmpsth($sth1,$sth2) {
	$tmpsth->execute($idx);
	while (my($paralogID,$species,$bootstrap) = $tmpsth->fetchrow_array()) {
	  next if ($species ne $spec || $paralogID eq $geneID); # discard orthologs and self
	  my $gh = $dbc->get_GeneHandle();
	  my $paralog = $gh->get_by_id($paralogID);
	  if ($paralog) {
	    push @{$self->{'paralogs'}},$paralog;
	    $paralog->bootstrap($bootstrap);
	  }
	}
      }
    }
  }
  if (defined($self->{'paralogs'})) {
    return @{$self->{'paralogs'}};
  }
  else {
    return ();
  }

}

=head2 get_domains

 Arg: (optional) e-value cut-off
 Description: Gets this gene's protein domains with e-value
              below given cut-off, default is 1e-2. Returned
              list is non-redundant.
              See also get_domain_sequence
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
    my $geneID = $self->sequence_id;
    my %seen;
    my $query = qq(SELECT DISTINCT PFAM_ID
                   FROM pfam
                   WHERE ID= ?
                   AND EVALUE<$cutoff
                   ORDER BY SEQ_START,EVALUE);
    my $sth = $dbh->prepare($query);
    $sth->execute($geneID);
    while (my ($pfamid) = $sth->fetchrow_array()) {
      push @{$self->{'domains'}},$pfamid if ($pfamid && !$seen{$pfamid}++);
    }
  }
  # remove undefined values
  @{$self->{'domains'}} = grep { $_ } @{$self->{'domains'}};

  return @{$self->{'domains'}} if (defined($self->{'domains'}));

}

=head2 get_domain_sequence

 Arg: (optional) e-value cut-off
 Description: Gets this gene's sequence of protein domains
              with e-value below given cut-off (default is
              1e-2).
 Returntype: list of strings (PFAM domain IDs)

=cut

sub get_domain_sequence {

  my $self = shift;
  my $cutoff = shift if @_;
  if(!$cutoff) {
    $cutoff = 1e-2;
  }
  if (!$self->{'domains'}) {
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->{'database_handle'};
    my $geneID = $self->sequence_id;
    my %seen;
    my $query = qq(SELECT DISTINCT PFAM_ID,SEQ_START
                   FROM pfam
                   WHERE ID= ?
                   AND EVALUE<$cutoff
                   ORDER BY SEQ_START,EVALUE);
    my $sth = $dbh->prepare($query);
    $sth->execute($geneID);
    while (my ($pfamid,$seq_start) = $sth->fetchrow_array()) {
      push @{$self->{'domains'}},$pfamid;
    }
  }
  # remove undefined values
  @{$self->{'domains'}} = grep { $_ } @{$self->{'domains'}};

  return @{$self->{'domains'}} if (defined($self->{'domains'}));

}

=head2 position

 Arg: optional, integer position of gene (coding start) on
      chromosome
 Description: Gets/sets position of gene (coding start) on
              chromosome. Returns 5'-most position in case
              of multiple starts.
 Returntype: integer

=cut

sub position {

  my ($self,$pos) = @_;
  $self->{'position'} = $pos if $pos;
  if (!defined $self->{'position'}) {
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->{'database_handle'};
    my @transcripts = $self->transcripts;
    my $position;
    my $sth = $dbh->prepare("SELECT m.C_START, m.C_STOP, m.STRAND
                             FROM map m, genes g
                             WHERE g.TID = ?
                             AND g.ID = m.ID");
    foreach my $transcript(@transcripts) {
      $sth->execute($transcript);
      while ((my @array) = $sth->fetchrow_array) {
	# Note that if a transcript is on the - strand, eg. if it is from
	# 1000-2000 on the - strand, then start=1000 stop=2000.
	my $start = $array[0] + 1; # The coordinates are counted from 0.
	my $end = $array[1] + 1;
	my $strand = $array[2];
	if ($strand eq '-') {
	  ($start,$end) = ($end,$start);
	}
	if (!(defined($position))) {
	  $position =$start;
	}
	else {
	  # If the gene is on the - strand, we want the coding start,
	  # eg. 2000 for a gene on the - strand from 1000-2000.
	  if    ($strand eq '+') {
	    if ($start < $position) {
	      $position = $start;
	    }
	  }
	  elsif ($strand eq '-') {
	    if ($start > $position) {
	      $position = $start;
	    }
	  }
	}
      }
      $sth->finish();
      $self->{'position'} = $position;
    }
  }
  return $self->{'position'};
}


=head2 chromosome

   Arg: optional, chromosome as string
   Description: Gets/sets chromosome that gene is on
   Returntype: string

=cut

sub chromosome {

   my ($self,$chrom) = @_;
   $self->{'chromosome'} = $chrom if $chrom;
   if (!defined $self->{'chromosome'}) {
    my $seqID = $self->sequence_id;
    my $dbc = $self->{'DBConnection'};
    my $dbh = $dbc->{'database_handle'};
    my $sth = $dbh->prepare("SELECT m.TARGET
                             FROM map m WHERE m.ID = ?");
    $sth->execute($seqID);
    ($self->{'chromosome'})=$sth->fetchrow_array();
    $sth->finish();
  }
  return $self->{'chromosome'};
}

=head2 AUTOLOAD

   Arg: any type
   Description: Useful for setting then getting a
                new attribute.
                e.g. $gene->color('red');
   Returntype: as Arg

=cut

sub AUTOLOAD {

  my $self = shift;
  my $attribute = our $AUTOLOAD;
  $attribute =~s/.*:://;
  $self->{$attribute} = shift if @_;
  return $self->{$attribute};
}


1;
