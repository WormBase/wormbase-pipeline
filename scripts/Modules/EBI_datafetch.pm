
package EBI_datafetch;

use strict;
use LWP::UserAgent;
use Getopt::Long;

my $UNIPROT_BASE = "http://www.uniprot.org/uniprot/?";
my $ENA_BASE     = "http://www.ebi.ac.uk/ena/data/view/";

sub new {
  my ($class) = @_;

  my $self = {
    _uniprot_base => $UNIPROT_BASE,
    _ena_base     => $ENA_BASE,
  };
  bless $self, $class;

  return $self;
}


####################################
# doc: to do
###################################
sub get_uniprot_info_for_species {
  my ($self, $species) = @_;

  if (! defined $species) {die("SRS: Full species name not given\n")};
  $species =~ s/\s+/\%20/g;

  my %entries;

  #my $ua       = LWP::UserAgent->new;
  my $query    = "query=organism:${species}&format=tab&columns=entry%20name,id,length,protein%20names,genes,ec,database(WORMBASE)";
  #my $tmp_file = "/tmp/srs_results.$$.txt";

  my $chall = "${UNIPROT_BASE}${query}";
  print STDERR "$chall\n";

  #my $qa2 = $ua->get($chall, ':content_file' => $tmp_file);
  #if (not $qa2->is_success) {
  #  die($qa2->status_line) 
  #}
  
  #open(my $f, $tmp_file);

  open(my $f, "wget -O - '$chall' |") or die("Could not wget $chall\n");

  while(<$f>) {
    next if /^Entry/;
    chomp;

    my ($entry_name,
        $entry_acc,
        $len,
        $description,
        $gene_names,
        $ec_number, 
        $wb_cds_names) = split(/\t/, $_);
        
    my @gene_names = split(/\s+/, $gene_names);
    my @wb_cds = split(/\s*;\s*/, $wb_cds_names);
    
    print "ID $entry_name\ndesc $description\nGN $gene_names\nWB $wb_cds_names\nEC $ec_number\nLN $len\n";

    $entries{$entry_acc} = {
      id => $entry_name,
      description => $description,
      gene_names => \@gene_names,
      wormbase_names => \@wb_cds,
      ec_number => $ec_number,
      length => $len,
    };
  }
  #unlink $tmp_file;

  return \%entries;
}


####################################
# doc: to do
###################################
sub get_ena_sequence {
  my ($self, $acc) = @_;

  my $url = "${ENA_BASE}/${acc}&display=fasta";
  
  open(my $entry, "wget -O - 'http://www.ebi.ac.uk/ena/data/view/$acc&display=fasta' |")
      or die("Could not open EBI sequence fetch command for $acc\n");

  my ($seq, $ver, $desc);
  
  while (<$entry>) {
    chomp;
    
    /^\>ENA\|\S+\|\S+\.(\d+)\s+(.+)$/ and do { 
      $ver  = $1;
      $desc = $2;
      next;
    };

    /^(\S+)$/ and $seq .= uc($1);
  }
  close($entry) or die("Could not successfully close EBI sequence command for $acc\n");

  return ($seq, $ver, $desc);
}

##################################
=head2 uniprot_base_url

 Title   : uniprot_base_url
 Usage   :
 Function:
 Returns : 
 Args    :

=cut

sub uniprot_base_url{
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_uniprot_base} = $value;
  }
  return $self->{_uniprot_base};
}


##################################
=head2 ena_base_url

 Title   : uniprot_base_url
 Usage   :
 Function:
 Returns : 
 Args    :

=cut

sub ena_base_url {
  my ($self,$value) = @_;
  
  if (defined $value) {
    $self->{_ena_base} = $value;
  }
  return $self->{_ena_base};
}


1;
