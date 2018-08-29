
package AGR;

use strict;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(get_bgi_genes);


############################
sub get_bgi_genes {
  my $bgi_json = shift;

  my (%bgi_genes, $json_string);
  open(my $json_fh, $bgi_json) or die "Could not open $bgi_json for reading\n";
  while(<$json_fh>) {
    $json_string .= $_;
  }
  
  my $json_reader = JSON->new;
  my $json = $json_reader->decode($json_string);
  
  foreach my $entry (@{$json->{data}}) {
    my $id = $entry->{primaryId};
    $bgi_genes{$id} = $entry;
  }

  return \%bgi_genes;

}

1;
