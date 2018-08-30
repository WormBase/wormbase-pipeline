
package AGR;

use strict;
use Exporter;
use JSON;

our @ISA = qw(Exporter);
our @EXPORT = qw(get_bgi_genes get_rfc_date get_file_metadata_json);


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


##############################################
sub get_rfc_date {

  my $date;
  
  open(my $date_fh, "date --rfc-3339=seconds |");
  while(<$date_fh>) {
    if (/^(\S+)\s+(\S+)/) {
      $date = "${1}T${2}";
    }
  }
  
  return $date;
}


sub get_file_metadata_json {
  my ($wb_rel, $date) = @_;

  $date = &get_rfc_date() if not defined $date;

  my $meta_data = {
    dateProduced => $date,
    dataProvider => { 
      crossReference => {
        id => "WB",
        pages => ["homepage"],
      },
      type => "curated",
    },
    release      => $wb_rel,
  };

  return $meta_data;
}


1;
