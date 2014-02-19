
package GAF;

use strict;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(print_wormbase_GAF_line get_GAF_date get_gene_info);

###################################
sub print_wormbase_GAF_line {
  my ($fh, $gene, $gene_public_name, $qualifier, $ont_id, $refs, $evi_code, $with_from, $aspect, $gene_syn, $tax_id, $date) = @_;
            
  printf($fh "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
         "WB", 
         $gene,
         $gene_public_name,
         $qualifier,
         $ont_id,
         $refs,
         $evi_code,
         $with_from,
         $aspect,
         "",
         $gene_syn eq $gene_public_name ? "" : $gene_syn,
         "gene", 
         "taxon:$tax_id",
         $date,
         "WB",
         "",
         "");
}

####################################
sub get_GAF_date {
  
  my $year  = (1900 + (localtime)[5]);
  my $month = (1 + (localtime)[4]);
  my $day   = (localtime)[3];
  
  my $date = sprintf("%04d%02d%02d", $year, $month, $day);
  
  return $date;
}


####################################
sub get_gene_info {
  my ($db, $wormbase, $species_name) = @_;

  my $def = _write_name_status_tm_def($species_name);
  
  my %data;

  my $tmq = $wormbase->table_maker_query($db, $def);
  while(<$tmq>) {
    chomp;
    s/\"//g;

    my @l = split(/\s+/, $_);
    
    next if $l[0] !~ /^WBGene/;

    $data{$l[0]} = {
      public_name   => $l[1],
      status        => $l[2],
      taxid         => $l[3],
      sequence_name => ($l[4]) ? $l[4] : "",
    };
  }

  unlink $def;
  return \%data;
}



sub _write_name_status_tm_def {
  my ($species) = @_;

  my $tm_def_file = "/tmp/gene_name_stat.TM.$$.def";
  open(my $def_fh, ">$tm_def_file") or die "Could not open $tm_def_file\n";

  my $species_condition = ($species) ? "Condition \"$species\"" : "";

  my $def = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Gene 
From 1 
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class Gene_name 
From 1 
Tag Public_name 

Colonne 3 
Width 12 
Mandatory 
Visible 
Next_Tag 
From 1 
Tag Status 
 
Colonne 4 
Width 12 
Mandatory 
Hidden 
Class 
Class Species 
From 1 
Tag Species 
$species_condition
 
Colonne 5 
Width 12 
Mandatory 
Visible 
Integer 
From 4 
Tag NCBITaxonomyID 

Colonne 6 
Width 12 
Optional 
Visible 
Class 
Class Gene_name 
From 1 
Tag Sequence_name 
 


EOF

  print $def_fh $def;
  close($def_fh);
  return $tm_def_file;
}
