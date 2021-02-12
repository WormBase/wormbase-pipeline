package GAF;

use strict;
use DateTime;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(print_wormbase_GAF_line print_wormbase_GAF_header get_GAF_date get_gene_info make_species_files);

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
sub print_wormbase_GAF_header {
    my ($fh, $release, $type, $version) = @_;
    $version||='2.0';

    my $date = DateTime->now;

    print $fh "\!gaf-version: $version\n" unless $type eq 'RNAi';
    print $fh "\!generated-by: WormBase\n";
    print $fh "\!date-generated: ".$date->ymd."\n";
    print $fh "\!project-URL: https://wormbase.org\n";
    if ($type eq 'GO') {
	print $fh "\!specification-URL: https://github.com/geneontology/geneontology.github.io/blob/issue-go-annotation-2917-gaf-2_2-doc/_docs/go-annotation-file-gaf-format-22.md\n";
    }
    else {
	print $fh "\!specification-URL: https://wiki.wormbase.org/index.php/WormBase_gene_association_file\n";
    }
    print $fh "\!project-release: $release\n";
    print $fh "\!Contact Email: help\@wormbase.org\n";
    
}

####################################
sub get_GAF_date {
  my ($in_date) = @_;

  my ($year, $month, $day);

  if (defined $in_date) {
    ($year, $month, $day) = split(/\-/, $in_date);
  } else {
    $year  = (1900 + (localtime)[5]);
    $month = (1 + (localtime)[4]);
    $day   = (localtime)[3];
  }

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

    my ($wbgene, $public_name, $sequence_name, $other_name, $status, $species) = split(/\t/, $_);
    next if $wbgene !~ /^WBGene/;

    if (not exists $data{$wbgene}) {
      $data{$wbgene} = {
        public_name   => $public_name,
        sequence_name => ($sequence_name) ? $sequence_name : "",
        other_name    => ($other_name) ? [$other_name] : [],
        status        => $status,
        species       => $species,
      };      
    } elsif ($other_name) {
      push @{$data{$wbgene}->{other_name}}, $other_name;
    }
  }

  unlink $def;
  return \%data;
}


###################################
sub make_species_files {
    my ($wormbase, $file, $daf) = @_;

    my %species_lines;
    my @headers;
    open (COLLATED, '<', $file) or die "Could not open $file for writing\n";
    while (<COLLATED>) {
	chomp;
	if ($_ =~ /^!/) {
	    push @headers, $_;
	    next;
	}
	my $taxon_id;
	if ($daf) {
	    ($taxon_id) = $_ =~ /^(\d+)\s/;
	}
	else {
	    ($taxon_id) = $_ =~ /taxon:(\d+)\s/;
	}
	push @{$species_lines{$taxon_id}}, $_;
    }
    close (COLLATED);

    my %accessors = ($wormbase->species_accessors);
    $accessors{$wormbase->species} = $wormbase;


    for my $wb(values %accessors) {
	my $species_file = $file . '.' . $wb->full_name('-g_species' => 1);
	open (SPECIES, '>', $species_file) or die "Could not open $species_file for writing\n";
	print SPECIES join("\n", @headers) . "\n";
	print SPECIES join("\n", @{$species_lines{$wb->ncbi_tax_id}}) . "\n" if exists $species_lines{$wb->ncbi_tax_id};;
	close (SPECIES);
    }

    return;
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
Condition WBGene* 

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
Optional 
Visible 
Class 
Class Gene_name 
From 1 
Tag Sequence_name 
 
Colonne 4
Width 12
Optional
Visible
Class
Class Gene_name
From 1
Tag Other_name

Colonne 5 
Width 12 
Mandatory 
Visible 
Next_Tag 
From 1 
Tag Status 
 
Colonne 6 
Width 12 
Mandatory 
Visible 
Class 
Class Species 
From 1 
Tag Species 
$species_condition


EOF

  print $def_fh $def;
  close($def_fh);
  return $tm_def_file;
}
