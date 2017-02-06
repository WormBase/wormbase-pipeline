#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;

use lib $ENV{CVS_DIR};
use Wormbase;
#use Log_files;

my %XREF_MAP = (
  "NCBI"       => "NCBIGene",
  "SwissProt"  => "UniProtKB",
  "TrEMBL"     => "UniProtKB",
  "RNAcentral" => "RNAcentral",
);



my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($outfile, $acedbpath, $ws_version, $out_fh);

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "test"        => \$test,
	    "verbose"     => \$verbose,
	    "store:s"     => \$store,
	    "database:s"  => \$acedbpath,
	    "outfile:s"   => \$outfile,
            "wsversion=s" => \$ws_version,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
      );
}

my $tace = $wormbase->tace;
#my $log  = Log_files->make_build_log($wormbase);

my $taxid = $wormbase->ncbi_tax_id;
my $full_name = $wormbase->full_name;

$acedbpath = $wormbase->autoace unless $acedbpath;

if (defined $outfile) {
  open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";
} else {
  $out_fh = \*STDOUT;
}

my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or die("Connection failure: ". Ace->error);

my $base_query = 'FIND Gene WHERE Sequence AND Live AND Species = "Caenorhabditis elegans"';

my $meta_data = {
  dateProduced => &get_rfc_date(),
  dataProvider => "WB", 
  release      => (defined $ws_version) ? $ws_version : $wormbase->get_wormbase_version_name(),
};


my @genes;

foreach my $sub_query (
  " AND Corresponding_CDS",
  " AND Corresponding_pseudogene",
  " AND Corresponding_transcript AND NOT Corresponding_CDS",
    ) {

  my $it = $db->fetch_many(-query=> $base_query . $sub_query);
  my $i = 0;
  
  while (my $obj=$it->next) {
    $i++;

    print STDERR $obj->name, "\n";
    
    next unless $obj->isObject();
    next unless $obj->Species;
    next unless $obj->Species->name eq $full_name;
    
    my $biotype = $obj->Biotype->name;
    
    my ($symbol, @synonyms);
    my $seq_name = $obj->Sequence_name->name;
    if ($obj->CGC_name) {
      $symbol = $obj->CGC_name->name;
      push @synonyms, $seq_name;
    } else {
      $symbol = $seq_name;
    }
    push @synonyms, map { $_->name } $obj->Other_name;
    
    
    my ($concise_desc, $auto_desc);
    if ($obj->Concise_description) {
      $concise_desc = $obj->Concise_description->name;
    }
    if ($obj->Automated_description) {
      $auto_desc = $obj->Automated_description->name; 
    }
    
    my @xrefs;
    foreach my $dblink ($obj->Database) {
      if (exists $XREF_MAP{$dblink}) {
        push @xrefs, {
          dataProvider => $XREF_MAP{$dblink},
          id =>  $dblink->right->right->name,
        };
      }
    }
    push @xrefs, {
      dataProvider => "Ensembl",
      id => $obj->name,
    };
    
    my ($gene_class_name, $brief_id_name);
    if ($obj->Gene_class and $obj->Gene_class->Description) {
      $gene_class_name = $obj->Gene_class->Description->name;
    }
    if ($obj->Corresponding_CDS and $obj->Corresponding_CDS->Brief_identification) {
      $brief_id_name = $obj->Corresponding_CDS->Brief_identification->name;
    }
    my $name;
    if ($gene_class_name) {
      my ($suff) = $symbol =~ /-(\S+)/; 
      $name = "$gene_class_name $suff";
    } elsif ($brief_id_name) {
      $name = $brief_id_name;
    } else {
      $name = undef;
    }
    
    my @secondary_ids;
    if ($obj->Acquires_merge) {
      foreach my $g ($obj->Acquires_merge) {
        push @secondary_ids, $g->name;
      }
    }
    
    my @pmids;
    foreach my $paper ($obj->Reference) {
      foreach my $pref ($paper->Database) {
        if ($pref->name eq 'MEDLINE') {
          push @pmids, "PMID:" . $pref->right->right->name;
        }
      }
    }

    my $json_gene = {
      primaryId          => $obj->name,
      symbol             => $symbol,
      name               => $name,
      taxonId            => $taxid,
      geneSynopsis       => defined($concise_desc) ? $concise_desc : $auto_desc,
      soTermId           => $biotype,
      synonyms           => \@synonyms,
      secondaryIds       => \@secondary_ids,
      crossReferences    => \@xrefs,
      geneLiteratureUrl  => "http://www.wormbase.org/species/c_elegans/gene/" . $obj->name ."-e-10",
      #references         => \@pmids,
    };
    
    #print "$obj $seq_name $cgc_name @other_names\n";
    
    push @genes, $json_gene;
  }
}

my $data = {
  metaData => $meta_data,
  data => \@genes,
};

my $json_obj = JSON->new;
my $string = $json_obj->allow_nonref->canonical->pretty->encode($data);
print $out_fh $string;

$db->close;

exit(0);


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
