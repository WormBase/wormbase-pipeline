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
  "NCBI"       => "NCBI_Gene",
  "SwissProt"  => "UniProtKB",
  "TrEMBL"     => "UniProtKB",
  "RNAcentral" => "RNAcentral",
);



my ($help, $debug, $test, $verbose, $store, $wormbase, $species);
my ($outfile, $acedbpath, $ws_version, $gtf_file, $out_fh, $locs);

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "test"        => \$test,
	    "verbose"     => \$verbose,
	    "store:s"     => \$store,
	    "database:s"  => \$acedbpath,
	    "outfile:s"   => \$outfile,
            "gtf=s"       => \$gtf_file,
            "species=s"   => \$species,
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

if (not defined $gtf_file) {
  my $suffix = $wormbase->species . ".gtf";
  $gtf_file = join("/", $wormbase->sequences, $suffix);
}

my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or die("Connection failure: ". Ace->error);

$locs = &get_location_data($db, $gtf_file);

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

    #last if $i % 10 == 0;

    print STDERR $obj->name, "\n";
    
    next unless $obj->isObject();
    next unless $obj->Species;
    next unless $obj->Species->name eq $full_name;
    
    my $biotype = $obj->Biotype->name;
    # temp hack to deal with data problem in WS258 (where some genes have "ncRNA" rather than "ncRNA_gene")
    $biotype = "SO:0001263" if $biotype eq "SO:0000655";


    my ($symbol, %synonyms);
    my $seq_name = $obj->Sequence_name->name;
    if ($obj->CGC_name) {
      $symbol = $obj->CGC_name->name;
      $synonyms{$seq_name} = 1;
    } else {
      $symbol = $seq_name;
    }
    map { $synonyms{$_->name} = 1 } $obj->Other_name;
        
    my $desc = "";
    if ($obj->Concise_description) {
      $desc = $obj->Concise_description->name;
    } elsif ($obj->Automated_description) {
      $desc = $obj->Automated_description->name; 
    }
    
    my @xrefs;
    foreach my $dblink ($obj->Database) {
      if (exists $XREF_MAP{$dblink}) {
        my $prefix = $XREF_MAP{$dblink};
        my $suffix = $dblink->right->right->name; 
        push @xrefs, "$prefix:$suffix";
      }
    }
    push @xrefs, "ENSEMBL:" . $obj->name;
    
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
      $name = "";
    }
    
    my @secondary_ids;
    if ($obj->Acquires_merge) {
      foreach my $g ($obj->Acquires_merge) {
        push @secondary_ids, "WB:" . $g->name;
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
      primaryId          => "WB:" . $obj->name,
      symbol             => $symbol,
      soTermId           => $biotype,
      taxonId            => "NCBITaxon:" . $taxid,
      geneLiteratureUrl  => "http://www.wormbase.org/species/c_elegans/gene/" . $obj->name ."#0e--10",
    };

    $json_gene->{name}              =  $name if $name;
    $json_gene->{geneSynopsis}      =  $desc if $desc;
    $json_gene->{synonyms}          =  [sort keys %synonyms] if keys %synonyms;
    $json_gene->{secondaryIds}      =  \@secondary_ids if @secondary_ids;
    $json_gene->{crossReferenceIds} =  \@xrefs if @xrefs;

    if (defined $locs) {
      if (exists $locs->{$obj->name}) {
        $json_gene->{genomeLocations} = [$locs->{$obj->name}];
      } else {
        # do not include genes without a location for now
        next;
      }
    }

    push @genes, $json_gene;
  }
}

my $data = {
  metaData => $meta_data,
  data => \@genes,
};


if (defined $outfile) {
  open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";
} else {
  $out_fh = \*STDOUT;
}

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

##############################################

sub get_location_data {
  my ($acedb, $gtf) = @_;

  #
  # get the assembly name for the canonical bioproject
  #

  my ($assembly_name, $fh, %locs);
  
  my $species_obj = $acedb->fetch(-class => 'Species', -name => $wormbase->full_name);
  my @seq_col = $species_obj->at('Assembly');
        
  foreach my $seq_col_name (@seq_col) {
    my ($bioproj);

    my $seq_col = $seq_col_name->fetch;
    my $this_assembly_name = $seq_col->Name;

    my @db = $seq_col->at('Origin.DB_info.Database');
    foreach my $db (@db) {
      if ($db->name eq 'NCBI_BioProject') {
        $bioproj = $db->right->right->name;
      }
    }
    
    if (defined $bioproj and $wormbase->ncbi_bioproject eq $bioproj) {
      $assembly_name = $this_assembly_name->name;
      last;
    }
  }
    
  if (not defined $assembly_name) {
    die "Could not find name of current assembly for " . $wormbase->species() . "\n";
  }

  #
  # parse the GTF
  #
  
  if ($gtf =~ /\.gz$/) {
    open( $fh, "gunzip -c $gtf |") 
        or die "Could not open gunzip stream GTF file $gtf\n";
  } else {
    open( $fh, $gtf) 
        or die "Could not open GTF file $gtf\n";
  }

  while(<$fh>) {
    /^\#/ and next;
    
    my @l = split(/\t/, $_);
    next if $l[2] ne "gene";

    my ($wbgid) = $l[8] =~ /gene_id\s+\"(\S+)\"/; 

    $locs{$wbgid} = {
      assembly       => $assembly_name,
      chromosome     => $l[0],
      startPosition  => $l[3] + 0, 
      endPosition    => $l[4] + 0,
      strand         => $l[6],
    };
  }
      
  return \%locs;
}
