#!/usr/bin/env perl
#===============================================================================
#         FILE:  worm_lite.pl
#===============================================================================

use strict;
use POSIX;
use YAML;
use Getopt::Long;
use Storable;

use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use FindBin;
use lib "$FindBin::Bin/../lib";
use WormBase2Ensembl;

my ( $debug, @species, @notspecies, $allspecies, $setup, $dna, $genes, $rules, $inputids, $meta, $pipeline_setup, $test, $yfile );

GetOptions(
  'species=s@'    => \@species,
  'notspecies=s@' => \@notspecies,
  'allspecies'    => \$allspecies,
  'setup'         => \$setup,
  'load_meta'     => \$meta,
  'load_dna'      => \$dna,
  'load_genes'    => \$genes,
  'load_pipeline' => \$pipeline_setup,
  'load_rules'    => \$rules,
  'load_iids'     => \$inputids,
  'debug'         => \$debug,
  'test'          => \$test,
  'yfile=s'       => \$yfile,

) || die("bad commandline parameter\n");


die "You must supply a valid YAML config file\n" if not defined $yfile or not -e $yfile;

my $global_config = YAML::LoadFile($yfile);
my $generic_config = $global_config->{generics};

my $cvsDIR = $generic_config->{cvsdir};

my ($prod_db_host, $prod_db_port, $prod_db_name) = 
    ($generic_config->{ensprod_host},
     $generic_config->{ensprod_port},
     $generic_config->{ensprod_dbname});

my ($tax_db_host, $tax_db_port, $tax_db_name) = 
    ($generic_config->{taxonomy_host},
     $generic_config->{taxonomy_port},
     $generic_config->{taxonomy_dbname});

if ($allspecies) {
  die "You cannot supply both -species and -allspecies!\n" if @species;

  @species = sort grep { $_ ne 'generics' } keys %$global_config;
  
  if (@notspecies) {
    my %not_species;
    map { $not_species{$_} = 1 } @notspecies;

    @species = grep { not exists $not_species{$_} } @species;

  }
} elsif (not @species) {
  die "You must supply either -species or -allspecies\n";
}

#
# check that all species are defined in the config
#
foreach my $species (@species) {
  if (not exists $global_config->{$species}) {
    die "Could not find config entry for species $species\n";
  }

  if ($test and $species !~ /_test/) {
    die "In test mode, all species entries in config should have a _test suffix (safety net)\n";
  }
}
foreach my $myspecies (@species) {
  my $myconfig = $global_config->{$myspecies};

  foreach my $db_property ('host', 'port', 'user', 'password') {
    if (not $myconfig->{database}->{$db_property}) {
      $myconfig->{database}->{$db_property} = $generic_config->{database}->{$db_property};
    }
  }

  &setupdb($myspecies, $myconfig)         if $setup;
  &load_meta_table($myspecies, $myconfig) if $setup or $meta;
  &load_assembly($myspecies, $myconfig)   if $dna;
  &load_genes($myspecies, $myconfig)      if $genes;
  &load_rules($myspecies, $myconfig)      if $rules or $pipeline_setup;
  &load_input_ids($myspecies, $myconfig)  if $inputids or $pipeline_setup;

}  

exit(0);


#########################################
sub setupdb {
  my ($species, $config) = @_;

  my $db = $config->{database};

  print ">>creating new database $db->{dbname} on $db->{host}\n";

  my $mysql = "mysql -h $db->{host} -P $db->{port} -u $db->{user} --password=$db->{password}";
  
  eval {
    print "Recreating database from scratch...\n";
    system("$mysql -e \"DROP DATABASE IF EXISTS $db->{dbname};\"") && die;
    system("$mysql -e \"create database $db->{dbname};\"")         && die;

    print "loading table.sql from ensembl...\n";
    system("$mysql $db->{dbname} < " . $cvsDIR . "/ensembl/sql/table.sql" ) && die;
    
    print "loading table.sql from ensembl-pipeline...\n";
    system("$mysql $db->{dbname} < " . $cvsDIR . "/ensembl-pipeline/sql/table.sql" ) && die("Could not load pipeline tables\n");
    
    my $species_lookup_params;
    if (exists $config->{taxon_id}) {
      $species_lookup_params = "-taxon_id $config->{taxon_id}";
    } elsif (exists $config->{species}) {
      $species_lookup_params = "-name \"$config->{species}\"";
    } else {
      die "Either taxon_id or species must be defined for $species to load the taxonomy\n";
    }

    print "Loading taxonomy...N";
    my $cmd = "perl $cvsDIR/ensembl-pipeline/scripts/load_taxonomy.pl $species_lookup_params "
        . "-taxondbhost $tax_db_host " 
        . "-taxondbport $tax_db_port "
        . "-taxondbname $tax_db_name "
        . "-lcdbhost $db->{host} "
        . "-lcdbport $db->{port} "
        . "-lcdbname $db->{dbname} "
        . "-lcdbuser $db->{user} "
        . "-lcdbpass $db->{password}";
    print "$cmd\n";        
    system($cmd) and die "Could not load taxonomy\n";
    
    print "Loading production table...\n";
    $cmd = "perl $cvsDIR/ensembl-production/scripts/production_database/populate_production_db_tables.pl "
        . "--host $db->{host} "
        . "--user $db->{user} "
        . "--pass $db->{password} "
        . "--port $db->{port} "
        . "--database $db->{dbname} "
        . "--mhost $prod_db_host "
        . "--mport $prod_db_port "
        . "--mdatabase $prod_db_name "
	. "--dropbaks "
	. "--dumppath /tmp/ ";
    print "$cmd\n";
    system($cmd) and die "Could not populate production tables\n";

    my $db_opt_string = sprintf("-dbhost %s -dbport %s -dbuser %s -dbpass %s -dbname %s", 
                                $db->{host},
                                $db->{port},
                                $db->{user},
                                $db->{password},
                                $db->{dbname});

    my @ana_conf_files;

    if ($generic_config->{analysisconf}) {
      push @ana_conf_files, $generic_config->{analysisconf};
    }
    if ($config->{analysisconf}) {
      push @ana_conf_files, $config->{analysisconf};
    }
    
    foreach my $cfile (@ana_conf_files) {
      if (-e $cfile) {
        print "Loading analyses...\n";
        $cmd = "perl $FindBin::Bin/analysis_setup.pl $db_opt_string -read -file $cfile";
        print "Running: $cmd\n";
        system($cmd) and die "Could not load analyses from $cfile\n";
      } else {
        die "Could not find analysis config file $cfile\n";
      }
    }


  };
  $@ and die("Error while building the database: $@");
}

##############################################
sub load_assembly {
  my ($species, $config) = @_;

  my $db = $config->{database};
  
  my $seq_level_coord_sys = $config->{seqlevel};
  my $top_level_coord_sys = $config->{toplevel};
  my $coord_sys_ver = $config->{assembly_version};

  my $seq_level_rank = ($seq_level_coord_sys eq $top_level_coord_sys) ? 1 : 2;
  my $top_level_rank = 1;

  if ($config->{agp}) {
    foreach my $glb (split(/,/, $config->{agp})) {
      foreach my $agp (glob("$glb")) {
        my $cmd = "perl $cvsDIR/ensembl-pipeline/scripts/load_seq_region.pl "
            . "-dbhost $db->{host} "
            . "-dbuser $db->{user} "
            . "-dbpass $db->{password} "
            . "-dbname $db->{dbname} "
            . "-dbport $db->{port} "
            . "-coord_system_name $top_level_coord_sys "
            . "-coord_system_version $coord_sys_ver "
            . "-rank $top_level_rank "
            . "-default_version "
            . "-agp_file $agp";
        print "Running: $cmd\n";
        system($cmd) and die "Could not load seq_regions from agp file\n";
      }
    }
  }

  foreach my $glb (split(/,/, $config->{fasta})) {
    foreach my $fasta (glob("$glb")) {
      my $cmd = "perl $cvsDIR/ensembl-pipeline/scripts/load_seq_region.pl "
          . "-dbhost $db->{host} "
          . "-dbuser $db->{user} "
          . "-dbpass $db->{password} "
          . "-dbname $db->{dbname} "
          . "-dbport $db->{port} "
          . "-coord_system_name $seq_level_coord_sys "
          . "-coord_system_version $coord_sys_ver "
          . "-rank $seq_level_rank "
          . "-default_version -sequence_level "
          . "-fasta_file $fasta";
      print "Running: $cmd\n";
      system($cmd) and die "Could not load seq_regions fasta file\n";
    }
  }

  if ($config->{agp}) {
    foreach my $glb (split(/,/, $config->{agp})) {
      foreach my $agp (glob("$glb")) {
        my $cmd = "perl $cvsDIR/ensembl-pipeline/scripts/load_agp.pl "
            . "-dbhost $db->{host} "
            . "-dbuser $db->{user} "
            . "-dbpass $db->{password} "
            . "-dbname $db->{dbname} "
            . "-dbport $db->{port} "
            . "-assembled_name $top_level_coord_sys "
            . "-assembled_version $coord_sys_ver "
            . "-component_name $seq_level_coord_sys "
            . "-component_version $coord_sys_ver "
            . "-agp_file $agp";
        print "Running: $cmd\n";
        system($cmd) and die "Could not load the assembly table from agp file\n";
      }
    }
  }
  
  my $cmd = "perl $cvsDIR/ensembl-pipeline/scripts/set_toplevel.pl "
      . "-dbhost $db->{host} "
      . "-dbport $db->{port} "
      . "-dbuser $db->{user} "
      . "-dbpass $db->{password} "
      . "-dbname $db->{dbname}" ;
  print "Running: $cmd\n";
  system($cmd) and die "Could not set toplevel\n";
  
  if ($config->{mitochondrial}) {
    my @mito_seqs = split(/,/, $config->{mitochondrial});

    $cmd = "perl $FindBin::Bin/set_codon_table.pl "
        . "-dbhost $db->{host} "
        . "-dbport $db->{port} "
        . "-dbuser $db->{user} "
        . "-dbpass $db->{password} "
        . "-dbname $db->{dbname} "
        . "-codontable 5 "
        . "@mito_seqs";
    print "Running: $cmd\n";
    system($cmd) and die "Could not set the mitochrondrial table";
  }

  # Finally, for elegans and briggsae, append the chromosome prefices (yuk, but it has to 
  # done to make BLAST dumping etc work properly
  if ($species eq 'elegans' or $species eq 'briggsae') {
    my $prefix = ($species eq 'elegans') ? 'CHROMOSOME_' : 'chr';

    my $mysql = "mysql -h $db->{host} -P $db->{port} -u $db->{user} --password=$db->{password} -D $db->{dbname}";
    my $sql = "UPDATE seq_region, coord_system "
        . "SET seq_region.name = CONCAT(\"$prefix\", seq_region.name) "
        . "WHERE seq_region.coord_system_id = coord_system.coord_system_id "
        . "AND coord_system.name = \"chromosome\"";
    print "Running: $mysql -e '$sql'\n";
    system("$mysql -e '$sql'") and die "Could not add chromosome prefixes to chromosome names\n";
  }


  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $db->{host},
    -user   => $db->{user},
    -dbname => $db->{dbname},
    -pass   => $db->{password},
    -port   => $db->{port},
      );
  
  foreach my $k ("assembly.default", "assembly.name") {
    $dba->dbc->do("DELETE FROM meta WHERE meta_key = \"$k\"");
    $dba->dbc->do("INSERT INTO meta (meta_key,meta_value) VALUES (\"$k\",\"$coord_sys_ver\")");
  }  

}

#############################################
sub load_genes {
  my ($species, $config) = @_;

  my $db = $config->{database};

  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $db->{host},
    -user   => $db->{user},
    -dbname => $db->{dbname},
    -pass   => $db->{password},
    -port   => $db->{port},
      );

  my (%ana_hash);

  foreach my $ana (@{$dba->get_AnalysisAdaptor->fetch_all}) {
    $ana_hash{$ana->logic_name} = $ana;
  }
  foreach my $logic ("wormbase", "wormbase_non_coding", "wormbase_pseudogene") {
    if (not exists $ana_hash{$logic}) {
      my $ana = Bio::EnsEMBL::Analysis->new(-logic_name => $logic,
                                            -gff_source => "WormBase",
                                            -module     => "WormBase");
      $dba->get_AnalysisAdaptor->store($ana);
      $ana_hash{$logic} = $ana;
    }
  }
  my $cod_analysis    = $ana_hash{wormbase};
  my $nc_analysis     = $ana_hash{wormbase_non_coding};
  my $pseudo_analysis = $ana_hash{wormbase_pseudogene};

  my (@path_globs, @gff2_files, @gff3_files); 

  if ($config->{gff}) {
    @path_globs = split(/,/, $config->{gff});
    foreach my $fglob (@path_globs) {
      my @files = glob($fglob);
      if (grep { not -e $_ }  @files) {
        die "GFF2 file '$fglob' could not be resolved to a real file name\n";
      }
      push @gff2_files, @files;
    }
  } elsif ($config->{gff3}) {
    @path_globs = split(/,/, $config->{gff3});
    foreach my $fglob (@path_globs) {
      my @files = glob($fglob);
      if (grep { not -e $_ }  @files) {
        die "GFF3 file '$fglob' could not be resolved to a real file name\n";
      }
      push @gff3_files, @files;
    }
  } else {
    die "No gff or gff3 stanza in config - death\n";
  }

  my $wb2ens = WormBase2Ensembl->new(-species => $species,
                                     -dbh     => $dba,
                                     -debug   => ($debug) ? 1 : 0,
                                     -verbose => 1,
      );

  my @genes;
  if (@gff2_files) {
    my $gff_fh;
    #
    # coding genes
    #
    open($gff_fh, "cat @gff2_files |") or die "Could not create GFF stream\n";
    push @genes, @{ $wb2ens->parse_protein_coding_gff2_fh( $gff_fh, $cod_analysis)};

    #
    # non-coding genes
    #
    foreach my $source_biotype_pair (['rRNA', 'rRNA'],
                                     ['ncRNA', 'ncRNA'],
                                     ['snRNA', 'snRNA'],
                                     ['snoRNA', 'snoRNA'],
                                     ['scRNA', 'scRNA'],
                                     ['piRNA', 'piRNA'],
                                     ['lincRNA', 'lincRNA'],
                                     ['asRNA', 'asRNA'],
                                     ['miRNA_mature', 'miRNA']) {
      my ($source, $biotype) = @$source_biotype_pair;

      open($gff_fh, "cat @gff2_files |") or die "Could not create GFF stream\n";
      push @genes, @{ $wb2ens->parse_non_coding_genes_gff2_fh( $gff_fh, $nc_analysis, $source, $biotype ) };
    }

    #
    # pseudogenes
    #
    open($gff_fh, "cat @gff2_files |") or die "Could not create GFF stream\n";
    push @genes, @{ $wb2ens->parse_non_coding_genes_gff2_fh( $gff_fh, $pseudo_analysis, 'Pseudogene', 'pseudogene')};

  } elsif (@gff3_files) {
    if (scalar(@gff3_files) == 1) {
      @genes = @{$wb2ens->parse_genes_gff3( $gff3_files[0], $cod_analysis, $nc_analysis, $pseudo_analysis)};
    } else {
      open(my $gff_fh, "cat @gff3_files |") or die "Could not create GFF stream\n";
      @genes = @{$wb2ens->parse_genes_gff3_fh( $gff_fh, $cod_analysis, $nc_analysis, $pseudo_analysis )};
    }
    if (scalar(@genes) == 0) {
      die "Could not extract any genes from GFF3 file. Exiting\n";
    }
  } else {
    die "No gff or gff3 files found - death\n";
  }

  $wb2ens->write_genes( \@genes );
  
  my $set_canon_cmd = "perl $cvsDIR/ensembl/misc-scripts/canonical_transcripts/select_canonical_transcripts.pl "
      . "-dbhost $db->{host} "
      . "-dbuser $db->{user} "
      . "-dbpass $db->{password} "
      . "-dbname $db->{dbname} "
      . "-dbport $db->{port} "
      . "-coord toplevel "
      . "-write";
  print "Running: $set_canon_cmd\n";
  system($set_canon_cmd) and die "Could not set canonical transcripts\n";

  my $timestamp = strftime("%Y-%m", localtime(time));
  my $versionstamp = "${timestamp}-WormBase";

  $dba->dbc->do('DELETE FROM meta WHERE meta_key = "genebuild.start_date"');
  $dba->dbc->do("INSERT INTO meta (meta_key,meta_value) VALUES (\"genebuild.start_date\",\"$versionstamp\")");
  
  $dba->dbc->do('DELETE FROM meta WHERE meta_key = "genebuild.version"');
  $dba->dbc->do("INSERT INTO meta (meta_key,meta_value) VALUES (\"genebuild.version\",\"$versionstamp\")");

}


#############################################
sub load_rules {
  my ($species, $config) = @_;
  
  my $db = $config->{database};

  my @conf_files;
  foreach my $path ($generic_config->{ruleconf}, $config->{ruleconf}) {
    if ($path) {
      foreach my $file (split(/,/, $path)) {
        if (-e $file) {
          push @conf_files, $file;
        } else {
          die "Rule config file $file could not be found\n";
        } 
      }
    }
  }
    
  my $load_rule_base = "perl $cvsDIR/ensembl-pipeline/scripts/rule_setup.pl "
      . "-dbhost $db->{host} "
      . "-dbuser $db->{user} "
      . "-dbpass $db->{password} "
      . "-dbname $db->{dbname} "
      . "-dbport $db->{port} "
      . "-read "
      . "-file ";
    
  foreach my $cfile (@conf_files) {
    if (-e $cfile) {
      print "Loading rules from $cfile...\n";
      my $cmd = "$load_rule_base $cfile";
      print "Running: $cmd\n";
      system($cmd) and die "Could not load analyses from $cfile\n";
    } else {
      die "Could not find analysis config file $cfile\n";
    }
  }
}



#############################################
sub load_input_ids {
  my ($species, $config) = @_;
  
  my $db = $config->{database};

  my $load_input_ids_base =  "perl $cvsDIR/ensembl-pipeline/scripts/make_input_ids "
      . "-dbhost $db->{host} "
      . "-dbuser $db->{user} "
      . "-dbpass $db->{password} "
      . "-dbname $db->{dbname} "
      . "-dbport $db->{port} ";

  my $slice_cmd = "$load_input_ids_base -slice -slice_size 75000 -coord_system toplevel -logic_name submitslice75k -input_id_type SLICE75K";
  my $trids_cmd = "$load_input_ids_base -translation_id -logic submittranslation";

  foreach my $cmd ($slice_cmd, $trids_cmd) {
    print "Running: $cmd\n";
    system($cmd) and die "Could not successfully make input ids\n";
  }

}

###############################################
sub load_meta_table {
  my ($species, $config) = @_;
  
  my $db = $config->{database};
  
  my $mysql = "mysql -h $db->{host} -P $db->{port} -u $db->{user} --password=$db->{password} -D $db->{dbname}";
  
  print "Populating meta table for $db->{dbname}...\n";
  foreach my $key (keys %$config) {
    if ($key =~ /^meta\.(\S+)/) {
      my $db_key = $1;
      my $val = $config->{$key};
      
      system("$mysql -e 'DELETE FROM meta WHERE meta_key = \"$db_key\"'") 
          and die "Could not delete $db_key from meta in $db->{dbname}\n";
      system("$mysql -e 'INSERT INTO meta (meta_key,meta_value) VALUES (\"$db_key\",\"$val\");'") 
          and die "Could not insert value for $db_key into meta\n";
    }
  }
}
