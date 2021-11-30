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
use Bio::EnsEMBL::Production::Utils::ProductionDbUpdater;
use FindBin;
use lib "$FindBin::Bin/../lib";
use WormBase2Ensembl;

my ( $debug, @species, @notspecies, $allspecies, $setup, $dna, $genes, $rules, $inputids, $meta, $pipeline_setup, $test, $yfile, $recognise_sources );

GetOptions(
  'species=s@'       => \@species,
  'notspecies=s@'    => \@notspecies,
  'allspecies'       => \$allspecies,
  'setup'            => \$setup,
  'load_meta'        => \$meta,
  'load_dna'         => \$dna,
  'load_genes'       => \$genes,
  'load_pipeline'    => \$pipeline_setup,
  'load_rules'       => \$rules,
  'load_iids'        => \$inputids,
  'debug'            => \$debug,
  'test'             => \$test,
  'yfile=s'          => \$yfile,
  'recognise_sources' => \$recognise_sources

) || die("bad commandline parameter\n");


die "You must supply a valid YAML config file\n" if not defined $yfile or (not -e $yfile and $yfile ne "-" ) ;

my $global_config = $yfile eq "-" ? YAML::Load(do {local $/; <>}) : YAML::LoadFile($yfile);
my $generic_config = $global_config->{generics};

my $cvsDIR;
$cvsDIR = $ENV{ENSEMBL_CVS_ROOT_DIR} if exists $ENV{ENSEMBL_CVS_ROOT_DIR};
# explicit definition in the YAML file over-rides environment definition
$cvsDIR = $generic_config->{cvsdir} if exists  $generic_config->{cvsdir};

die "You must define a location of the Ensembl code " .
    "either with cvsdir (in the YAML file) " .
    "or ENSEMBL_CVS_ROOT_DIR in the environment\n" if not defined $cvsDIR; 

my $ensembl_version;
$ensembl_version = $ENV{ENSEMBL_VERSION} if exists $ENV{ENSEMBL_VERSION};
$ensembl_version = $generic_config->{ensembl_version} if exists $generic_config->{ensembl_version};

die "You must define the version of the Ensembl code that you're using " .
    "either with ensembl_version (in the YAML file) " .
    "or ENSEMBL_VERSION in the environment\n" if not defined $ensembl_version;

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

  foreach my $key (keys %$generic_config) {
    if ($key =~ /_database/) {
      foreach my $key2 (keys %{$generic_config->{$key}}) {
        if (not $myconfig->{$key} or not $myconfig->{$key}->{$key2}) {
          $myconfig->{$key}->{$key2} = $generic_config->{$key}->{$key2};
        }
      }
    } else {
      if (not $myconfig->{$key}) {
        $myconfig->{$key} = $generic_config->{$key};
      } else {
        $myconfig->{$key} = join(",", $generic_config->{$key}, $myconfig->{$key});
      }
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

  my $db = $config->{core_database};
  my $tax_db = $config->{taxonomy_database};
  my $prod_db = $config->{production_database};

  print STDERR ">>creating new database $db->{dbname} on $db->{host}\n";

  my $mysql = "mysql -h $db->{host} -P $db->{port} -u $db->{user} --password=$db->{password}";
  
  eval {
    print STDERR "Recreating database from scratch...\n";
    system("$mysql -e \"DROP DATABASE IF EXISTS $db->{dbname};\"") && die;
    system("$mysql -e \"create database $db->{dbname};\"")         && die;

    print STDERR "loading table.sql from ensembl...\n";
    system("$mysql $db->{dbname} < " . $cvsDIR . "/ensembl/sql/table.sql" ) && die;
    
    print STDERR "loading table.sql from ensembl-pipeline...\n";
    system("$mysql $db->{dbname} < " . $cvsDIR . "/ensembl-pipeline/sql/table.sql" ) && die("Could not load pipeline tables\n");
    
    my $species_lookup_params;
    if (exists $config->{taxon_id}) {
      $species_lookup_params = "-taxon_id $config->{taxon_id}";
    } elsif (exists $config->{species}) {
      $species_lookup_params = "-name \"$config->{species}\"";
    } else {
      die "Either taxon_id or species must be defined for $species to load the taxonomy\n";
    }

    print STDERR "Loading taxonomy...\n";
    my $cmd = "perl $cvsDIR/ensembl-pipeline/scripts/load_taxonomy.pl $species_lookup_params "
        . "-taxondbhost $tax_db->{host} " 
        . "-taxondbport $tax_db->{port} "
        . "-taxondbname $tax_db->{dbname} "
        . "-lcdbhost $db->{host} "
        . "-lcdbport $db->{port} "
        . "-lcdbname $db->{dbname} "
        . "-lcdbuser $db->{user} "
        . "-lcdbpass $db->{password}";
    print STDERR "$cmd\n";        
    system($cmd) and die "Could not load taxonomy\n";
    
    print STDERR "Loading production table...\n";

    if ($ensembl_version < 100){
    	$cmd = "perl $cvsDIR/ensembl-production/scripts/production_database/populate_production_db_tables.pl "
            . "--host $db->{host} "
            . "--user $db->{user} "
            . "--pass $db->{password} "
            . "--port $db->{port} "
            . "--database $db->{dbname} "
            . "--mhost $prod_db->{host} "
            . "--mport $prod_db->{port} "
            . "--muser $prod_db->{user} "
            . "--mdatabase $prod_db->{dbname} "
	    . "--dropbaks "
	    . "--dumppath /tmp/ ";

       print STDERR "$cmd\n";
       system($cmd) and die "Could not populate production tables\n";
    }

   # populate_production_db_tables.pl is no longer supported in E! 100 and above.
   # now use the ProductionDbUpdater module instead.
 
   else{
       my $prod_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
           -user    => $prod_db->{user},	
           -dbname  => $prod_db->{dbname},
           -host    => $prod_db->{host},
           -port    => $prod_db->{port} 
       );

       my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
           -host   => $db->{host},
           -user   => $db->{user},
           -dbname => $db->{dbname},
           -pass   => $db->{password},
           -port   => $db->{port},
       );

       &populate_prod_tables($prod_dba, $dba);
   }

    my $db_opt_string = sprintf("-dbhost %s -dbport %s -dbuser %s -dbpass %s -dbname %s", 
                                $db->{host},
                                $db->{port},
                                $db->{user},
                                $db->{password},
                                $db->{dbname});

    my @ana_conf_files = split(/,/, $config->{analysisconf});
    
    foreach my $cfile (@ana_conf_files) {
      if (-e $cfile) {
        print STDERR "Loading analyses...\n";
        $cmd = "perl $FindBin::Bin/analysis_setup.pl $db_opt_string -read -file $cfile";
        print STDERR "Running: $cmd\n";
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

  my $db = $config->{core_database};
  
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
        print STDERR "Running: $cmd\n";
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
      print STDERR "Running: $cmd\n";
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
        print STDERR "Running: $cmd\n";
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
  print STDERR "Running: $cmd\n";
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
    print STDERR "Running: $cmd\n";
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
    print STDERR "Running: $mysql -e '$sql'\n";
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

  my $db = $config->{core_database};

  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $db->{host},
    -user   => $db->{user},
    -dbname => $db->{dbname},
    -pass   => $db->{password},
    -port   => $db->{port},
      );

  my (%ana_hash);

  &empty_out_old_gene_set($dba);

  foreach my $ana (@{$dba->get_AnalysisAdaptor->fetch_all}) {
    $ana_hash{$ana->logic_name} = $ana;
  }

  my $cod_logic     = (exists $config->{gff_codinggene_logic}) 
      ? $config->{gff_codinggene_logic} 
      : "wormbase";
  my $nc_logic      = (exists $config->{gff_noncodinggene_logic}) 
      ? $config->{gff_noncodinggene_logic} 
      : "wormbase_non_coding";
  my $pseudo_logic  = (exists $config->{gff_pseudogene_logic}) 
      ? $config->{gff_pseudogene_logic} 
  : "wormbase_pseudogene";

  foreach my $logic ($cod_logic, $nc_logic, $pseudo_logic) {
    if (not exists $ana_hash{$logic}) {
      my $ana = Bio::EnsEMBL::Analysis->new(-logic_name => $logic,
                                            -gff_source => "WormBase_imported",
                                            -module     => "WormBase");
      $dba->get_AnalysisAdaptor->store($ana);
      $ana_hash{$logic} = $ana;
    }
  }
  my $cod_analysis    = $ana_hash{$cod_logic};
  my $nc_analysis     = $ana_hash{$nc_logic};
  my $pseudo_analysis = $ana_hash{$pseudo_logic};

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
				     -recognise_sources => ($recognise_sources) ? 1 : 0,
                                     -verbose => 1,
                                     -ignoregffphases => (exists $config->{gff_ignore_phases}) 
                                       ? $config->{gff_ignore_phases} 
                                       : 0);

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
                                     ['tRNA_pseudogene', 'tRNA_pseudogene'],
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
    my $source_hash;
    if ($config->{gff_sources}) {
      $source_hash = {};
      foreach my $source (split(/,/, $config->{gff_sources})) {
        $source_hash->{$source} = 1;
      }
    } else {
      $source_hash = { WormBase => 1, WormBase_imported => 1 };
    }

    if (scalar(@gff3_files) == 1) {
      @genes = @{$wb2ens->parse_genes_gff3( $gff3_files[0], $cod_analysis, $nc_analysis, $pseudo_analysis, $source_hash )};
    } else {
      open(my $gff_fh, "cat @gff3_files |") or die "Could not create GFF stream\n";
      @genes = @{$wb2ens->parse_genes_gff3_fh( $gff_fh, $cod_analysis, $nc_analysis, $pseudo_analysis, $source_hash )};
    }
    if (scalar(@genes) == 0) {
      die "Could not extract any genes from GFF3 file. Exiting\n";
    }
  } else {
    die "No gff or gff3 files found - death\n";
  }

  if (exists $config->{fix_phases} and $config->{fix_phases}) {
    $wb2ens->phase_fix(\@genes );
  }

  $wb2ens->write_genes( \@genes );

  # correct seleno proteins
  if ($config->{seleno}) {
    my %seleno_genes = map { $_ => 1 } split(/,/, $config->{seleno});

    foreach my $gene (@genes) {
      if (exists $seleno_genes{$gene->stable_id}) {
        print STDERR "Fixing selenoprotein " . $gene->stable_id . "\n";
        foreach my $t (@{$gene->get_all_Transcripts}) {
          $wb2ens->translation_fix($t);
        }
      }
    }
  }
    
  my $set_canon_cmd = "perl $cvsDIR/ensembl/misc-scripts/canonical_transcripts/select_canonical_transcripts.pl "
      . "-dbhost $db->{host} "
      . "-dbuser $db->{user} "
      . "-dbpass $db->{password} "
      . "-dbname $db->{dbname} "
      . "-dbport $db->{port} "
      . "-coord toplevel "
      . "-write";
  print STDERR "Running: $set_canon_cmd\n";
  system($set_canon_cmd) and die "Could not set canonical transcripts\n";

  my $timestamp = strftime("%Y-%m", localtime(time))."-WormBase";
  my $versionstamp = $timestamp;
  $versionstamp = $config->{"meta.genebuild.version"} if exists $config->{"meta.genebuild.version"};
  $versionstamp = $config->{"meta"}{"genebuild.version"} if exists $config->{"meta"}{"genebuild.version"}; 

  $dba->dbc->do('DELETE FROM meta WHERE meta_key = "genebuild.start_date"');
  $dba->dbc->do("INSERT INTO meta (meta_key,meta_value) VALUES (\"genebuild.start_date\",\"$timestamp\")");
  
  $dba->dbc->do('DELETE FROM meta WHERE meta_key = "genebuild.version"');
  $dba->dbc->do("INSERT INTO meta (meta_key,meta_value) VALUES (\"genebuild.version\",\"$versionstamp\")");

}


#############################################
sub load_rules {
  my ($species, $config) = @_;
  
  my $db = $config->{core_database};

  my @conf_files = split(/,/, $config->{ruleconf});
    
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
      print STDERR "Loading rules from $cfile...\n";
      my $cmd = "$load_rule_base $cfile";
      print STDERR "Running: $cmd\n";
      system($cmd) and die "Could not load analyses from $cfile\n";
    } else {
      die "Could not find analysis config file $cfile\n";
    }
  }
}



#############################################
sub load_input_ids {
  my ($species, $config) = @_;
  
  my $db = $config->{core_database};

  my $load_input_ids_base =  "perl $cvsDIR/ensembl-pipeline/scripts/make_input_ids "
      . "-dbhost $db->{host} "
      . "-dbuser $db->{user} "
      . "-dbpass $db->{password} "
      . "-dbname $db->{dbname} "
      . "-dbport $db->{port} ";

  my $slice_cmd = "$load_input_ids_base -slice -slice_size 75000 -coord_system toplevel -logic_name submitslice75k -input_id_type SLICE75K";
  my $trids_cmd = "$load_input_ids_base -translation_id -logic submittranslation";

  foreach my $cmd ($slice_cmd, $trids_cmd) {
    print STDERR "Running: $cmd\n";
    system($cmd) and die "Could not successfully make input ids\n";
  }

}

###############################################
sub load_meta_table {
  my ($species, $config) = @_;
  
  my $db = $config->{core_database};

  my $meta = $config->{meta};
  if ($meta){
    while (my ($k,$v)=each %$meta){
      &update_meta($db,$k,$v);
    }
  }
  else{
   foreach my $key (keys %$config) {
    if ($key =~ /^meta\.(\S+)/) {
      my $db_key = $1;
      my $val = $config->{$key};
      &update_meta($db,$db_key,$val);      
    }
   }
  }
}

sub update_meta{
   my ($db,$db_key,$db_val)=@_;

   my $mysql = "mysql -h $db->{host} -P $db->{port} -u $db->{user} --password=$db->{password} -D $db->{dbname}"; 
 
   print STDERR "Populating meta table for $db->{dbname}...\n";

   system("$mysql -e 'DELETE FROM meta WHERE meta_key = \"$db_key\"'") 
       and die "Could not delete $db_key from meta in $db->{dbname}\n";
   system("$mysql -e 'INSERT INTO meta (meta_key,meta_value) VALUES (\"$db_key\",\"$db_val\");'") 
       and die "Could not insert value $db_val for $db_key into meta in $db->{dbname}\n";
}

#################################################
sub empty_out_old_gene_set {
  my ($dba) = @_;

  print STDERR "Emptying out old gene set and associated tables...\n";

  my $dbc = $dba->dbc;

  $dbc->do('TRUNCATE TABLE gene')  or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE gene_attrib')  or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE transcript') or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE transcript_attrib')  or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE exon')  or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE exon_transcript') or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE translation')  or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE translation_attrib') or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE protein_feature')  or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE interpro')  or die $dbc->errstr;
  
  $dbc->do('TRUNCATE TABLE xref') or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE object_xref') or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE ontology_xref') or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE identity_xref') or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE external_synonym') or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE density_type') or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE density_feature') or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE dependent_xref') or die $dbc->errstr;
  $dbc->do('TRUNCATE TABLE genome_statistics') or die $dbc->errstr;
  
}

######################################################
sub populate_prod_tables {

    my ($prod_dba, $dba) = @_;
	
    my $updater = Bio::EnsEMBL::Production::Utils::ProductionDbUpdater->new(
        -PRODUCTION_DBA => $prod_dba
    );

    my $tables = [
      'attrib_type',
      'biotype',
      'external_db',
      'misc_set',
      'unmapped_reason'
     ];
     
     print STDERR "Populating prod tables with ProductionDbUpdater...\n";
     $updater->update_controlled_tables($dba->dbc, $tables);
	
}
