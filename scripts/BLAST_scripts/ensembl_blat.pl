#!/usr/bin/env perl
#
# DESCRIPTION:
#   setting up the BLAT pipeline
#
# Last edited by: $Author: klh $
# Last edited on: $Date: 2012-01-13 16:16:57 $

use lib $ENV{'CVS_DIR'};

use constant USAGE => <<HERE;
ensembl_blat.pl options:
            -debug USER_NAME    sets email address and debug mode
            -store FILE_NAME    use a Storable wormbase configuration file
            -species SPECIES_NAME species name a.e. elegans
            -user NAME          database user name
            -password PASSWORD  database password
	    -host DBHOSTNAME    database host
	    -port DBPORT        database port
HERE

use Getopt::Long;
use Storable;

use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;

# some magic to turn off the deprecation warnings
use Bio::EnsEMBL::Utils::Exception;
verbose('OFF');

use Wormbase;
use strict;

my $CONFIG = {
  query_file_location => '/nfs/wormpub/BUILD/cDNA',
  logic_names => {
    'blat_brugia_ests'    => 'brugia/EST.masked*',
    'blat_brugia_cdnas'   => 'brugia/mRNA.masked*',

    'blat_elegans_ests'   => 'elegans/EST.masked*',
    'blat_elegans_cdnas'  => 'elegans/mRNA.masked*',
    'blat_elegans_osts'   => 'elegans/OST.masked*',
    'blat_elegans_rsts'   => 'elegans/RST.masked*',
    'blat_elegans_ncrnas' => 'elegans/ncRNA.masked*',
    'blat_elegans_tc1s'   => 'elegans/tc1.masked*',
    
    'blat_briggsae_ests'  => 'briggsae/EST.masked*',
    'blat_briggsae_cdnas' => 'briggsae/mRNA.masked*',
    
    'blat_brenneri_ests'  => 'brenneri/EST.masked*',
    'blat_brenneri_cdnas' => 'brenneri/mRNA.masked*',
    
    'blat_remanei_ests'   => 'remanei/EST.masked*',
    'blat_remanei_cdnas'  => 'remanei/mRNA.masked*',
    
    'blat_japonica_ests'  => 'japonica/EST.masked*',
    'blat_japonica_cdnas' => 'japonica/mRNA.masked*',
    
    'blat_heterorhabditis_ests'  => 'heterorhabditis/EST.masked*',
    'blat_heterorhabditis_cdnas' => 'heterorhabditis/mRNA.masked*',
    
    'blat_pristionchus_ests'  => 'pristionchus/EST.masked*',
    'blat_pristionchus_cdnas' => 'pristionchus/mRNA.masked*',
    
    'blat_nembase_ests'   => 'nembase/EST.masked*',
    'blat_washu_ests'     => 'washu/EST.masked*',
    'blat_nematode_ests'  => 'nematode/EST.masked*',
  },
};

my($debug,$store,$port,$user,$password,$species,$host,$dbname);
GetOptions(
 'debug=s'    => \$debug,
 'store=s'    => \$store,
 'user=s'     => \$user,
 'password=s' => \$password,
 'species=s'  => \$species,
 'host=s'     => \$host,
 'port=s'     => \$port,
 'dbname=s'   => \$dbname,
)||die(USAGE);

# WormBase setup
my $wormbase;
 if ($store) {
    $wormbase = Storable::retrieve($store)
      or croak("Can't restore wormbase from $store\n");
} else {
    $wormbase = Wormbase->new(
        -debug    => $debug,
        -organism  => $species,
    );
}

my $database = sprintf('worm_ensembl_%s',lc(ref $wormbase));

# more setup
my $log = Log_files->make_build_log($wormbase);
$log->write_to("Updating BLAT input files for $database\n");

$host||='farmdb1';
$port||=3306;

# MYSQL setup
my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
        -host     => $host,
        -user     => $user,
        -dbname   => $dbname || $database,
        -port     => $port,
        -driver   => 'mysql', 
    );
$db->password($password);

foreach my $ana (&get_all_blat_analyses()) {
  my $logic_name = $ana->logic_name();

  # a.) clean out all dna_align_features
  &clean_dna_align_features($ana);

  my $dependent_ana = &get_dependent_analysis($ana);
  &clean_input_ids($ana, $dependent_ana);
  &make_input_ids($ana, $dependent_ana);
}

$log->write_to("Finished.\n");
$log->mail();

###########################
sub get_dependent_analysis {
  my ($ana) = @_;

  my $logic_name = $ana->logic_name;

  my $rule = $db->get_RuleAdaptor->fetch_by_goal($ana);
  die("No rule for $logic_name\n") if not defined $rule;

  my $cond_list = $rule->list_conditions;
  die("Expecting exactly one condition for $logic_name\n")
      if not defined $cond_list or scalar(@$cond_list) != 1;
  
  my $dependent_ana = $db->get_AnalysisAdaptor->fetch_by_logic_name($cond_list->[0]);
  die("Could not find dependent analysis for $logic_name\n")
      if not defined $dependent_ana;
  
  die(sprintf("Input_id type for $logic_name (%s) and its dependent analysis (%s) do not agree\n", 
              $ana->input_id_type, $dependent_ana->input_id_type))
      if $ana->input_id_type ne $dependent_ana->input_id_type;

  return $dependent_ana;

}

#######################
sub make_input_ids {
  my ($ana, $dep_ana) = @_;

  my $sth=$db->prepare(
    'INSERT INTO input_id_analysis (input_id,input_id_type,analysis_id,created,result) VALUES (?,?,?,NOW(),0)'
      );

  my $file_base = $CONFIG->{query_file_location};

  if (exists $CONFIG->{logic_names}->{$ana->logic_name}) {
    my $file_pattern = sprintf("%s/%s", $file_base, $CONFIG->{logic_names}->{$ana->logic_name});

    my @files = glob($file_pattern);
    foreach my $file (@files){
      $sth->execute($file,$dep_ana->input_id_type,$dep_ana->dbID);
    }
  } else {
    print STDERR "Nothing in config for logic_name " . $ana->logic_name .  ", so skipping\n";
  }
}

##########################
sub clean_dna_align_features {
  my @ana = @_;

  foreach my $ana (@ana) {
    my $dbID=$ana->dbID;
    $db->dbc->do("DELETE FROM dna_align_feature WHERE analysis_id = $dbID");
  }
}


##########################
sub clean_input_ids {
  my @ana = @_;

  foreach my $ana (@ana) {
    my $dbID=$ana->dbID;
    $db->dbc->do("DELETE FROM input_id_analysis WHERE analysis_id = $dbID");
  }
}

##########################
sub get_all_blat_analyses {
  my @blat_ana;
  foreach my $ana (@{$db->get_AnalysisAdaptor->fetch_all}) {
    if ($ana->program_file =~ /blat$/) {
      push @blat_ana, $ana;
    } 
  }

  return @blat_ana;
}


