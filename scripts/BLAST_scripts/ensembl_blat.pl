#!/usr/bin/env perl
#
# DESCRIPTION:
#   setting up the BLAT pipeline
#
# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2015-04-27 11:24:44 $

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
  logic_names => {
    'bmalayi_est'    => 'brugia/EST.masked*',
    'bmalayi_mrna'   => 'brugia/mRNA.masked*',
    
    'celegans_est'   => 'elegans/EST.masked*',
    'celegans_mrna'  => 'elegans/mRNA.masked*',
    'celegans_ost'   => 'elegans/OST.masked*',
    'celegans_rst'   => 'elegans/RST.masked*',
    'celegans_ncrna' => 'elegans/ncRNA.masked*',
    'celegans_tc1'   => 'elegans/tc1.masked*',
    
    'cbriggsae_est'  => 'briggsae/EST.masked*',
    'cbriggsae_mrna' => 'briggsae/mRNA.masked*',
    
    'cbrenneri_est'  => 'brenneri/EST.masked*',
    'cbrenneri_mrna' => 'brenneri/mRNA.masked*',
    
    'cremanei_est'   => 'remanei/EST.masked*',
    'cremanei_mrna'  => 'remanei/mRNA.masked*',
    
    'cjaponica_est'  => 'japonica/EST.masked*',
    'cjaponica_mrna' => 'japonica/mRNA.masked*',
    
    'ppacificus_est'  => 'pristionchus/EST.masked*',
    'ppacificus_mrna' => 'pristionchus/mRNA.masked*',
    
    'ovolvulus_est'  => 'ovolvulus/EST.masked*',
    'ovolvulus_mrna' => 'ovolvulus/mRNA.masked*',
 
    'sratti_est'     => 'sratti/EST.masked*',
    'sratti_mrna'    => 'sratti/mRNA.masked*',
   
    'nembase_nematode_est'   => 'nembase/EST.masked*',
    'washu_nematode_est'     => 'washu/EST.masked*',
    'embl_nematode_est'  => 'nematode/EST.masked*',
  }
};

my($debug,$store,$port,$user,$password,$species,$host,$dbname,$WS_version);
GetOptions(
  'debug=s'    => \$debug,
  'store=s'    => \$store,
  'user=s'     => \$user,
  'password=s' => \$password,
  'species=s'  => \$species,
  'host=s'     => \$host,
  'port=s'     => \$port,
  'dbname=s'   => \$dbname,
  'version=s'  => \$WS_version,
)||die(USAGE);

# WormBase setup
my $wormbase;
 if ($store) {
    $wormbase = Storable::retrieve($store)
      or croak("Can't restore wormbase from $store\n");
} else {
    $wormbase = Wormbase->new(
      -debug    => $debug,
      -organism => $species,
      -version  => $WS_version,
    );
}

$WS_version = $wormbase->get_wormbase_version_name;
my $database = "worm_ensembl_${species}";

# more setup
my $log = Log_files->make_build_log($wormbase);
$log->write_to("Updating BLAT input files for $database\n");

$host||=$ENV{'WORM_DBHOST'};
$port||=$ENV{'WORM_DBPORT'};

# MYSQL setup
my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
        -host     => $host,
        -user     => $user,
        -dbname   => $dbname || $database,
        -port     => $port,
    );
$db->password($password);

foreach my $ana (&get_all_blat_analyses()) {
  my $logic_name = $ana->logic_name();

  &update_analysis_entry($ana);
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
sub update_analysis_entry {
  my ($ana) = @_;

  # set species-specific genome fasta file to search against
  # set parameters depending on self-hits or non-self hits

  my $genome_fa = $wormbase->genome_seq;
  if (not -e $genome_fa) {
    die "Could not find genome fasta file $genome_fa\n";
  }
  $ana->db_file($genome_fa);
  $ana->db_version($WS_version);

  my $parameters = "";
  if (lc($ana->db) eq $species) {
    if ($ana->logic_name =~ /mrna$/) {
      $parameters = "-fine";
    }
  } else {
    $parameters = "-q=dnax -t=dnax";
  }

  $ana->parameters($parameters);
  $ana->adaptor->update($ana);
}


#######################
sub make_input_ids {
  my ($ana, $dep_ana) = @_;

  my $sth=$db->prepare(
    'INSERT INTO input_id_analysis (input_id,input_id_type,analysis_id,created,result) VALUES (?,?,?,NOW(),0)'
      );

  my $file_base = $wormbase->basedir . "/cDNA/";

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
  foreach my $logic_name (keys %{$CONFIG->{logic_names}}) {
    my $ana = $db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
    if (defined $ana) {
      push @blat_ana, $ana;
    } 
  }

  return @blat_ana;
}


