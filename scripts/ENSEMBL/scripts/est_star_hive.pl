#!/software/bin/perl -w
#
# est_star_hive.pl
# 
# by Gary Williams            
#

# This sets up the environement for the EST STAR hive aligner pipeline
# for each of the logic_names ('embl_star', 'washu_star',
# 'nembase_star') and runs the beekeeper scripts.

#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2014-06-13 13:01:38 $      

use strict;                                      
use Getopt::Long;
use Carp;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


######################################
# variables and command-line options # 
######################################

my ($help, $SPECIES);
my (@logic_names, %logic_names);

GetOptions ("help"           => \$help,
	    "species:s"      => \$SPECIES,
	    );

# Display help if required
&usage("Help") if ($help);

if (!defined $SPECIES) {die "-species is not set\n"}

##########################
# MAIN BODY OF SCRIPT
##########################

# set up the locations of things
my $project_dir = "/nfs/panda/ensemblgenomes/wormbase/BUILD_DATA/cDNA/$SPECIES";
my $WORM_PACKAGES = "/nfs/panda/ensemblgenomes/wormbase/software/packages";
my $PIPELINE_CODE_DIR = "/nfs/panda/ensemblgenomes/wormbase/software/packages/ensembl_genomes";
my $PIPELINE_DATA_DIR = "/nfs/nobackup/ensemblgenomes/wormbase/BUILD/EST-data/$SPECIES";
my $CORE_DB_NAME = "worm_ensembl_$SPECIES";
my $DB_HOST = "mysql-wormbase-pipelines";
my $DB_PORT = 4331;
my $DB_USER = "wormadmin";
my $DB_PASS = "worms";
my $HIVE_DB_HOST = "mysql-wormbase-pipelines";
my $HIVE_DB_PORT = 4331;
my $HIVE_DB_USER = "wormadmin";
my $HIVE_DB_PASS = "worms";
my $FORCE = 0;
if (! -d $PIPELINE_DATA_DIR) {
  mkdir mkdir $PIPELINE_DATA_DIR, 0777;
} else{
  # tidy up any existing files
  system("rm -rf $PIPELINE_DATA_DIR/*");
}
chdir $PIPELINE_DATA_DIR;

my $HIVEDB = "$ENV{USER}_EST_Pipeline_${SPECIES}";
my $HIVE_URL = "mysql://${HIVE_DB_USER}:${HIVE_DB_PASS}\@${HIVE_DB_HOST}:${HIVE_DB_PORT}/${HIVEDB}";

$ENV{'PATH'} = "/nfs/panda/ensemblgenomes/external/STAR:$ENV{PATH}";
$ENV{'PERL5LIB'} = "/nfs/panda/ensemblgenomes/wormbase/software/packages/ensembl_genomes/eg-estalignfeature/lib/:$ENV{'PERL5LIB'}";
$ENV{'PERL5LIB'} = "/nfs/panda/ensemblgenomes/wormbase/software/packages/ensembl_genomes/eg-pipelines/modules/:$ENV{'PERL5LIB'}";
$ENV{'PERL5LIB'} = "/nfs/panda/ensemblgenomes/wormbase/software/packages/bioperl/bioperl-live:$ENV{'PERL5LIB'}";
$ENV{'PERL5LIB'} = "/nfs/panda/ensemblgenomes/wormbase/software/packages/bioperl/bioperl-run/lib:$ENV{'PERL5LIB'}";
$ENV{'PERL5LIB'} = "/nfs/panda/ensemblgenomes/wormbase/software/packages/ensembl-compara/modules:$ENV{'PERL5LIB'}";
$ENV{'PERL5LIB'} = "/nfs/panda/ensemblgenomes/wormbase/software/packages/ensembl-variation/modules:$ENV{'PERL5LIB'}";
$ENV{'PERL5LIB'} = "/nfs/panda/ensemblgenomes/wormbase/software/packages/ensembl-pipeline/modules:$ENV{'PERL5LIB'}";
$ENV{'PERL5LIB'} = "/nfs/panda/ensemblgenomes/wormbase/software/packages/ensembl-pipeline/scripts:$ENV{'PERL5LIB'}";
$ENV{'PERL5LIB'} = "/nfs/panda/ensemblgenomes/wormbase/software/packages/ensembl-analysis/modules:$ENV{'PERL5LIB'}";
$ENV{'PERL5LIB'} = "/nfs/panda/ensemblgenomes/wormbase/software/packages/ensembl-production/modules:$ENV{'PERL5LIB'}";
$ENV{'PERL5LIB'} = "/nfs/panda/ensemblgenomes/wormbase/software/packages/ensembl/modules:$ENV{'PERL5LIB'}";
$ENV{'PERL5LIB'} = "/nfs/panda/ensemblgenomes/wormbase/software/lib:$ENV{'PERL5LIB'}";
$ENV{'PERL5LIB'} = "$ENV{CVS_DIR}:$ENV{'PERL5LIB'}";
$ENV{'PERL5LIB'} = "$ENV{CVS_DIR}/ENSEMBL/lib:$ENV{'PERL5LIB'}";

$ENV{'ENSEMBL_REGISTRY'} = "$ENV{CVS_DIR}/ENSEMBL/etc/Registry.pm";

# Get the Registry
# This assumes the environment variable ENSEMBL_REGISTRY is set, this is
#     used as the name of the configuration file to read.
# Or, the file .ensembl_init exists in the home directory, it is
#     used as the configuration file.
Bio::EnsEMBL::Registry->load_all();

# get the maximum intron size (Arnaud may incorporate this into the hive script at some time)
my $MAX_INTRON_LENGTH = get_max_intron_size($SPECIES);



foreach my $LOGIC_NAME ('embl_star', 'washu_star', 'nembase_star'){

  print "\nRunning STAR hive pipeline for ${LOGIC_NAME}.\n\n";

# set up the EST_FILE based on the logic_name
  my %est_files = (
		   'embl_star' => 'EST',
		   'washu_star' => 'Nematode.net',
		   'nembase_star' => 'Nembase',
		  );
  my $EST_FILE = "${project_dir}/$est_files{$LOGIC_NAME}";
  
  if (! -e $EST_FILE) {
    print "ERROR: The $LOGIC_NAME sequence file ${EST_FILE} has not been set up yet.\n";
    exit 1;
  }

  if (-z $EST_FILE) {
    print "The $LOGIC_NAME sequence file ${EST_FILE} has no sequences in it. Skipping it ...\n";
    next;
  }

# Delete the old hive database and any old results
  system("mysql --host=${DB_HOST} --port=${DB_PORT} --user=wormadmin --password=worms -e \"DROP DATABASE IF EXISTS ${HIVEDB}\"") == 0 or die "mysql drop database failed: $?";
  system("mysql --host=${DB_HOST} --port=${DB_PORT} --user=wormadmin --password=worms -e \"USE ${CORE_DB_NAME}; DELETE FROM dna_align_feature WHERE analysis_id IN (SELECT analysis_id FROM analysis WHERE logic_name = '${LOGIC_NAME}')\"") == 0 or die "mysql delete old data failed: $?";

  system("perl ${PIPELINE_CODE_DIR}/ensembl-hive/scripts/init_pipeline.pl Bio::EnsEMBL::Hive::PipeConfig::EST_STAR_pipeline_conf -registry $ENV{ENSEMBL_REGISTRY} -species ${SPECIES} -ensembl_cvs_root_dir ${PIPELINE_CODE_DIR} -hive_password ${HIVE_DB_PASS} -hive_port ${HIVE_DB_PORT} -hive_host ${HIVE_DB_HOST} -hive_user ${HIVE_DB_USER} -hive_dbname ${HIVEDB} -pipeline_dir ${PIPELINE_DATA_DIR} -core_db_name ${CORE_DB_NAME} -est_file ${EST_FILE} -force 0 -logic_name ${LOGIC_NAME} -max_intron_length ${MAX_INTRON_LENGTH} -load_core 1") == 0 or die "$LOGIC_NAME init_pipeline failed: $?";

  system("${PIPELINE_CODE_DIR}/ensembl-hive/scripts/beekeeper.pl -reg_conf $ENV{ENSEMBL_REGISTRY} -url ${HIVE_URL} -sync") == 0 or die "system $LOGIC_NAME beekeeper -sync failed: $?";

  system("${PIPELINE_CODE_DIR}/ensembl-hive/scripts/beekeeper.pl -reg_conf $ENV{ENSEMBL_REGISTRY} -url ${HIVE_URL} -submit_workers_max 100 -sleep 2 -life_span 2000 -loop") == 0 or die "$LOGIC_NAME beekeeper -loop failed: $?";

    print "\nSTAR hive pipeline for $LOGIC_NAME has finished.\n\n";


}


print "Script has finished.\n";
exit(0);






##############################################################
#
# Subroutines
#
##############################################################



##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################
sub get_max_intron_size {

  my ($species) = @_;
  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor( "$species", 'Core', 'gene' );
  
  if (! defined $gene_adaptor) {
    die "can't get a gene adaptor for species, $species\n";
  }
  
  # Get all transcripts through all genes
  print "Finding maximum intron size for species $species\n";
  my $MAX_INTRON_SIZE = 0;
  my $genes_aref = $gene_adaptor->fetch_all_by_biotype('protein_coding');
  print STDERR "processing " . @$genes_aref . " genes\n";
  foreach my $gene (@$genes_aref) {
    my $transcripts_aref = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts_aref) {
      my $introns_aref = $transcript->get_all_Introns();
      foreach my $intron (@$introns_aref) {
	if ($intron->length() > $MAX_INTRON_SIZE) {
	  $MAX_INTRON_SIZE = $intron->length();
	}
      }
    }
  }
  print "Maximum intron size: $MAX_INTRON_SIZE\n";
  return $MAX_INTRON_SIZE;
}



# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - script_template.pl

=head1 USAGE

=over 4

=item script_template.pl  [-options]

=back

This script does...blah blah blah

script_template.pl MANDATORY arguments:

=over 4

=item None at present.

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item xxx (xxx@sanger.ac.uk)

=back

=cut
