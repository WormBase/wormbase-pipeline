#!/usr/bin/env perl
#
# DESCRIPTION:
#   setting up the BLAT pipeline
#
# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2009-06-01 10:24:58 $

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

###################
use lib '/software/worm/ensembl/ensembl/modules/';
use lib '/software/worm/lib/bioperl-live/';
###################

use Getopt::Long;
use Storable;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

# some magic to turn off the deprecation warnings
use Bio::EnsEMBL::Utils::Exception;
verbose('OFF');

use Wormbase;
use strict;

my($debug,$store,$database,$port,$user,$password,$species,$host);
GetOptions(
 'debug=s'    => \$debug,
 'store=s'    => \$store,
 'user=s'     => \$user,
 'password=s' => \$password,
 'species=s'  => \$species,
 'host=s'     => \$host,
 'port=s'     => \$port,
)||die(USAGE);

# WormBase setup
my $wormbase;
 if ($store) {
    $wormbase = Storable::retrieve($store)
      or croak("Can't restore wormbase from $store\n");
} else {
    $wormbase = Wormbase->new(
        -debug    => $debug,
        -species  => $species,
    );
}

my $database = sprintf('worm_ensembl_%s',lc(ref $wormbase));

# more setup
my $log = Log_files->make_build_log($wormbase);
$log->write_to("Updating BLAT input files for $database\n");

$host||='ia64d';
$port||=3306;

# MYSQL setup
my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host     => $host,
        -user     => $user,
        -dbname   => $database,
        -port     => $port,
        -driver   => 'mysql', 
    );
$db->password($password);

# a.) clean out all dna_align_features
$db->dbc->do('DELETE FROM dna_align_feature');

# b.) clean out all input_ids for the blat analysis
&clean_input_ids();

# c.) make new input_ids
&make_input_ids();

$log->write_to("Finished.\n");
$log->mail();

# slight misnomer, first cleans the input_id_analysis table and puts new input_ids back in.
sub make_input_ids {
    # locations for each fasta_file
    my %est_files = (
        'blat_brugia_ests'    => '/nfs/disk100/wormpub/BUILD/cDNA/brugia/EST.masked*',
        'blat_brugia_cdnas'   => '/nfs/disk100/wormpub/BUILD/cDNA/brugia/mRNA.masked*',

        'blat_elegans_ests'   => '/nfs/disk100/wormpub/BUILD/cDNA/elegans/EST.masked*',
        'blat_elegans_cdnas'  => '/nfs/disk100/wormpub/BUILD/cDNA/elegans/mRNA.masked*',
        'blat_elegans_osts'   => '/nfs/disk100/wormpub/BUILD/cDNA/elegans/OST.masked*',
        'blat_elegans_rsts'   => '/nfs/disk100/wormpub/BUILD/cDNA/elegans/RST.masked*',
        'blat_elegans_ncrnas' => '/nfs/disk100/wormpub/BUILD/cDNA/elegans/ncRNA.masked*',
        'blat_elegans_tc1s'   => '/nfs/disk100/wormpub/BUILD/cDNA/elegans/tc1.masked*',

        'blat_briggsae_ests'  => '/nfs/disk100/wormpub/BUILD/cDNA/briggsae/EST.masked*',
        'blat_briggsae_cdnas' => '/nfs/disk100/wormpub/BUILD/cDNA/briggsae/mRNA.masked*',
        
        'blat_brenneri_ests'  => '/nfs/disk100/wormpub/BUILD/cDNA/brenneri/EST.masked*',
        'blat_brenneri_cdnas' => '/nfs/disk100/wormpub/BUILD/cDNA/brenneri/mRNA.masked*',

        'blat_remanei_ests'   => '/nfs/disk100/wormpub/BUILD/cDNA/remanei/EST.masked*',
        'blat_remanei_cdnas'  => '/nfs/disk100/wormpub/BUILD/cDNA/remanei/mRNA.masked*',
        
        'blat_japonica_ests'  => '/nfs/disk100/wormpub/BUILD/cDNA/japonica/EST.masked*',
        'blat_japonica_cdnas' => '/nfs/disk100/wormpub/BUILD/cDNA/japonica/mRNA.masked*',
        
        'blat_heterorhabditis_ests'  => '/nfs/disk100/wormpub/BUILD/cDNA/heterorhabditis/EST.masked*',
        'blat_heterorhabditis_cdnas' => '/nfs/disk100/wormpub/BUILD/cDNA/heterorhabditis/mRNA.masked*',

        'blat_pristionchus_ests'  => '/nfs/disk100/wormpub/BUILD/cDNA/pristionchus/EST.masked*',
        'blat_pristionchus_cdnas' => '/nfs/disk100/wormpub/BUILD/cDNA/pristionchus/mRNA.masked*',
        
        'blat_nembase_ests'   => '/nfs/disk100/wormpub/BUILD/cDNA/nembase/EST.masked*',
        'blat_washu_ests'     => '/nfs/disk100/wormpub/BUILD/cDNA/washu/EST.masked*',
        'blat_nematode_ests'  => '/nfs/disk100/wormpub/BUILD/cDNA/nematode/EST.masked*',
        ); 

    # SQL handles
    my $sth=$db->prepare(
        'INSERT INTO input_id_analysis (input_id,input_id_type,analysis_id,created,result) VALUES (?,?,?,NOW(),0)'
        );
    my $del_handle=$db->prepare('DELETE FROM input_id_analysis WHERE analysis_id=?');

    # create_input_ids for each analysis
    foreach my $analysis(&get_all_blat_analysis()) {
        
        $del_handle->execute($analysis->dbID);
        $del_handle->execute($analysis->input_id_type_analysis->dbID); 

        if ($est_files{$analysis->logic_name}) {
            my @files = glob($est_files{$analysis->logic_name});
            foreach my $file (@files){
                $sth->execute($file,$analysis->input_id_type_analysis->logic_name,$analysis->input_id_type_analysis->dbID);
            }
        }else {
            printf STDERR "ERROR: cannot find files for %s\n",$analysis->logic_name;
        }
    }
}

# helper functions 
sub clean_input_ids {
    my @analysis_to_clean = &get_all_blat_analysis();
    foreach my $analysis (@analysis_to_clean){
        my $dbID=$analysis->dbID;
        $db->dbc->do("DELETE FROM input_id_analysis WHERE analysis_id = $dbID");
    }    
}

# get all analysis objects for blat
sub get_all_blat_analysis {
    my $adb = $db->get_AnalysisAdaptor();
    return $adb->fetch_by_program_name('/software/worm/bin/blat/blat');
}


# extensions to EnsEMBL classes
package Bio::EnsEMBL::Analysis;

sub input_id_type_analysis {
    my ($self)=@_;
    my $sth = $self->adaptor->prepare('SELECT a.analysis_id FROM input_id_type_analysis i JOIN analysis a on(i.input_id_type=a.logic_name) WHERE i.analysis_id = ?');
    $sth->execute($self->dbID);
    my $id_type;
    while (my ($tmp) = $sth->fetchrow_array){$id_type=$tmp}
    
    return $self->adaptor->fetch_by_dbID($id_type);
}

package Bio::EnsEMBL::DBSQL::AnalysisAdaptor;

# to fetch all analysis with a certain program_file
sub fetch_by_program_name {
    my ($self,$program_name) = @_;
    my @returned_analysis;
    my $sth = $self->prepare('SELECT analysis_id FROM analysis WHERE program_file = ?');
    $sth->execute($program_name);
    my $ids = $sth->fetchall_arrayref()||[]; # array_ref->array_refs
    foreach my $analysis(@$ids){
        push @returned_analysis , $self->fetch_by_dbID($$analysis[0]);
    }
    return @returned_analysis;
}

