#!/software/bin/perl
#
# parse the semi-fake nembase file

use lib '/software/worm/ensembl/ensembl/modules';
use lib '/software/worm/ensembl/ensembl-pipeline/modules';
use lib '/software/worm/lib/bioperl-live';

use YAML;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(verbose warning);
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Pipeline::Tools::ExonUtils;
use DBI qw(:sql_types);

verbose('OFF');

my $yfile="$ENV{CVS_DIR}/ENSEMBL/etc/ensembl_lite.conf";
my $config = ( YAML::LoadFile($yfile) )->{nembase};

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
      -host   => $config->{database}->{host},    # ia64d
      -user   => $config->{database}->{user},    # wormadmin
      -dbname => $config->{database}->{dbname},  # worm_nembase
      -pass   => $config->{database}->{password},# whateva
      -port   => $config->{database}->{port},    # 3306
  );

my $analysis = $db->get_AnalysisAdaptor()->fetch_by_logic_name('wormbase');
my $gene_adaptor = $db->get_GeneAdaptor;

# ACP00005_1      curated coding_exon     1       324     .       +       0       CDS "ACP00005_1"
#  ^^^^
#  seq_region                             start   stop            +       phase 1  gene/cds_name

while (<>){
    
    s/\"//g; # get rid of the flanking quotes
    my ($seq_region,$skip1,$skip2,$start,$stop,$skip3,$strand,$phase,$skip,$cds) = split;
    
    # get the seq_region slice
    my $slice = $db->get_SliceAdaptor->fetch_by_region( 'chromosome', $seq_region );
    
    ############################
    # make an exon

    my $exon = &_make_exon($start,$stop,$slice,$cds);
  
    ############################
    # make transcripts / translations
    
    my $translation = &_make_translation($exon,$cds);
    my $transcript  = &_make_transcript($exon,$translation);
    
    ##############################
    # make gene
    
    my $gene = &_make_gene($cds,$transcript);
        
    # save gene
    $gene_adaptor->save($gene);
}

$db->dbc->do('UPDATE gene SET biotype="protein_coding"');


#############################
# make an exon from GFF bits and bobs
#  le haque:
#      * we know that there is only one exon per gene on + starting at 1 in phase 1

sub _make_exon {
    my ($start,$stop,$seq_region,$cds)=@_;
    my $_exon  = new Bio::EnsEMBL::Exon;
    $_exon->start($start);
    $_exon->end($stop);
    $_exon->analysis($analysis);
    $_exon->slice($seq_region);
    $_exon->phase(1);
    $_exon->strand(1);
    $_exon->version(1);
    $_exon->stable_id($cds.'1');
    my $end_phase = ( $phase + ( $_exon->end - $_exon->start ) + 1 ) % 3;
    $_exon->end_phase($end_phase);

    return $_exon;
}

###########################
# make an translation from a exon and a cds_id

sub _make_translation {
    my($exon,$cds)   = @_;
    my $_translation = new Bio::EnsEMBL::Translation;
    $_translation->start_Exon($exon);
    $_translation->end_Exon($exon);
    $_translation->start(1);
    $_translation->end($exon->end);
    $_translation->stable_id($cds);
    $_translation->version(1);

    return $_translation;
}

#############################
# make a transcript from the translation and the exon

sub _make_transcript {
    my ($exon,$translation)=@_;
    my $_transcript = new Bio::EnsEMBL::Transcript;
    $_transcript->add_Exon($exon);    
    $_transcript->translation($translation);
    $_transcript->version(1);
    $_transcript->stable_id($translation->stable_id);

    return $_transcript;
}

##################################
# make a gene from the transcript

sub _make_gene {
    my ($cds,$transcript)=@_;
    my $_gene = new Bio::EnsEMBL::Gene;
    $_gene->biotype( $analysis->logic_name ); # maybe not needed
    $_gene->version(1);
    $_gene->stable_id($cds);
    $_gene->add_Transcript($transcript);
    $_gene->stable_id($cds);
    
    return $_gene;
}
