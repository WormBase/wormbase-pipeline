#!/usr/local/ensembl/bin/perl  -w

use WormBase;
use WormBaseConf;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Wormy;


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $WB_DBHOST,
    -user   => $WB_DBUSER,
    -dbname => $WB_DBNAME,
    -pass   => $WB_DBPASS,
    -port   => $WB_DBPORT,
);

# adding assembly type to meta table

my $analysis_adaptor = $db->get_AnalysisAdaptor();

my $analysis  = $analysis_adaptor->fetch_by_logic_name(Uniprot);
#my $analysis  = $analysis_adaptor->fetch_by_logic_name(WormPep);

foreach my $chromosome_info ( @{$WB_CHR_INFO} ) {
    last if $chromosome_info->{chr_name} eq 'MtDNA';
    print "handling " . $chromosome_info->{'chr_name'} . " with files " . $chromosome_info->{'agp_file'} . " and " . $chromosome_info->{'gff_file'} . "\n"
      if ($WB_DEBUG);

    my $chr = $db->get_SliceAdaptor->fetch_by_region( 'Chromosome', $chromosome_info->{'chr_name'},
        1, ( $chromosome_info->{'length'}, 1, $WB_NEW_COORD_SYSTEM_VERSION ) );

     my $blast = Wormy->new($chr,'/nfs/disk100/wormpub/DATABASES/WS158/CHROMOSOMES/'.$chromosome_info->{gff_file},$analysis);
     $blast->save($db);
     print "have " . scalar @{$blast->{hits}} . " blastxfeatures.\n"; #if ($WB_DEBUG);

}

