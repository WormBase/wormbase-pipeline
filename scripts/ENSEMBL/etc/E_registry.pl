use strict;
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ia64c',-user=>'wormro',-port=> 3306, -species => 'Caenorhabditis elegans',-group => 'core',-dbname => 'worm_ensembl_elegans');
new Bio::EnsEMBL::DBSQL::DBAdaptor(-host=>'ia64c',-user=>'wormro',-port=>'3306',-species=>'Caenorhabditis brenneri',-group => 'core',-dbname => 'worm_ensembl_pb2801');
new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ia64c',-user=>'wormro',-port=> 3306, -species => 'Caenorhabditis briggsae',-group => 'core',-dbname => 'worm_ensembl_briggsae');
new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ia64c',-user=>'wormro',-port=> 3306, -species => 'Caenorhabditis remanei',-group => 'core',-dbname => 'worm_ensembl_remanei');

new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(-host => 'ia64c',-user=>'wormro',-port=> 3306, -species => 'Compara',-dbname => 'worm_compara_lagan');

1;
