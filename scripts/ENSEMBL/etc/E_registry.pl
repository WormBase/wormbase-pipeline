use strict;
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(-host => 'farmdb1',
                                            -user => 'wormro',
                                            -port => 3306,
                                            -dbname => 'worm_compara',
                                            -group => 'compara',
                                            -species => 'Compara trees'
                                            );

new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(-host => 'farmdb1',
                                            -user => 'wormro',
                                            -port => 3306,
                                            -dbname => 'worm_compara_pecan',
                                            -group => 'compara',
                                            -species => 'Compara alignments'
                                            );

new Bio::EnsEMBL::DBSQL::DBAdaptor(-host    => 'farmdb1',
                                   -user    =>'wormro',
                                   -port    => 3306, 
                                   -species => 'Caenorhabditis elegans',
                                   -group   => 'core',
                                   -dbname  => 'worm_ensembl_elegans');

new Bio::EnsEMBL::DBSQL::DBAdaptor(-host    => 'farmdb1',
                                   -user    =>'wormro',
                                   -port    => 3306, 
                                   -species => 'Caenorhabditis briggsae',
                                   -group   => 'core',
                                   -dbname  => 'worm_ensembl_briggsae');


new Bio::EnsEMBL::DBSQL::DBAdaptor(-host    => 'farmdb1',
                                   -user    =>'wormro',
                                   -port    => 3306, 
                                   -species => 'Caenorhabditis remanei',
                                   -group   => 'core',
                                   -dbname  => 'worm_ensembl_remanei');


new Bio::EnsEMBL::DBSQL::DBAdaptor(-host    => 'farmdb1',
                                   -user    =>'wormro',
                                   -port    => 3306, 
                                   -species => 'Caenorhabditis brenneri',
                                   -group   => 'core',
                                   -dbname  => 'worm_ensembl_brenneri');


new Bio::EnsEMBL::DBSQL::DBAdaptor(-host    => 'farmdb1',
                                   -user    =>'wormro',
                                   -port    => 3306, 
                                   -species => 'Caenorhabditis japonica',
                                   -group   => 'core',
                                   -dbname  => 'worm_ensembl_japonica');


new Bio::EnsEMBL::DBSQL::DBAdaptor(-host    => 'farmdb1',
                                   -user    =>'wormro',
                                   -port    => 3306, 
                                   -species => 'Caenorhabditis angaria',
                                   -group   => 'core',
                                   -dbname  => 'worm_ensembl_cangaria');


new Bio::EnsEMBL::DBSQL::DBAdaptor(-host    => 'farmdb1',
                                   -user    =>'wormro',
                                   -port    => 3306, 
                                   -species => 'Caenorhabditis sp.7',
                                   -group   => 'core',
                                   -dbname  => 'worm_ensembl_csp7');


new Bio::EnsEMBL::DBSQL::DBAdaptor(-host    => 'farmdb1',
                                   -user    =>'wormro',
                                   -port    => 3306, 
                                   -species => 'Caenorhabditis sp.9',
                                   -group   => 'core',
                                   -dbname  => 'worm_ensembl_csp9');


new Bio::EnsEMBL::DBSQL::DBAdaptor(-host    => 'farmdb1',
                                   -user    =>'wormro',
                                   -port    => 3306, 
                                   -species => 'Caenorhabditis sp.11',
                                   -group   => 'core',
                                   -dbname  => 'worm_ensembl_csp11');


new Bio::EnsEMBL::DBSQL::DBAdaptor(-host    => 'farmdb1',
                                   -user    =>'wormro',
                                   -port    => 3306, 
                                   -species => 'Haemonchus contortus',
                                   -group   => 'core',
                                   -dbname  => 'worm_ensembl_hcontortus');


new Bio::EnsEMBL::DBSQL::DBAdaptor(-host    => 'farmdb1',
                                   -user    =>'wormro',
                                   -port    => 3306, 
                                   -species => 'Pristionchus pacificus',
                                   -group   => 'core',
                                   -dbname  => 'worm_ensembl_pristionchus');


new Bio::EnsEMBL::DBSQL::DBAdaptor(-host    => 'farmdb1',
                                   -user    =>'wormro',
                                   -port    => 3306, 
                                   -species => 'Strongyloides ratti',
                                   -group   => 'core',
                                   -dbname  => 'worm_ensembl_sratti');


new Bio::EnsEMBL::DBSQL::DBAdaptor(-host    => 'farmdb1',
                                   -user    =>'wormro',
                                   -port    => 3306, 
                                   -species => 'Meloidogyne hapla',
                                   -group   => 'core',
                                   -dbname  => 'worm_ensembl_mhapla');


new Bio::EnsEMBL::DBSQL::DBAdaptor(-host    => 'farmdb1',
                                   -user    =>'wormro',
                                   -port    => 3306, 
                                   -species => 'Meloidogyne incognita',
                                   -group   => 'core',
                                   -dbname  => 'worm_ensembl_mincognita');


new Bio::EnsEMBL::DBSQL::DBAdaptor(-host    => 'farmdb1',
                                   -user    =>'wormro',
                                   -port    => 3306, 
                                   -species => 'Brugia malayi',
                                   -group   => 'core',
                                   -dbname  => 'worm_ensembl_brugia');


new Bio::EnsEMBL::DBSQL::DBAdaptor(-host    => 'farmdb1',
                                   -user    =>'wormro',
                                   -port    => 3306, 
                                   -species => 'Trichinella spiralis',
                                   -group   => 'core',
                                   -dbname  => 'worm_ensembl_trspiralis');

######################

Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Compara trees', -alias => ['compara_trees']);
Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Compara alignments', -alias => ['compara_alignments']);

Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Caenorhabditis elegans', -alias => ['elegans']);
Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Caenorhabditis briggsae', -alias => ['briggsae']);
Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Caenorhabditis remanei', -alias => ['remanei']);
Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Caenorhabditis brenneri', -alias => ['brenneri']);
Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Caenorhabditis japonica', -alias => ['japonica']);
Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Caenorhabditis angaria', -alias => ['cangaria']);
Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Caenorhabditis sp.7', -alias => ['csp7']);
Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Caenorhabditis sp.9', -alias => ['csp9']);
Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Caenorhabditis sp.11', -alias => ['csp11']);
Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Haemonchus contortus', -alias => ['hcontortus']);
Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Pristionchus pacificus', -alias => ['pristionchus']);
Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Strongyloides ratti', -alias => ['sratti']);
Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Meloidogyne hapla', -alias => ['mhapla']);
Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Meloidogyne incognita', -alias => ['mincognita']);
Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Brugia malayi', -alias => ['brugia']);
Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(-species => 'Trichinella spiralis', -alias => ['tspiralis']);


1;
