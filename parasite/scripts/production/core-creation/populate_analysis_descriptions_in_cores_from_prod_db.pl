#!/usr/bin/env perl

=pod

=head1 NAME

  populate_analysis_descriptions_in_cores_from_prod_db.pl

=head1 SYNOPSIS

  Populate a core database with analysis descriptions from a production database

=head1 ARGUMENTS

  perl populate_analysis_descriptions_in_cores_from_prod_db.pl
         -prod_host
         -prod_port
         -prod_user
         -prod_db
         -host
         -port
         -user
	 -pass
	 -db

=cut


use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Production::Utils::ProductionDbUpdater;

my ($prod_host, $prod_port, $prod_user, $prod_db, $host, $port, $user, $pass, $db);

GetOptions(
  'prod_host=s'    => \$prod_host,
  'prod_port=s'    => \$prod_port,
  'prod_user=s'    => \$prod_user,
  'prod_db=s'      => \$prod_db,
  'host=s'         => \$host,
  'port=s'         => \$port,
  'user=s'         => \$user,
  'pass=s'         => \$pass, 
  'db=s'           => \$db
) || die("bad commandline parameter\n");


my $prod_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user    => $prod_user,	
    -dbname  => $prod_db,
    -host    => $prod_host,
    -port    => $prod_port 
);

my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host   => $host,
    -user   => $user,
    -dbname => $db,
    -pass   => $pass,
    -port   => $port,
);

my $updater = Bio::EnsEMBL::Production::Utils::ProductionDbUpdater->new(
    -PRODUCTION_DBA => $prod_dba
);

print STDERR "Updating analysis descriptions for $db\n";
$updater->update_analysis_description($dba->dbc);
