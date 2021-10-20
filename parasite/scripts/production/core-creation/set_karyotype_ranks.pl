#!/usr/bin/env perl

use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my ($dbhost, $dbport, $dbuser, $dbpass, $dbname, $karyotypes_file);

GetOptions(
  'host=s'            => \$dbhost,
  'port=s'            => \$dbport,
  'user=s'            => \$dbuser,
  'pass=s'            => \$dbpass,
  'dbname=s'          => \$dbname,
  'karyotypes_file=s' => \$karyotypes_file
) or die "Couldn't get options" ;


## connect to db


my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-dbname  => $dbname,
	-host    => $dbhost,
	-port    => $dbport,
	-user    => $dbuser,
	-pass    => $dbpass,
);

## read file

open(KFILE, "<", $karyotypes_file) or die "Can not read file: $karyotypes_file";

my %ranks;

while (<KFILE>){
	# name	karyotype_rank
	# LG1	1

	chomp;
	my @fields = split(/\t/, $_);
	die "expected 2 columns: name	karyotype_rank" if scalar(@fields) != 2;

	$ranks{$fields[0]} = $fields[1];
}

## prepare SQL

my $check_for_karyotypes = $dba->dbc->prepare("SELECT COUNT(*) FROM 
					seq_region_attrib LEFT JOIN attrib_type
					ON seq_region_attrib.attrib_type_id = attrib_type.attrib_type_id
					WHERE attrib_type.code = 'karyotype_rank' ");

my $fetch_seq_region_id = $dba->dbc->prepare("SELECT seq_region_id FROM
					seq_region WHERE
					name = ?");


my $fetch_attrib_type_id = $dba->dbc->prepare("SELECT attrib_type_id FROM
					attrib_type WHERE
					code = 'karyotype_rank' ");

my $insert_karyotype_rank = $dba->dbc->prepare("INSERT INTO 
					seq_region_attrib (seq_region_id, attrib_type_id, value)
					VALUES (?,?,?)");


## do the business

# check we don't already have any karyotype_ranks for this genome

$check_for_karyotypes->execute();
my $check = $check_for_karyotypes->fetchrow();
die "Genome already has a karyotype" unless $check == 0;
$check_for_karyotypes->finish();

# fetch attrib_type_id

$fetch_attrib_type_id->execute();
my $attrib_type_id = $fetch_attrib_type_id->fetchrow();
$fetch_attrib_type_id->finish();

# insert rank values

foreach my $name (keys %ranks){

    # fetch seq_region_id

    $fetch_seq_region_id->execute($name);
    my $seq_region_id = $fetch_seq_region_id->fetchrow();
    $fetch_seq_region_id->finish();

    # insert
    $insert_karyotype_rank->execute($seq_region_id, $attrib_type_id, $ranks{$name});
    print "added karyotype_rank $ranks{$name} for $name\n";
    $insert_karyotype_rank->finish();
}


