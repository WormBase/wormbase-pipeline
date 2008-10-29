#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  p2msql.pl
#
#        USAGE:  ./p2msql.pl 
#
#  DESCRIPTION:  filter to change PostgreSQL to mySQL dumps for the clustal database
#                works best with "pg_dump -d -x -O -t clustal" dumps.
#      CREATED:  10/29/08 10:21:01 GMT by mh6@sanger.ac.uk
#        USAGE:  perl p2msql.pl < postgresdump.sql > mysqldump.sql
#===============================================================================

use strict;

my $header =1 ;
while (<>){

	# skip the header
	$header = 0 if /^CREATE/;
	next if $header;

	# schema
	s/peptide_id character varying(25)/peptide_id VARCHAR(25)/; # varchar bit

	# index
        next if /ALTER TABLE ONLY clustal/;
        if (/ADD CONSTRAINT/){
		print  "ALTER TABLE clustal ADD UNIQUE (peptide_id);";
	}
	else {
         print $_;
        }
}

# the single inserts are needed because mysql can't bulk load from stdin
# the pg_dump also uses some iffy syntax in the schema instead of SQL :-(
