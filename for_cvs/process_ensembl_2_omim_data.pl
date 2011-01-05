#!/usr/bin/perl
# from Norie de la Cruz
# Version 2nd February 2009

# usage: process_ensembl_2_omim_data.pl
#    process file and generates a list of ensemble_id to omim_id relationships

# needed changes:
#    change aceserver port/name

use strict;
use Ace;

my $DB = Ace->connect(
    -host => 'localhost',
    -port => 23100
);

my $it = $DB->fetch_many( -query => 'find Protein Species="Homo sapiens"' );

open OUTFILE, ">./ensembl_id2omim_id.txt" or die "Cannot open outfile\n";

while ( my $object = $it->next ) {

    print STDERR "processing $object (${\$object->Species})\n" if $ENV{DEBUG};
    next unless $object->Species eq 'Homo sapiens';
    my $db_info = $object->DB_info;
    unless ($db_info) {
        print STDERR "ERROR: cannot find DB_info for $object\n";
        next;
    }

    my @data = $db_info->col;

    foreach my $db_data (@data) {

        if ( $db_data =~ m/OMIM/ ) {

            my @db_data = $db_data->col;
            foreach my $omim_data (@db_data) {
                if ( $omim_data =~ m/disease/ ) {
                    my $disease_id = $omim_data->right;
                    my ( $ensembl, $ensembl_id ) = split /:/, "$object";
                    print OUTFILE "$ensembl_id\=\>";    #OUTFILE
                    print OUTFILE "$disease_id\n";      #OUTFILE

                }

            }

        } ## end if ($db_data =~ m/OMIM/)

    }    # end foreach my $db_data (@data)

}

print "ok\n";
