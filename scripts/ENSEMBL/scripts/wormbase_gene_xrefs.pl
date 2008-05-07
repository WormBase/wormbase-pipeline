#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  wormbase_gene_xrefs.pl
#
#        USAGE:  ./wormbase_gene_xrefs.pl
#
#  DESCRIPTION:
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:   (), <>
#      COMPANY:
#      VERSION:  1.0
#      CREATED:  21/11/06 16:50:14 GMT
#     REVISION:  ---
#===============================================================================

use strict;
use Ace;
use Getopt::Long;

my ( $database, $insert, $update ,$fix);
GetOptions(
    'database=s'   => \$database,
    'insert_xref'  => \$insert,
    'update_xrefs' => \$update,
    'fix_xrefs'    => \$fix,
  )
  || die "bad commandline option \n";

if ($insert || $fix) {
    my $dbh = Ace->connect( -path => $database ) || die "cannot open connection to db", Ace->error;

    my @genes = $dbh->fetch( -class => 'Gene', -name => '*', -Species => 'Caenorhabditis elegans' );
    foreach my $gene (@genes) {
        next unless ( $gene->Sequence );
        my $primary_id  = $gene->Sequence_name;
        my $public_name = $gene->Public_name;
        my $description = $gene->Concise_description ? $gene->Concise_description . ' [Source: WormBase]' : '';
        $description =~ s/\;/\./g;
	if ($insert){
	        print
"INSERT IGNORE INTO xref (external_db_id,dbprimary_acc,display_label,info_type,info_text,description) VALUES (2400,\"$primary_id\",\"$public_name\",\"DIRECT\",\"Externally assigned relationship between $primary_id and $public_name\",\"$description\");\n";
	}
	else {
		print
"UPDATE xref SET display_label=\"$public_name\", info_type=\"DIRECT\", info_text=\"Externally assigned relationship between $primary_id and $public_name\", description=\"$description\" WHERE external_db_id=2400 AND dbprimary_acc=\"$primary_id\";\n";
		print
"UPDATE gene g,gene_stable_id s SET g.description=\"$description\" WHERE g.gene_id=s.gene_id AND s.stable_id=\"$primary_id\";\n" if length $description > 4;
	
	}
    }
}
if ($update){
	print "UPDATE gene g, gene_stable_id s, xref x SET display_xref_id=xref_id,g.description=x.description WHERE g.gene_id=s.gene_id AND s.stable_id=x.dbprimary_acc AND external_db_id=2400 AND length(x.description) > 4;\n";
}
