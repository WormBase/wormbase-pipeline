#!/software/bin/perl

# script to generate xrefs from the WormBase AceDB database

use strict;
use Ace;
use Getopt::Long;

my ( $database, $species);
GetOptions(
    'database=s'   => \$database,
    'species=s'    => \$species,
  )
  || die "bad commandline option \n";

$species ||='Caenorhabditis elegans';

print "-- xref data\n";

my $dbh = Ace->connect( -path => $database ) || die "cannot open connection to db", Ace->error;

my $genes = $dbh->fetch_many( -query => "find Gene Species=\"$species\"");

while(my $gene = $genes->next) {
    next unless ( $gene->Sequence );
    next unless ( "${\$gene->Status}" eq 'Live' );
    
    my $primary_id  = $gene->Sequence_name; # the gene stable id
    my $public_name = $gene->Public_name;
    my $description = $gene->Concise_description ? $gene->Concise_description . ' [Source: WormBase]' : ''; # fluf for the description bit
    $description =~ s/\;/\./g;
    
    my @cdses = $gene->Corresponding_CDS;
    
    # gene
    print 'INSERT IGNORE INTO xref (external_db_id,dbprimary_acc,display_label,info_type,info_text) VALUES (2400,"',
        $primary_id,'","',$gene,'","DIRECT","Externally assigned relationship between ',
        "$primary_id and $gene\");\n";

    # locus
    print 'INSERT IGNORE INTO xref (external_db_id,dbprimary_acc,display_label,info_type,info_text,description) VALUES (2440,"',
      $primary_id,'","',$public_name,'","DIRECT","Externally assigned relationship between ',
      "$primary_id and $public_name\",\"$description\");\n";
    
    # wormpep
    foreach my $cds(@cdses){
        print 'INSERT IGNORE INTO xref (external_db_id,dbprimary_acc,display_label,info_type,info_text) VALUES (2420,"',
          "$cds\",\"${\$cds->Corresponding_protein}\",\"DIRECT\",\"Externally assigned relationship between ",
          "$cds and ${\$cds->Corresponding_protein}\");\n";
    }        
}

print "-- update object xref links for genes\n";
print "INSERT IGNORE INTO object_xref (ensembl_id,ensembl_object_type,xref_id) SELECT gene_id,'Gene',xref_id FROM gene join gene_stable_id using (gene_id),xref WHERE xref.dbprimary_acc = gene_stable_id.stable_id AND xref.external_db_id = 2440;\n";
print "INSERT IGNORE INTO object_xref (ensembl_id,ensembl_object_type,xref_id) SELECT gene_id,'Gene',xref_id FROM gene join gene_stable_id using (gene_id),xref WHERE xref.dbprimary_acc = gene_stable_id.stable_id AND xref.external_db_id = 2400;\n";

print "-- update object xref links for proteins\n";
print "INSERT IGNORE INTO object_xref (ensembl_id,ensembl_object_type,xref_id) SELECT translation_id,'Translation',xref_id FROM translation join translation_stable_id using (translation_id),xref WHERE xref.dbprimary_acc = translation_stable_id.stable_id AND xref.external_db_id = 2420;\n";

print "-- update the display xrefs and gene descriptions\n";
print "UPDATE gene g, gene_stable_id s, xref x SET display_xref_id=xref_id,g.description=x.description WHERE g.gene_id=s.gene_id AND s.stable_id=x.dbprimary_acc AND external_db_id=2440;\n";
