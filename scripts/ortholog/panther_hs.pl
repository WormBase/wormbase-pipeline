#!/bin/env perl
# this script will 
# * translate panther human uniprot protein ids to ensembl protein ids
# * take a panther ortholog flatfile and create an ACE output
# caveat: might barf on a new schema as the stable_ids were moved around

use DBI;

my $db = DBI->connect( "DBI:mysql:database=homo_sapiens_core_64_37;host=ensembldb.ensembl.org;port=5306",'anonymous');

my $sth=$db->prepare('SELECT stable_id FROM translation_stable_id JOIN object_xref ox ON(ensembl_id=translation_id AND ensembl_object_type="Translation") JOIN xref x ON (ensembl_object_type="Translation" AND ox.xref_id=x.xref_id) WHERE dbprimary_acc=? AND external_db_id=2200');

while(<>){
    ($from)=$_=~/(WBGene\d+)/;
    ($to)=$_=~/(ENSG\d+)/;
    ($uniprot)=$_=~/UniProtKB=(\w+)/;
        
    $sth->execute($uniprot);
    my @proteins;
    while(my @ary = $sth->fetchrow_array){push @proteins,$ary[0]}

    next unless $proteins[0];
    print "Gene : $from\n";
    map{print "Ortholog_other ENSEMBL:$_ Panther\n"}@proteins;

    foreach my $prot(@proteins){
    print <<HERE

Protein : ENSEMBL:$prot
Species "Homo sapiens"
DB_info Database EnsEMBL ENSEMBL_geneID $to
DB_info Database EnsEMBL ENSEMBL_proteinID $prot
DB_info Database UniProt UniProtAcc $uniprot
DB_info Database UniProt UniProtID $uniprot

HERE
     }
}
