#!/usr/bin/env perl

# Loads a TSV file of gene synonyms matching a gene.stable_id
#
# adapted from https://github.com/Ensembl/plant_tools/blob/master/production/synonyms_gbk/load_synonyms.pl
#
# ..which was adapted from eg-pipelines/scripts/xrefs_pipeline/load_sgd_gene_names.pl
# 

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Utils::CliHelper;

use Carp;
use Data::Dumper;
use Log::Log4perl qw(:easy);

my $logger = get_logger();

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

my $optsd = [ @{ $cli_helper->get_dba_opts() } ];
push( @{$optsd}, "file:s" );
push( @{$optsd}, "verbose" );
push( @{$optsd}, "extdb:s" );

my $opts = $cli_helper->process_args( $optsd, \&pod2usage );

if( !$opts->{'file'} ){
  die "# need -file synonyms.edited.tsv  NOTE: make sure you edit it before loading!\n";
}

if( !$opts->{'extdb'} ){
  die "# need -extdb arg, example -extdb EntrezGene\n";
}

if( $opts->{'verbose'} ) {
  Log::Log4perl->easy_init($DEBUG);
}
else {
  Log::Log4perl->easy_init($INFO);
}


## 1) connect to core db in indicated server

$logger->info( "Loading " . $opts->{'dbname'} );
my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                     -USER   => $opts->{'user'},
                     -PASS   => $opts->{'pass'},
                     -HOST   => $opts->{'host'},
                     -PORT   => $opts->{'port'},
                     -DBNAME => $opts->{'dbname'} );


my $gene_adaptor = $dba->get_GeneAdaptor;
my $dea = $dba->get_DBEntryAdaptor;

## SQL statements

my $synonym_store_sth = $dba->dbc->prepare("INSERT INTO 
        external_synonym (xref_id, synonym)
        VALUES (?, ?)
        ");

my $add_external_db_sth = $dba->dbc->prepare("INSERT INTO
        external_db (db_name, status, priority, type)
        VALUES (? , 'KNOWNXREF', 50, 'MISC')
        ");

##

## 2) read TSV file with synonyms
my $file = $opts->{'file'};

my ($stableid, $synonym, $word);
my (%syns);

open(TSV,'<',$file) || die "# ERROR: cannot read $file: $!";

LINE: while ( my $line = <TSV> ) {

    #SSTP_0000076900	Ss_dcap_2
    #SSTP_0001285400	Ss_mcm_4

    next if($line =~ m/^#/ || $line !~ m/\t/);

    chomp($line); 
    ($stableid, $synonym) = split(/\t/, $line);

    # accumulate synonyms of the same stable_id
    push @{$syns{$stableid}}, $synonym;
}
close(TSV);

## 3) Check that the external_db exists, and add it if not

my $dbname = $opts->{'extdb'};
my $dbRefs = $dea->get_external_db_ids($dbname, 'NULL', 1);
unless (scalar(@$dbRefs)>0){
  $logger->info( "Adding new external_db entry for ", $opts->{'extdb'});
  $add_external_db_sth->execute($opts->{'extdb'});
  $add_external_db_sth->finish();
  $dbRefs = $dea->get_external_db_ids($dbname, 'NULL', 1);
}

## 4) create display_xrefs linked to synonyms

foreach $stableid (keys(%syns)) {
   
    # check target gene exists 
        my $gene = $gene_adaptor->fetch_by_stable_id($stableid);
        if ( !$gene ) {
        $logger->info( "Cannot find $stableid, skip it");
        next;
     }
 	
    # check whether gene already has display_xref  
      my $old_display_xref = $gene->display_xref();

      if( $old_display_xref ) {
        $logger->info( "$stableid has display_xref_id set to ", $old_display_xref->display_id() );
        # get its xref_id
	my $xref_id = $old_display_xref->dbID();
        $logger->info( "Its xref_id is $xref_id");
	# get existing synonyms for this xref
        my @existing_synonyms = @{$old_display_xref->get_all_synonyms };
        # hang the synonyms off the existing display xref
	SYN: foreach $synonym (@{$syns{$stableid}}){
	  # check the synonym isn't already added
	  foreach my $existing_synonym (@existing_synonyms){
            if ($synonym eq $existing_synonym){
	      $logger->info( "$synonym already associated with $stableid, skipping");
	      next SYN;
	    }
	  }
          $synonym_store_sth->execute($xref_id, $synonym);
	  $logger->info( "Added $synonym to $stableid");
        }
	$synonym_store_sth->finish();
    }

   # no existing display_xref so we make a new one
   # where there are multiple synonyms, the display xref will just be the first one
    else{
      my $new_display_xref = Bio::EnsEMBL::DBEntry -> new (
	     -PRIMARY_ID  => ${$syns{$stableid}}[0],
             -DBNAME      => $opts->{'extdb'},
             -DISPLAY_ID  => ${$syns{$stableid}}[0],
             -INFO_TYPE   => 'SEQUENCE_MATCH',
         );
      # add all synonyms to the new xref  
      foreach $synonym (@{$syns{$stableid}}){
        $new_display_xref->add_synonym($synonym);
	$logger->info( "Added $synonym to $stableid");
      } 
      $dbRefs = $dea->get_external_db_ids($dbname, 'NULL', 1);
      my $dbRef = shift(@$dbRefs);
      my $xref_id = $dea->_store_or_fetch_xref($new_display_xref,$dbRef);
      $new_display_xref->dbID($xref_id);

      # and update the gene
      $gene->display_xref($new_display_xref);
      $gene_adaptor->update($gene);
   }
}   
