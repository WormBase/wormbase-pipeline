#!/usr/local/bin/perl5.6.0 -w

# Keith Bradnam (krb@sanger.ac.uk)
# Finds and removes duplicate clones in worm01 database (duplicates arise 
# from new sequence versions of clones being added to the database)
# script removes redundant entries from clone, contig, and InputIdAnalysis


use strict;
use DBI;


####################################################################
# set some parameters
####################################################################
# mysql database
my $dbhost = "ecs1f";
my $dbuser = "wormadmin";
my $dbname = "worm01";
my $dbpass = "worms";


####################################################################
# connect to the Mysql database
####################################################################

my $dbh = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
    || die "cannot connect to db, $DBI::errstr";

my $clone_table     = $dbh->prepare (q{ SELECT internal_id, embl_id, embl_version FROM clone 
				       ORDER BY embl_id,embl_version} );

my $contig_table    = $dbh->prepare (q{ SELECT internal_id, id, dna FROM contig where clone = ? } );

my $delete_clone    = $dbh->prepare (q{ DELETE FROM clone           WHERE internal_id = ? } );
my $delete_contig   = $dbh->prepare (q{ DELETE FROM contig          WHERE internal_id = ? } );
my $delete_dna      = $dbh->prepare (q{ DELETE FROM dna             WHERE id = ? } );
my $delete_inputId  = $dbh->prepare (q{ DELETE FROM InputIdAnalysis WHERE inputId = ? } );


#grab clone information, sorted by Accession then version number
$clone_table->execute;
my @row = $clone_table->fetchrow_array;

#store first line of table to be able to compare to
my $last_id           = $row[0];
my $last_embl_id      = $row[1];
my $last_embl_version = $row[2];
my $counter = 0;

# loop through rest of clone table

while (@row = $clone_table->fetchrow_array) {

  # Do adjacent rows share the same accession?
  if ($row[1] eq $last_embl_id){
    print "$row[1].$row[2] ($row[0])  replaces   $last_embl_id.$last_embl_version ($last_id)\n";

    # now grab contig internal_id and dna id entries from contig table
    $contig_table->execute($last_id);
    my @new_row = $contig_table->fetchrow_array;
    my $contig_internal_id = $new_row[0];    
    my $contig_id = $new_row[1];
    my $dna_id = $new_row[2];

    print "\tDELETE FROM clone           WHERE internal_id = $last_id\n";
    print "\tDELETE FROM contig          WHERE internal_id = $contig_internal_id\n";
    print "\tDELETE FROM dna             WHERE id = $dna_id\n";
    print "\tDELETE FROM InputIdAnalysis WHERE inputId = $contig_id\n";
    print "\n\n";
    $delete_clone->execute($last_id);
    $delete_contig->execute($contig_internal_id);
    $delete_dna->execute($dna_id);
    $delete_inputId->execute($contig_id);

    

    $counter++;
  }
  $last_id = $row[0];
  $last_embl_id = $row[1];
  $last_embl_version = $row[2];
   
}

print "\n$counter entries deleted in worm01 database.\n";


# Close active database handles

$clone_table->finish;
$contig_table->finish;
$delete_clone->finish;
$delete_contig->finish;
$delete_dna->finish;
$delete_inputId->finish;

$dbh->disconnect;


exit 0;
