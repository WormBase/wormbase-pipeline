#!/usr/local/bin/perl -wT
#author ar2

use strict;
use SangerPaths qw(core celegans);
use SangerWeb;
use NameDB_handler;
use SangerWeb;
use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
$| = 1;

my $DB;
my $USER = param('user');
my $sw = SangerWeb->new();
if ($sw->is_dev()) {
    $DB = 'test_wbgene_id;utlt-db;3307';
} else {
    $DB = 'wbgene_id;shap;3303';
}

print <<END;
Content-Type: text; charset=iso-8859-1

END

croak ("no user or database") unless ($DB and $USER);
#print "USER:$USER<br>DB:$DB<br>";

my $db = get_db_connection($USER,$DB);
my $query =<<END;
SELECT primary_identifier.object_public_id, secondary_identifier.object_name, primary_identifier.object_live
    FROM primary_identifier,secondary_identifier 
    WHERE secondary_identifier.object_id = primary_identifier.object_id and primary_identifier.domain_id = 3
    ORDER by object_public_id
    ;
END
    my $sth = $db->dbh->prepare($query);
    $sth->execute or die "Unable to execute query: $db->errstr\n";
    my $row;
   #print "query executed<br>";
    while ($row = $sth->fetchrow_arrayref) {
      print "$row->[0]\t$row->[1]\t";
      if ($row->[2] eq "1") {print "Live\n";}
      elsif ($row->[2] eq "0") {print "Dead\n";}
      else {print "ERROR no status\n";}
    }
    $sth->finish;



sub get_db_connection {
    my $DOMAIN = 'Variation';
    my ($USER,$DB) = @_;

    my $db = NameDB->connect($DB,$USER,$USER,1); #1 is for web output
    
    $db || return(undef);
    $db->setDomain($DOMAIN);
    
    return($db);
}
