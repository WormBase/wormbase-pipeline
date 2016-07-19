#!/usr/local/bin/perl -T 
#########
## Author:        mh6
## Group:         wormbase
## Maintainer:    mh6
## Last modified: $Date: 2013-10-18 12:11:15 $
## Id:            $Id: nameserver_json.pl,v 1.1 2013-10-18 12:11:15 mh6 Exp $
## Source:        $Source: /cvsroot/CVSmaster/wormbase/scripts/NAMEDB/scripts/nameserver_json.pl,v $

use warnings;
use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use JSON;
use DBI;

my $q = CGI->new;

my %valid_domains = (Gene => 'Gene', Variation =>'Variation', Feature => 'Feature');

my $domain = $q->param('domain')||'Gene';

print $q->header('application/json');

my $dbh=DBI->connect('DBI:mysql:database=nameserver_live;host=web-wwwdb-core-02;port=3449','mh6','mh6')or die $DBI::errstr;
my $sth=$dbh->prepare('SELECT object_public_id,name_type_name,object_name,object_live FROM primary_identifier pi LEFT JOIN  secondary_identifier using (object_id) LEFT JOIN name_type nt USING(name_type_id) LEFT JOIN domain d ON (d.domain_id=pi.domain_id) WHERE domain_name=?');

$sth->execute($valid_domains{$domain});
my $array_ref = $sth->fetchall_arrayref();
my $json=JSON->new();
$json->pretty(1);
print $json->encode($array_ref);
$dbh->disconnect;
