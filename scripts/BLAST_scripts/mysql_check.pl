#!/usr/local/perl -w

use DBI;
use Getopt::Long;

my ($dbname, $debug);

GetOptions ( 
	    "dbname:s" => \$dbname,
	    "debug"    => \$debug
	   );


my $dbh = DBI->connect("dbi:mysql:database=$dbname;host=ecs1f","wormro") or die "Cant connect to $dbname\n$DBI::errstr";
print "connected\n------------\n\n";

# get table list

my $sth = $dbh->prepare("SHOW TABLES");
$sth->execute;

my @row;
my @tables;
my @corrupt;
while( @row = $sth->fetchrow_array ) {
  push( @tables,$row[0]);
}

#test tables
foreach my $table (@tables) {
  $sth = $dbh->prepare("CHECK TABLE $table");
  $sth->execute;
  while( @row = $sth->fetchrow_array ) {
    print "@row\n" if $debug;
    foreach ( @row ) {
      if ($_ =~ /Corrupt/) {
	print "CORRUPT TABLE $dbname.$table\n";
	push (@corrupt,$dbname);
      }
    }
  }
}

$dbh->disconnect;

print "\n------------\ndisconnected\n\n\n";

if ( $corrupt[0] ) {
  print "DATABASE $dbname CORRUPT TABLES @corrupt \nUse REPAIR TABLE  to fix table eg\n\nREPAIR TABLE $corrupt[0]\n\n";
}
else {
  print "DATABASE $dbname OK\n";
}
