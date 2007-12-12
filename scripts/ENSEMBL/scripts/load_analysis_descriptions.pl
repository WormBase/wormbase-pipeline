#!/usr/local/ensembl/bin/perl -w

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;

$! = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 3306;
my $dbname;

&GetOptions( 
            'dbhost=s'      => \$dbhost,
            'dbname=s'      => \$dbname,
            'dbuser=s'      => \$dbuser,
            'dbpass=s'      => \$dbpass,
            'dbport=s'      => \$dbport,
           );



my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor
  (
   -host   => $dbhost,
   -user   => $dbuser,
   -dbname => $dbname,
   -pass   => $dbpass,
   -port   => $dbport,
  );

my $aa = $db->get_AnalysisAdaptor;
my $analyses = $aa->fetch_all;
my %hash;
foreach my $analysis(@$analyses){
  $hash{lc($analysis->logic_name)} = $analysis;
}
LINE:while(<>){
  /^(\S+)\s+(.+)/ and do {
    my ($logic_name, $description) = ($1, $2);

    $description =~ s/^\s+//;
    $description =~ s/\s+$//;

    next if not $description;

    if (exists $hash{lc($logic_name)}) {
      my $analysis = $hash{lc($logic_name)};
      
      $analysis->description($description);

      my $display_label = $logic_name;
      $display_label =~ s/_//g;
      $analysis->display_label($display_label);

      $aa->update($analysis);

      delete $hash{lc($logic_name)};
    }
  }
}


if ( scalar(keys %hash)==0) { 
  print "\nAll analysis descriptions have been updated, every analysis has a description now\n" ; 
}else{
  foreach my $k (keys %hash) {
    warn "No description was found for logic name $k ( ".$hash{$k}->dbID." ) \n";
  }
}
