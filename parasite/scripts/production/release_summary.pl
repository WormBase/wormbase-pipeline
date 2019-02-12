#!/usr/bin/env perl
use ProductionMysql;
use YAML;
my %result;
for my $core_db (ProductionMysql->previous_staging->core_databases) {
  my ($spe,$cies,$bioproject) = split "_", $core_db;
  my $species = "${spe}_${cies}";
  $result{$species}{$bioproject}{previous_core_db} = $core_db;
  $result{$species}{$bioproject}{previous}=join("\t",
     ProductionMysql->previous_staging->meta_value($core_db, "assembly.name"),
     ProductionMysql->previous_staging->meta_value($core_db, "genebuild.version")
  );
 }
for my $core_db (ProductionMysql->staging->core_databases) {
  my ($spe,$cies,$bioproject) = split "_", $core_db;
  my $species = "${spe}_${cies}";
  $result{$species}{$bioproject}{current_core_db} = $core_db;
  $result{$species}{$bioproject}{current}=join("\t",
     ProductionMysql->staging->meta_value($core_db, "assembly.name"),
     ProductionMysql->staging->meta_value($core_db, "genebuild.version")
  );
}
for my $species (keys %result){
   for my $bioproject (keys %{$result{$species}}){
      delete $result{$species}{$bioproject} if $result{$species}{$bioproject}{previous} eq $result{$species}{$bioproject}{current};
   }
   delete $result{$species} unless %{$result{$species}};
}
print Dump(\%result);
print "\n";
 
my @links;
for my $species (keys %result){
   for my $bioproject (keys %{$result{$species}}){
      my $r = $result{$species}{$bioproject};
      my $Species=ucfirst($species);
      $Species =~ s/_/ /g;
      push @links, sprintf("<a href=\"http://parasite.wormbase.org/${species}_$bioproject\">%s (%s)</a>", $Species, uc($bioproject));
   }
}
print join ", ", @links;
print "\n";
unless (grep {$_ =~ /elegans/} keys %result ){
   print "\n!!!\nWHERE IS C ELEGANS YO\n";
}