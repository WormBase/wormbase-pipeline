ncbi_assemblies(){
  cat \
     <( curl -s 'https://www.ncbi.nlm.nih.gov/assembly/organism/browse-organism/6157/all/?p$l=Ajax&page=1&pageSize=1000&sortCol=4&sortDir=-1') \
     <( curl -s 'https://www.ncbi.nlm.nih.gov/assembly/organism/browse-organism/6231/all/?p$l=Ajax&page=1&pageSize=1000&sortCol=4&sortDir=-1') \
     | grep -B 2 href | perl -e '$/="--"; while(<>){ /<td>(.*?) \(/ ; print "$1-"; /<a href="(.*?)">(.*?)</; print "$2\t$1\n"; }' 
}

our_assemblies(){
  perl -MProductionMysql -e 'for my $core_db (ProductionMysql->previous_staging->core_databases){
       next unless $core_db =~ /core_$ENV{PREVIOUS_PARASITE_VERSION}/;
       my ($spe, $cies, $bp) = split "_", $core_db;
       my $species = join " ", ucfirst($spe), $cies;
       my $url = "http://parasite.wormbase.org/${spe}_${cies}_${bp}";
       printf("%s-%s\t%s?accession=%s\n", $species, ProductionMysql->previous_staging->meta_value($core_db,"assembly.name"), $url, ProductionMysql->previous_staging->meta_value($core_db,"assembly.accession"));}
'
}

join -a1 -a2 -11 -21 -t $'\t' \
  <( ncbi_assemblies | sort -k1,1 -t $'\t' ) <(our_assemblies | sort -k1,1 -t $'\t' ) \
  | perl -e 'my $species; while(<>){my ($species_here, $l) = split ("-", $_, 2); if ($species_here ne $species){print "\n$species_here\n"; $species=$species_here}; print $l}'  \
  | less
