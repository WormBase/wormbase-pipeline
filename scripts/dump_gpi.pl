#!/usr/bin/env perl
#
# exporter to dump gene / transcript / protein info as GPI file
#   specs: http://www.geneontology.org/page/gene-product-information-gpi-format
#
# usage:
#   perl dump_gpi.pl -species elegans


use strict;
use lib $ENV{CVS_DIR};
use Wormbase;
use Log_files;
use Getopt::Long;
use Storable;
use Ace;
use IO::File;


my ($debug,$test,$species,$store,$output,$database);

GetOptions(
   'debug=s'   => \$debug,   # send log mails only to one person
   'species=s' => \$species, # specify the species to run on
   'test'      => \$test,    # use the test database instead of the live one
   'store=s'   => \$store,   # pass a storable (for the build)
   'output=s'  => \$output,  # write somewhere else and not to REPORTS/
   'database=s'=> \$database,# specify a different database than BUILD/$species
)||die(@!);


my $wormbase;
if ($store) {
  $wormbase = retrieve($store) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

$output||=$wormbase->reports . "/" . $wormbase->species . ".gene_product_info.gpi";
my $outfile = IO::File->new($output,'w')||$log->log_and_die(@!);

$log->write_to("creating a GPI file at $output for ${\$wormbase->long_name}\n");

my $db = Ace->connect(-path => ($database ||$wormbase->autoace));

my $genes = $db->fetch_many(-query => 'Find Gene;Species="'.$wormbase->long_name.
            '";Live;Corresponding_transcript OR Corresponding_CDS OR Corresponding_pseudogene')
            or $log->log_and_die(Ace->error);

print $outfile "!gpi-version: 1.2\n";

while (my $g = $genes->next){
   # Gene block
   printf $outfile "WB\t%s\t%s\t%s\t%s\t%s\ttaxon:%i\t\t%s\n", 
       $g,
       $g->Public_name,
       ($g->Gene_class?$g->Gene_class->Description:''),
       join('|',$g->Other_name),
       $g->Biotype->SO_name,
       $g->Species->NCBITaxonomyID,
       uniprot_xref($g)
   ;

   foreach my $t($g->Corresponding_CDS){
     # Transcript/CDS block
     printf $outfile "WB\t%s\t%s\t%s\t%s\ttranscript\ttaxon:%i\t\%s\t%s\n", 
       $t,
       $g->Public_name,
       ($g->Gene_class?$g->Gene_class->Description:''),
       join('|',$g->Other_name),
       $t->Species->NCBITaxonomyID,
       "WB:$g",
       ''
     ;
     foreach my $p($t->Corresponding_protein){
       # Protein block
       printf $outfile "WB\t%s\t%s\t%s\t%s\tprotein\ttaxon:%i\t\%s\t%s\n", 
       $p,
       uc($g->Public_name),
       ($g->Gene_class?$g->Gene_class->Description:''),
       join('|',$g->Other_name),
       $p->Species->NCBITaxonomyID,
       "WB:$t",
       uniprot_xref($p->fetch)
     ;
     }
   }

   foreach my $t($g->Corresponding_transcript){
       next if "${\$t->Method}" eq 'Coding_transcript';
       # ncRNA transcript block
       printf $outfile "WB\t%s\t%s\t%s\t%s\t%s\ttaxon:%i\t\%s\t%s\n", 
       $t,
       $g->Public_name,
       ($g->Gene_class?$g->Gene_class->Description:''),
       '',
       $t->Method->GFF_SO->SO_name,
       $t->Species->NCBITaxonomyID,
       "WB:$g",
       uniprot_xref($t->fetch)
     ;
   }
}


# return the UNiProtIsoformAcc if it exists, else return the UniProtAcc if there is only *one*
sub uniprot_xref{
  my ($object)=@_;
  
  my @result;
  my @rnacentral;

  # case for Genes
  map {my $type="$_";
       map{
         my $subtype="$_";
         map{
          $result[0]='UniProtKB:'.$_ if $type eq 'SwissProt';
          push @result,'UniProtKB:'.$_ if $type eq 'TrEMBL' && ! $result[0];
          push @rnacentral,'RNAcentral:'.$_ if $type eq 'RNAcentral';
         }$_->col 
       }$_->col
  }$object->at('Identity.DB_info.Database');

  # case for Proteins/CDS
  map {my $type=$_;
       map{
        my $subtype="$_";
        $result[0]= "UniProtKB:".$_->right(1) if $subtype eq 'UniProtIsoformAcc';
        push @result,"UniProtKB:".$_->right(1) if $subtype eq 'UniProtAcc' && ! $result[0];
        push @rnacentral,'UniProtKB_GCRP:'.$_->right(1) if $type eq 'UniProtGCRP';
        push @rnacentral,'RNAcentral:'.$_->right(1) if $type eq 'RNAcentral';
       }$_->col
  }$object->at('DB_info.Database');
 
  push (@rnacentral,$result[0]) if scalar(@result)==1;
  return join('|',@rnacentral);
}

$log->mail;

# returns status 13 when finishing, so try to clean up
$outfile->close;
$db->close;

exit(0);
