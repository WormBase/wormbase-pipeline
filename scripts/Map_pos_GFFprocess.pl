#!/software/bin/perl -w

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use Ace;
use File::Copy;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species);
my( $dbpath);
GetOptions (
            "help"      => \$help,
            "debug=s"	=> \$debug,
	    "test"	=> \$test,
	    "store:s"	=> \$store,
	    "species:s" => \$species,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

$dbpath = $wormbase->orgdb;
my $db = Ace->connect('-path' => $dbpath) or $log->log_and_die("cant open Ace connection to $dbpath\n".Ace->error."\n");


# get the interpolated physical mapping position of the Variations
print "Reading interpolated Variation mapping data\n";
my %variation;
my $query = "find Variation where Interpolated_map_position";
my $vars = $db->fetch_many('-query' => $query);
while (my $var = $vars->next){
  $variation{$var->name} = $var->Interpolated_map_position(2);
}


# get the interpolated physical mapping position of the Genes
print "Reading interpolated Gene mapping data\n";
my %gene;
$query = "find Gene where Interpolated_map_position";
my $genes = $db->fetch_many('-query' => $query);
while (my $gene = $genes->next){
  $gene{$gene->name} = $gene->Interpolated_map_position(2);
}

# get the exact physical mapping position of the Genes
print "Reading exact Gene mapping data\n";
my %gene_exact;
$query = "find Gene where Map AND NEXT AND NEXT = Position";
$genes = $db->fetch_many('-query' => $query);
while (my $gene = $genes->next){
  $gene_exact{$gene->name} = $gene->Map(3);
}

# get the CGC/WGN name of the Genes
print "Reading WGN names of genes\n";
my %locus;
$query = "find Gene where CGC_name";
$genes = $db->fetch_many('-query' => $query);
while (my $gene = $genes->next){
  $locus{$gene->name} = $gene->CGC_name;
}

undef $db;

#read GFF files
my @chroms = $wormbase->get_chromosome_names(-prefix => 1);

foreach my $chroms (@chroms) {
  print "Processing chromosome $chroms\n";

  my $gff_file = $wormbase->GFF_file_name($chroms);
  
  my $INF = $wormbase->open_GFF_file($chroms,undef,$log);
  
  open (NEW,">>$gff_file.new") or $log->log_and_die("cant write new GFF file: $!\n");

  while (my $line = <$INF>) {
    chomp $line;
    if ($line =~ /^#/) {next;}

    my @f = split(/\t/, $line);

    if (defined $f[8] && $f[1] eq 'Allele') {
      my ($allele) = $f[8] =~ /Variation\s+\"(\S+)\"/;
      #if (!defined $allele) {die "no allele ID found Line: $line";}
      $f[8] .= " ; Interpolated_map_position \"$variation{$allele}\"" if (exists $variation{$allele});
	
      print NEW join("\t",@f) . "\n";
      #print join("\t",@f) . "\n";

    } elsif ($f[1] eq 'gene' && $f[2] eq 'gene') {
      my ($gene) = $f[8] =~ /Gene\s+\"(\S+)\"/;
      #if (!defined $gene) {die "no gene ID found Line: $line";}
      $f[8] .= " ; Interpolated_map_position \"$gene{$gene}\""if (exists $gene{$gene});
      $f[8] .= " ; Position \"$gene_exact{$gene}\""if (exists $gene_exact{$gene});
      $f[8] .= " ; Locus \"$locus{$gene}\""if (exists $locus{$gene});
	
      print NEW join("\t",@f) . "\n";
      #print join("\t",@f) . "\n";

    } else {
      print NEW "$line\n";
    }
  }
  close $INF;
  close NEW;
  move("$gff_file.new",$gff_file) unless $wormbase->assembly_type eq 'contig' ;
}


$log->mail();
exit(0);
