#!/usr/bin/perl -w

use lib (-e "/wormsrv2/scripts") ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;

my $database;
my $version;
GetOptions ( 'database:s' => \$database,
	     'version:s'  => \$version);

print STDERR "\nStart $0 at ",&runtime,"\n=======================\n\n";

$database = "/wormsrv2/autoace" unless $database;
my $GFF_dir = glob("$database/CHROMOSOMES");

#table maker defs used
my $peptides = "wormpep_lengths.def";
my $gene_details = "orthologs.def";


my %chrom_lengths = (
		     'I'   => '15080552',
		     'II'  => '15279311',
		     'III' => '13783317',
		     'IV'  => '17493785',
		     'V'   => '20922231',
		     'X'   => '17718851'
		    );

my $tace = &tace;

# Get protein lengths
print STDERR "Getting peptide lengths from $database\n";
open ( PEP ,"echo table-maker -p $database/wquery/$peptides | $tace $database |") or die "cant access $database\t$!\n";
my %pep_length;
while( <PEP> ) {
  my @data = split;
  next unless /WP/;
  $data[0] =~ s/WP://;
  $data[0] =~ s/\"//g;
  $pep_length{$data[0]} = $data[1];
}
close PEP;

# Get gene info 
print STDERR "Getting gene details from $database\n";
open ( GENE ,"echo table-maker -p $database/wquery/$gene_details | $tace $database |") or die "cant access $database\t$!\n";
my %genes;
while( <GENE> ) {
  #"D1007.5a"      "D1007" "I"     -1.0361 "curated"       "WP:CE29034"
  s/\"//g;
  my @d = split;
  next unless $d[3];
  $genes{$d[0]}->{'clone'}  = $d[1];
  $genes{$d[0]}->{'chrom'}  = $d[2];
  $genes{$d[0]}->{'map'}    = $d[3];
  $genes{$d[0]}->{'method'} = $d[4];
  $d[5] =~ s/WP://;
  $genes{$d[0]}->{'peptide'}= $d[5];
}

close GENE;

my $genome_base = 0;
print STDERR "Parsing data from GFF files . . . \n";
foreach my $chrom ( qw(I II III IV V X ) ) {
  print STDERR "\tchromosome_$chrom \n";
  open( GFF, "<$GFF_dir/CHROMOSOME_${chrom}.gff") or die "cant open gff $GFF_dir/CHROMOSOME_${chrom}.gff\t$!\n";
  while( <GFF> ) {
    #CHROMOSOME_X    curated CDS     39220   41753   .       +       .       CDS "Y73B3A.21"
    my @d = split(/\s+/,$_);
    next unless ( $d[1] eq "curated" and $d[2] eq "CDS" );
    $d[9] =~ s/\"//g;
    my $start;
    my $end;
    my $strand;
    if( $d[6] eq "+" ) {
      $start = $d[3];
      $end = $d[4];
      $strand = 1;
    }
    else {
      $start = $d[4];
      $end = $d[3];
      $strand = -1;
    }

    $genes{$d[9]}->{'start'} = $start;
    $genes{$d[9]}->{'end'} = $end;
    $genes{$d[9]}->{'strand'} = $strand;
    $genes{$d[9]}->{'genome_start'} = $genome_base + $start;
    $genes{$d[9]}->{'genome_end'} = $genome_base + $end;
  }
  close GFF;
  $genome_base += $chrom_lengths{"$chrom"};
}

$version = &get_wormbase_version unless $version;
print STDERR "Writing output to \n";
foreach my $gene ( keys %genes ) {
  print "$gene\t",uc "$gene","\t",                      #gene protein
    $genes{$gene}->{'clone'},"\t",                       #clone
    $genes{$gene}->{'chrom'},"\t",                       #chromosome
    $genes{$gene}->{'strand'},"\t",                      #strand
    $pep_length{"$genes{$gene}->{'peptide'}"},"\t",      #peptide length
    "test\t",                                            #spliced length
    "1\t1\t",  
    $genes{$gene}->{'start'},"\t",                       #chrom start
    $genes{$gene}->{'end'},"\t",                         #chrom end
    $genes{$gene}->{'genome_start'},"\t",                #genome start
    $genes{$gene}->{'genome_end'},"\t",                  #genome end
    $genes{$gene}->{'map'},"\t",                         #map position
    "source\t",
    "WS$version\t",                                          #version
    $genes{$gene}->{'method'},"\n";                      #method
  }

print STDERR "\n=======================\nFinished at ",&runtime,"\n";
