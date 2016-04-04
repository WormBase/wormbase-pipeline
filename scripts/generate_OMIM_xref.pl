#!/software/bin/perl -w
#
# generate_OMIM_xref.pl
#
# Usage : generate_OMIM_xref.pl [-options]
#
# Target file unless specified - misc_output/WS<wormbase_version>_OMIMXREF.dat
# Soirce database unless specified - autoace
#
use strict;
use Getopt::Long;
use Storable;
use lib $ENV{CVS_DIR};
use Wormbase;

my ($test,
    $debug,
    $store,
    $species,
    $database,
    $wormbase,
    $outfile,
    );

GetOptions (
  "test"            => \$test,
  "debug=s"         => \$debug,
  "store:s"         => \$store,
  "species:s"       => \$species,
  "database:s"      => \$database,
  "outfile:s"       => \$outfile,
    );


if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species,
      );
}

#################################
# Set up some useful stuff      #
#################################
# establish log file.
my $log = Log_files->make_build_log($wormbase);

$species = $wormbase->species;
my $tace = $wormbase->tace;
my $full_species_name = $wormbase->full_name;
my $wormbase_version = $wormbase->get_wormbase_version_name;
my $dbdir;

if (defined $database) {
  $dbdir = $database;
}
unless (defined $database) {
  $dbdir = $wormbase->autoace;
}

if (not defined $outfile) {
  $outfile = $wormbase->misc_output . "/${wormbase_version}_OMIMXREF.dat";
}

##########################
# MAIN BODY OF SCRIPT
##########################

my (%genes,
    %geneomims,);

my @Unique_genes;
my $out_fh;
$log->write_to("Generating OMIM_XREF Query\n");
my $query = &generate_OMIMXREF_query();
$log->write_to("Retrieving OMIM data from $dbdir\n");
my $command = "Table-maker -p $query\nquit\n";
my $count2;
open(my $tacefh,  "echo '$command' | $tace $dbdir |");
## Gene OMIM_ID ##
print "Table_maker_finished\n";
while(<$tacefh>) {
  $count2++;

  unless (/\"/) {next;}
  chomp; s/\"//g;
  #print "$_\n";
  my $F;
  my @F = split"\t";

  #Store WBID :: Human readables in associated arrays.
  if (!grep (/$F[0]/, @Unique_genes))        {push @Unique_genes, $F[0];}
  if ((exists $F[1])&&!grep (/$F[1]/, @{$geneomims{$F[0]}}))      {push @{$geneomims{$F[0]}}, $F[1];}
}
$log->write_to("Finished Retrieving OMIM_XREF data\n");


open($out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");
print "Output file opened $outfile\n";

my $count;
#WBGene00026070=>111100%136836%111100%136836%111100%136836%111100%136836
print $out_fh "# WBGene OMIM gene/disaeas XREF data from $wormbase_version release of WormBase.\n";
print $out_fh "# WBGene_ID=>OMIM_IDs[%]\n";
foreach my $gene (@Unique_genes) {
  $count++;
  if ($test) {print "$count Processing $gene\n";}
  print $out_fh "$gene=>";

  #Variations genetic modification of organism if modified, 
  if (exists $geneomims{$gene}) {
    if (defined $geneomims{$gene}) {
      print $out_fh "$geneomims{$gene}[0]";
      foreach my $i ( 1 .. $#{ $geneomims{$gene} } ) {
	print $out_fh "%$geneomims{$gene}[$i]";
      }
      print $out_fh "\t";
    }
  }
  else {
    print $out_fh "EMPTY\t";
  }
  print $out_fh "\n";
}
close($out_fh);
$log->write_to("\nOMIM_XREF data generated for $count genes........ ;) \n\n");
$log->mail();
exit(0);


sub generate_OMIMXREF_query {
  my ($full_species) = @_;

  my $tmdef = "/tmp/gene_tmquery.$$.def";
  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  

  my $condition = "";

  my $tablemaker_template = <<"EOF";
 
Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Genes_elegans 
From 1 
 
Colonne 2 
Width 12 
Mandatory 
Hidden 
Class 
Class Database 
From 1 
Tag Database 
Condition "OMIM"
 
Colonne 3 
Width 12 
Optional 
Hidden 
Class 
Class Database_field 
Right_of 2 
Tag  HERE
Condition "gene"  
 
Colonne 4 
Width 12 
Optional 
Visible 
Class 
Class Text 
Right_of 3 
Tag  HERE  

 
EOF

  print $qfh $tablemaker_template;
  return $tmdef;


}
__END__

=pod

=head2 NAME - generate_OMIM_xref.pl

=head1 USAGE

=over 4

=item generate_OMIM_xref.pl  [-options]

=back

This script generates an xref file between WBGene_ID and OMIM gene and disease IDs.

generate_OMIM_xref.pl

MANDATORY arguments:

=over 4

=item None at present.

=back

generate_OMIM_xref.pl  OPTIONAL arguments:

=over 4

=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.

=back

=over 4

=item -test, Test mode.

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Paul Davis (paul.davis@ebi.ac.uk)

=back

=cut
