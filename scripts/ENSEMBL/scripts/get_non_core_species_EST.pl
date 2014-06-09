#!/software/bin/perl -w
#
# get_non_core_EST.pl
# 
# by Gary Williams            
#

# This gets the EST sequences in fasta format from the ENA and formats
# them so that just their accession number of on the title line for
# use by the EST-STAR hive aligner

#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2014-06-09 12:25:42 $      

use strict;                                      
use Getopt::Long;
use Carp;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


######################################
# variables and command-line options # 
######################################

my ($help, $project_dir);
my (@only_species, %only_species);

GetOptions ("help"           => \$help,
	    "project_dir:s"  => \$project_dir, # usually $wormbase or $parasite
	    "species=s@"     => \@only_species,
	    );

# Display help if required
&usage("Help") if ($help);

if (!defined $project_dir) {$project_dir = "/nfs/panda/ensemblgenomes/wormbase/BUILD_DATA/cDNA/"}

##########################
# MAIN BODY OF SCRIPT
##########################

my %nembase_species_decode = (
                   
			      "Angiostrongylus cantonensis" => "AAC",
			      "Ancylostoma braziliense" => "ABC",
			      "Ancylostoma caninum" => "ACC",
			      "Anisakis simplex" => "AIC",
			      "Ascaris lumbricoides" => "ALC",
			      "Ascaris suum" => "ASC",
			      "Ancylostoma ceylanicum" => "AYC",
			      "Brugia malayi" => "BMC",
			      "Brugia pahangi" => "BPC",
			      "Bursaphelenchus mucronatus" => "BUC",
			      "Bursaphelenchus xylophilus" => "BXC",
			      "Caenorhabditis brenneri" => "CBC",
			      "Caenorhabditis elegans" => "CEC",
			      "Caenorhabditis briggsae" => "CGC",
			      "Caenorhabditis japonica" => "CJC",
			      "Caenorhabditis remanei" => "CRC",
			      "Caenorhabditis sp. 5" => "CSC",
			      "Ditylenchus africanus" => "DAC",
			      "Dirofilaria immitis" => "DIC",
			      "Dictyocaulus viviparus" => "DVC",
			      "Globodera mexicana" => "GMC",
			      "Globodera pallida" => "GPC",
			      "Globodera rostochiensis" => "GRC",
			      "Heterorhabditis bacteriophora" => "HBC",
			      "Haemonchus contortus" => "HCC",
			      "Heterodera glycines" => "HGC",
			      "Heterodera schachtii" => "HSC",
			      "Loa loa" => "LLC",
			      "Litomosoides sigmodontis" => "LSC",
			      "Meloidogyne arenaria" => "MAC",
			      "Meloidogyne chitwoodi" => "MCC",
			      "Meloidogyne hapla" => "MHC",
			      "Meloidogyne incognita" => "MIC",
			      "Meloidogyne javanica" => "MJC",
			      "Meloidogyne paranaensis" => "MPC",
			      "Necator americanus" => "NAC",
			      "Nippostrongylus brasiliensis" => "NBC",
			      "Onchocerca ochengi" => "OCC",
			      "Oesophagostomum dentatum" => "ODC",
			      "Onchocerca flexuosa" => "OFC",
			      "Ostertagia ostertagi" => "OOC",
			      "Onchocerca volvulus" => "OVC",
			      "Parelaphostrongylus tenuis" => "PAC",
			      "Pratylenchus penetrans" => "PEC",
			      "Pristionchus pacificus" => "PPC",
			      "Panagrolaimus superbus" => "PSC",
			      "Parastrongyloides trichosuri" => "PTC",
			      "Pratylenchus vulnus" => "PVC",
			      "Radopholus similis" => "RSC",
			      "Steinernema carpocapsae" => "SCC",
			      "Steinernema feltiae" => "SFC",
			      "Strongyloides ratti" => "SRC",
			      "Strongyloides stercoralis" => "SSC",
			      "Toxocara canis" => "TCC",
			      "Teladorsagia circumcincta" => "TDC",
			      "Trichostrongylus vitrinus" => "TIC",
			      "Toxascaris leonina" => "TLC",
			      "Trichuris muris" => "TMC",
			      "Trichinella spiralis" => "TSC",
			      "Trichuris vulpis" => "TVC",
			      "Wuchereria bancrofti" => "WBC",
			      "Xiphinema index" => "XIC",
			      "Zeldia punctata" => "ZPC",
		     );

my %nematode_net_species_decode = (

				   "Ancylostoma caninum" => "030219.ancylostoma_caninum",
				   "Ancylostoma ceylanicum" => "030325.ancylostoma_ceylanicum",
				   "Ascaris suum" => "020516.ascaris_suum",
				   "Brugia malayi" => "030818.brugia_malayi",
				   "Caenorhabditis remanei" => "060510.c.remanei",
				   "Dirofilaria immitis" => "030709.dirofilaria_immitis",
				   "Ditylenchus africanus" => "071107.ditylenchus_africanus",
				   "Globodera pallida" => "060810.globodera_pallida",
				   "Globodera rostochiensis" => "071219.globodera_rostochiensis",
				   "Heterodera glycines" => "071019.heterodera_glycines",
				   "Heterodera schachtii" => "051230.heterodera_schachtii",
				   "Haemonchus contortus" => "040702_haemonchus_contortus",
				   "Meloidogyne arenaria" => "051212_meloidogyne_arenaria",
				   "Meloidogyne chitwoodi" => "050421.meloidogyne_chitwoodi",
				   "Meloidogyne hapla" => "050523.meloidogyne_hapla",
				   "Meloidogyne incognita" => "050621.meloidogyne_incognita",
				   "Meloidogyne javanica" => "060103_meloidogyne_javanica",
				   "Meloidogyne parananesis" => "051214.meloidogyne_parananesis",
				   "Nippostrongylus brasiliensis" => "060608.nippostrongylus_brasiliensis",
				   "Onchocerca flexuosa" => "071212.onchocerca_flexuosa",
				   "Ostertagia ostertagi" => "060103.ostertagia_ostertagi",
				   "Parastrongyloides trichosuri" => "060103.parastrongyloides_trichosuri",
				   "Pratylenchus penetrans" => "030402.pratylenchus_penetrans",
				   "Pristionchus pacificus" => "010706.pristionchus_pacificus",
				   "Radopholus similis" => "071114.radopholus_similis",
				   "Strongyloides ratti" => "031211.strongyloides_ratti",
				   "Strongyloides stercoralis" => "011025.strongyloides_stercoralis",
				   "Toxocara canis" => "060103.toxocara_canis",
				   "Trichinella spiralis" => "021105.trichinella_spiralis",
				   "Trichuris muris" => "071220.trichuris_muris",
				   "Xiphinema index" => "050325_xiphinema_index",
				   "Zeldia punctata" => "010618.zeldia_punctata",
				 );

# Get the Registry
# This assumes the environment variable ENSEMBL_REGISTRY is set, this is
#     used as the name of the configuration file to read.
# Or, the file .ensembl_init exists in the home directory, it is
#     used as the configuration file.
Bio::EnsEMBL::Registry->load_all();

map { $only_species{$_} = 1 } @only_species;

foreach my $species (keys %only_species){

  my $meta_container = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'MetaContainer');
  my $taxon_id = $meta_container->get_taxonomy_id();
  my $binomial_species_name = $meta_container->get_scientific_name();

  print "$species taxon ID: $taxon_id\n";

  # the EST files go here
  my $target_dir = "$project_dir/$species";

  if (!-d $target_dir) {mkdir $target_dir, 0777}
  my $outfile = "${target_dir}/EST";
  $outfile =~ s/\s/_/g;

  my $old_size = undef;
  if (-e $outfile) {
    $old_size = -s $outfile;
  }

  open (OUT, ">$outfile") || die("Can't open the ouput file $outfile");

  # now get the EST/mRNA sequences
  my $query = 'http://www.ebi.ac.uk/ena/data/warehouse/search?query=%22mol_type=%22mRNA%22%20AND%20tax_tree(TAXON)%22&result=sequence_release&display=fasta';
  $query =~ s/TAXON/$taxon_id/;
  open (DATA, "wget -q -O - '$query' |") || die("Can't get mRNA data for $species\n");

  my $entry_count=0;
  while(my $line = <DATA>) {
    if ($line =~ /^>ENA\|\w+\|(\S+)\s/) {
      $line = ">$1\n";
      $entry_count++;
    }
    print OUT $line;
  }

  print "$entry_count EST/mRNA entries found\n";

  close(DATA);
  close (OUT);

  if (defined $old_size) {
    my $new_size = -s $outfile;
    if ($new_size < $old_size) {
      print "WARNING: the new ${species} file is smaller than the old file (old size: $old_size bytes, new_size: $new_size bytes).\n";
    } elsif ($new_size > $old_size) {
      print "The ${species} file has been updated (old size: $old_size bytes, new_size: $new_size bytes).\n";
    } else {
      print "The ${species} EST file has not changed - no updates to do.\n";
    }
  } else {
      print "The ${species} EST file has been created.\n";
  }
  print "\n";


# NEMBASE4
  my $nembase = $nembase_species_decode{$binomial_species_name};
  if (defined $nembase) {
    my $query = "http://www.nematodes.org/downloads/databases/NEMBASE4/${nembase}_nuc.fsa";
    system("wget -q -O ${target_dir}/Nembase '$query' ");
  } else {
    # make an empty file
    open(OUT, ">${target_dir}/Nembase");
    close OUT;
  }
  
# Nematode.net
  my $nematode_net = $nematode_net_species_decode{$binomial_species_name};

  if (defined $nematode_net) {
    my $query = "http://nematode.net/Data/cluster_ftp/SUMMARY/CLUSTER_SUMMARY.${nematode_net}";
    system("wget -q -O ${target_dir}/Nematode.net '$query' ");
  } else {
    # make an empty file
    open(OUT, ">${target_dir}/Nematode.net");
    close OUT;
  }

}

print "Finished.\n";
exit(0);






##############################################################
#
# Subroutines
#
##############################################################



##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################




# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - script_template.pl

=head1 USAGE

=over 4

=item script_template.pl  [-options]

=back

This script does...blah blah blah

script_template.pl MANDATORY arguments:

=over 4

=item None at present.

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item xxx (xxx@sanger.ac.uk)

=back

=cut
