#!/usr/local/bin/perl -w
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2015-03-03 12:04:16 $      

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use File::Copy "cp";


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($homol, $misc, $species, $blat);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "homol"      => \$homol,
	    "misc"       => \$misc,
	    "species:s"   => \$species,
	    "blat"       => \$blat
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species
			     );
}
$species = $wormbase->species;
# establish log file.
my $log = Log_files->make_build_log($wormbase);

my @core_organisms = $wormbase->core_species;
if( $species eq 'elegans') {
  &parse_misc_elegans_files	if $misc;
  &parse_homol_data    		if $homol;
  &parse_elegans_homol_data     if $homol;
} else {
  if(grep(/$species/, map(lc $_,  @core_organisms))){  #other core (tierII) species)
    &parse_homol_data           if $homol;
    &parse_briggsae_data        if ($misc && $species eq 'briggsae');
    &parse_remanei_data         if ($misc && $species eq 'remanei');
    &parse_brugia_data          if ($misc && $species eq 'brugia');
    &parse_nematode_seqs        if $misc;
    &parse_genBlastG            if $homol;
  }
}
$log->mail();
exit(0);

sub parse_misc_elegans_files {
  my %files_to_load = (
		       $wormbase->misc_dynamic."/misc_genefinder.ace"           => "genefinder_predictions",
		       $wormbase->misc_dynamic."/misc_twinscan.ace"             => "twinscan_predictions"  ,
		       $wormbase->misc_dynamic."/misc_mgene.ace"                => "mgene_predictions"     ,
		       $wormbase->misc_dynamic."/misc_RNASEQ_CDS.ace"           => "RNASEQ_CDS_predictions",
		       $wormbase->misc_dynamic."/misc_modENCODE_aggregate_transcripts.ace" => "modENCODE_aggregate_transcripts",
		       $wormbase->misc_dynamic."/misc_modENCODE_Tiling_array_TARs.ace" => "modENCODE_Tiling_array_TARs",
		       $wormbase->misc_dynamic."/misc_jigsaw.ace"               => "jigsaw_predictions"    ,
		       $wormbase->misc_dynamic."/misc_TEC_RED_homol_data.ace"   => "TEC_RED"               ,
		       $wormbase->misc_dynamic."/misc_TEC_RED_homol.ace"        => "TEC_RED"               ,
		       $wormbase->misc_static."/orthologs.ensembl.ace"          => "ensembl_orthologs"     ,
		       $wormbase->misc_static."/orthologs.genes.ace"            => "ortholog_stubs"        ,
		       $wormbase->misc_static."/orthologs.inparanoid.ace"       => "inparanoid_orthologs"  ,
		       $wormbase->misc_static."/orthologs.panther.ace"          => "panther_orthologs"     ,
		       #$wormbase->misc_static."/orthologs.treefam.ace"          => "treefam_orthologs"     ,
		       $wormbase->misc_static."/misc_TEC_RED_sequence_data.ace" => "TEC_RED"               ,
		       $wormbase->misc_static."/nembase_nematode_contigs.ace"   => "nembase_ace"           ,
		       $wormbase->misc_static."/other_nematode_ESTs.ace"        => "other_nematode_ace"    ,
		       $wormbase->misc_static."/washu_nematode_contigs.ace"     => "washu_nem_ace"         ,
		       #$wormbase->misc_dynamic."/misc_mass_spec_GenniferMerrihew.ace"  => "mass_spec"      ,
		       $wormbase->misc_static."/GINUMBERS/GI_numbers.ace"       => "gi_number"             ,
		       $wormbase->misc_static.'/misc_mtce_protein_IDs.ace'      => 'mtce_protein_IDs'      ,
		       $wormbase->misc_dynamic.'/Caenorhabditae_sequence_data_to_load.ace'     => 'Caenorhabditae_seq_data',
		       $wormbase->misc_dynamic.'/fosmids.ace'                   => 'vancouver_fosmids'     ,
		       $wormbase->misc_dynamic.'/misc_21urna_homol.ace'         => '21uRNAs'               ,
		       $wormbase->misc_dynamic.'/misc_Expression_pattern_homol.ace'  => 'Expression_patterns',
		       $wormbase->misc_dynamic.'/misc_Tijsterman_G4.ace'        => 'Tijsterman_G4',
		       $wormbase->misc_dynamic.'/misc_Lamm_polysomes.ace'       => 'Lamm_polysomes',
		       #$wormbase->misc_dynamic.'/misc_RNASeq_hits_elegans.ace'  => 'RNASeq_hits',
		       #$wormbase->misc_dynamic.'/RNASeq_splice_elegans.ace'     => 'RNASeq_splice',
                       $wormbase->misc_static.'/homology_groups.ace' 	        => 'homol_groups',
                       $wormbase->misc_static.'/eggNOG.ace' 	                => 'eggNOG_groups',
		       $wormbase->misc_static.'/TF-modENCODE2014.ace'           => 'TF-modENCODE2014',
		      );

  $log->write_to("Loading files to ".$wormbase->autoace."\n==================================\n");
  foreach my $file (sort keys %files_to_load) {
    $log->write_to("\tloading $file -tsuser $files_to_load{$file}\n");
    $wormbase->load_to_database($wormbase->autoace,$file, $files_to_load{$file},$log);
  }
}

# Brugia malayi specific files
sub parse_brugia_data {
  my %files_to_load = (
                    $wormbase->misc_dynamic.'/misc_TIGR_brugia.ace' => 'TIGR_gene_models',
                    $wormbase->misc_static.'/brugia_bacs.ace.gz'    => 'Brugia BACS',
  );
  $log->write_to("Loading files to ".$wormbase->autoace."\n==================================\n");
  foreach my $file (sort keys %files_to_load) {
    $log->write_to("\tloading $file -tsuser $files_to_load{$file}\n");
    $wormbase->load_to_database($wormbase->autoace,$file, $files_to_load{$file},$log);
  }

}

sub parse_nematode_seqs {
  my %files2load = (
		    "nembase_nematode_contigs.ace"   => "nembase_ace"           ,
		    "other_nematode_ESTs.ace"        => "other_nematode_ace"    ,
		    "washu_nematode_contigs.ace"     => "washu_nem_ace"         ,
		   );
  $log->write_to("loading nematode sequences to $species database\n");
  foreach my $file ( sort keys %files2load ) {
    my $tsuser = $files2load{"$file"};
    $log->write_to("\tloading $file -tsuser -$tsuser\n");
    $wormbase->load_to_database($wormbase->autoace,$wormbase->misc_static."/$file",$tsuser, $log);
  }
}

sub parse_homol_data {

  my @files2Load = (
		    #BLAST data
		    "${species}_blastp.ace",
		    "${species}_blastx.ace",
		    #motif info
		    "worm_ensembl_${species}_motif_info.ace",
		    'interpro_motifs.ace',
		    #protein info
		    "worm_ensembl_${species}_interpro_motif_info.ace",
		    #other data
		    'repeat_homologies.ace',
		    'inverted_repeats.ace',
		   );

  $log->write_to("\nLoading homol data\n==============================\n");
  
  foreach my $file ( @files2Load ) {
    my $tsuser = substr($file,0,-4); #file name without ace
    $log->write_to("\tloading $file -tsuser -$tsuser\n");
    $wormbase->load_to_database($wormbase->autoace,$wormbase->acefiles."/$file",$tsuser, $log);
  }
}

sub parse_elegans_homol_data {

  my @files2Load = (
		    "ensembl_protein_info.ace",
		   );

  $log->write_to("\nLoading homol data\n==============================\n");
  
  foreach my $file ( @files2Load ) {
    my $tsuser = substr($file,0,-4); #file name without ace
    $log->write_to("\tloading $file -tsuser -$tsuser\n");
    $wormbase->load_to_database($wormbase->autoace,$wormbase->acefiles."/$file",$tsuser, $log);
  }


}


sub parse_remanei_data {

  # RNASeq_splice data
  #$wormbase->load_to_database($wormbase->autoace, $wormbase->misc_dynamic.'/RNASeq_splice_remanei.ace', 'RNASeq_splice_remanei', $log);
}


sub parse_briggsae_data {
  #
  # briggsae BAC end data
  #
  my $brigdb = $wormbase->database('briggsae');
  my $bac_end_dir = "BAC_ENDS";
  if (-d  "$brigdb/$bac_end_dir") {
    $log->write_to("\nLoading briggsae BAC ends from $brigdb/$bac_end_dir\n===========================\n");
    
    foreach my $be_file ("briggsae_BAC_ends.fasta",
                         "briggsae_BAC_end_homols.ace") {
      $log->write_to("\tLoading $be_file\n");
      $wormbase->load_to_database($wormbase->autoace,"$brigdb/$bac_end_dir/$be_file","BAC_ends", $log);
    }
  } else {
    $log->write_to("\tWARNING : $brigdb/$bac_end_dir dir not found; will NOT load BAC_END data\n");
  }


  # this could be in the parse_homol_data() routine but it is restricted to briggsae
  # load the briggsae TEC-RED acefiles/ homol data
  my @files = (
	    "misc_briggsae_TEC_RED_homol_data.ace",
	    "misc_briggsae_TEC_RED_homol.ace"
	   );
  foreach my $file (@files){
    $log->write_to("\tload $file\n");
    my $tsuser = substr($file,0,-4); #file name without ace
    $wormbase->load_to_database($wormbase->autoace,$wormbase->acefiles."/$file", $tsuser, $log);
  }

  # RNASeq_splice data
  #$wormbase->load_to_database($wormbase->autoace, $wormbase->misc_dynamic.'/RNASeq_splice_briggsae.ace', 'RNASeq_splice_briggsae', $log);

}

# genBlastG projections of elegans proteins onto other species' genomes
sub parse_genBlastG {

  if ($species ne 'elegans' and $species ne 'brugia' and $species ne 'ovolvulus') {
        $wormbase->load_to_database($wormbase->autoace, $wormbase->acefiles."/genblast.ace", "genBlastG", $log);
  }

}



__END__

=pod

=head2 NAME - script_template.pl

=head1 USAGE

=over 4

=item load_data_sets.pl  [-options]

=back

Loads lots of fairly static files that need to go in to each release.

script_template.pl MANDATORY arguments:

=over 4

=item None at present.

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -homol 

* load results of farm analyses

=item -misc

*loads static datasets like TEC-RED, gene_predictions, nematode ESTs etc.

=back

=item -brig

*loads data from briggsae eg BAC end and proteins

=back

=item -blat

*loads all of the BLAT data ( this should be done by BLAT but here just in case you need it eg after database corruption)

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

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
