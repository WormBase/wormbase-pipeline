#!/usr/local/bin/perl -w
# Last updated by: $Author: ar2 $     
# Last updated on: $Date: 2006-02-27 21:59:46 $      

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($homol, $misc, $brig, $blat);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "homol"      => \$homol,
	    "misc"       => \$misc,
	    "briggsae"   => \$brig,
	    "blat"       => \$blat
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

&parse_misc_files    if $misc;
&parse_homol_data    if $homol;
&parse_briggsae_data if $brig;
&parse_blat_data     if $blat;

$log->mail();
exit(0);

sub parse_misc_files {
  my @nem_files = ("other_nematode_ESTs.ace", "nembase_nematode_contigs.ace", "washu_nematode_contigs.ace");
  foreach (@nem_files) {
    $wormbase->run_command("cp ".$wormbase->wormpub."/analysis/ESTs/$_ ".$wormbase->acefiles."/") unless -e ($wormbase->acefiles."/$_");
  }
  my %files_to_load = (
#		       $wormbase->misc_dynamic."/misc_genefinder.ace"           => "genefinder_predictions",
#		       $wormbase->misc_dynamic."/misc_twinscan.ace"             => "twinscan_predictions"  ,
#		       $wormbase->misc_dynamic."/misc_TEC_RED_homol_data.ace"   => "TEC_RED"               ,
#		       $wormbase->misc_dynamic."/misc_TEC_RED_homol.ace"        => "TEC_RED"               ,
#		       $wormbase->misc_dynamic."/WS145_refseq.ace"              =>  'refseq_IDs'           ,
#		       $wormbase->misc_dynamic."/WS145_aceview.ace"             =>  'aceview_IDs'          ,
#		       $wormbase->misc_static."/ortholog_WS131.ace"             => "briggsae_orthologs"    ,
		       $wormbase->acefiles."/other_nematode_ESTs.ace"           => 'nematode_ESTs'         ,
		       $wormbase->acefiles."/nembase_nematode_contigs.ace"      => 'nembase_ESTs'          ,
		       $wormbase->acefiles."/washu_nematode_contigs.ace"        => 'washu_ESTs'            ,
#		       $wormbase->wormpub."/analysis/GI_numbers/GI_numbers.ace" => "gi_number"             ,
		      );

  $log->write_to("Loading files to ".$wormbase->autoace."\n==================================\n");
  foreach my $file (keys %files_to_load) {
    $log->write_to("\tloading $file -tsuser $files_to_load{$file}\n");
    $wormbase->load_to_database($wormbase->autoace,$file, $files_to_load{$file},$log);
  }
}

sub parse_homol_data {

  my @files2Load = (
		    #BLAST data
		    "worm_pep_blastp.ace",
		    "worm_brigpep_blastp.ace",
		    "worm_dna_blastx.ace",
		    #motif info
		    "worm_pep_motif_info.ace",
		    "worm_brigpep_motif_info.ace",
		    #protein info
		    "ensembl_protein_info.ace",
		    "worm_pep_interpro_motif_info.ace",
		    "worm_brigpep_interpro_motif_info.ace",
		    #other data
		    "repeat_homologies.ace",
		    "waba.ace",
		    "TRF.ace",
		    "inverted_repeats.ace"
		   );


  $log->write_to("\nLoading homol data\n==============================\n");
  foreach my $file ( @files2Load ) {
    my $tsuser = substr($file,0,-4); #file name without ace
    $log->write_to("\tloading $file -tsuser -$tsuser\n");
    $wormbase->load_to_database($wormbase->autoace,$wormbase->acefiles."/$file",$tsuser, $log);
  }
}

sub parse_briggsae_data {

  # briggsae BAC end data
  my @files = ("briggsae_BAC_ends.fasta",
	       "briggsae_homol_data.ace",
	       "briggsae_BAC_ends_data.ace",
	       "briggsae_bac_clone_ends.ace",
	       "bac_ends_unique.ace"
	      );

  my $brig_dir = $wormbase->database('brigace')."/BAC_ENDS";
  $log->write_to("\nLoading briggsae BAC ends from $brig_dir\n===========================\n");
  foreach my $file (@files){
    $log->write_to("\tload $file\n");
    $wormbase->load_to_database($wormbase->autoace,"$brig_dir/$file","BAC_ends", $log);
  }
  # and the brigpep file
  $wormbase->load_to_database($wormbase->autoace, $wormbase->database('brigace')."/brigpep.ace","brigpep", $log);
}

sub parse_blat_data {
  $log->write_to("loading BLAT data\n");
  my @files = (
	       'autoace.blat.embl.ace',	       'autoace.blat.est.ace',
	       'autoace.blat.mrna.ace',	       'autoace.blat.ncrna.ace',
	       'autoace.blat.nematode.ace',    'autoace.blat.nembase.ace',
	       'autoace.blat.ost.ace',	       'autoace.blat.tc1.ace',
	       'autoace.blat.washu.ace',       'autoace.ci.est.ace',
	       'autoace.ci.mrna.ace',	       'autoace.ci.ost.ace',
	       'autoace.good_introns.est.ace',
	       'autoace.good_introns.mrna.ace',
	       'autoace.good_introns.ost.ace',
	       'virtual_objects.autoace.blat.embl.ace',
	       'virtual_objects.autoace.blat.est.ace',
	       'virtual_objects.autoace.blat.mrna.ace',
	       'virtual_objects.autoace.blat.ncrna.ace',
	       'virtual_objects.autoace.blat.nematode.ace',
	       'virtual_objects.autoace.blat.nembase.ace',
	       'virtual_objects.autoace.blat.ost.ace',
	       'virtual_objects.autoace.blat.tc1.ace',
	       'virtual_objects.autoace.blat.washu.ace',
	       'virtual_objects.autoace.ci.est.ace',
	       'virtual_objects.autoace.ci.mrna.ace',
	       'virtual_objects.autoace.ci.ost.ace'
	      );

 foreach my $file (@files){
    $log->write_to("\tload $file\n");
    my $db = $wormbase->autoace;
    $wormbase->load_to_database($db,$wormbase->blat."/$file",undef, $log);
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
