#!/usr/local/bin/perl5.8.0 -w
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2007-06-01 10:07:05 $


use lib $ENV{'CVS_DIR'};

use strict;
use Wormbase;
use Getopt::Long;
use File::Copy;
use File::Path;
use Log_files;
use Storable;

my ($test, $database, $debug);
my ($mask, $dump, $run, $postprocess, $ace, $load, $process, $virtual);
my @types;
my $all;
my $store;
my ($species, $qspecies, $nematode);

GetOptions (
	    'debug:s'     => \$debug,
	    'test'        => \$test,
	    'database:s'  => \$database,
	    'store:s'     => \$store,
	    'species:s'   => \$species,  #target species (ie genome seq)
	    'mask'        => \$mask,
	    'dump'        => \$dump,
	    'process'     => \$process,
	    'virtual'     => \$virtual,
	    'run'         => \$run,
	    'postprocess' => \$postprocess,
	    'ace'         => \$ace,
	    'load'        => \$load,
	    'types:s'     => \@types,
	    'all'         => \$all,
	    'qspecies:s'  => \$qspecies,    #query species (ie cDNA seq)
	    'nematode'    => \$nematode
	   );

my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test     => $test,
			     -organism => $species
			   );
}

$species = $wormbase->species;#for load
my $log = Log_files->make_build_log($wormbase);

my $wormpub = $wormbase->wormpub;
$database = $wormbase->autoace unless $database;
my $blat_dir = $wormbase->blat;

#The mol_types available for each species is different
#defaults lists - can be overridden by -types
my %mol_types = ( 'elegans'   => [qw(ESTs mRNA ncRNA OSTs tc1 )],
				  'briggsae'  => [qw( mRNA ESTs )],
				  'remanei'   => [qw( mRNA ESTs )],
				  'brenneri'  => [qw( mRNA ESTs )],
				  'japonica'  => [qw( mRNA ESTs )]
				);

my @nematodes = qw(nematode washu nembase);

#remove other species if single one specified
if( $qspecies ){
	if( grep(/$qspecies/, keys %mol_types)) {
		foreach (keys %mol_types){
			delete $mol_types{$_} unless ($qspecies eq $_);
		}
	}
	else {
		$log->log_and_die("we only deal in Caenorhabditidae species!\n");
	}
}
	
#set specific mol_types if specified.
if(@types) {
	foreach (keys %mol_types){
		($mol_types{$_}) = @types;
	}
	@nematodes = ();
}

#only do the "other nematode" stuff
if($nematode) {
  foreach (keys %mol_types){
    delete $mol_types{$_};
  }
}

# mask the sequences based on Feature_data within the species database (or autoace for elegans.)
if( $mask ) {
	foreach my $species ( keys %mol_types ) {
		foreach my $moltype (@{$mol_types{$species}}) {
			$wormbase->bsub_script("transcriptmasker.pl -species $species -mol_type $moltype", $species, $log);
		}
	}
	
	#copy the nematode ESTs from BUILD_DATA
	foreach (@nematodes) {
	  mkdir ($wormbase->basedir."/cDNA/$_") unless  -e ($wormbase->basedir."/cDNA/$_");
	  copy($wormbase->build_data."/cDNA/$_/ESTs", $wormbase->basedir."/cDNA/$_/ESTs.masked");
	}
}

# Now make blat target database using autoace (will be needed for all possible blat jobs)

&dump_dna if $dump;


# run all blat jobs on cluster ( eg cbi4 )
if ( $run ) {
	#create Wormbase.pm objects for each species to access cdna directores
	my %accessors = $wormbase->species_accessors;
	
	# run other species
	foreach my $species (keys %accessors) { #doesn't include own species.
		my $cdna_dir = $accessors{$species}->maskedcdna;
		foreach my $moltype ( @{$mol_types{$accessors{$species}->species}} ) {
			my $split_count = 1;
			
			&check_and_shatter($accessors{$species}->maskedcdna, "$moltype.masked");
			foreach my $seq_file (glob ($accessors{$species}->maskedcdna."/$moltype.masked*")) {
				my $cmd = "bsub -J ".$accessors{$species}->pepdir_prefix."_$moltype \"/software/worm/bin/blat/blat -noHead -t=dnax -q=dnax ";
				$cmd .= $wormbase->genome_seq ." $seq_file ";
				$cmd .= $wormbase->blat."/${species}_${moltype}_${split_count}.psl\"";
				$wormbase->run_command($cmd, $log);
				$split_count++;
			}
		}
	}			
	
	#run own species.
	foreach my $moltype (@{$mol_types{$wormbase->species}} ){
		my $split_count = 1;
		my $seq_dir = $wormbase->maskedcdna;
		&check_and_shatter($wormbase->maskedcdna, "$moltype.masked");
		foreach my $seq_file (glob ($seq_dir."/$moltype.masked*")) {
			my $cmd = "bsub -J ".$wormbase->pepdir_prefix."_$moltype \"/software/worm/bin/blat/blat -noHead ";
			$cmd .= $wormbase->genome_seq ." $seq_file ";
			$cmd .= $wormbase->blat."/".$wormbase->species."_${moltype}_${split_count}.psl\"";
			$wormbase->run_command($cmd, $log);	
			$split_count++;	
		}		
	}
	#run other nematodes 
	foreach my $moltype (@nematodes ){
		my $split_count = 1;
		my $seq_dir = $wormbase->basedir."/cDNA/$moltype";
		&check_and_shatter($seq_dir, "ESTs.masked");
		foreach my $seq_file (glob ($seq_dir."/EST*")) {
			my $cmd = "bsub -J ".$wormbase->pepdir_prefix."_$moltype \"/software/worm/bin/blat/blat -noHead -q=dnax -t=dnax ";
			$cmd .= $wormbase->genome_seq ." $seq_file ";
			$cmd .= $wormbase->blat."/${moltype}_${split_count}.psl\"";
			$wormbase->run_command($cmd, $log);	
			$split_count++;	
		}		
	}	
	
}

if( $postprocess ) {
  # merge psl files and convert to ace format
  $log->write_to("merging PSL files \n");
  my $blat_dir = $wormbase->blat;
  my $species = $wormbase->species;
  foreach my $species (keys %mol_types) {
  	foreach my $moltype ( @{$mol_types{$species}}){
 	 $wormbase->run_command("cat $blat_dir/${species}_${moltype}_*   > $blat_dir/${species}_${moltype}_out.psl", $log); # /d causes compiler warning (?)
    }
  }
}

if ( $process or $virtual ) {
	foreach my $species (keys %mol_types) {
	  foreach my $type (@{$mol_types{$species}} ) {
    	#create virtual objects
    	$wormbase->run_script("blat2ace.pl -virtual -type $type -qspecies $species", $log) if $virtual;
    	$wormbase->run_script("blat2ace.pl -type $type -qspecies $species", $log) if $process;
     }
   }
}

if( $load ) {
  foreach my $type (@types){
    $log->write_to("loading BLAT data - $type\n");

    # virtual objs
    my $file =  "$blat_dir/virtual_objects.autoace.blat.$type.ace";
    $wormbase->load_to_database( $database, $file,"virtual_objects_$type");

    # Don't need to add confirmed introns from nematode data (because there are none!)
    unless ( ($type eq "nematode") || ($type eq "washu") || ($type eq "nembase") || ($type eq "tc1") || ($type eq "embl")|| ($type eq "ncrna") ) {
      $file = "$blat_dir/virtual_objects.autoace.ci.$type.ace"; 
      $wormbase->load_to_database($database, $file, "blat_confirmed_introns_$type");

      $file = "$blat_dir/autoace.good_introns.$type.ace";
      $wormbase->load_to_database($database, $file, "blat_good_introns_$type");
    }

    # BLAT results
    $file = "$blat_dir/autoace.blat.${species}_$type.ace";
    $wormbase->load_to_database($database, $file, "blat_${species}_${type}_data");
  }
}

$log->mail;
exit(0);


###############################################################################################################

sub check_and_shatter {
	my $dir = shift;
	my $file = shift;
	
	my $seq_count = qx(grep -c '>' $dir/$file);
	if( $seq_count > 10000) {
		$wormbase->run_script("shatter $dir/$file 5000 $dir/$file", $log);
		$wormbase->run_command("rm -f $dir/$file", $log);
	}
}

#############################################################################
# dump_dna                                                                  #
# gets data out of autoace/camace, runs tace query for chromosome DNA files #
# and chromosome link files.                                                #
#############################################################################

sub dump_dna {
  # this really just makes sure the list of sequence files to BLAT against is written. Used the seq files under the organism database.

  my %accessors = $wormbase->species_accessors;
  $accessors{$wormbase->species} = $wormbase;
  foreach my $species ( keys %accessors ) {
    # genome sequence dna files are either .dna or .fa
    my @files = glob($accessors{$species}->chromosomes."/*.dna");
    push(@files,glob($accessors{$species}->chromosomes."/*.fa"));

    open(GENOME,">".$accessors{$species}->autoace."/genome_seq") or $log->log_and_die("cant open genome sequence file".$accessors{$species}->autoace."/genome_seq: $!\n");
    foreach (@files){
      print GENOME "$_\n";
    }
    close GENOME;
  }
}
__END__

=pod

=head1 NAME - BLAT_controller.pl

=head2 USAGE

BLAT_controller.pl is a wrapper for the several scripts involved in running the various BLAT jobs in the WormBase build

BLAT_controller.pl  arguments:

=over 4

=item 

-debug user  =  debug mode to specify who gets log mail

=item 

-test        =  use the test build

=item 

-store string  =  load a previously serialised Wormbase object ( from Wormbase.pm ) to maintain test and debug setting when called from autoace_builder.pl

=item 

-database string = run on database other than build


 * script action options

=item 

-mask        runs transcriptmasker.pl  to mask ESTs polyA TSLs etc

=item 

-dump        dumps the target DNA sequences to

=item 

-process     runs blat_them_all.pl -process which runs blat2ace.pl to convert psl -> ace files

=item 

-virtual     runs blat_them_all.pl -virtual to create the virtual object the BLAT results hang off

=item 

-run         runs batch_BLAT.pl which shatters ESTs and nematode and submits bsub jobs to cbi1

=item 

-postprocess merges the psl files from the shattered EST and nematode files

=item 

-load       loads the BLAT acefiles in to the specified database

=item 

-types string  allows the specification of certain types to run BLAT. Valid types : est mrna ncrna ost nematode embl tc1 washu nembase

=item 

-all        runs all of the above tpyes.

=cut
