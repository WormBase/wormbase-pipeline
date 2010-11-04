#!/usr/local/bin/perl5.8.0 -w
#
# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2010-11-04 15:29:55 $


use lib $ENV{'CVS_DIR'};

use strict;
use Wormbase;
use Getopt::Long;
use File::Copy;
use File::Path;
use Log_files;
use Storable;
use Sequence_extract;
use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;

my ($test, $database, $debug);
my ($mask, $dump, $run, $postprocess, $load, $process, $virtual, $intron);
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
	    'load'        => \$load,
	    'types:s'     => \@types,
	    'all'         => \$all,
	    'qspecies:s'  => \$qspecies,    #query species (ie cDNA seq)
	    'nematode'    => \$nematode,
	    'intron'      => \$intron
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
my $seq_obj      = Sequence_extract->invoke($database, undef, $wormbase) if $intron;

#The mol_types available for each species is different
#defaults lists - can be overridden by -types

my %mol_types = ( 'elegans'          => [qw( EST mRNA ncRNA OST tc1 RST )],
		  'briggsae'         => [qw( mRNA EST )],
		  'remanei'          => [qw( mRNA EST )],
		  'brenneri'         => [qw( mRNA EST )],
		  'japonica'         => [qw( mRNA EST )],
		  'heterorhabditis'  => [qw( mRNA EST )],
		  'brugia'           => [qw( mRNA EST )],
		  'pristionchus'     => [qw( mRNA EST )],
		  'nematode'         => [qw( EST )],
		  'nembase'          => [qw( EST )],
		  'washu'            => [qw( EST )],

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
		$log->log_and_die("we only deal in certain species!\n");
	}
	@nematodes = ();
}
	
#set specific mol_types if specified.
if(@types) {
	foreach (keys %mol_types){
		($mol_types{$_}) = [(@types)];
	}
	@nematodes = ();
}

#only do the "other nematode" stuff
if($nematode) {
  foreach my $mt (keys %mol_types){
  	next if (grep /$mt/ ,@nematodes);
    delete $mol_types{$mt};
  }
}

# mask the sequences based on Feature_data within the species database (or autoace for elegans.)
if( $mask ) {
	foreach my $qspecies ( keys %mol_types ) {
		next if (grep /$qspecies/, @nematodes);
		foreach my $moltype (@{$mol_types{$qspecies}}) {
			#transcriptmasker is designed to be run on a sinlge species at a time.
			#therefore this uses the query species ($qspecies) as the -species parameter so that 
			#transcriptmasker creates a different Species opbject and gets the correct paths! 
			$wormbase->bsub_script("transcriptmasker.pl -species $qspecies -mol_type $moltype", $qspecies, $log);
		}
	}
	
	#copy the nematode ESTs from BUILD_DATA
	foreach (@nematodes) {
	  mkpath ($wormbase->basedir."/cDNA/$_") unless  -e ($wormbase->basedir."/cDNA/$_");
	  if (! -e $wormbase->build_data."/cDNA/$_/EST") {$wormbase->log_and_die("Can't find file " . $wormbase->build_data."/cDNA/$_/EST\n");};
	  copy($wormbase->build_data."/cDNA/$_/EST", $wormbase->basedir."/cDNA/$_/EST.masked") 
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
				my $cmd = "bsub -E \"ls /software/worm\" -J ".$accessors{$species}->pepdir_prefix."_$moltype \"/software/worm/bin/blat/blat -noHead -t=dnax -q=dnax ";
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
			my $cmd = "bsub -E \"ls /software/worm\" -o /dev/null -e /dev/null -J ".$wormbase->pepdir_prefix."_$moltype \"/software/worm/bin/blat/blat -noHead ";
			$cmd .= $wormbase->genome_seq ." $seq_file ";
			$cmd .= $wormbase->blat."/".$wormbase->species."_${moltype}_${split_count}.psl\"";
			$wormbase->run_command($cmd, $log);	
			$split_count++;	
		}		
	}
	#run other nematodes - these dont have Species.pm so cant use maskedcdna method.
	foreach my $moltype (@nematodes ){
		my $split_count = 1;
		my $seq_dir = $wormbase->basedir."/cDNA/$moltype";
		&check_and_shatter($seq_dir, "EST.masked");
		foreach my $seq_file (glob ($seq_dir."/EST*")) {
			my $cmd = "bsub -E \"ls /software/worm\" -J ".$wormbase->pepdir_prefix."_$moltype \"/software/worm/bin/blat/blat -noHead -q=dnax -t=dnax ";
			$cmd .= $wormbase->genome_seq ." $seq_file ";
			$cmd .= $wormbase->blat."/${moltype}_EST_${split_count}.psl\"";
			$wormbase->run_command($cmd, $log);	
			$split_count++;				
		}		
	}			
}

if( $postprocess ) {
  # merge psl files and convert to ace format
  $log->write_to("merging PSL files \n");
  my $blat_dir = $wormbase->blat;
  foreach my $species (keys %mol_types) {
  	foreach my $moltype ( @{$mol_types{$species}}){
 	 $wormbase->run_command("cat $blat_dir/${species}_${moltype}_* |sort -u  > $blat_dir/${species}_${moltype}_out.psl", $log); # /d causes compiler warning (?)
    }
  }
}

if ( $process or $virtual ) {

  # first run all the blat2ace.pl -virtual jobs
  my $lsf1 = LSF::JobManager->new();
  foreach my $species (keys %mol_types) {
    foreach my $type (@{$mol_types{$species}} ) {
      #create virtual objects
      $log->write_to("Submitting $species $type for virtual procesing\n");

      my $job_name = "worm_".$wormbase->species."_blat";
      my $cmd1 = $wormbase->build_cmd("blat2ace.pl -virtual -type $type -qspecies $species");
      # ask for a file size limit of 2 Gb and a memory limit of 4 Gb
      my @bsub_options = (-F => "2000000", 
			  -M => "4000000", 
			  -R => "\"select[mem>4000] rusage[mem=4000]\"",
			  -J => $job_name);
      $lsf1->submit(@bsub_options, $cmd1) if $virtual;
    	
    }
  }
  $lsf1->wait_all_children( history => 1 );
  $log->write_to("All blat2ace runs have completed!\n");
  for my $job ( $lsf1->jobs ) {    # much quicker if history is pre-cached
    $log->error("$job exited non zero\n") if $job->history->exit_status != 0;
  }
  $lsf1->clear;   

  # now run all the blat2ace.pl -intron jobs
  my $lsf2 = LSF::JobManager->new();
  foreach my $species (keys %mol_types) {
    foreach my $type (@{$mol_types{$species}} ) {
      #create virtual objects
      $log->write_to("Submitting $species $type for intron processing\n");
    	
      my $job_name = "worm_".$wormbase->species."_blat";
      my $cmd2 = $wormbase->build_cmd("blat2ace.pl -type $type -qspecies $species -intron");
      my @bsub_options = (-F => "2000000", 
			  -M => "4000000", 
			  -R => "\"select[mem>4000] rusage[mem=4000]\"",
			  -J => $job_name);
      $lsf2->submit(@bsub_options, $cmd2) if $process;

    }
  }
  $lsf2->wait_all_children( history => 1 );
  $log->write_to("All blat2ace runs have completed!\n");
  for my $job ( $lsf2->jobs ) {    # much quicker if history is pre-cached
    $log->error("$job exited non zero\n") if $job->history->exit_status != 0;
  }
  $lsf2->clear;   

  if($intron) {
    $log->write_to("confirming introns . . \n");
    #only do for self species matches
    foreach my $type (@{$mol_types{$wormbase->species}} ) {
      $log->write_to("\t$type\n");
      &confirm_introns($type);
    }
  }
}


if( $load ) {
  foreach my $species (keys %mol_types){
    $log->write_to("Loading $species BLAT data\n");
    foreach my $type (@{$mol_types{$species}}){
      $log->write_to("\tloading BLAT data - $type\n"); 
      # virtual objs
      my $file =  "$blat_dir/virtual_objects." . $wormbase->species . ".blat.$type.$species.ace";
      # skip the virtual_objects ace files that are not created
      if ($wormbase->species eq $species || -e $file) {
	$wormbase->load_to_database( $database, $file,"virtual_objects_$type", $log);
      } else {
	$log->write_to("\tskipping $file\n"); 
      }

      # Don't need to add confirmed introns from nematode data (because there are none!)
      unless ( ($type eq "nematode") || ($type eq "washu") || ($type eq "nembase") || ($type eq "tc1") || ($type eq "embl")|| ($type eq "ncrna") ) {
	$file = "$blat_dir/virtual_objects.".$wormbase->species.".ci.$type.$species.ace"; 
	$wormbase->load_to_database($database, $file, "blat_confirmed_introns_$type", $log);

	# only load in the good_intron files that have been created in confirm_introns()
        if ($species eq $wormbase->species) {
	  $file = "$blat_dir/".$wormbase->species.".good_introns.$type.ace";
	  $wormbase->load_to_database($database, $file, "blat_good_introns_$type", $log);
	}
      }

      # BLAT results
      $file = "$blat_dir/".$wormbase->species.".blat.${species}_$type.ace";
      $wormbase->load_to_database($database, $file, "blat_${species}_${type}_data", $log,1);
    }
  }
}


#confirm introns
sub confirm_introns {

  my $type = shift;
  local (*GOOD,*BAD);
  
  # open the output files
  open (GOOD, ">$blat_dir/".$wormbase->species.".good_introns.$type.ace") or die "$!";
  open (BAD,  ">$blat_dir/".$wormbase->species.".bad_introns.$type.ace")  or die "$!";
  
  my ($link,@introns,$dna,$switch);
 
  
  $/ = "";
  	
  # set qspecies to be just 'elegans' for now
  # is this meant to iterate over all species?
  # if so, then we need separate $blat_dir/".$wormbase->species.".{good,bad}_introns.$type.ace files for each species
  my $qspecies = $wormbase->species;
  	
  open (CI, "<$blat_dir/".$wormbase->species.".ci.${qspecies}_${type}.ace")  or $log->log_and_die("Cannot open $blat_dir/".$wormbase->species.".ci.${qspecies}_${type}.ace $!\n");
  while (<CI>) {
    next unless /^\S/;
    if (/Sequence : \"(\S+)\"/) {
      $link = $1;
      print "Sequence : $link\n" if $debug;
      @introns = split /\n/, $_;
       
      # evaluate introns #
      $/ = "";
      foreach my $test (@introns) {
	if ($test =~ /Confirmed_intron/) {
	    my @f = split / /, $test;
	  
	  #######################################
	  # get the donor and acceptor sequence #
	  #######################################
	  
	    my ($first,$last,$start,$end,$pastfirst,$prelast);
	    if ($f[1] < $f[2]) {
		($first,$last,$pastfirst,$prelast) = ($f[1]-1,$f[2]-1,$f[1],$f[2]-2);
	    }
	    else {
		($first,$last,$pastfirst,$prelast) = ($f[2]-1,$f[1]-1,$f[2],$f[1]-2);
	    }	
	    
	    $start = $seq_obj->Sub_sequence($link,$first,2);
	    $end   = $seq_obj->Sub_sequence($link,$prelast,2);
	  
	    print "Coords start $f[1] => $start, end $f[2] => $end\n" if $debug;
	  
	  ##################
	  # map to S_child #
	  ##################
	  
	  my $lastvirt = int((length $seq_obj->Sub_sequence($link)) /100000) + 1;
	  my ($startvirt,$endvirt,$virtual);
	  if ((int($first/100000) + 1 ) > $lastvirt) {
	    $startvirt = $lastvirt;
	  }
	  else {
	    $startvirt = int($first/100000) + 1;
	  }
	  if ((int($last/100000) + 1 ) > $lastvirt) {
	    $endvirt = $lastvirt;
	  }
	  else {
	    $endvirt = int($first/100000) + 1;
	  }
	  
	  if ($startvirt == $endvirt) { 
	    $virtual = "Confirmed_intron_$type:" .$link."_".$startvirt;
	   }
	  elsif (($startvirt == ($endvirt - 1)) && (($last%100000) <= 50000)) {
	    $virtual = "Confirmed_intron_$type:" .$link."_".$startvirt;
	  }
	  
	  #################
	  # check introns #
	  #################
	  
	    my $firstcalc = int($f[1]/100000);
	    my $seccalc   = int($f[2]/100000);
	    print STDERR "Problem with $test\n" unless (defined $firstcalc && defined $seccalc); 
	    my ($one,$two);
	    if ($firstcalc == $seccalc) {
		$one = $f[1]%100000;
		$two = $f[2]%100000;
	    }
	    elsif ($firstcalc == ($seccalc-1)) {
		$one = $f[1]%100000;
		$two = $f[2]%100000 + 100000;
		print STDERR "$virtual: $one $two\n";
	    }
	    elsif (($firstcalc-1) == $seccalc) {
		$one = $f[1]%100000 + 100000;
		$two = $f[2]%100000;
		print STDERR "$virtual: $one $two\n";
	    } 
	    print STDERR "Problem with $test\n" unless (defined $one && defined $two); 
	    
	    if ( ( (($start eq 'gt') || ($start eq 'gc')) && ($end eq 'ag')) ||
		 (  ($start eq 'ct') && (($end eq 'ac') || ($end eq 'gc')) ) ) {	 
		print GOOD "Feature_data : \"$virtual\"\n";
		
		# check to see intron length. If less than 25 bp then mark up as False
		# dl 040414
		
		if (abs ($one - $two) <= 25) {
		    print GOOD "Confirmed_intron $one $two False $f[4]\n\n";
		}
		else {
		  if ($type eq "mRNA"){
		    print GOOD "Confirmed_intron $one $two cDNA $f[4]\n\n";
		  }
		  else {
		    print GOOD "Confirmed_intron $one $two EST $f[4]\n\n";
		  }
		}
	      }
	    else {
	      if ($type eq "mRNA"){
		print BAD "Feature_data : \"$virtual\"\n";
		print BAD "Confirmed_intron $one $two cDNA $f[4]\n\n";		
	      }
	      else {
		print BAD "Feature_data : \"$virtual\"\n";
		print BAD "Confirmed_intron $one $two EST $f[4]\n\n";
	      }
	    }
	  }
      }
    }
  }
  close CI;
  
  close GOOD;
  close BAD;
  
}


$log->mail;
exit(0);


###############################################################################################################

sub check_and_shatter {
	my $dir = shift;
	my $file = shift;
	
	unless( -e "$dir/$file" ) {
		$log->write_to("$file doesnt exist - hopefully already shattered for other species\n");
		my @shatteredfiles = glob("$dir/$file*");
		if(scalar @shatteredfiles == 0){
			$log->log_and_die("shattered files also missing - not good");
		}
	}else {		
		my $seq_count = qx(grep -c '>' $dir/$file);
		if( $seq_count > 10000) {
			$wormbase->run_script("shatter $dir/$file 10000 $dir/$file", $log);
			$wormbase->run_command("rm -f $dir/$file", $log);
		}
	}
}

#############################################################################
# dump_dna                                                                  #
# gets data out of autoace/camace, runs tace query for chromosome DNA files #
# and chromosome link files.                                                #
#############################################################################

sub dump_dna {
    my @files = glob($wormbase->chromosomes."/*.dna");
    push(@files,glob($wormbase->chromosomes."/*.fa"));
    
    open(GENOME,">".$wormbase->autoace."/genome_seq") or $log->log_and_die("cant open genome sequence file".$wormbase->autoace."/genome_seq: $!\n");
    foreach (@files){
	next if (/supercontig/ && scalar @files>1); # don't touch this
	print GENOME "$_\n";
    }
    close GENOME;
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
