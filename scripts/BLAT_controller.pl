#!/usr/local/bin/perl5.8.0 -w
#
# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2015-04-27 13:29:35 $


use lib $ENV{'CVS_DIR'};

use strict;
use Wormbase;
use Getopt::Long;
use File::Copy;
use File::Path;
use Log_files;
use Storable;
use Sequence_extract;
use Modules::WormSlurm;

my ($test, $database, $debug);
my ($mask, $dump_dna, $run, $postprocess, $load, $process, $intron, $blat_exe);
my @types;
my $store;
my ($species, $qspecies, $nematode, $no_backup_on_load, $min_coverage);

GetOptions (
  'debug:s'        => \$debug,
  'test'           => \$test,
  'database:s'     => \$database,
  'store:s'        => \$store,
  'species:s'      => \$species,  #target species (ie genome seq)
  'mask'           => \$mask,
  'dump'           => \$dump_dna,
  'process'        => \$process,
  'run'            => \$run,
  'postprocess'    => \$postprocess,
  'load'           => \$load,
  'nobackuponload' => \$no_backup_on_load,
  'types:s'        => \@types,
  'qspecies:s'     => \$qspecies,    #query species (ie cDNA seq)
  'intron'         => \$intron,
  'blatexe=s'      => \$blat_exe,
  'mincoverage=s'  => \$min_coverage,
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

$blat_exe = 'blat' unless $blat_exe;
$database    = $wormbase->autoace unless $database;
$species  = $wormbase->species;

my $log      = Log_files->make_build_log($wormbase);
my $wormpub  = $wormbase->wormpub;
my $blat_dir = $wormbase->blat;
my $lsfdir   = $wormbase->build_lsfout;
my $seq_obj  = Sequence_extract->invoke($database, undef, $wormbase) if $intron;

#The mol_types available for each species is different
#defaults lists - can be overridden by -types

my %mol_types = ( 'elegans'          => [qw( EST mRNA ncRNA OST tc1 RST Trinity Nanopore)],
		  'briggsae'         => [qw( mRNA EST Trinity)],
		  'remanei'          => [qw( mRNA EST)],
		  'brenneri'         => [qw( mRNA EST)],
		  'japonica'         => [qw( mRNA EST Trinity)],
		  'brugia'           => [qw( mRNA EST Trinity IsoSeq)],
		  'pristionchus'     => [qw( mRNA EST)],
                  'ovolvulus'        => [qw( mRNA EST Trinity)],
                  'sratti'           => [qw( mRNA EST)],
		  'tmuris'           => [qw( mRNA EST Trinity IsoSeq)],
		  'nematode'         => [qw( EST )],
		  'nembase'          => [qw( EST )],
		  'washu'            => [qw( EST )],

				);

my @other_nematodes = ('nematode',
                       'washu',
                       'nembase');

#remove other species if single one specified
if( $qspecies ){
  if( grep(/$qspecies/, keys %mol_types)) {
    foreach (keys %mol_types){
      delete $mol_types{$_} unless ($qspecies eq $_);
    }
  } else {
    $log->log_and_die("we only deal in certain species!\n");
  }
}
	
#set specific mol_types if specified.
if(@types) {
  foreach (keys %mol_types){
    ($mol_types{$_}) = [(@types)];
  }
}


# mask the sequences based on Feature_data within the species database (or autoace for elegans.)
my %slurm_jobs;
if( $mask ) {
  
    foreach my $moltype (@{$mol_types{$species}}) {
	my $cmd;
	if ($store) {
	    $cmd = $wormbase->build_cmd_line("transcriptmasker.pl -mol_type $moltype", $store);
	} else {
	    $cmd = $wormbase->build_cmd("transcriptmasker.pl -mol_type $moltype");
	}
	
	$cmd .= " -test" if $test;
	$cmd .= " -debug $debug" if $debug;
	
	my $job_id = WormSlurm::submit_job_with_name($cmd, 'production', '1500m', '01:00:00', "$lsfdir/BLAT_mask_${species}_${moltype}.slurmout", "$lsfdir/BLAT_mask_${species}_${moltype}.slurmerr", "BLAT_mask_${species}_${moltype}");
	$slurm_jobs{$job_id} = $cmd;
    }
    WormSlurm::wait_for_jobs(keys %slurm_jobs);

    $log->write_to("All transcriptmasker runs have completed!\n");
    for my $job_id ( keys %slurm_jobs ) {
	$log->error("Slurm job $job_id (" . $slurm_jobs{$job_id} . ") exited non zero\n") if WormSlurm::get_exit_code($job_id) != 0;
    }


    if ($log->report_errors) {
	$log->log_and_die("Some transcriptmasker jobs failed; bailing out now\n");
    }

    for my $moltype ( @{$mol_types{$species}} ) {
	&check_and_shatter( $wormbase->maskedcdna, "$moltype.masked" );
    }

    # finally, update the other_nematode collection, only for elegans
    if ($species eq 'elegans') {
	foreach my $on (@other_nematodes) {
	    my $source_path = $wormbase->build_data . "/cDNA/$on";
	    my $target_path = $wormbase->basedir .  "/cDNA/$on";
      
	    foreach my $mt (@{$mol_types{$on}}) {
		foreach my $file (glob("$target_path/$mt.masked*")) {
		    unlink $file;
		}
        
		my $sfile = "$source_path/$mt";
		my $tfile = "$target_path/${mt}.masked";
        
		if (not -e $sfile) {
		    $log->error("Could not find file $sfile in BLAT preparation\n");
		}
		$wormbase->run_command("cp $sfile $tfile", $log);
        
		&check_and_shatter( $target_path, "${mt}.masked" );
	    }
	}
    }
}

# Now make blat target database using autoace (will be needed for all possible blat jobs)

&dump_dna if $dump_dna;


# run all blat jobs on cluster ( eg cbi4 )
if ( $run ) {
    
    my %slurm_blat_jobs;
    # run other species
    foreach my $qs (keys %mol_types) {
	# skip own species
	next if $qs eq $species;
	
	my $seq_dir = join("/", $wormbase->basedir, "cDNA", $qs );
	
	foreach my $moltype( @{$mol_types{$qs}} ) {
	    
	    &check_and_shatter( $seq_dir, "$moltype.masked" );
	    
	    foreach my $seq_file (glob("${seq_dir}/${moltype}.masked*")) {
		my $chunk_num = 1; 
		if ($seq_file =~ /_(\d+)$/) {
		    $chunk_num = $1;
		}
		
		my $cmd = "$blat_exe -noHead -t=dnax -q=dnax ";
		$cmd .= $wormbase->genome_seq ." $seq_file ";
		$cmd .= $wormbase->blat."/${qs}_${moltype}_${chunk_num}.psl";
		
		my $job_id = WormSlurm::submit_job_with_name($cmd, 'production', '1g', '1:00:00', "$lsfdir/BLAT_run_${species}_${qs}_${moltype}_${chunk_num}.slurmout", "$lsfdir/BLAT_run_${species}_${qs}_${moltype}_${chunk_num}.slurmerr", "BLAT_run_${species}_${qs}_${moltype}_${chunk_num}");
		$slurm_blat_jobs{$job_id} = $cmd;
	    }
	}
    }
    

    # run own species.
    foreach my $moltype (@{$mol_types{$species}} ){
	my $seq_dir = $wormbase->maskedcdna;
	&check_and_shatter($wormbase->maskedcdna, "$moltype.masked");
	foreach my $seq_file (glob ($seq_dir."/$moltype.masked*")) {
	    my $chunk_num = 1; 
	    if ($seq_file =~ /_(\d+)$/) {
		$chunk_num = $1;
	    }
	    my $cmd = "$blat_exe -noHead ";
	    $cmd .= $wormbase->genome_seq ." $seq_file ";
	    $cmd .= $wormbase->blat."/${species}_${moltype}_${chunk_num}.psl";


	    my $job_id = WormSlurm::submit_job_with_name($cmd, 'production', '1g', '1:00:00', "$lsfdir/BLAT_run_${species}_${moltype}_${chunk_num}.slurmout", "$lsfdir/BLAT_run_${species}_${moltype}_${chunk_num}.lsferr", "BLAT_run_${species}_${species}_${moltype}_${chunk_num}");
	    $slurm_blat_jobs{$job_id} = $cmd;
	}
    }

    WormSlurm::wait_for_jobs(keys %slurm_blat_jobs);

    $log->write_to("All BLAT runs have completed!\n");
    for my $job_id ( keys %slurm_blat_jobs ) {    # much quicker if history is pre-cached
	$log->error("Slurm job $job_id (" . $slurm_blat_jobs{$job_id} . ") exited non zero\n") if WormSlurm::get_exit_code($job_id) != 0;
    }

    if ($log->report_errors) {
	$log->log_and_die("Some BLAT run jobs failed - bailing\n");
    }
}

if( $postprocess ) {
    # merge psl files and convert to ace format
    $log->write_to("merging PSL files \n");
    my $blat_dir = $wormbase->blat;
    foreach my $species (keys %mol_types) {
	foreach my $moltype ( @{$mol_types{$species}}){
	    $wormbase->run_command("cat $blat_dir/${species}_${moltype}_* |sort -u  > $blat_dir/${species}_${moltype}_out.psl", $log);
	}
    }
}

if ( $process ) {

    my %slurm_process_jobs;
    foreach my $qspecies (keys %mol_types) {
	foreach my $type (@{$mol_types{$qspecies}} ) {
	    #create virtual objects
	    $log->write_to("Submitting $qspecies $type for virtual procesing\n");
      
	    my $cmd = "blat2ace.pl -groupaligns -type $type -qspecies $qspecies";
	    if ($qspecies eq $species) {
		if ($intron) {
		    $cmd .= " -intron";
		} 
	    } else {
		if ($min_coverage) {
		    $cmd .= " -mincoverage $min_coverage";
		}
	    }
	    $cmd = $wormbase->build_cmd($cmd);

	    if ($test) {
		$cmd .= " -test";
	    }
	    if ($debug) {
		$cmd .= " -debug $debug";
	    }

	    my $job_name = "BLAT_blat2ace_${species}_${qspecies}_${type}";
	    my $mem = '8g';
	    if ($type eq 'Nanopore' && $qspecies ne $species) {
		$mem = '24g';
		$cmd .= ' -mincoverage 90' unless $min_coverage;
	    }
	    
	    my $job_id = WormSlurm::submit_jobs_with_name($cmd, 'production', $mem, '1:00:00', "$lsfdir/$job_name.slurmout", "$lsfdir/$job_name.slurmerr", $job_name);
	    $slurm_process_jobs{$job_id} = $cmd;
	}
    }
    WormSlurm::wait_for_jobs(keys %slurm_process_jobs);
    $log->write_to("All blat2ace runs have completed!\n");
    for my $job_id (keys %slurm_process_jobs) {
	$log->error("Slurm job $job_id (" . $slurm_process_jobs{$job_id} . ") exited non zero\n") if WormSlurm::get_exit_code($job_id) != 0;
    }
  
    if ($log->report_errors) {
	$log->log_and_die("Some blat2ace jobs died - bailing\n");
    }

    if($intron) {
	$log->write_to("confirming introns . . \n");
	#only do for self species matches
	foreach my $type (@{$mol_types{$species}} ) {
	    my $virt_hash = &confirm_introns($wormbase->species, $type);
	    my $vfile = "$blat_dir/virtual_objects.${species}.ci.${type}.${species}.ace";
	    open(my $vfh, ">$vfile") or $log->log_and_die("Could not open $vfile for writing\n");

	    foreach my $tname (keys %$virt_hash) {
		my @schild = sort { 
		    my ($na) = ($a =~ /_(\d+)$/ or 0); my ($nb) = ($b =~ /_(\d+)$/ or 0); $na <=> $nb;
		} keys %{$virt_hash->{$tname}};
		
		print $vfh "\nSequence : \"$tname\"\n";
		foreach my $child (@schild) {
		    printf $vfh "S_Child Feature_data %s %d %d\n", $child, @{$virt_hash->{$tname}->{$child}};
		}
	    }
	    close($vfh);
	}
    }
}


if( $load and not $log->report_errors) {
    foreach my $qspecies (keys %mol_types){
	$log->write_to("Loading $qspecies BLAT data\n");
	foreach my $type (@{$mol_types{$qspecies}}){
	    $log->write_to("\tloading BLAT data - $type\n"); 
	    
	    # virtual objs
	    my $file =  "$blat_dir/virtual_objects.$species.blat.$type.$qspecies.ace";
	    if (-e $file) {
		$wormbase->load_to_database( $database, $file,"virtual_objects_$type", $log);
	    } else {
		$log->write_to("\tskipping $file\n"); 
	    }
	    
	    # confirmed introns - will only be any for within-species alignments
	    if ($qspecies eq $species) {
		
		my ($in_v_tag, $in_v_file) = ("blat_introns_virtual$type",  
					      "$blat_dir/virtual_objects.".$wormbase->species.".ci.$type.$qspecies.ace");
		my ($in_tag, $in_file) = ("blat_good_introns_${type}", "$blat_dir/$species.good_introns.$type.ace");
		
		if (-e $in_v_file and -e $in_file) {
		    $wormbase->load_to_database($database, $in_v_file, $in_v_tag, $log);
		    $wormbase->load_to_database($database, $in_file, $in_tag, $log);
		} else {
		    $log->write_to("\tskipping ci file(s) for $type because one or both is missing\n"); 
		}
	    }
	    
	    # BLAT results
	    $file = "$blat_dir/$species.blat.${qspecies}_$type.ace";
	    $wormbase->load_to_database($database, 
					$file, 
					"blat_${qspecies}_${type}_data", 
					$log, 
					($test or $no_backup_on_load) ? 0 : 1);
	}
    }
}


#confirm introns
sub confirm_introns {
    my ($qspecies, $type) = @_;
  
    # open the output files
    open (my $good_fh, ">$blat_dir/$species.good_introns.$type.ace") or die "$!";
    open (my $bad_fh,  ">$blat_dir/$species.bad_introns.$type.ace")  or die "$!";
  
    my ($link,@introns, %seqlength, %virtuals);
   
    $/ = "";
  	
    open (CI, "<$blat_dir/$species.ci.${qspecies}_${type}.ace")  
	or $log->log_and_die("Cannot open $blat_dir/$species.ci.${qspecies}_${type}.ace $!\n");

    while (<CI>) {
	next unless /^\S/;
	if (/Sequence : \"(\S+)\"/) {
	    $link = $1;
	    
	    if (not exists $seqlength{$link}) {
		$seqlength{$link} = length($seq_obj->Sub_sequence($link));
	    }
	    
	    @introns = split /\n/, $_;
	    
	    # evaluate introns #
	    $/ = "";
	    foreach my $test (@introns) {
		if ($test =~ /Confirmed_intron/) {
		    my @f = split /\s+/, $test;
		    
		    #######################################
		    # get the donor and acceptor sequence #
		    #######################################
	            
		    my ($tstart, $tend, $strand) = ($f[1], $f[2], 1);
		    if ($tend < $tstart) {
			($tstart, $tend, $strand) = ($tend, $tstart, -1);
		    }
		    
		    
		    my $start_splice = $seq_obj->Sub_sequence($link,$tstart - 1, 2);
		    my $end_splice   = $seq_obj->Sub_sequence($link,$tend - 2,2);
		    
		    print "Coords start => $tstart, end $tend\n" if $debug;
		    
		    ##################
		    # map to S_child #
		    ##################
		    my $binsize = 150000;
		    
		    my $bin = 1 +  int( $tstart / $binsize );
		    my $bin_start = ($bin - 1) * $binsize + 1;
		    my $bin_end   = $bin_start + $binsize - 1;
		    
		    if ($bin_end > $seqlength{$link}) {
			$bin_end = $seqlength{$link};
		    }
		    my $bin_of_end = 1 +  int( $tend / $binsize );
		    
		    # propagate rule from old code: if feature spans more than 2 bins, junk it
		    if ($bin != $bin_of_end and $bin != $bin_of_end - 1) {
			next;
		    }
		    
		    if ($tend > $bin_end) {
			$bin = 0;
			$bin_start = 1;
			$bin_end   = $seqlength{$link};
		    } else {
			$tstart = $tstart - $bin_start + 1;
			$tend   = $tend - $bin_start + 1;
		    }
		    
		    if ($strand < 0) {
			($tstart, $tend) = ($tend, $tstart);
		    }
		    
		    my $virtual = sprintf("Confirmed_intron_%s:%s%s", 
					  $type, 
					  $link, 
					  $bin ? "_$bin" : "");
		    
		    if (not exists $virtuals{$link}->{$virtual}) {
			$virtuals{$link}->{$virtual} = [$bin_start, $bin_end];
		    }
		    
		    if ( ( (($start_splice eq 'gt') || ($start_splice eq 'gc')) && ($end_splice eq 'ag')) ||
			 (  ($start_splice eq 'ct') && (($end_splice eq 'ac') || ($end_splice eq 'gc')) ) ) {	 
			
			print $good_fh "Feature_data : \"$virtual\"\n";
			
			# check to see intron length. If less than 25 bp then mark up as False
			# dl 040414
			
			if (abs($tend - $tstart) <= 25) {
			    print $good_fh "Confirmed_intron $tstart $tend False $f[4]\n\n";
			} else {
			    if ($type eq "mRNA"){
				print $good_fh "Confirmed_intron $tstart $tend cDNA $f[4]\n\n";
			    }
			    else {
				print $good_fh "Confirmed_intron $tstart $tend EST $f[4]\n\n";
			    }
			}
		    }
		    else {
			if ($type eq "mRNA"){
			    print $bad_fh "Feature_data : \"$virtual\"\n";
			    print $bad_fh "Confirmed_intron $tstart $tend cDNA $f[4]\n\n";		
			}
			else {
			    print $bad_fh "Feature_data : \"$virtual\"\n";
			    print $bad_fh "Confirmed_intron $tstart $tend EST $f[4]\n\n";
			}
		    }
		}
	    }
	}
    }
    close(CI);
    close($good_fh);
    close($bad_fh);
    
    return \%virtuals;
}


$log->mail;
exit(0);


###############################################################################################################

sub check_and_shatter {
    my $dir = shift;
    my $file = shift;
    
    unless( -e "$dir/$file" ) {
	my @shatteredfiles = glob("$dir/$file*");
	if(scalar @shatteredfiles == 0){
	    $log->log_and_die("shattered files $dir/$file * also missing - not good");
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

