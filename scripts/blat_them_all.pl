#!/usr/local/bin/perl5.8.0 -w
#
# blat_them_all.pl
# 
# by Kerstin Jekosch
#
# Gets sequences ready for blatting, blats sequences, processes blat output, makes confirmed introns
# and virtual objects to hang the data onto
#
# Last edited by: $Author: dl1 $
# Last edited on: $Date: 2004-04-14 14:50:03 $


use strict;
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};           
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Carp;


##############################
# Misc variables and paths   #
##############################

my ($help, $debug, $verbose, $est, $mrna, $ost, $nematode, $embl, 
    $blat, $tc1, $process, $virtual, $dump, $camace, $fine);
my $maintainers = "All";
my $errors      = 0;
my $bin         = "/wormsrv2/scripts";
our $log;
our $blat_dir   = "/wormsrv2/autoace/BLAT";    # default BLAT directory, can get changed if -camace used
our $dbpath     = "/wormsrv2/autoace";         # default database location
our %homedb;                                   # for storing superlink->lab connections
our $blatex     = '/nfs/disk100/wormpub/bin.ALPHA/blat';
our $giface     = &giface;
our %word = (
	     est      => 'BLAT_EST',
	     mrna     => 'BLAT_mRNA',
	     embl     => 'BLAT_EMBL',
	     nematode => 'BLAT_NEMATODE',
	     ost      => 'BLAT_OST',
	     tc1      => 'BLAT_TC1',
	     );


#########################
# Command line options  #
#########################

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "verbose"    => \$verbose,
	    "est"        => \$est,
	    "mrna"       => \$mrna,
	    "ost"        => \$ost,
	    "nematode"   => \$nematode,
	    "embl"       => \$embl,
	    "tc1"        => \$tc1,
	    "dump"       => \$dump,
	    "blat"       => \$blat,
	    "process"    => \$process,
	    "virtual"    => \$virtual,
	    "camace"     => \$camace,
	    "fine"       => \$fine
	    );


########################################
# command-line options & ramifications #
########################################


# set to camace or autoace
$blat_dir = "/wormsrv1/camace/BLAT" if $camace;
$dbpath   = "/wormsrv1/camace"      if $camace;
our $seq  = "$blat_dir/autoace.fa";               

# Help pod documentation
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}


# Exit if no dump/process/blat/virtual process option chosen
&usage(1) unless ($dump || $process || $blat || $virtual); 

# Exit if no data type choosen [EST|mRNA|EMBL|NEMATODE|OST] (or -dump not chosen)
&usage(2) unless ($est || $mrna || $embl || $tc1 || $nematode || $ost || $dump); 

# Exit if multiple data types choosen [EST|mRNA|EMBL|NEMATODE|OST]
# ignore if -dump is being run
unless($dump){
  my $flags = 0;
  $flags++ if $est;
  $flags++ if $ost;
  $flags++ if $mrna;
  $flags++ if $embl;
  $flags++ if $tc1;
  $flags++ if $nematode;
  &usage(3) if ($flags > 1);
}

# Exit if -fine option is being used without -mrna
&usage(4) if ($fine && !$mrna);


# assign type variable
my $data;
($data = 'est')      if ($est);
($data = 'ost')      if ($ost);
($data = 'mrna')     if ($mrna);
($data = 'embl')     if ($embl);
($data = 'tc1')      if ($tc1);
($data = 'nematode') if ($nematode);


# Select the correct set of query sequences for blat
my $query = "/nfs/disk100/wormpub/analysis/ESTs/";
$query   .= 'elegans_ESTs.masked'  if ($est);      # EST data set
$query   .= 'elegans_OSTs'         if ($ost);      # OST data set
$query   .= 'elegans_TC1s'         if ($tc1);      # TC1 data set
$query   .= 'elegans_mRNAs.masked' if ($mrna);     # mRNA data set
$query   .= 'other_nematode_ESTs'  if ($nematode); # ParaNem EST data set
$query   .= 'elegans_embl_cds'     if ($embl);     # Other CDS data set, DNA not peptide!

&create_log_files;




###########################################################################
#                                                                         #
#                        MAIN PART OF PROGRAM                             #
#                                                                         #
###########################################################################



# Write sequence data (chromsomes) from autoace/camace
&dump_dna if ($dump);
    

# BLAT the query data type 
if ($blat) {
  
  &usage(5) unless (-e "$seq");   # Check that autoace.fa exists
  
  # Check that autoace.fa file was created prior to start of (re)build 
  &usage(6) if ( (-M "/wormsrv2/autoace/logs/A1:Build_in_progress" < -M "${blat_dir}/autoace.fa") && (!$camace) );

  my $runtime = &runtime;   
  print LOG "$runtime: Running blat and putting the results in $blat_dir/${data}_out.psl\n";
  

  # BLAT system call
  # nematode ESTs handled differently
  # use -fine for mrna sequences is specified
  if ($nematode) {
    system("$blatex $seq $query -t=dnax -q=dnax $blat_dir/${data}_out.psl") && croak "Blat failed $!\n";
  }
  elsif($fine) {
    system("$blatex -fine $seq $query $blat_dir/${data}_out.psl") && croak "Blat failed $!\n";
  }
  else{
    system("$blatex $seq $query $blat_dir/${data}_out.psl") && croak "Blat failed $!\n";
  }
}



####################################################################
# process blat output (*.psl) file and convert results to acefile
# by running blat2ace
# Then produce confirmed introns files
#####################################################################

if ($process) {

    my $runtime = &runtime;
    print "Mapping blat data to autoace\n" if ($verbose);      
    print LOG "$runtime: Processing blat ouput file, running blat2ace.pl\n";
    
    # treat slightly different for nematode data (no confirmed introns needed)
    if ( ($nematode) || ($tc1) ) {
	&run_command("$bin/blat2ace.pl -$data"); 
    }
    elsif($camace){
	&run_command("/nfs/team71/worm/dl1/wormbase/wormbase/scripts/blat2ace.pl -$data -intron -camace"); 	
    }
    else {
	&run_command("$bin/blat2ace.pl -$data -intron"); 
    }

    $runtime = &runtime;
    print "Producing confirmed introns in databases\n\n" if $verbose;
    print LOG "$runtime: Producing confirmed introns in databases\n";

    # produce confirmed introns for all but nematode and tc1 data
    unless ( ($nematode) || ($tc1) ) {
	print "Producing confirmed introns using $data data\n" if $verbose;
	&confirm_introns('autoace',"$data");
	&confirm_introns('camace', "$data");
	&confirm_introns('stlace', "$data");
    }
}


#########################################
# produce files for the virtual objects #
#########################################
if ($virtual) {
    my $runtime = &runtime;
    print LOG "$runtime: Producing $data files for the virtual objects\n";
    
    print "// Assign laboratories to superlinks*\n";
    # First assign laboratories to each superlink object (stores in %homedb)
    &sequence_to_lab;
    &virtual_objects_blat($data);
}


##############################
# Clean up and say goodbye   #
##############################

close(LOG);

# send log
# warn about errors in subject line if there were any
#if($errors == 0){
#  &mail_maintainer("BUILD SCRIPT: blat_them_all",$maintainers,$log);
#}
#else{
#  &mail_maintainer("BUILD SCRIPT: blat_them_all: $errors ERROR!",$maintainers,$log);
#}


exit(0);





#################################################################################
#                                                                               #
#                     T H E    S U B R O U T I N E S                            #
#                                                                               #
#################################################################################

sub sequence_to_lab {
  # Connect superlink objects to their corresponding laboratory object
  # store in global %homedb
  local (*LINK);
  my $name;
  
  print LOG "Assign LINK* objects to laboratory\n\n";
  # deal with the superlink objects
  open (LINK, "<$blat_dir/superlinks.ace") || croak "Couldn't open superlinks.ace $!";
  while (<LINK>) {
      
      if (/^Sequence\s+\:\s+\"(\S+)\"/) {
	  $name = $1;
#	  print "// New sequence is $name\n";
	  next;
      }
      if (/^From_laboratory\s+\"(\S+)\"/) {
	  $homedb{$name} = $1;
	  print LOG "assigning $1 to $name\n";
#	  print "assigning $1 to $name\n";
	  undef ($name);
	  next;
      }
  }
  close(LINK);
  
  print LOG "\n";
  
}


#############################################################################
# dump_dna                                                                  #
# gets data out of autoace/camace, runs tace query for chromosome DNA files #
# and chromosome link files.                                                #
#############################################################################

sub dump_dna {

  local (*CHANGE,*NEW);

  my $command;

  unless ($camace) {
    $command  = "query find Sequence \"CHROMOSOME*\"\n";
    $command .= "show -a -f /wormsrv2/autoace/BLAT/chromosome.ace\n";
    $command .= "follow Subsequence\n";
    $command .= "show -a -f /wormsrv2/autoace/BLAT/superlinks.ace\n";
    $command .= "dna -f /wormsrv2/autoace/BLAT/autoace.first\nquit\n";
  }
  else {
    $command  = "query find Sequence \"SUPERLINK*\"\n";
    $command .= "show -a -f /wormsrv1/camace/BLAT/superlinks.ace\n";
    $command .= "dna -f /wormsrv1/camace/BLAT/autoace.first\nquit\n";
  }
  
  # tace dump chromosomal DNA and superlinks file
  &run_command("echo '$command' | $giface $dbpath");

  # Check that superlinks file created ok
  &usage(11) unless (-e "${blat_dir}/superlinks.ace");

  # Change '-'s in chromosome sequences into 'n's because blat excludes '-'
  # Not strictly needed anymore but left in for safety
  my $sequence;

  open (CHANGE, "<$blat_dir/autoace.first");
  open (NEW, ">$seq");
  while (<CHANGE>) {
    chomp;
    $sequence = $_;
    $sequence =~ tr/-/n/;
    print NEW "$sequence\n";
  }
  close(CHANGE);
  close(NEW);
  
  # remove intermediary sequence file
  unlink ("${blat_dir}/autoace.first") if (-e "${blat_dir}/autoace.first");



}

#########################################################################################################



sub confirm_introns {

  my ($db,$data) = @_;
  local (*GOOD,*BAD,*SEQ);
  
  # open the output files
  open (GOOD, ">$blat_dir/$db.good_introns.$data.ace") or die "$!";
  open (BAD,  ">$blat_dir/$db.bad_introns.$data.ace")  or die "$!";
  
  my ($link,@introns,$dna,$switch,$tag);
  
  ($tag = "cDNA") if ($mrna || $embl);
  ($tag = "EST")  if ($est || $ost); 
  
  
  $/ = "";
  open (CI, "<$blat_dir/${db}.ci.${data}.ace")      or die "Cannot open $blat_dir/$db.ci.$data.ace $!\n";
  while (<CI>) {
    next unless /^\S/;
    if (/Sequence : \"(\S+)\"/) {
      $link = $1;
      print "Sequence : $link\n";
      @introns = split /\n/, $_;
      
      # get the link sequence #
      print "Extracting DNA sequence for $link\n";
      undef ($dna);
      
      open(SEQ, "<$blat_dir/autoace.fa") || &usage(5);
      $switch = 0;
      $/ = "\n";
      
      # added shortcuts next & last to speed this section
      
      while (<SEQ>) {
	if (/^\>$link$/) {
	  $switch = 1;
	  next;
	}
	elsif (/^(\w+)$/) {
	  if ($switch == 1) {
	    chomp;
	    $dna .= $1;
	  }			
	}
	elsif ($switch == 1) {
	  $switch = 0;
	  last;
	}
	else { 
	  $switch = 0;
	}
      }
      close SEQ;
      
      print "DNA sequence is " . length($dna) . " bp long.\n";
      
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
	    
	    $start = substr($dna,$first,2);
	    $end   = substr($dna,$prelast,2);
	  
#	    print "Coords start $f[1] => $start, end $f[2] => $end\n";
	  
	  ##################
	  # map to S_child #
	  ##################
	  
	  my $lastvirt = int((length $dna) /100000) + 1;
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
	    $virtual = "Confirmed_intron_EST:" .$link."_".$startvirt     if ($est);
	    $virtual = "Confirmed_intron_OST:" .$link."_".$startvirt     if ($ost);
	    $virtual = "Confirmed_intron_mRNA:".$link."_".$startvirt     if ($mrna);
	    $virtual = "Confirmed_intron_EMBL:".$link."_".$startvirt     if ($embl);
	  }
	  elsif (($startvirt == ($endvirt - 1)) && (($last%100000) <= 50000)) {
	    $virtual = "Confirmed_intron_EST:" .$link."_".$startvirt     if ($est);
	    $virtual = "Confirmed_intron_OST:" .$link."_".$startvirt     if ($ost);
	    $virtual = "Confirmed_intron_mRNA:".$link."_".$startvirt     if ($mrna);
	    $virtual = "Confirmed_intron_EMBL:".$link."_".$startvirt     if ($embl);
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
		    print GOOD "Confirmed_intron $one $two $tag $f[4]\n\n";
		}
	    }
	    else {
		print BAD "Feature_data : \"$virtual\"\n";
		print BAD "Confirmed_intron $one $two $tag $f[4]\n\n";		
	    }
	}
    }
  }
}
  close CI;
  
  close GOOD;
  close BAD;
  
}


#############################
# virtual object generation #
#############################

sub virtual_objects_blat {
    
  my ($data) = shift;
  local (*OUT_autoace_homol,*OUT_camace_homol,*OUT_stlace_homol);
  local (*OUT_autoace_feat,*OUT_camace_feat,*OUT_stlace_feat);
  my ($name,$length,$total,$first,$second,$m,$n);
  
  # autoace
  open (OUT_autoace_homol, ">$blat_dir/virtual_objects.autoace.blat.$data.ace") or die "$!";
  open (OUT_autoace_feat,  ">$blat_dir/virtual_objects.autoace.ci.$data.ace")     or die "$!";
  # camace
  open (OUT_camace_homol,  ">$blat_dir/virtual_objects.camace.blat.$data.ace")  or die "$!";
  open (OUT_camace_feat,   ">$blat_dir/virtual_objects.camace.ci.$data.ace")      or die "$!";
  # stlace
  open (OUT_stlace_homol,  ">$blat_dir/virtual_objects.stlace.blat.$data.ace")  or die "$!";
  open (OUT_stlace_feat,   ">$blat_dir/virtual_objects.stlace.ci.$data.ace")      or die "$!";
  
  open (ACE, "<$blat_dir/chromosome.ace") || die &usage(11);
  while (<ACE>) {
    if (/Subsequence\s+\"(\S+)\" (\d+) (\d+)/) {
      $name   = $1;
      $length = $3 - $2 + 1;
      $total = int($length/100000) +1;
      # autoace
      print OUT_autoace_homol "Sequence : \"$name\"\n";
      print OUT_autoace_feat  "Sequence : \"$name\"\n";

      # Have to ignore MTCE sequence as there is no lab (RW or HX) associated with it
      unless($name eq "MTCE"){
	# camace
	print OUT_camace_homol  "Sequence : \"$name\"\n" if ($homedb{$name} eq "HX");
	print OUT_camace_feat   "Sequence : \"$name\"\n" if ($homedb{$name} eq "HX");
	# stlace
	print OUT_stlace_homol  "Sequence : \"$name\"\n" if ($homedb{$name} eq "RW");
	print OUT_stlace_feat   "Sequence : \"$name\"\n" if ($homedb{$name} eq "RW");
      }
      for ($n = 0; $n <= $total; $n++) {
	$m      = $n + 1;
	$first  = ($n*100000) + 1;
	$second = $first + 149999;
	if (($length - $first) < 100000) {
	  $second = $length;
	  # autoace
	  print OUT_autoace_homol "S_child Homol_data $word{$data}:$name"."_$m $first $second\n";
	  print OUT_autoace_feat  "S_child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n";

	  unless($name eq "MTCE"){
	    # camace
	    print OUT_camace_homol  "S_child Homol_data $word{$data}:$name"."_$m $first $second\n"             if ($homedb{$name} eq "HX");
	    print OUT_camace_feat   "S_child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n" if ($homedb{$name} eq "HX");
	    # stlace
	    print OUT_stlace_homol  "S_child Homol_data $word{$data}:$name"."_$m $first $second\n"             if ($homedb{$name} eq "RW");
	    print OUT_stlace_feat   "S_child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n" if ($homedb{$name} eq "RW");
	  }
	  last;
	}					
	else {
	  ($second = $length) if ($second >  $length);
	  # autoace
	  print OUT_autoace_homol "S_child Homol_data $word{$data}:$name"."_$m $first $second\n";
	  print OUT_autoace_feat  "S_child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n";
	  # camace
	  print OUT_camace_homol  "S_child Homol_data $word{$data}:$name"."_$m $first $second\n"             if ($homedb{$name} eq "HX");
	  print OUT_camace_feat   "S_child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n" if ($homedb{$name} eq "HX");
	  # stlace
	  print OUT_stlace_homol  "S_child Homol_data $word{$data}:$name"."_$m $first $second\n"             if ($homedb{$name} eq "RW");
	  print OUT_stlace_feat   "S_child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n" if ($homedb{$name} eq "RW");
	}
      }
      print OUT_autoace_homol "\n";
      print OUT_autoace_feat  "\n";

      unless($name eq "MTCE"){
	print OUT_camace_homol  "\n" if ($homedb{$name} eq "HX");
	print OUT_camace_feat   "\n" if ($homedb{$name} eq "HX");
	print OUT_stlace_homol  "\n" if ($homedb{$name} eq "RW");
	print OUT_stlace_feat   "\n" if ($homedb{$name} eq "RW");
      }
    }
  }
  close ACE;
  close OUT_autoace_homol;
  close OUT_autoace_feat;
  close OUT_camace_homol;
  close OUT_camace_feat;
  close OUT_stlace_homol;
  close OUT_stlace_feat;
  
  # clean up if you are dealing with parasitic nematode conensus or TC1 insertion data
  # dl 040315 - this is crazy. we make all of the files and then delete the ones we don't want.
  #             don't rock the boat...

  if ( ($data eq "nematode") || ($data eq "tc1") ) {
    unlink ("$blat_dir/virtual_objects.autoace.ci.$data.ace");
    unlink ("$blat_dir/virtual_objects.camace.ci.$data.ace");
    unlink ("$blat_dir/virtual_objects.stlace.ci.$data.ace");
  }

}

##################################################################################################


#################################################################################################
#
# Usage / Help subroutine
#
##################################################################################################


sub usage {
    my $error = shift;

    if ($error eq "Help") {
      # Normal help menu
      system ('perldoc',$0);
      exit (0);
    }

    if ($error == 1) {
      # no option supplied
      print "\nNo process option choosen [-dump|-blat|-process|virtual]\n";
      print "Run with one of the above options\n\n";
      exit(0);
    }
    elsif ($error == 2) {
      # No data-type choosen
      print "\nNo data option choosen [-est|-mrna|-ost|-embl|-nematode]\n";
      print "Run with one of the above options\n\n";
      exit(0);
    }
    elsif ($error == 3) {
      # 'Multiple data-types choosen
      print "\nMultiple data option choosen [-est|-mrna|-ost|-nematode|-embl]\n";
      print "Run with one of the above options\n\n";
      exit(0);
    }
    elsif ($error == 4) {
      # -fine used without -mrna
      print "-fine can only be specified if you are using -mrna\n\n";
      exit(0);
    }
    elsif ($error == 5) {
      # 'autoace.fa' file is not there or unreadable
      print "\nThe WormBase 'autoace.fa' file you does not exist or is non-readable.\n";
      print "Check File: '${blat_dir}/autoace.fa'\n\n";
      exit(0);
    }
    elsif ($error == 6) {
      # BLAT failure
      print "BLAT failure.\n";
      print "Whoops! you're going to have to start again.\n\n";
      exit(0);
    }
    elsif ($error == 11) {
      # 'superlinks.ace' file is not there or unreadable
      print "\nThe WormBase 'superlinks.ace' file you does not exist or is non-readable.\n";
      print "Check File: '${blat_dir}/superlinks.ace'\n\n";
      exit(0);
    }
    elsif ($error == 0) {
      # Normal help menu
      exec ('perldoc',$0);
    }
}

######################################################################################################


sub run_command{
  my $command = shift;
  print LOG &runtime, ": started running $command\n";
  my $status = system($command);
  if($status != 0){
    $errors++;
    print LOG "ERROR: $command failed\n";
  }
  print LOG &runtime, ": finished running $command\n";

  # for optional further testing by calling subroutine
  return($status);
}

################################################################
sub create_log_files{

  my $WS_version = &get_wormbase_version_name;


  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate     = `date +%y%m%d`; chomp $rundate;
  $log        = "/wormsrv2/logs/${script_name}.${WS_version}.${rundate}.$$";

  open (LOG, ">$log") or die "cant open $log";
  print LOG "$script_name\n";
  print LOG "started at ",`date`,"\n";
  print LOG "WormBase version : ${WS_version}\n";
  print LOG "======================================================================\n";
  print LOG "blat_them_all.pl query sequence options:\n";
  print LOG " -est      : perform blat for ESTs\n"                              if ($est);
  print LOG " -mrna     : perform blat for mRNAs\n"                             if ($mrna);
  print LOG " -fine     : used with -mrna option\n"                             if ($fine);
  print LOG " -ost      : perform blat for OSTs\n"                              if ($ost);
  print LOG " -embl     : perform blat for misc. non-WormBase CDS from EMBL\n"  if ($embl);
  print LOG " -nematode : perform blatx for parasitic nematode ESTs\n"          if ($nematode);
  print LOG "\n";
  print LOG "blat_them_all.pl process options:\n";
  print LOG " -dump      : dump chromosome sequences from autoace\n"            if ($dump);
  print LOG " -blat      : run blat against dumped chromosome sequences\n"      if ($blat);
  print LOG " -process   : sort and process raw blat output\n"                  if ($process);
  print LOG " -virtual   : create virtual sequence objects\n"                   if ($virtual);
  print LOG "\n";
}

##########################################################################################################

# Old comments from header part of script
#
# 16.10.01 Kerstin Jekosch
# 17.10.01 kj: modified to get everything onto wormsrv2 and to include an mRNA and parasitic nematode blatx option
# 26.11.01 kj: runs everything for autoace AND camace now
# 13.11.01 kj: added some file removers to tidy up at the end
# 14.11.01 kj: added option to just produce virtual objects
# 01.02.02 dl: added option to search miscPep file
# 01.02.02 dl: uncommented report logging & mail routines
# 01.02.02 dl: routine to convert '-' -> 'N' needs to be within the same BLOCK as the tace command
#            : else you get a zero length fasta file each time and the confirm intron routines fail
# 02.02.21 dl: typos in the naming of the confirmed_intron virtual objects
# 02.04.08 dl: old style logging for autoace.fa check, prevented complete run of subs
#



__END__

=pod

=head2   NAME - blat_them_all.pl

=head1 USAGE

=over 4

=item  blat_them_all.pl -options

=back

A wrapper script to generate blat data by: 

1) getting target sequence out of autoace

2) blatting it against all ESTs/OSTs, mRNAs, and CDSs from EMBL entries

3) Processing blat output files, and mapping hits back to autoace/camace 

4) Producing confirmed introns

5) Producing virtual objects to 'hang' the data onto

All output is stored in /wormsrv2/autoace/BLAT/ (or /wormsrv1/camace/BLAT/ if -camace
is specified)

blat_them_all mandatory arguments:

=over 4

=item -est

run everything for ESTs

=back

or

=item -mrna

run everything for mRNAs

=back

or

=item -embl

run everything for the CDSs of non-WormBase gene predictions in EMBL

=back

or

=item -nematode   

run everything for non-C.elegans nematode ESTs

=back

or 

=item -ost   

run everything for OST data

=back



blat_them_all optional arguments:

=item -dump      

start by first dumping out target chromosome sequences and acefiles from autoace/camace

=back

=item -blat      

start with blatting (i.e. autoace.fa & chromosome.ace already present)

=back

=item -process   

start later by processing (and sorting/mapping) existing *.psl file

=back

=item -camace    

Use /wormsrv1/camace rather than the default /wormsrv2/autoace

=back

=item -fine      

Forces use of new -fine option in blat, only in conjunction with -mrna

=back

=item -verbose   

Show more output to screen (useful when running on command line)

=back

=item -debug <user> 

Send output only to user and not to everyone in group

=back

=item -help      

This help

=back

Script written by Kerstin Jekosch with heavy rewrite by Keith Bradnam

=cut



