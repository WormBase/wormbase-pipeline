#!/usr/local/bin/perl5.8.0 -w
#
# blat_them_all.pl
# 
# by Kerstin Jekosch
#
# Gets sequences ready for blatting, blats sequences, processes blat output, makes confirmed introns
# and virtual objects to hang the data onto
#
# Last edited by: $Author: krb $
# Last edited on: $Date: 2003-09-04 09:30:31 $

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Carp;


##############################
# Misc variables and paths   #
##############################

my ($help, $debug, $verbose, $est, $mrna, $ost, $nematode, $embl, 
    $blat, $process, $virtual, $dump, $camace, $fine);
my $maintainers = "All";
our $log;
our $blat_dir = "/wormsrv2/autoace/BLAT"; # default BLAT directory, can get changed if -camace used
our $dbpath = "/wormsrv2/autoace"; # default database location
my $bin     = "/wormsrv2/scripts";
our %homedb; # for storing superlink->lab connections
our $blatex  = '/nfs/disk100/wormpub/bin.ALPHA/blat';
our $giface  = &giface;
our %word = (
	     EST      => 'BLAT_EST',
	     mRNA     => 'BLAT_mRNA',
	     EMBL     => 'BLAT_EMBL',
	     NEMATODE => 'BLATX_NEMATODE',
	     OST      => 'BLAT_OST',
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
$dbpath = "/wormsrv1/camace"     if $camace;
our $seq   = "$blat_dir/autoace.fa";               

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
&usage(2) unless ($est || $mrna || $embl || $nematode || $ost || $dump); 

# Exit if multiple data types choosen [EST|mRNA|EMBL|NEMATODE|OST]
&usage(3) if (($est + $mrna + $ost + $embl + $nematode) > 1);

# Exit if -fine option is being used without -mrna
&usage(4) if ($fine && !$mrna);


# assign type variable
my $data;
($data = 'EST')      if ($est);
($data = 'OST')      if ($ost);
($data = 'mRNA')     if ($mrna);
($data = 'EMBL')     if ($embl);
($data = 'NEMATODE') if ($nematode);


# Select the correct set of query sequences for blat
my $query = "/nfs/disk100/wormpub/analysis/ESTs/";
$query   .= 'C.elegans_nematode_ESTs'     if ($est); # EST data set
$query   .= 'C.elegans_nematode_OSTs'     if ($ost); # OST data set
$query   .= 'C.elegans_nematode_mRNAs'    if ($mrna); # mRNA data set
$query   .= 'non_C.elegans_nematode_ESTs' if ($nematode); # ParaNem EST data set
$query   .= 'C.elegans_nematode_miscPep'  if ($embl); # Other CDS data set, DNA not peptide!


&create_log_files;




###########################################################################
#                                                                         #
#                        MAIN PART OF PROGRAM                             #
#                                                                         #
###########################################################################



# Write sequence data (chromsomes) from autoace/camace
# Also assigns laboratories to each superlink object (stores in %homedb)
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
  

  if ($est) {
    unless ($camace) {
      system("$bin/blat2ace.pl -est -intron") && croak "Mapping failed\n"; 
    }
    else {
      system("$bin/blat2ace.pl -est -intron -camace") && croak "Mapping failed\n"; 	
    }
  }

  if ($ost) {
    unless ($camace) {
      system("$bin/blat2ace.pl -ost -intron") && croak "Mapping failed\n"; 
    }
    else {
      system("$bin/blat2ace.pl -ost -intron -camace") && croak "Mapping failed\n"; 	
    }
  }

  if ($mrna) {
    unless ($camace) {
      system("$bin/blat2ace.pl -mrna -intron") && croak "Mapping failed\n"; 
    }
    else {
      system("$bin/blat2ace.pl -mrna -intron -camace") && croak "Mapping failed\n"; 
    }
  }

  if ($embl) {  
    unless ($camace) {
      system("$bin/blat2ace.pl -embl -intron") && croak "Mapping failed\n"; 
    }
    else {
      system("$bin/blat2ace.pl -embl -intron -camace") && croak "Mapping failed\n"; 
    }
  }
  
  if ($nematode) {
    system("$bin/blat2ace.pl -nematode") && croak "Mapping failed\n"; 
  }


  $runtime = &runtime;
  print "Producing confirmed introns in databases\n\n" if $verbose;
  print LOG "$runtime: Producing confirmed introns in databases\n";

  # produce confirmed introns #
  if ($est) {    
    print "Producing confirmed introns using EST data\n" if $verbose;
    &confirm_introns('autoace','EST');
    &confirm_introns('camace', 'EST');
    &confirm_introns('stlace', 'EST');
  }
  if ($ost) {
    print "Producing confirmed introns using OST data\n" if $verbose;
    &confirm_introns('autoace','OST');
    &confirm_introns('camace', 'OST');
    &confirm_introns('stlace', 'OST');
  }
  if ($mrna) {
    print "Producing confirmed introns using mRNA data\n" if $verbose;
    &confirm_introns('autoace','mRNA');
    &confirm_introns('camace', 'mRNA');
    &confirm_introns('stlace', 'mRNA');
  }
  if ($embl) {
    print "Producing confirmed introns using EMBL CDS data\n" if $verbose;
    &confirm_introns('autoace','EMBL');
    &confirm_introns('camace', 'EMBL');
    &confirm_introns('stlace', 'EMBL');
  }
}


#########################################
# produce files for the virtual objects #
#########################################
if ($virtual) {
  my $runtime = &runtime;
  print LOG "$runtime: Producing $data files for the virtual objects\n";
  &virtual_objects_blat($data);
}


##############################
# Clean up and say goodbye   #
##############################

close(LOG);

# mail log to maintainer
&mail_maintainer("WormBase Report: blat_them_all ",$maintainers,$log);

exit(0);





#################################################################################
#                                                                               #
#                     T H E    S U B R O U T I N E S                            #
#                                                                               #
#################################################################################





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
    $command .= "follow subsequence\n";
    $command .= "show -a -f /wormsrv2/autoace/BLAT/superlinks.ace\n";
    $command .= "dna -f /wormsrv2/autoace/BLAT/autoace.first\nquit\n";
  }
  else {
    $command  = "query find Sequence \"SUPERLINK*\"\n";
    $command .= "show -a -f /wormsrv1/camace/BLAT/superlinks.ace\n";
    $command .= "dna -f /wormsrv1/camace/BLAT/autoace.first\nquit\n";
  }
  
  # tace dump chromosomal DNA and superlinks file
  system("echo '$command' | $giface $dbpath") && &usage(5);

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



  # Now connect superlink objects to their corresponding laboratory object
  # store in global %homedb
  local (*LINK);
  my $name;
  
  print LOG "Assign LINK* objects to laboratory\n\n";
  # deal with the superlink objects
  open (LINK, "<$blat_dir/superlinks.ace") || croak "Couldn't open superlinks.ace $!";
  while (<LINK>) {
    if (/^Sequence\s+\:\s+\"(\S+)\"/) {
      $name = $1;
      next;
    }
    if (/^From_Laboratory\s+\"(\S+)\"/) {
      $homedb{$name} = $1;
      print LOG "assigning $1 to $name\n";
      undef ($name);
      next;
    }
  }
  close(LINK);
  
  print LOG "\n";
  
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
	  
#		    print "Coords start $f[1] => $start, end $f[2] => $end\n";
	  
	  ##################
	  # map to S_Child #
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
	    print GOOD "Confirmed_intron $one $two $tag\n\n";
	  }  	
	  else {
	    print BAD "Feature_data : \"$virtual\"\n";
	    print BAD "Confirmed_intron $one $two $tag\n\n";		
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
  open (OUT_autoace_homol, ">$blat_dir/virtual_objects.autoace.$word{$data}.ace") or die "$!";
  open (OUT_autoace_feat,  ">$blat_dir/virtual_objects.autoace.ci.$data.ace")     or die "$!";
  # camace
  open (OUT_camace_homol,  ">$blat_dir/virtual_objects.camace.$word{$data}.ace")  or die "$!";
  open (OUT_camace_feat,   ">$blat_dir/virtual_objects.camace.ci.$data.ace")      or die "$!";
  # stlace
  open (OUT_stlace_homol,  ">$blat_dir/virtual_objects.stlace.$word{$data}.ace")  or die "$!";
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
      # camace
      print OUT_camace_homol  "Sequence : \"$name\"\n" if ($homedb{$name} eq "HX");
      print OUT_camace_feat   "Sequence : \"$name\"\n" if ($homedb{$name} eq "HX");
      # stlace
      print OUT_stlace_homol  "Sequence : \"$name\"\n" if ($homedb{$name} eq "RW");
      print OUT_stlace_feat   "Sequence : \"$name\"\n" if ($homedb{$name} eq "RW");
      
      for ($n = 0; $n <= $total; $n++) {
	$m      = $n + 1;
	$first  = ($n*100000) + 1;
	$second = $first + 149999;
	if (($length - $first) < 100000) {
	  $second = $length;
	  # autoace
	  print OUT_autoace_homol "S_Child Homol_data $word{$data}:$name"."_$m $first $second\n";
	  print OUT_autoace_feat  "S_Child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n";
	  # camace
	  print OUT_camace_homol  "S_Child Homol_data $word{$data}:$name"."_$m $first $second\n"             if ($homedb{$name} eq "HX");
	  print OUT_camace_feat   "S_Child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n" if ($homedb{$name} eq "HX");
	  # stlace
	  print OUT_stlace_homol  "S_Child Homol_data $word{$data}:$name"."_$m $first $second\n"             if ($homedb{$name} eq "RW");
	  print OUT_stlace_feat   "S_Child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n" if ($homedb{$name} eq "RW");
	  last;
	}					
	else {
	  ($second = $length) if ($second >  $length);
	  # autoace
	  print OUT_autoace_homol "S_Child Homol_data $word{$data}:$name"."_$m $first $second\n";
	  print OUT_autoace_feat  "S_Child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n";
	  # camace
	  print OUT_camace_homol  "S_Child Homol_data $word{$data}:$name"."_$m $first $second\n"             if ($homedb{$name} eq "HX");
	  print OUT_camace_feat   "S_Child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n" if ($homedb{$name} eq "HX");
	  # stlace
	  print OUT_stlace_homol  "S_Child Homol_data $word{$data}:$name"."_$m $first $second\n"             if ($homedb{$name} eq "RW");
	  print OUT_stlace_feat   "S_Child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n" if ($homedb{$name} eq "RW");
	}
      }
      print OUT_autoace_homol "\n";
      print OUT_autoace_feat  "\n";
      print OUT_camace_homol  "\n" if ($homedb{$name} eq "HX");
      print OUT_camace_feat   "\n" if ($homedb{$name} eq "HX");
      print OUT_stlace_homol  "\n" if ($homedb{$name} eq "RW");
      print OUT_stlace_feat   "\n" if ($homedb{$name} eq "RW");
    }
  }
  close ACE;
  close OUT_autoace_homol;
  close OUT_autoace_feat;
  close OUT_camace_homol;
  close OUT_camace_feat;
  close OUT_stlace_homol;
  close OUT_stlace_feat;
  
  # clean up if you are dealing with parasitic nematode conensus data
  if ($data eq "NEMATODE") {
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



