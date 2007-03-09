#!/usr/local/bin/perl5.8.0 -w
#
# blat2ace.pl
# 
# by Kerstin Jekosch
#
# Exporter to map blat data to genome and to find the best match for each EST, mRNA, OST, etc.
#
# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2007-03-09 15:17:10 $

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

#########################
# Command line options  #
#########################

my ($est, $mrna, $ncrna, $ost, $tc1, $nematode, $embl, $camace, $intron, $washu, $nembase);
my ($help, $debug, $test, $verbose, $store, $wormbase);
my $database;

GetOptions (
	    "help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "database:s" => \$database,
            "est"        => \$est,
            "mrna"       => \$mrna,
            "ncrna"      => \$ncrna,
            "ost"        => \$ost,
            "tc1"        => \$tc1,
            "nematode"   => \$nematode,
	    "nembase"    => \$nembase,
	    "washu"      => \$washu,
            "embl"       => \$embl,
	    "intron"     => \$intron,
	   );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage(0) if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


#############################
# variables and directories #
#############################

# set database paths, default to autoace unless -camace
my $ace_dir         = $wormbase->autoace;     # AUTOACE DATABASE DIR
my $blat_dir  = $wormbase->blat;
my $canonical = $wormbase->database('camace');

if ($camace) {
    $blat_dir  = "$canonical/BLAT";
}

#############################
# CommonData hash retrieval #
#############################

my %NDBaccession2est = $wormbase->FetchData('NDBaccession2est');     # EST accession => name
my %estorientation   = $wormbase->FetchData('estorientation');       # EST accession => orientation [5|3]

my %hash;
my (%best,%other,%bestclone,%match,%ci);

our %camace;
our %stlace;

our $type = "";
our %word = (
	     est      => 'BLAT_EST',
	     ost      => 'BLAT_OST',
	     mrna     => 'BLAT_mRNA',
	     ncrna    => 'BLAT_ncRNA',
	     embl     => 'BLAT_EMBL',
	     tc1      => 'BLAT_TC1',
	     nematode => 'BLAT_NEMATODE',
	     nembase    => 'BLAT_NEMBASE',
	     washu    => 'BLAT_WASHU'
	     );


########################################
# command-line options & ramifications #
########################################

# Exit if no data type choosen [EST|mRNA|EMBL|NEMATODE|OST|WASHU|NEMBASE]
# or if multiple data types are chosen


$log->log_and_die("no type specified\n") unless 
    ($est || $mrna || $ost || $ncrna || $tc1 || $nematode || $embl || $washu || $nembase);


my $flags = 0;
$flags++ if $est;
$flags++ if $ost;
$flags++ if $mrna;
$flags++ if $ncrna;
$flags++ if $embl;
$flags++ if $tc1;
$flags++ if $nematode;
$flags++ if $nembase;
$flags++ if $washu;
&usage(2) if ($flags > 1);

# assign type variable
($type = 'est')      if ($est);
($type = 'ost')      if ($ost);
($type = 'mrna')     if ($mrna);
($type = 'ncrna')    if ($ncrna);
($type = 'embl')     if ($embl);
($type = 'tc1')      if ($tc1);
($type = 'nematode') if ($nematode);
($type = 'nembase')  if ($nembase);
($type = 'washu')    if ($washu);

#########################################
# get links for database                #
#########################################

# parse links for camace
my @camclones = qw(SUPERLINK_CB_I SUPERLINK_CB_II SUPERLINK_CB_IIIL SUPERLINK_CB_IIIR SUPERLINK_CB_IR SUPERLINK_CB_IV SUPERLINK_CB_V SUPERLINK_CB_X); 
foreach my $camclone (@camclones) {
  $camace{$camclone} = 1;
}

# parse links for stlace
my @stlclones = qw(SUPERLINK_RW1 SUPERLINK_RW1R SUPERLINK_RW2 SUPERLINK_RW3A SUPERLINK_RW3B SUPERLINK_RW4 SUPERLINK_RW5 SUPERLINK_RWXL SUPERLINK_RWXR);
foreach my $stlclone (@stlclones) {
  $stlace{$stlclone} = 1;
}


##########################################################################################
# map the blat hits to ace - i.e. process blat output (*.psl) file into set of ace files #
##########################################################################################
$log->write_to($wormbase->runtime.": Start mapping\n\n");

# open input and output filehandles
open(ACE,  ">$blat_dir/autoace.$type.ace_uncompressed")  or die "Cannot open $blat_dir/autoace.${type}.ace_uncompressed $!\n";
open(BLAT, "<$blat_dir/${type}_out.psl")    or die "Cannot open $blat_dir/${type}_out.psl $!\n";

my $number_of_replacements = 0;
my %reported_this_query_before;

# loop through each blat hit
while (<BLAT>) {
  next unless (/^\d/);
  my @f            = split "\t";

  my $match        = $f[0];                    # number of bases matched by blat
  my $strand       = $f[8];                    # strand that match is on
  my $query        = $f[9];                    # query sequence name
  my $query_size   = $f[10];                   # query sequence length
  my $superlink    = $f[13];                   # name of superlink that was used as blat target sequence
  my $slsize       = $f[14];                   # superlink size
  my $lastvirt     = int($slsize/100000) + 1;  # for tracking how many virtual sequences have been created???
  my $matchstart   = $f[15];                   # target (superlink) start coordinate...
  my $matchend     = $f[16];                   # ...and end coordinate
  my $block_count  = $f[17];                   # block count
  my @lengths      = split (/,/, $f[18]);      # sizes of each blat 'block' in any individual blat match
  my @query_starts = split (/,/, $f[19]);      # start coordinates of each query block
  my @slink_starts = split (/,/, $f[20]);      # start coordinates of each target (superlink) block

  # replace EST name (usually accession number) by yk... name 
  if ( ($est || $ost) && (exists $NDBaccession2est{$query}) ) {
    my $estname  = $NDBaccession2est{$query};
    if ($query ne $estname) {
#	  $log->write_to("EST name $query was replaced by $estname\n\n");
      $number_of_replacements++;
      $query = $estname;
    }
  }

  ###############################
  # find virtual superlink part #
  ###############################
	
  my ($virtual,$startvirtual,$endvirtual);

  if ((int($matchstart/100000) + 1) > $lastvirt) { 
      $startvirtual = $lastvirt;
  }
  else {
      $startvirtual = int($matchstart/100000) +1;
  }  
    
  if ( (int($matchend/100000) +1) > $lastvirt) { 
      $endvirtual = $lastvirt;
  }
  else {
      $endvirtual = int($matchend/100000) +1;
  }  
  
  if ($startvirtual == $endvirtual) {
      $virtual = "$word{$type}:${superlink}_${startvirtual}";
  }	
  elsif (($startvirtual == ($endvirtual - 1)) && (($matchend%100000) <= 50000)) {
      $virtual = "$word{$type}:${superlink}_${startvirtual}";
  }
  else {
    if (! exists $reported_this_query_before{$query}) {
      $reported_this_query_before{$query} = 1; # don't want to report this one again
      $log->write_to("$query wasn't assigned to a virtual object as match size was too big\n");
      $log->write_to("Start is $matchstart, end is $matchend on $superlink\n\n");
    }
    next;
  }

  # calculate (acedb) score for each blat match
  # new way of calculating score, divide by query size rather than sum of matching blocks, 
  my $score = ($match/$query_size)*100;
  
  #########################
  # calculate coordinates #
  #########################
    
  # need to allow for est exons in the next virtual object, otherwise they get remapped to the start 
  # of the virtual by performing %100000
  my @exons = ();  
  my $calc = int(($slink_starts[0]+1)/100000);
  
  for (my $x = 0;$x < $block_count; $x++) {
      my $newcalc = int(($slink_starts[$x]+1)/100000);
      my $virtualstart;

      if ($calc == $newcalc) {	
	  $virtualstart =  ($slink_starts[$x] +1)%100000;
      }
      elsif ($calc == ($newcalc - 1)) {
	  $virtualstart = (($slink_starts[$x] +1)%100000) + 100000;
      }
      

      my $virtualend = $virtualstart + $lengths[$x] -1;
      
      if ($calc != $newcalc) {
	  $log->write_to("// MISMATCH: $query [$strand] $virtualstart $virtualend :: [virtual slice $calc -> $newcalc, offset ".($matchend%100000)."]\n\n");
      }

      if (!defined $virtualstart) {
	  $log->write_to("$query will be discarded as the match is too long\n");
	  $log->write_to("$query [$strand] $virtualstart $virtualend  [virtual slice $calc -> $newcalc, offset ".($matchend%100000)."]\n\n");
	  next;
      }

      my ($query_start,$query_end);
      
        # blatx 6-frame translation v 6-frame translation
      if ($nematode || $washu || $nembase) {
	  my $temp;
	  if (($strand eq '++') || ($strand eq '-+')) {
	      $query_start = $query_starts[$x] +1;
	      $query_end   = $query_start + $lengths[$x] -1;
	      if ($strand eq '-+') {
		  $temp        = $query_end;
		  $query_end   = $query_start;
		  $query_start = $temp; 
	      }
	  }
	  elsif (($strand eq '--') || ($strand eq '+-')) {
	      $temp         = $virtualstart;
	      $virtualstart = $virtualend;
	      $virtualend   = $temp;

	      $query_start  = $query_size  - $query_starts[$x];
	      $query_end    = $query_start - $lengths[$x] +1;

	      if ($strand eq '--') {
		  $temp        = $query_end;
		  $query_end   = $query_start;
		  $query_start = $temp; 
	      }
	  }			
      }
      else {
	  if ($strand eq '+'){
	      $query_start   = $query_starts[$x] +1;
	      $query_end     = $query_start + $lengths[$x] -1;
	  }
	  elsif ($strand eq '-') {
	      $query_start   = $query_size - $query_starts[$x];
	      $query_end     = $query_start - $lengths[$x] +1;
	  }		
      }		
      
      # 1 out of 8000,000 cDNAs ends up with 0 as a start and breaks the mapping, this hack is to stop it
      $virtualstart++ if($virtualstart == 0);
      $query_start++ if ($query_start == 0);

      # write to output file
      print ACE "Homol_data : \"$virtual\"\n";
      if ($type eq "nematode" || $type eq "washu" || $type eq "nembase") {
	  printf ACE "DNA_homol\t\"%s\"\t\"$word{$type}\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$query,$score,$virtualstart,$virtualend,$query_start,$query_end;
	  
#      print "// ERROR: $query [$strand] $virtualstart $virtualend $query_start $query_end ::: [$debug_start,$debug_end]  $newcalc - $calc {$slink_starts[$x]}\n" unless ((defined $virtualstart) && (defined $virtualend));
	  
      }
      else {
	  printf ACE "DNA_homol\t\"%s\"\t\"$word{$type}_OTHER\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$query,$score,$virtualstart,$virtualend,$query_start,$query_end;
      }
      push @exons, [$virtualstart,$virtualend,$query_start,$query_end];				
  }
  
    
  # collect best hits for each query sequence 
  # Choose hit with highest score (% of query length which are matching bases) 
  # If multiple hits have same scores (also meaning that $match must be same) store 
  # details of extra hits against same primary key in %best
  if (exists $best{$query}) {
      if (($score > $best{$query}->{'score'})) { 
	  # Add all new details if score is better...
	  $best{$query}->{'score'} = $score;
	  $best{$query}->{'match'} = $match;
	  @{$best{$query}->{'entry'}} = ({'clone' => $virtual,'link' => $superlink,'exons' => \@exons});
      }
      elsif($score == $best{$query}->{'score'}){
	  #...only add details (name and coordinates) of extra hits if scores are same
	  push @{$best{$query}->{'entry'}}, {'clone' => $virtual,'link' => $superlink,'exons' => \@exons};
      }
  }
  else {
      $best{$query}->{'match'} = $match;
      $best{$query}->{'score'} = $score;
      @{$best{$query}->{'entry'}} = ({'clone' => $virtual,'link' => $superlink,'exons' => \@exons});
  }
}
close(BLAT);
close(ACE);
#close (OUTBLAT);

# concise report
if ($est || $ost) {
  $log->write_to("\nThere were $number_of_replacements replacements of EST names\n\n");
}

####################################
# produce outfile for best matches #
####################################
if ($nematode || $washu || $nembase) {
  $wormbase->run_command("mv $blat_dir/autoace.$type.ace_uncompressed $blat_dir/autoace.blat.$type.ace_uncompressed", $log);
} else {
  open (AUTBEST, ">$blat_dir/autoace.best.$type.ace_uncompressed");
  open (STLBEST, ">$blat_dir/stlace.best.$type.ace_uncompressed");
  open (CAMBEST, ">$blat_dir/camace.best.$type.ace_uncompressed");

  foreach my $found (sort keys %best) {
    if (exists $best{$found}) {
      foreach my $entry (@{$best{$found}->{'entry'}}) {
	if (@{$best{$found}->{'entry'}} < 2) {
	  my $virtual   = $entry->{'clone'};
	  my $superlink = $entry->{'link'};
	  foreach my $ex (@{$entry->{'exons'}}) {
	    my $score        = $best{$found}->{'score'};
	    my $virtualstart = $ex->[0];
	    my $virtualend   = $ex->[1];
	    my $query_start  = $ex->[2];
	    my $query_end    = $ex->[3];
		    
	    # print output for autoace, camace, and stlace
	    print  AUTBEST "Homol_data : \"$virtual\"\n";
	    printf AUTBEST "DNA_homol\t\"%s\"\t\"$word{$type}_BEST\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$found,$score,$virtualstart,$virtualend,$query_start,$query_end;
	    if ($camace{$superlink}) {
	      print  CAMBEST "Homol_data : \"$virtual\"\n";
	      printf CAMBEST "DNA_homol\t\"%s\"\t\"$word{$type}_BEST\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$found,$score,$virtualstart,$virtualend,$query_start,$query_end;
	    } elsif ($stlace{$superlink}) {
	      print  STLBEST "Homol_data : \"$virtual\"\n";
	      printf STLBEST "DNA_homol\t\"%s\"\t\"$word{$type}_BEST\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$found,$score,$virtualstart,$virtualend,$query_start,$query_end;
	    }	  
	  }

		
	  #############################
	  # produce confirmed introns #
	  #############################
	  if ($intron) {
	    #$log->write_to("Producing confirmed introns\n");
	    my ($n) = ($virtual =~ /\S+_(\d+)$/);
	    for (my $y = 1; $y < @{$entry->{'exons'}}; $y++) {
	      my $last   = $y - 1;
	      my $first  =  (${$entry->{"exons"}}[$last][1] + 1) + (($n-1)*100000);
	      my $second =  (${$entry->{'exons'}}[$y][0]    - 1) + (($n-1)*100000);
	      $estorientation{$found} = 5 if ($mrna || $embl);
	      if (${$entry->{'exons'}}[0][2] < ${$entry->{'exons'}}[0][3]) {
		if ((${$entry->{'exons'}}[$y][2] == ${$entry->{'exons'}}[$last][3] + 1) && (($second - $first) > 2)) {
		  if (exists $estorientation{$found} && $estorientation{$found} eq '3') {
		    push @{$ci{$superlink}}, [$second,$first,$found];
		  } elsif (exists $estorientation{$found} && $estorientation{$found} eq '5') {
		    push @{$ci{$superlink}}, [$first,$second,$found];
		  } else {
		    $log->write_to("WARNING: Direction not found for $found\n\n");
		  }
		}
	      } elsif (${$entry->{'exons'}}[0][2] > ${$entry->{'exons'}}[0][3]) {
		if ((${$entry->{'exons'}}[$last][3] == ${$entry->{'exons'}}[$y][2] + 1) && (($second - $first) > 2)) {
		  if (exists $estorientation{$found} && $estorientation{$found} eq '3') {
		    push @{$ci{$superlink}}, [$first,$second,$found];
		  } elsif (exists $estorientation{$found} && $estorientation{$found} eq '5') {
		    push @{$ci{$superlink}}, [$second,$first,$found]; 
		  } else {
		    $log->write_to("WARNING: Direction not found for $found\n\n");
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  close(AUTBEST);
  close(CAMBEST);
  close(STLBEST);
}
########################################################
# produce final BLAT output (including BEST and OTHER) #
########################################################

unless ($nematode || $washu || $nembase) {
  # Open new (final) output files for autoace, camace, and stlace
  open (OUT_autoace, ">$blat_dir/autoace.blat.$type.ace_uncompressed") or die "$!";
  open (OUT_camace,  ">$blat_dir/camace.blat.$type.ace_uncompressed")  or die "$!";
  open (OUT_stlace,  ">$blat_dir/stlace.blat.$type.ace_uncompressed")  or die "$!";

  # Change input separator to paragraph mode, but store what it old mode in $oldlinesep
  my $oldlinesep = $/;
  $/ = "";

  my (%line);
  my $superlink = "";

  # assign 
  open(ABEST,  "<$blat_dir/autoace.best.$type.ace_uncompressed");
  while (<ABEST>) {
    if ($_ =~ /^Homol_data/) {
      # flag each blat hit which is best (all of them) - set $line{$_} to 1
      # %line thus stores keys which are combos of virtual object name + blat hit details
      $line{$_} = 1;
      ($superlink) = (/\"$word{$type}\:(\S+)\_\d+\"/);

      # Print blat best hits to final output file
      print OUT_autoace "// Source $superlink\n\n";
      print OUT_autoace $_;
    
      # camace
      if ($camace{$superlink}) {
	print OUT_camace $_;
      }
      # and stlace
      elsif ($stlace{$superlink}) {
	print OUT_stlace $_;
      }
    }
  }
  close ABEST;


  # Now look through original output file (where everything is set to BLAT_OTHER) to
  # output those blat OTHER hits which are not flagged as BLAT_BEST in the .best.ace file
  # Does this by comparing entries in %line hash

  open(AOTHER, "<$blat_dir/autoace.$type.ace_uncompressed");
  while (<AOTHER>) {
    if ($_ =~ /^Homol_data/) {
      my $line = $_;
      # for comparison to %line hash, need to change OTHER to BEST in $_
      s/BLAT_EST_OTHER/BLAT_EST_BEST/g     if ($est);
      s/BLAT_OST_OTHER/BLAT_OST_BEST/g     if ($ost); 
      s/BLAT_mRNA_OTHER/BLAT_mRNA_BEST/g   if ($mrna);
      s/BLAT_ncRNA_OTHER/BLAT_ncRNA_BEST/g if ($ncrna);
      s/BLAT_EMBL_OTHER/BLAT_EMBL_BEST/g   if ($embl);
      s/BLAT_TC1_OTHER/BLAT_TC1_BEST/g     if ($tc1);
    
      # Only output BLAT_OTHER hits in first output file which we now know NOT to
      # really be BEST hits
      unless (exists $line{$_}) {
	print OUT_autoace $line;
      
	# camace
	if ($camace{$superlink}) {
	  print OUT_camace $line;
	}
	# and stlace
	elsif ($stlace{$superlink}) {
	  print OUT_stlace $line;
	}
      
      }	
    }
  }
  close AOTHER;

  # reset input line separator
  $/= $oldlinesep;

  ###################################
  # produce confirmed intron output #
  ###################################

  if ($intron) {
    
    open(CI_auto, ">$blat_dir/autoace.ci.${type}.ace");
    open(CI_cam,  ">$blat_dir/camace.ci.${type}.ace");
    open(CI_stl,  ">$blat_dir/stlace.ci.${type}.ace");
  
    foreach my $link (sort keys %ci) {
      my %double;
	
      print CI_auto "\nSequence : \"$link\"\n";
      print CI_stl  "\nSequence : \"$link\"\n" if ($stlace{$link});
      print CI_cam  "\nSequence : \"$link\"\n" if ($camace{$link});
	
      for (my $i = 0; $i < @{$ci{$link}}; $i++) {
	my $merge = $ci{$link}->[$i][0].":".$ci{$link}->[$i][1];
	if (!exists $double{$merge}) {
	  if ($mrna) {
	    printf CI_auto "Confirmed_intron %d %d mRNA $ci{$link}->[$i][2]\n",  $ci{$link}->[$i][0], $ci{$link}->[$i][1];
	    (printf CI_cam "Confirmed_intron %d %d mRNA $ci{$link}->[$i][2]\n",  $ci{$link}->[$i][0], $ci{$link}->[$i][1]) if ($camace{$link});
	    (printf CI_stl "Confirmed_intron %d %d mRNA $ci{$link}->[$i][2]\n",  $ci{$link}->[$i][0], $ci{$link}->[$i][1]) if ($stlace{$link});
	  }
	  if ($embl) {
	    printf CI_auto "Confirmed_intron %d %d Homol $ci{$link}->[$i][2]\n",  $ci{$link}->[$i][0], $ci{$link}->[$i][1];
	    (printf CI_cam "Confirmed_intron %d %d Homol $ci{$link}->[$i][2]\n",  $ci{$link}->[$i][0], $ci{$link}->[$i][1]) if ($camace{$link});
	    (printf CI_stl "Confirmed_intron %d %d Homol $ci{$link}->[$i][2]\n",  $ci{$link}->[$i][0], $ci{$link}->[$i][1]) if ($stlace{$link});
	  }
	  if ($est) {
	    printf CI_auto "Confirmed_intron %d %d EST $ci{$link}->[$i][2]\n",  $ci{$link}->[$i][0], $ci{$link}->[$i][1];
	    (printf CI_cam "Confirmed_intron %d %d EST $ci{$link}->[$i][2]\n",  $ci{$link}->[$i][0], $ci{$link}->[$i][1]) if ($camace{$link});
	    (printf CI_stl "Confirmed_intron %d %d EST $ci{$link}->[$i][2]\n",  $ci{$link}->[$i][0], $ci{$link}->[$i][1]) if ($stlace{$link});
	  }
	  if ($ost) {
	    printf CI_auto "Confirmed_intron %d %d EST $ci{$link}->[$i][2]\n",  $ci{$link}->[$i][0], $ci{$link}->[$i][1];
	    (printf CI_cam "Confirmed_intron %d %d EST $ci{$link}->[$i][2]\n",  $ci{$link}->[$i][0], $ci{$link}->[$i][1]) if ($camace{$link});
	    (printf CI_stl "Confirmed_intron %d %d EST $ci{$link}->[$i][2]\n",  $ci{$link}->[$i][0], $ci{$link}->[$i][1]) if ($stlace{$link});
	  }
	  $double{$merge} = 1;
	}
      }
    }
    
    close CI_auto;
    close CI_cam;
    close CI_stl;
  }
}


#compress acefiles so that all object data is loaded together, expanded to run on all homol
my @filenames = ("autoace.$type.ace", "autoace.best.$type.ace", "autoace.blat.$type.ace", "camace.best.$type.ace", "camace.blat.$type.ace", "stlace.best.$type.ace", "stlace.blat.$type.ace");
my $filename;
$log->write_to("\n#########################################\nCompressing DNA_homolo acefiles\n#########################################\n");
foreach $filename (@filenames) {
#  if (-e ("$blat_dir/${filename}"."_uncompressed") ) {
    $log->write_to("Compressing ${filename}"."_uncompressed\n");
    $wormbase->run_script("acecompress.pl -file $blat_dir/${filename}_uncompressed -homol -build", $log);
    $log->write_to("Compressed........\n");
#  }
}

#$wormbase->run_script("acecompress.pl -file $blat_dir/autoace.best.$type.ace -homol", $log);
#$wormbase->run_script("acecompress.pl -file $blat_dir/autoace.blat.$type.ace -homol", $log);
#$wormbase->run_script("acecompress.pl -file $blat_dir/autoace.$type.ace -homol", $log);

$log->mail;
print "Finished.\n" if ($verbose);
exit(0);




#################################################################################
#                                                                               #
#                          Subroutines                                          #
#                                                                               #
#################################################################################

sub usage {
    my $error = shift;
    
    if ($error == 1) {
	# No data-type choosen
	print "\nNo data option choosen [-est|-mrna|-ost|-nematode|-ost|-washu|-nembase]\n";
	print "Run with one of the above options\n\n";
	exit(0);
    }
    if ($error == 2) {
	# 'Multiple data-types choosen
	print "\nMultiple data option choosen [-est|-mrna|-ost|-nematode|-embl|-washu|-nembase]\n";
	print "Run with one of the above options\n\n";
	exit(0);
    }
    if ($error == 3) {
	# 'chromosome.ace' file is not there or unreadable
	print "\nThe WormBase 'chromosome.ace' file you specified does not exist or is non-readable.\n";
	print "Check File: ''\n\n";
	exit(0);
    }
    if ($error == 20) {
	# 
	print "\nDon't want to do this for the -nematode or -washu or -nembase options.\n";
	print "hasta luego\n\n";
	exit(0);
    }
    elsif ($error == 0) {
	# Normal help menu
	exec ('perldoc',$0);
    }
}



__END__

=pod

=head1 NAME - blat2ace.pl

=head2 USAGE

blat2ace.pl maps blat output to acedb. Thereby, it produces output for autoace and camace
(autoace.ace and camace,ace). In addition, it produces files assigning the ESTs to one place 
in the genome (autoace.best.ace and camace.best.ace). ESTs that have more than one best 
match are reported in morethan1match.txt. 

blat2ace.pl  arguments:

=over 4

=item 

-camace => produce output for camace (camace.blat.ace, /helpfiles/camace.best.ace, /helpfiles/camace.ace)

=item 

-intron => produce output for confirmed introns (autoace.ci.ace, camace.ci.ace)

=item 

-mrna => perform everything for mRNAs

=item 

-ncrna => perform everything for ncRNAs

=item

-est => perform everything for ESTs

=item

-ost => perform everything for OSTs

=item

-nematode => perform everything for non-C. elegans ESTs

=item

-nembase => perform everything for NemBase contigs

=item

-washu => perform everything for WashU Nematode.net - John Martin <jmartin@watson.wustl.edu> contigs

=item

-embl => perform everything for non-WormBase CDSs in EMBL

=back

=head1 AUTHOR

Kerstin Jekosch (kj2@sanger.ac.uk)

=cut

