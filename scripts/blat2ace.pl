#/software/bin/perl -w
#
# blat2ace.pl
# 
# by Kerstin Jekosch
#
# Exporter to map blat data to genome and to find the best match for each EST, mRNA, OST, etc.
#
# Last edited by: $Author: gw3 $
# Last edited on: $Date: 2009-05-06 11:13:22 $

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use File::Copy;

#########################
# Command line options  #
#########################

my ($est, $mrna, $ncrna, $ost, $tc1, $nematode, $embl, $camace, $intron, $washu, $nembase);
my ($help, $debug, $test, $verbose, $store, $wormbase, $species);
my ($database, $virtualobjs, $type, $qspecies);
my ($acefile,$pslfile);

GetOptions (
	    "help"       => \$help,
	    "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
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
	    "virtual"    => \$virtualobjs,
	    "type:s"     => \$type,
	    "qspecies:s" => \$qspecies, #query species
	    "ace:s"      => \$acefile,
	    "psl:s"      => \$pslfile
	   );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species
			     );
}
# establish log file.
my $log = Log_files->make_build_log($wormbase);


#############################
# variables and directories #
#############################

# set database paths, default to autoace unless -camace
my $ace_dir   = $wormbase->orgdb;
my $blat_dir = (defined $database ? $database."/BLAT" : $wormbase->blat);
my @nematodes = qw(nematode washu nembase);

#############################
# CommonData hash retrieval #
#############################
my %NDBaccession2est;
my %estorientation;
if ($qspecies eq 'elegans'){ #this should change if other species need this too.
# if we're processing elegans cDNAs we need to use the elegans common data to convert acc -> yk
	my $accessor = $wormbase;
	unless ($wormbase->species eq 'elegans') {
		$accessor = Wormbase->new( 	-debug   => $wormbase->debug,
                             		-test    => $wormbase->test,
                            		-organism => 'elegans'
                             );
    }
	%NDBaccession2est = $accessor->FetchData('NDBaccession2est');     # EST accession => name
	%estorientation   = $accessor->FetchData('estorientation');       # EST accession => orientation [5|3]
}
my %hash;
my (%best,%other,%bestclone,%match,%ci);

$log->log_and_die("no type specified\n") unless $type;

##########################################################################################
# map the blat hits to ace - i.e. process blat output (*.psl) file into set of ace files #
##########################################################################################
$log->write_to($wormbase->runtime.": Start mapping $qspecies $type\n\n");

# open input and output filehandles
$acefile = "$blat_dir/".$wormbase->species.".${qspecies}_$type.ace" unless $acefile;
$pslfile = "$blat_dir/${qspecies}_${type}_out.psl" unless $pslfile;
open(ACE,  ">$acefile")  or die "Cannot open $acefile $!\n";
open(BLAT, "<$pslfile")    or die "Cannot open $pslfile $!\n";

my $number_of_replacements = 0;
my %reported_this_query_before;
my %make_virt_obj;
# loop through each blat hit
my $method = $qspecies eq $wormbase->species ? "BLAT_${type}_OTHER" : "BLAT_Caen_${type}_OTHER";
$method = "BLAT_".uc($qspecies) if grep(/$qspecies/, @nematodes);

while (<BLAT>) {

  next unless (/^\d/);
  my @f            = split "\t";

  my $match        = $f[0];                    # number of bases matched by blat
  my $strand       = $f[8];                    # strand that match is on
  my $query        = $f[9];                    # query sequence name
  my $query_size   = $f[10];                   # query sequence length
  my $superlink    = $f[13];                   # name of superlink that was used as blat target sequence
  my $slsize       = $f[14];                   # target seq size (used to be superlink hence sl)
  my $lastvirt     = int($slsize/100000) + 1;  # for tracking how many virtual sequences have been created???
  my $matchstart   = $f[15];                   # target (superlink) start coordinate...
  my $matchend     = $f[16];                   # ...and end coordinate
  my $block_count  = $f[17];                   # block count
  my @lengths      = split (/,/, $f[18]);      # sizes of each blat 'block' in any individual blat match
  my @query_starts = split (/,/, $f[19]);      # start coordinates of each query block
  my @slink_starts = split (/,/, $f[20]);      # start coordinates of each target (superlink) block

  next if ($strand eq '--');
  $make_virt_obj{$superlink} = $slsize  if($virtualobjs);# store which sequence to make virtual objects for 
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
      $virtual = "BLAT_${type}:${superlink}_${startvirtual}";
  }	
  elsif (($startvirtual == ($endvirtual - 1)) && (($matchend%100000) <= 50000)) {
      $virtual = "BLAT_${type}:${superlink}_${startvirtual}";
  }
  else {
    if (! exists $reported_this_query_before{$query}) {
      $reported_this_query_before{$query} = 1; # don't want to report this one again
    # $log->write_to("$query wasn't assigned to a virtual object as match size was too big\n") if $wormbase->debug;
    #  $log->write_to("Start is $matchstart, end is $matchend on $superlink\n\n") if $wormbase->debug;
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
	  #$log->write_to("// MISMATCH: $query [$strand] $virtualstart $virtualend :: [virtual slice $calc -> $newcalc, offset ".($matchend%100000)."]\n\n");
      }

      if (!defined $virtualstart) {
	  #$log->write_to("$query will be discarded as the match is too long\n");
	  #$log->write_to("$query [$strand] $virtualstart $virtualend  [virtual slice $calc -> $newcalc, offset ".($matchend%100000)."]\n\n");
	  next;
      }

      my ($query_start,$query_end);
      
        # blatx 6-frame translation v 6-frame translation
#      if ($nematode || $washu || $nembase) {
	  if( $wormbase->species ne $qspecies) {
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
	  elsif ($strand eq '+-') {
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
	  printf ACE "DNA_homol\t\"%s\"\t\"BLAT_${type}\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$query,$score,$virtualstart,$virtualend,$query_start,$query_end;
      }
      else {
	printf ACE "DNA_homol\t\"%s\"\t$method\t%.1f\t%d\t%d\t%d\t%d\n\n",$query,$score,$virtualstart,$virtualend,$query_start,$query_end;
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
&make_virt_objs($type) if $virtualobjs;

# concise report
if ($est || $ost) {
  $log->write_to("\nThere were $number_of_replacements replacements of EST names\n\n");
}

####################################
# produce outfile for best matches #
####################################
if ($nematode || $washu || $nembase) {
  #$wormbase->run_command("mv $blat_dir/".$wormbase->species.".$type.ace $blat_dir/".$wormbase->species.".blat.$type.ace", $log);
  my $fromace = "$blat_dir/".$wormbase->species.".$type.ace";
  my $toace = "$blat_dir/".$wormbase->species.".blat.$type.ace";
  move($fromace, $toace) || do {$log->write_to("ERROR: move of $fromace to $toace failed: $!\n"); $log->error};
} else {
  my $no_direction = 0;		# count of transcripts with no specified orientation

  open (AUTBEST, ">$blat_dir/".$wormbase->species.".best.${qspecies}_$type.ace");
  my $method = $qspecies eq $wormbase->species ? "BLAT_${type}_BEST" : "BLAT_Caen_${type}_BEST";
  $method = "BLAT_".uc($qspecies) if grep(/$qspecies/, @nematodes);
  foreach my $found (sort keys %best) {
    if (exists $best{$found}) {
      foreach my $entry (@{$best{$found}->{'entry'}}) {
# comment this out to allow two or more $entry's with equal scores to all be 'BEST'
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
	    printf AUTBEST "DNA_homol\t\"%s\"\t$method\t%.1f\t%d\t%d\t%d\t%d\n\n",$found,$score,$virtualstart,$virtualend,$query_start,$query_end;
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
	      if ($type eq 'mRNA' && !exists $estorientation{$found}) {	# default orientation of mRNA is forwards
		$estorientation{$found} = 5;
	      }
	      if (${$entry->{'exons'}}[0][2] < ${$entry->{'exons'}}[0][3]) {
		if ((${$entry->{'exons'}}[$y][2] == ${$entry->{'exons'}}[$last][3] + 1) && (($second - $first) > 2)) {
		  if (exists $estorientation{$found} && $estorientation{$found} eq '3') {
		    push @{$ci{$superlink}}, [$second,$first,$found];
		  } elsif (exists $estorientation{$found} && $estorientation{$found} eq '5') {
		    push @{$ci{$superlink}}, [$first,$second,$found];
		  } else {
		    $no_direction++
		  }
		}
	      } elsif (${$entry->{'exons'}}[0][2] > ${$entry->{'exons'}}[0][3]) {
		if ((${$entry->{'exons'}}[$last][3] == ${$entry->{'exons'}}[$y][2] + 1) && (($second - $first) > 2)) {
		  if (exists $estorientation{$found} && $estorientation{$found} eq '3') {
		    push @{$ci{$superlink}}, [$first,$second,$found];
		  } elsif (exists $estorientation{$found} && $estorientation{$found} eq '5') {
		    push @{$ci{$superlink}}, [$second,$first,$found]; 
		  } else {
		    $no_direction++
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  $log->write_to("WARNING: Direction not found for $no_direction transcripts\n\n") if ($no_direction);
  close(AUTBEST);
}
########################################################
# produce final BLAT output (including BEST and OTHER) #
########################################################

unless ($nematode || $washu || $nembase) {
  # Open new (final) output files for autoace, camace, and stlace

  open (OUT_autoace, ">$blat_dir/".$wormbase->species.".blat.${qspecies}_$type.ace") or $log->log_and_die("cant open $blat_dir/".$wormbase->species.".blat.${qspecies}_$type.ace :$!");

  # Change input separator to paragraph mode, but store what it old mode in $oldlinesep
  my $oldlinesep = $/;
  $/ = "";

  my (%line);
  my $superlink = "";

  # assign 
  open(ABEST,  "<$blat_dir/".$wormbase->species.".best.${qspecies}_$type.ace") or $log->log_and_die("cant open $blat_dir/".$wormbase->species.".best.${qspecies}_$type.ace $!\n");

  while (<ABEST>) {
    if ($_ =~ /^Homol_data/) {
      # flag each blat hit which is best (all of them) - set $line{$_} to 1
      # %line thus stores keys which are combos of virtual object name + blat hit details
      $line{$_} = 1;
      ($superlink) = (/\"BLAT_${type}\:(\S+)\_\d+\"/);

      # Print blat best hits to final output file
      print OUT_autoace "// Source $superlink\n\n";
      print OUT_autoace $_;
    }
  }
  close ABEST;

  # we expect to have BEST hits in elegans
  if ($qspecies eq 'elegans' && scalar (keys %line) == 0) {
    $log->error;
    $log->write_to("ERROR: we don't have any BEST $type in:\n$blat_dir/".$wormbase->species.".best.${qspecies}_$type.ace\n\n");
  }
  $log->write_to("Have " . scalar (keys %line) . " BEST $type in:\n$blat_dir/".$wormbase->species.".best.${qspecies}_$type.ace\n");

  # Now look through original output file (where everything is set to BLAT_OTHER) to
  # output those blat OTHER hits which are not flagged as BLAT_BEST in the .best.ace file
  # Does this by comparing entries in %line hash

  my $count_other_hits = 0;
  my $count_total_hits = 0;
  open(AOTHER, "<$blat_dir/".$wormbase->species.".${qspecies}_$type.ace");
  while (<AOTHER>) {
    if ($_ =~ /^Homol_data/) {
      my $line = $_;
      if ($line =~ /OTHER/) {$count_total_hits++};
      # for comparison to %line hash, need to change OTHER to BEST in $_
      s/OTHER/BEST/g;
      # Only output BLAT_OTHER hits in first output file which we now know NOT to
      # really be BEST hits
      unless (exists $line{$_}) {
	print OUT_autoace $line;
	$count_other_hits++;
      }	
    }
  }
  close AOTHER;

  # we expect to have OTHER hits in elegans
  if ($qspecies eq 'elegans' && $count_other_hits == 0) {
    $log->error;
    $log->write_to("ERROR: we don't have any OTHER $type in:\n$blat_dir/".$wormbase->species.".${qspecies}_$type.ace\n");
  }
  $log->write_to("Have $count_other_hits OTHER $type in:\n$blat_dir/".$wormbase->species.".${qspecies}_$type.ace\n");
  $log->write_to("Have a total of $count_total_hits (BEST & OTHER) $type in:\n$blat_dir/".$wormbase->species.".${qspecies}_$type.ace\n");

  # reset input line separator
  $/= $oldlinesep;

  ###################################
  # produce confirmed intron output #
  ###################################

  if ($intron) {
    open(CI_auto, ">$blat_dir/".$wormbase->species.".ci.${qspecies}_${type}.ace");

    foreach my $link (sort keys %ci) {
      my %double;
      
      print CI_auto "\nSequence : \"$link\"\n";
      
      for (my $i = 0; $i < @{$ci{$link}}; $i++) {
	my $merge = $ci{$link}->[$i][0].":".$ci{$link}->[$i][1];
	if (!exists $double{$merge}) {
	  # If RST or OST modify output.
	  if (($type eq "RST") or ($type eq "OST")) {
	    printf CI_auto "Confirmed_intron %d %d EST $ci{$link}->[$i][2]\n",  $ci{$link}->[$i][0], $ci{$link}->[$i][1];
	  }
	  elsif ($type eq "mRNA") {
	    printf CI_auto "Confirmed_intron %d %d cDNA $ci{$link}->[$i][2]\n",  $ci{$link}->[$i][0], $ci{$link}->[$i][1];
	  }
	  else {
	    printf CI_auto "Confirmed_intron %d %d $type $ci{$link}->[$i][2]\n",  $ci{$link}->[$i][0], $ci{$link}->[$i][1];
	  }
	}
	$double{$merge} = 1;
      }
    }
    
    close CI_auto;
  }
}


#compress acefiles so that all object data is loaded together, expanded to run on all homol

my @filenames = ($wormbase->species.".${qspecies}_$type.ace", $wormbase->species.".best.${qspecies}_$type.ace", $wormbase->species.".blat.${qspecies}_$type.ace", );
my $filename;
$log->write_to("\n#########################################\nCompressing DNA_homolo acefiles\n#########################################\n");
foreach $filename (@filenames) {
  # $wormbase->run_command("/bin/mv $blat_dir/$filename $blat_dir/${filename}"."_uncompressed",$log);
  my $fromace = "$blat_dir/$filename";
  my $toace = "$blat_dir/${filename}_uncompressed";
  move($fromace, $toace) || $log->write_to("WARNING: move of $fromace to $toace failed: $!\n");

  if (-e $toace ) {
    $log->write_to("Compressing $toace\n");
    $wormbase->run_script("acezip.pl -file $toace -build", $log);
    $log->write_to("Compressed........\n");
  }
}

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

#############################
# virtual object generation #
#############################

sub make_virt_objs {
    
  my $data = shift;
  local (*OUT_autoace_homol);
  local (*OUT_autoace_feat);
  my ($name,$length,$total,$first,$second,$m,$n);
  
  # autoace
  open (OUT_autoace_homol, ">$blat_dir/virtual_objects.".$wormbase->species.".blat.$data.$qspecies.ace") or die "$!";
  open (OUT_autoace_feat,  ">$blat_dir/virtual_objects.".$wormbase->species.".ci.$data.$qspecies.ace")   or die "$!";
  
  foreach my $name (keys %make_virt_obj) {
	$length = $make_virt_obj{$name};
    $total = int($length/100000) +1;
      # autoace
    print OUT_autoace_homol "Sequence : \"$name\"\n";
    print OUT_autoace_feat  "Sequence : \"$name\"\n";

    for ($n = 0; $n <= $total; $n++) {
		$m      = $n + 1;
		$first  = ($n*100000) + 1;
		$second = $first + 149999;
		if (($length - $first) < 100000) {
			$second = $length;
		  # autoace
	  		print OUT_autoace_homol "S_child Homol_data BLAT_$data:$name"."_$m $first $second\n";
	  		print OUT_autoace_feat  "S_child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n";

	  		last;
		}					
		else {
	  		($second = $length) if ($second >  $length);
	  		# autoace
	  		print OUT_autoace_homol "S_child Homol_data BLAT_$data:$name"."_$m $first $second\n";
	  		print OUT_autoace_feat  "S_child Feature_data Confirmed_intron_$data:$name"."_$m $first $second\n";
		}
	}
    print OUT_autoace_homol "\n";
    print OUT_autoace_feat  "\n";
  }
  close OUT_autoace_homol;
  close OUT_autoace_feat;

}


__END__

=pod

=head1 NAME - blat2ace.pl

=head2 USAGE

blat2ace.pl maps blat output to acedb. Thereby, it produces output for autoace and camace
(species.ace and camace,ace). In addition, it produces files assigning the ESTs to one place 
in the genome (species.best.ace and camace.best.ace). ESTs that have more than one best 
match are reported in morethan1match.txt. 

blat2ace.pl  arguments:

=over 4

=item 

-camace => produce output for camace (camace.blat.ace, /helpfiles/camace.best.ace, /helpfiles/camace.ace)

=item 

-intron => produce output for confirmed introns (species.ci.ace, camace.ci.ace)

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

