#!/usr/local/bin/perl5.8.0 -w
#
# blat2ace.pl
# 
# by Kerstin Jekosch
#
# Exporter to map blat data to genome and to find the best match for each EST, mRNA, OST, etc.
#
# Last edited by: $Author: krb $
# Last edited on: $Date: 2003-09-04 15:32:11 $


use strict;
use Data::Dumper;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Getopt::Long;

#########################
# Command line options  #
#########################

my ($help, $est, $mrna, $ost, $nematode, $embl, $camace, $intron);

GetOptions ("help"       => \$help,
            "est"        => \$est,
            "mrna"       => \$mrna,
            "ost"        => \$ost,
            "nematode"   => \$nematode,
            "embl"       => \$embl,
            "camace"     => \$camace,
	    "intron"     => \$intron
);

#############################
# variables and directories #
#############################

our $log;


# set database paths, default to autoace unless -camace
my $blat_dir  = "/wormsrv2/autoace/BLAT";
my $tace      = &tace ." /wormsrv2/autoace";
if ($camace) {
    $blat_dir  = "/wormsrv1/camace/BLAT";
    $tace      = &tace." /wormsrv1/camace";
}

my %EST_name;    # EST accession => name
my %EST_dir;     # EST accession => orientation [5|3]

my %hash;
my (%best,%other,%bestclone,%match,%ci);

our %camace;
our %stlace;

our $type = "";
our %word = (
	     EST      => 'BLAT_EST',
	     OST      => 'BLAT_OST',
	     mRNA     => 'BLAT_mRNA',
	     EMBL     => 'BLAT_EMBL',
	     NEMATODE => 'BLATX_NEMATODE',
	     );



&create_log_files;


########################################
# command-line options & ramifications #
########################################

# Help pod documentation
&usage(0) if ($help);

# Exit if no data type choosen [EST|mRNA|EMBL|NEMATODE|OST]
# or if multiple data types are chosen

&usage(1) unless ($est || $mrna || $ost || $nematode || $embl); 
&usage(2) if (($est + $mrna + $ost + $nematode || $embl) > 1);

# assign type variable
($type = 'EST')      if ($est);
($type = 'OST')      if ($ost);
($type = 'mRNA')     if ($mrna);
($type = 'EMBL')     if ($embl);
($type = 'NEMATODE') if ($nematode);


############################################
# EST data from autoace (name,orientation) #
############################################

# check to see if EST hash data exists, make it via tablemaker queries if absent
# else read it into memory
unless (-e "$blat_dir/EST.dat") {
    (%EST_name,%EST_dir) = &make_EST_hash;
}
else {
  open (FH, "<$blat_dir/EST.dat") or die "EST.dat : $!\n";
  undef $/;
  my $data = <FH>;
  eval $data;
  die if $@;
  $/ = "\n";
  close FH;
}

#########################################
# get links for database                #
#########################################

# parse links for camace
my @camclones = qw(cTel3X cTel4X cTel7X cTel33B cTel54X 6R55 SUPERLINK_CB_I SUPERLINK_CB_II SUPERLINK_CB_IIIL SUPERLINK_CB_IIIR SUPERLINK_CB_IR SUPERLINK_CB_IV SUPERLINK_CB_V SUPERLINK_CB_X); 
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

my $runtime = &runtime;
print LOG "$runtime: Start mapping\n\n";

# output filehandle
open(ACE,  ">$blat_dir/autoace.$type.ace")  or die "Cannot open $blat_dir/autoace.${type}.ace $!\n";

# input filehandle
open(BLAT, "<$blat_dir/${type}_out.psl")  or die "Cannot open $blat_dir/${type}_out.psl $!\n";

while (<BLAT>) {
  next unless (/^\d/);
  my @f           = split "\t";
  my $match       = $f[0];                    # number of bases matched by blat
  my $query       = $f[9];                    # query sequence name
  my $query_size  = $f[10];                   # query sequence length
  my $superlink   = $f[13];                   # name of superlink that was used as blat target sequence
  my $slsize      = $f[14];                   # superlink size
  my $lastvirt    = int($slsize/100000) + 1;  # for tracking how many virtual sequences have been created???
  my $matchstart  = $f[15];                   # target (superlink) start coordinate...
  my $matchend    = $f[16];                   # ...and end coordinate
  my @lengths     = split (/,/, $f[18]);      # sizes of each blat 'block' in any individual blat match
  my @querystarts = split (/,/, $f[19]);      # start coordinates of each query block
  my @slinkstarts = split (/,/, $f[20]);      # start coordinates of each target (superlink) block



  # replace EST name (usually accession number) by yk... name 
  if (( $est || $ost )  && (exists $EST_name{$query})) {
    my $estname  = $EST_name{$query};
    if ($query ne $estname) {
      print LOG "EST name '$query' was replaced by '$estname'\n\n";
      $query = $estname;
    }
  }


  ###############################
  # find virtual superlink part #
  ###############################
	
  my ($virtual,$startvirtual,$endvirtual);
  if ((int($matchstart/100000) +1) > $lastvirt) { $startvirtual = $lastvirt;}
  else {$startvirtual = int($matchstart/100000) +1;}  
    
  if ((int($matchend/100000) +1) > $lastvirt) { $endvirtual = $lastvirt;}
  else {$endvirtual = int($matchend/100000) +1;}  
  
  if ($startvirtual == $endvirtual) {
    $virtual = "$word{$type}:${superlink}_${startvirtual}";
  }	
  elsif (($startvirtual == ($endvirtual - 1)) && (($matchend%100000) <= 50000)) {
    $virtual = "$word{$type}:${superlink}_${startvirtual}";
  }
  else {
    print LOG "$query wasn't assigned to a virtual object as match size was too big\n";
    print LOG "Start is $matchstart, end is $matchend on $superlink\n\n";
    next;
  }
    
  ###################
  # calculate score #
  ###################
  
  my $sum   = 0;
  foreach my $length (@lengths) {
    $sum = $sum + $length;
  }	

  # new way of calculating score, divide by query size rather than sum of matching blocks, 
  # so short 30 bp polyA chunks will get 30/2000 rather than 30/30, should improve things
  # that are currently assigned OTHER but should be BEST
  my $score = ($match/$query_size)*100;
  
  my @exons = ();
  
  #########################
  # calculate coordinates #
  #########################
    
  # need to allow for est exons in the next virtual object, otherwise they get remapped to the start 
  # of the virtual by performing %100000
  
  my $calc = int(($slinkstarts[0]+1)/100000);
  
  for (my $x = 0;$x < $f[17]; $x++) {
    my $newcalc      = int(($slinkstarts[$x]+1)/100000);
    my $virtualstart;
    if ($calc == $newcalc) {	
      $virtualstart = ($slinkstarts[$x] +1)%100000;
    }
    elsif ($calc == ($newcalc-1)) {
      $virtualstart = (($slinkstarts[$x] +1)%100000) + 100000;
    }
    my $virtualend   = $virtualstart + $lengths[$x] -1;
    my ($query_start,$query_end);
    
    # blatx 6-frame translation v 6-frame translation
    if ($nematode) {
      my $temp;
      if (($f[8] eq '++') || ($f[8] eq '-+')) {
	$query_start   = $querystarts[$x] +1;
	$query_end     = $query_start + $lengths[$x] -1;
	if ($f[8] eq '-+') {
	  $temp     = $query_end;
	  $query_end   = $query_start;
	  $query_start = $temp; 
	}
      }
      elsif (($f[8] eq '--') || ($f[8] eq '+-')) {
	$temp         = $virtualstart;
	$virtualstart = $virtualend;
	$virtualend   = $temp;
	$query_start     = $f[10] - $querystarts[$x];
	$query_end       = $query_start - $lengths[$x] +1;
	if ($f[8] eq '--') {
	  $temp     = $query_end;
	  $query_end   = $query_start;
	  $query_start = $temp; 
	}
      }			
    }
    else {
      if ($f[8] eq '+'){
	$query_start   = $querystarts[$x] +1;
	$query_end     = $query_start + $lengths[$x] -1;
      }
      elsif ($f[8] eq '-') {
	$query_start   = $f[10] - $querystarts[$x];
	$query_end     = $query_start - $lengths[$x] +1;
      }		
    }		
    print LOG "$query was mapped to $virtual\n\n";
    
    # write to output file
    print  ACE "Homol_data : \"$virtual\"\n";
    if ($type eq "NEMATODE") {
      printf ACE "DNA_homol\t\"%s\"\t\"$word{$type}\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$query,$score,$virtualstart,$virtualend,$query_start,$query_end;
    }
    else {
      printf ACE "DNA_homol\t\"%s\"\t\"$word{$type}_OTHER\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$query,$score,$virtualstart,$virtualend,$query_start,$query_end;
    }
    push @exons, [$virtualstart,$virtualend,$query_start,$query_end];				
  }
    
  ########################
  # collect best matches #
  ########################
  
  if (exists $best{$query}) {
    if ($score >= $best{$query}->{'score'}) {
      if ( ($score > $best{$query}->{'score'}) || ($match > $best{$query}->{'match'})) { 
	$best{$query}->{'score'} = $score;
	$best{$query}->{'match'} = $match;
	@{$best{$query}->{'entry'}} = ({'clone' => $virtual,'link' => $superlink,'exons' => \@exons});
      }
      elsif ($match == $best{$query}->{'match'}) {
	$best{$query}->{'score'} = $score;
	push @{$best{$query}->{'entry'}}, {'clone' => $virtual,'link' => $superlink,'exons' => \@exons};
      }
    }
  }
  else {
    $best{$query}->{'match'} = $match;
    $best{$query}->{'score'} = $score;
    @{$best{$query}->{'entry'}} = ({'clone' => $virtual,'link' => $superlink,'exons' => \@exons});
  }
}
close BLAT;
close ACE;

#########################################
# 
#########################################


####################################
# produce outfile for best matches #
####################################

&usage(20) if ($nematode);

open (AUTBEST, ">$blat_dir/autoace.best.$type.ace");
open (STLBEST, ">$blat_dir/stlace.best.$type.ace");
open (CAMBEST, ">$blat_dir/camace.best.$type.ace");

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
	  
	  # print line for autoace
	  print  AUTBEST "Homol_data : \"$virtual\"\n";
	  printf AUTBEST "DNA_homol\t\"%s\"\t\"$word{$type}_BEST\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$found,$score,$virtualstart,$virtualend,$query_start,$query_end;
	  # camace
	  if ($camace{$superlink}) {
	    print  CAMBEST "Homol_data : \"$virtual\"\n";
	    printf CAMBEST "DNA_homol\t\"%s\"\t\"$word{$type}_BEST\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$found,$score,$virtualstart,$virtualend,$query_start,$query_end;
	  }
	  # and stlace
	  elsif ($stlace{$superlink}) {
	    print  STLBEST "Homol_data : \"$virtual\"\n";
	    printf STLBEST "DNA_homol\t\"%s\"\t\"$word{$type}_BEST\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$found,$score,$virtualstart,$virtualend,$query_start,$query_end;
	  }
	  
	}
	
	#############################
	# produce confirmed introns #
	#############################
	
	# this section 'produce confirmed introns'
	if ($intron) {
	  print LOG "Producing confirmed introns\n";
	  my ($n) = ($virtual =~ /\S+_(\d+)$/);
	  for (my $y = 1; $y < @{$entry->{'exons'}}; $y++) {
	    my $last   = $y - 1;
	    my $first  =  (${$entry->{"exons"}}[$last][1] + 1) + (($n-1)*100000);
	    my $second =  (${$entry->{'exons'}}[$y][0]    - 1) + (($n-1)*100000);
	    $EST_dir{$found} = 5 if ($mrna || $embl);
	    if (${$entry->{'exons'}}[0][2] < ${$entry->{'exons'}}[0][3]) {
	      if ((${$entry->{'exons'}}[$y][2] == ${$entry->{'exons'}}[$last][3] + 1) && (($second - $first) > 2)) {
		if (exists $EST_dir{$found} && $EST_dir{$found} eq '3') {
		  push @{$ci{$superlink}}, [$second,$first];
		}
		elsif (exists $EST_dir{$found} && $EST_dir{$found} eq '5') {
		  push @{$ci{$superlink}}, [$first,$second];
		}
		else {
		  print LOG "WARNING: Direction not found for $found\n\n";
		}
	      }
	    }
	    elsif (${$entry->{'exons'}}[0][2] > ${$entry->{'exons'}}[0][3]) {
	      if ((${$entry->{'exons'}}[$last][3] == ${$entry->{'exons'}}[$y][2] + 1) && (($second - $first) > 2)) {
		if (exists $EST_dir{$found} && $EST_dir{$found} eq '3') {
		  push @{$ci{$superlink}}, [$first,$second];
		}
		elsif (exists $EST_dir{$found} && $EST_dir{$found} eq '5') {
		  push @{$ci{$superlink}}, [$second,$first]; 
		}
		else {
		  print LOG "WARNING: Direction not found for $found\n\n";
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

########################################################
# produce final BLAT output (including BEST and OTHER) #
########################################################

&usage(20) if ($nematode);

# autoace
open (OUT_autoace, ">$blat_dir/autoace.blat.$type.ace") or die "$!";
# camace
open (OUT_camace,  ">$blat_dir/camace.blat.$type.ace")  or die "$!";
# stlace
open (OUT_stlace,  ">$blat_dir/stlace.blat.$type.ace")  or die "$!";


my (%line);
my $temp = $/;
$/ = "";

my $superlink = "";

# assign 
open(ABEST,  "<$blat_dir/autoace.best.$type.ace");
while (<ABEST>) {
    if ($_ =~ /^Homol_data/) {
	$line{$_} = 1;
	($superlink) = (/\"BLAT\_$type\:(\S+)\_\d+\"/);

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


open(AOTHER, "<$blat_dir/autoace.$type.ace");
while (<AOTHER>) {
  if ($_ =~ /^Homol_data/) {
    my $line = $_;
    s/BLAT_EST_OTHER/BLAT_EST_BEST/g unless ($mrna || $embl || $nematode || $ost);
    s/BLAT_OST_OTHER/BLAT_OST_BEST/g     if ($ost); 
    s/BLAT_mRNA_OTHER/BLAT_mRNA_BEST/g   if ($mrna);
    s/BLAT_EMBL_OTHER/BLAT_EMBL_BEST/g   if ($embl);
    
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

$/= $temp;

###################################
# produce confirmed intron output #
###################################

if ($intron) {
  
  open(CI_auto, ">$blat_dir/autoace.ci.${type}.ace");
  open(CI_cam,  ">$blat_dir/camace.ci.${type}.ace");
  open(CI_stl,  ">$blat_dir/stlace.ci.${type}.ace");
  
  foreach my $superlink (sort keys %ci) {
    my %double;
    
    print CI_auto "\nSequence : \"$superlink\"\n";
    print CI_stl  "\nSequence : \"$superlink\"\n" if ($stlace{$superlink});
    print CI_cam  "\nSequence : \"$superlink\"\n" if ($camace{$superlink});
    
    for (my $i = 0; $i < @{$ci{$superlink}}; $i++) {
      my $merge = $ci{$superlink}->[$i][0].":".$ci{$superlink}->[$i][1];
      if (!exists $double{$merge}) {
	if ($mrna) {
	  printf CI_auto "Confirmed_intron %d %d mRNA\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1];
	  (printf CI_cam "Confirmed_intron %d %d mRNA\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1]) if ($camace{$superlink});
	  (printf CI_stl "Confirmed_intron %d %d mRNA\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1]) if ($stlace{$superlink});
	}
	if ($embl) {
	  printf CI_auto "Confirmed_intron %d %d Homol\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1];
	  (printf CI_cam "Confirmed_intron %d %d Homol\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1]) if ($camace{$superlink});
	  (printf CI_stl "Confirmed_intron %d %d Homol\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1]) if ($stlace{$superlink});
	}
	if ($est) {
	  printf CI_auto "Confirmed_intron %d %d EST\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1];
	  (printf CI_cam "Confirmed_intron %d %d EST\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1]) if ($camace{$superlink});
	  (printf CI_stl "Confirmed_intron %d %d EST\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1]) if ($stlace{$superlink});
	}
	if ($ost) {
	  printf CI_auto "Confirmed_intron %d %d EST\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1];
	  (printf CI_cam "Confirmed_intron %d %d EST\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1]) if ($camace{$superlink});
	  (printf CI_stl "Confirmed_intron %d %d EST\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1]) if ($stlace{$superlink});
	}
	$double{$merge} = 1;
      }
    }
  }
  
  close CI_auto;
  close CI_cam;
  close CI_stl;

}

##############################
# hasta luego                #
##############################

close(LOG);
exit(0);




#################################################################################
#                                                                               #
#                          Subroutines                                          #
#                                                                               #
#################################################################################


#########################################
# get EST names  (-e option only)       #
#########################################

sub commands {
    
my $command1=<<EOF;
Table-maker -p "/wormsrv2/autoace/wquery/ESTacc2names.def"
quit
EOF

my $command2=<<EOF;
Table-maker -p "/wormsrv2/autoace/wquery/ESTorient.def"
quit
EOF

    return($command1,$command2);

}

sub make_EST_hash {
    
  my ($command1,$command2) = &commands;
  my ($acc,$name,$orient);
  
  my %EST_name = ();
  my %EST_dir  = ();
  
  # get EST names  (-e option only)       #
  open (TACE, "echo '$command1' | $tace | ");
  while (<TACE>) {
    chomp;
    next if ($_ eq "");
    next if (/\/\//);
    s/acedb\>\s//g;
    s/\"//g;
    s/EMBL://g;
    ($acc,$name) = ($_ =~ /^(\S+)\s(\S+)/);
    $name = $acc unless ($name);
    $EST_name{$acc} = $name;
  }
  close TACE;
  
  # get EST orientation (5' or 3')    #
  open (TACE, "echo '$command2' | $tace | ");
  while (<TACE>) {
    chomp;
    next if ($_ eq "");
    next if (/\/\//);
    s/acedb\>\s//g;
    s/\"//g;
    ($name,$orient) = ($_ =~ /^(\S+)\s+EST_(\d)/);
    $EST_dir{$name} = $orient if ($orient);
  }
  close TACE;
  
  # Data::Dumper write hash to /wormsrv2/autoace/BLAT/EST.dat
  open (OUT, ">/wormsrv2/autoace/BLAT/EST.dat") or die "EST.dat : $!";
  print OUT Data::Dumper->Dump([\%EST_name],['*EST_name']);
  print OUT Data::Dumper->Dump([\%EST_dir],['*EST_dir']);
  close OUT;
  
  return (%EST_name,%EST_dir);
  
  
}





###################################

sub usage {
    my $error = shift;
    
    if ($error == 1) {
	# No data-type choosen
	print "\nNo data option choosen [-est|-mrna|-ost|-nematode|-ost]\n";
	print "Run with one of the above options\n\n";
	exit(0);
    }
    if ($error == 2) {
	# 'Multiple data-types choosen
	print "\nMultiple data option choosen [-est|-mrna|-ost|-nematode|-embl]\n";
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
	print "\nDon't want to do this for the -nematode option.\n";
	print "hasta luego\n\n";
	exit(0);
    }
    elsif ($error == 0) {
	# Normal help menu
	exec ('perldoc',$0);
    }
}


##############################################################

sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`
");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $WS_version = &get_wormbase_version_name;
  my $rundate = &runtime;
  $log        = "/wormsrv2/logs/$script_name.${WS_version}.$rundate.$$";

  open (LOG, ">$log") or die "cant open $log";
  print LOG "$script_name\n";
  print LOG "started at ",`date`,"\n";
  print LOG "=============================================\n";
  print LOG "\n";

}


###################################

#
# -i  : get confirmed introns
#
# -c  : get best matches for camace
# -s  : get best matches for stlace
#
# -e  : create output for ESTs 
# -y  : create output for OSTs 
# -m  : create output for mRNAs 
# -x  : create output for parasitic nematode ESTs (blatx)
# -o  : create output for other CDS
#
# -h  : print help



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

=back

=item

-est => perform everything for ESTs

=back

=item

-ost => perform everything for OSTs

=back

=item

-nematode => perform everything for non-C. elegans ESTs

=back

=item

-embl => perform everything for non-WormBase CDSs in EMBL

=back

=head1 AUTHOR

Kerstin Jekosch (kj2@sanger.ac.uk)

=cut

