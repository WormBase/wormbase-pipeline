#!/usr/local/bin/perl5.8.0 -w
#
# confirm_genes.pl
#
# by Kerstin Jekosch
#
# checks gene models for EST and cDNA evidence and adds a tag to the sequence object 
# if the whole model is confirmed
#
# Last updated by: $Author: ar2 $     
# Last updated on: $Date: 2004-07-06 15:36:58 $      


#use strict;
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
use Wormbase;
use Carp;
use IO::Handle;
use Getopt::Long;


##############################
# command-line options       #
##############################

my ($est, $mrna, $help, $verbose, $test, $quicktest); # the main options

GetOptions (
            "est"          => \$est,
            "mrna"         => \$mrna,
	    "help"         => \$help,
            "verbose"      => \$verbose,
	    "test"         => \$test,
	    "quicktest"    => \$quicktest
            );

##############################
# Other script variables     #
##############################

my @chrom  = qw(I II III IV V X);
@chrom = qw( III ) if $quicktest;
my (@todo,%unconf);


&printhelp() unless ($est || $mrna);
&printhelp() if ($help);


# what analysis to do?
push(@todo, 'BLAT_TRANSCRIPT_BEST') if ($est); # this file includes OST data
push(@todo, 'BLAT_mRNA_BEST')       if ($mrna); 


my $exon   = 'exon';
my $intron = 'intron';

# for debugging with $verbose
my @trygene = qw() if ($verbose);
my @tryest  = qw() if ($verbose);

###############
# directories #
###############

my $basedir = "/wormsrv2";
$basedir = glob("~wormpub/TEST_BUILD") if ($test or $quicktest);

my $gffdir = "$basedir/autoace/GFF_SPLITS/GFF_SPLITS";
my $outdir = "$basedir/wormbase/misc_dynamic";



####################################################################################################################



#############
# get genes #
#############

foreach my $todo (@todo) {

  %unconf = ();
  
  if  ($todo =~ /BLAT_TRANSCRIPT_BEST/) {
    unlink("$outdir/misc_confirmed_by_EST.ace") if (-e "$outdir/misc_confirmed_by_EST.ace");
    open(ACE,">>$outdir/misc_confirmed_by_EST.ace") || die "Cannot open $outdir/misc_confirmed_by_EST.ace\n";
  }
  elsif ($todo =~ /BLAT_mRNA_BEST/) {
    unlink("$outdir/misc_confirmed_by_mRNA.ace") if (-e "$outdir/misc_confirmed_by_mRNA.ace");
    open(ACE,">>$outdir/misc_confirmed_by_mRNA.ace") || die "Cannot open $outdir/misc_confirmed_by_mRNA.ace\n";
  }
  
  foreach my $chrom (@chrom) {

    #############
    # get genes #
    #############
    
    print &runtime, " : Getting coordinates of 'coding_exon' features for chromosome $chrom\n" if ($verbose);
    my @refgen  = &read_gff('coding_exon',$chrom);
    my %gen     = %{$refgen[0]}; 
    my @genlist = @{$refgen[1]};
    print &runtime, " : Finished getting coordinates\n\n" if ($verbose);    

    ########################
    # get EST and/or cDNAs #
    ########################

    print &runtime, " : Getting $todo feature coordinates for chromosome $chrom\n"  if ($verbose);
    my @refest  = &read_gff($todo,$chrom);
    my %est     = %{$refest[0]};
    my @estlist = @{$refest[1]};
    print &runtime, " : Finished getting coordinates\n\n"  if ($verbose);

    ###############
    # check exons #
    ###############

    print &runtime, " : Checking exons\n" if ($verbose);
    my %confirm = %{&find_match($exon,\%est,\@estlist,\%gen,\@genlist)};

    #check for exon confirmation
    foreach my $gene (@genlist) {
      my $exons = scalar @{$gen{$gene}};
      # check exons
      my $ex = 1;
      for (my $c = 0; $c < $exons; $c++) {
	my $first = $gen{$gene}->[$c][0];
	my $last  = $gen{$gene}->[$c][1];
	for (my $d = $first; $d <= $last; $d++) {
	  if (!exists $confirm{$gene}->{$d}) {
	    $ex = 0;
	    last;
	  }
	}
      }
      if ($ex == 0) {
	$unconf{$gene} = 1;
	foreach my $trygene (@trygene) {
	  print "$gene unconfirmed\n" if (($gene =~ $trygene) && $verbose);
	}
      }
    }
    print &runtime, " : Finished checking\n\n" if ($verbose);

    # save space
    foreach my $old (keys %confirm) {
      delete $confirm{$old};
    }

    #################
    # check introns #
    #################

    my $refig = &make_intron(\%gen);
    my %ingen = %{$refig};
    my $refie = &make_intron(\%est);
    my %inest = %{$refie};

    #save space
    undef %gen;
    undef %est;

    print &runtime, " : Checking introns\n" if ($verbose);
    %confirm = %{&find_match($intron,\%inest,\@estlist,\%ingen,\@genlist)};

    if  (($verbose) && ($todo =~ /BLAT_TRANSCRIPT_BEST/)) {
      print "Checking intron confirmation and putting output for chrom $chrom into $outdir/misc_confirmed_by_EST.ace\n";
    }
    elsif (($verbose) && ($todo =~ /BLAT_mRNA_BEST/)) {
      print "Checking intron confirmation and putting output for chrom $chrom into $outdir/misc_confirmed_by_mRNA.ace\n";
    }

    #check for intron confirmation
    foreach my $gene (@genlist) {
      next if $unconf{$gene};
      unless ( $gen{$gene} ) {
	print "$gene not in gen hash\n" if $verbose;
      }
      my $introns = scalar @{$gen{$gene}}-1;
      # check introns
      my $in = 1;
      for (my $c = 0; $c < $introns; $c++) {
	my $first = $ingen{$gene}->[$c][0];
	my $last  = $ingen{$gene}->[$c][1];
	for (my $d = $first; $d <= $last; $d++) {
	  if (!exists $confirm{$gene}->{$d}) {
	    $in = 0;
	    last;
	  }
	}
      }
      if ($in == 1) {
	print ACE "CDS : \"$gene\"\n";
	print ACE "Confirmed_by EST\n\n"  if ($todo =~ /BLAT_TRANSCRIPT_BEST/);
	print ACE "Confirmed_by cDNA\n\n" if ($todo =~ /BLAT_mRNA_BEST/);
      }
    }
  }
  close(ACE);
}


exit(0);

####################################################################################################################



################
#     subs     #
################

sub read_gff {
  my $filename = shift;
  my $chrom    = shift; 
  my $name     = "";
  my %hash     = ();
  open (GFF, "$gffdir/CHROMOSOME_$chrom.$filename.gff") or die "Cannot open $gffdir/CHROMOSOME_$chrom.$filename.gff $!\n";
  while (<GFF>) {
    next if /\#/;
    my @f = split /\t/;
    my $name;
    ($name) = ($f[8] =~ /\"(\S+)\"/)          if ($filename eq 'coding_exon');
    ($name) = ($f[8] =~ /\"Sequence:(\S+)\"/) if ($filename ne 'coding_exon');
    push @{$hash{$name}} , [$f[3],$f[4]]; 
    # push @{$hash{$name}}, { start => $f[3],end => $f[4] };
    
  }
  
  # need to fix shifts in the alignment
  if ($filename =~ /BLAT/) {
    my %trash;
    foreach my $gene (keys %hash) {
      for (my $g = 0; $g < (@{$hash{$gene}} - 1); $g++) {
	my $exonstart     = $hash{$gene}->[$g][0]; 
	my $exonend       = $hash{$gene}->[$g][1]; 
	my $nextexonstart = $hash{$gene}->[$g+1][0];
	my $nextexonend   = $hash{$gene}->[$g+1][1];
	# extend exons if gap is small enough
	if (abs($nextexonstart - $exonend) < 5) {
	  $hash{$gene}->[$g+1][0] = $exonstart;
	  push @{$trash{$gene}}, $g;
	}  
      }
    }
    
    # delete redundant exons
    foreach my $snod (sort keys %trash) {
      for (my $m = (@{$trash{$snod}})-1; $m >= 0; $m--) {
	splice (@{$hash{$snod}},$trash{$snod}->[$m], 1);
      }
    }
  }
  my @list = sort { $hash{$a}->[0][0] <=> $hash{$b}->[0][0] || $a cmp $b } keys %hash;  

  close(GFF);
  return (\%hash,\@list);
}

################

sub make_intron {
  my $ref   = shift;
  my %struc = %{$ref};
  my %intron;
  foreach my $gene (sort keys %struc) {
    for (my $k = 0; $k < ((scalar @{$struc{$gene}})-1); $k++) {
      my $next        = $k + 1;
      my $chromistart = $struc{$gene}->[$k][1] +1;
      my $chromiend   = $struc{$gene}->[$next][0]-1;
      push @{$intron{$gene}}, [$chromistart,$chromiend];  
    }
  } 
  return (\%intron);
}

################

sub find_match {    
  
  my $status   =   shift;
  my %ESTs     = %{shift;};
  my @ESTlist  = @{shift;};
  my %genes    = %{shift;};
  my @genelist = @{shift;};
  my $lastfail = 0;
  my %store_match = ();
  print &runtime, " : Start searching, modus $status\n" if ($verbose);
  
  
  ##################
  # loop over ESTs #
  ##################

 GENE:for (my $y = 0; $y < @genelist; $y++) {

    my $testgene   = $genelist[$y];
    next if $unconf{$testgene};

    my $geneblocks = (scalar @{$genes{$testgene}})-1; # -1 for array pos last exon
    my $genestart = $genes{$testgene}->[0][0];
    my $geneend   = $genes{$testgene}->[$geneblocks][1];

  EST:for (my $x = $lastfail; $x < @ESTlist; $x++) {
      my $testest   = $ESTlist[$x];
      next unless $ESTs{$testest};
      my $estblocks = (scalar @{$ESTs{$testest}})-1; # -1 for array pos last exon
      my $eststart  = $ESTs{$testest}->[0][0];
      my $estend    = $ESTs{$testest}->[$estblocks][1];

      if ($verbose) {
	foreach my $try (@trygene) { 
	  print "Testing $testgene against $testest\n" if ($testgene =~ $try);
	}
	foreach my $tryest (@tryest) {
	  if ($testest =~ $tryest) {
	    print "\nTesting EST $testest now\n"; 
	    for (my $z = 0; $z <= $estblocks; $z++) {
	      print "Exon $z: $ESTs{$testest}->[$z][0]:$ESTs{$testest}->[$z][1]\n";
	    }
	  }
	}
      }

      ##################
      # loop over ESTs #
      ##################

   # next gene if geneend is left of ESTstart 
   if ($eststart > $geneend) {
     if ($verbose) {
       foreach my $try (@trygene) { 
	 print "With $testest jumped from $testgene to next gene: $genelist[$y+1]\n\n" if $testgene =~ $try; 
       }
     }
     next GENE;
   }
   
   # next EST if genestart is right of ESTend
   elsif ($estend < $genestart) {
     $lastfail = $x;
     if ($verbose) {
       foreach my $try (@trygene) { 
	 print "For $testgene jumped from $testest to next one: $ESTlist[$x+1]\n\n"  if $testgene =~ $try;
       }
     }    
     next EST;
   }
   
   ############################
   # EST matches gene somehow #
   ############################
   
   else {
     
     if ($verbose) {
       foreach my $tryest (@tryest) {
	 print "\nFound match with EST $testest\n" if $testest =~ $tryest; 
       }
     }
     
     # loop over gene exons/introns
   GENEEXON:for (my $w = 0; $w <= $geneblocks; $w++) {
     
     my $CDS_exon_start = $genes{$testgene}->[$w][0];
     my $CDS_exon_end   = $genes{$testgene}->[$w][1];
     
     # loop over EST exons/introns
   ESTEXON:for (my $v = 0; $v <= $estblocks; $v++) {
     
     my $est_exon_start = $ESTs{$testest}->[$v][0];
     my $est_exon_end   = $ESTs{$testest}->[$v][1];
     
     if ($verbose) {
       foreach my $tryest (@tryest) {
	 print "Testing exon $v, $est_exon_start:$est_exon_end\n" if $testest =~ $tryest;
       }
     }    
     
     ####################
     # WHAT WE ALL WANT #
     ####################
     
     if (($CDS_exon_start == $est_exon_start) && ($CDS_exon_end == $est_exon_end)) {
#                            &confirm($testgene,$est_exon_start,$est_exon_end,\%store_match);
       map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$est_exon_end);
       if ($verbose) {
	 foreach my $try (@trygene) { 
	   print "EST $testest\n" if $testgene =~ $try;
	   print "allconfirm:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n" if $testgene =~ $try;
	 }
       }
       next ESTEXON;
     }
     
     #################
     # SPECIAL CASES #
     #################
     
     # okay, something might be redundant, but better twice than never...
     
     # $v = 0 (eststart)
     elsif ($v == 0) {
       
       if (($est_exon_start < $CDS_exon_end) && ($est_exon_end == $CDS_exon_end)) { 
	 if ($est_exon_start >= $CDS_exon_start) {
	   map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$est_exon_end);
	   if ($verbose) {
	     foreach my $try (@trygene) { 
	       print "EST $testest\n" if $testgene =~ $try;
	       print "ESTstart:\tEST $est_exon_start\tgene $CDS_exon_start\tEST $est_exon_end\tgene $CDS_exon_end\tgene $w\tEST $v\n" if $testgene =~ $try;
	     }
	   }    
	   next GENEEXON;
	 }
	 elsif ($est_exon_start < $CDS_exon_start) {
	   map {$store_match{$testgene}->{$_}=1;} ($CDS_exon_start..$est_exon_end);
	   if ($verbose) {
	     foreach my $try (@trygene) { 
	       print "EST $testest\n" if $testgene =~ $try;
	       print "ESTstart:\tEST $est_exon_start\tgene $CDS_exon_start\tEST $est_exon_end\tgene $CDS_exon_end\tgene $w\tEST $v\n" if $testgene =~ $try;
	     }
	   }    
	   next GENEEXON;
	 }
	 else {
	   if ($verbose) {
	     foreach my $tryest (@tryest) {
	       if ($testest =~ $tryest) {
		 print "Missed and overlap for eststart:\tgene $testgene\tEST $testest\n";
		 print "Coordinates:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n\n";
	       }
	     }
	   }             
	 }
	 
       }
       
       # $v = 0 AND $v = $estblocks (single exon EST)
       elsif ($v == $estblocks) {
	 
	 # single exon EST within geneexon
	 if (($est_exon_start >= $CDS_exon_start) && ($est_exon_end <= $CDS_exon_end)) {
	   map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$est_exon_end);
	   if ($verbose) {
	     foreach my $try (@trygene) { 
	       print "EST $testest\n" if $testgene =~ $try;
	       print "singleEST:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n" if $testgene =~ $try;
	     }
	   }    
	   next EST;
	 }
	 
	 # single exon EST overlapping gene start
	 elsif (($est_exon_start <= $CDS_exon_start) && ($est_exon_end <= $CDS_exon_end) && ($w == 0)) {
	   map {$store_match{$testgene}->{$_}=1;} ($CDS_exon_start..$est_exon_end);
	   if ($verbose) {
	     foreach my $try (@trygene) { 
	       print "EST $testest\n" if $testgene =~ $try;
	       print "singleEST:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n" if $testgene =~ $try;
	     }
	   }    
	   next EST;
	 }
	 
	 # single exon EST overlapping gene end
	 elsif (($est_exon_start >= $CDS_exon_start) && ($est_exon_end >= $CDS_exon_end) && ($w == $geneblocks)) {
	   map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$CDS_exon_end);
	   if ($verbose) {
	     foreach my $try (@trygene) { 
	       print "EST $testest\n" if $testgene =~ $try;
	       print "singleEST:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n" if $testgene =~ $try;
	     }
	   }    
	   next EST;
	 }
	 
	 # single exon EST covering whole single exon gene
	 elsif (($est_exon_start < $CDS_exon_start) && ($est_exon_end > $CDS_exon_end) && ($geneblocks == 0)) {
	   map {$store_match{$testgene}->{$_}=1;} ($CDS_exon_start..$CDS_exon_end);
	   if ($verbose) {
	     foreach my $try (@trygene) { 
	       print "EST $testest\n" if $testgene =~ $try;
	       print "singleEST:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n" if $testgene =~ $try;
	     }
	   }    
	   next GENE;
	 }
	 else {
	   if ($verbose) {
	     foreach my $tryest (@tryest) {
	       if ($testest =~ $tryest) {
		 print "Missed and overlap for single exon EST:\tgene $testgene\tEST $testest\n";
		 print "Coordinates:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n\n";
	       }
	     } 
	   }           
	 }
       }
     }
     # $v = $estblocks (estend)
     elsif ($v == $estblocks) {
       if (($est_exon_start == $CDS_exon_start) && ($est_exon_end > $CDS_exon_start)) {
	 if ($est_exon_end <= $CDS_exon_end) {
	   map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$est_exon_end);
	   if ($verbose) {
	     foreach my $try (@trygene) { 
	       print "EST $testest\n" if $testgene =~ $try;
	       print "ESTend:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n" if $testgene =~ $try;
	     }
	   }    
	   next EST;
	 }
	 elsif ($est_exon_end > $CDS_exon_end) {
	   map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$CDS_exon_end);
	   if ($verbose) {
	     foreach my $try (@trygene) { 
	       print "EST $testest\n" if $testgene =~ $try;
	       print "ESTend:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n" if $testgene =~ $try;
	     }
	   }    
	   next EST;
	 }  
	 else {
	   if ($verbose) {
	     foreach my $tryest (@tryest) {
	       if ($testest =~ $tryest) {
		 print "Missed and overlap for estend:\tgene $testgene\tEST $testest\n";
		 print "Coordinates:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n\n";
	       }
	     } 
	   }            
	 }
       }
     }
     # $w = 0 (genestart)
     elsif ($w == 0) {
       if (($CDS_exon_start < $est_exon_end) && ($CDS_exon_end == $est_exon_end)) {
	 if ($CDS_exon_start >= $est_exon_start) {
	   map {$store_match{$testgene}->{$_}=1;} ($CDS_exon_start..$est_exon_end);
	   if ($verbose) {
	     foreach my $try (@trygene) { 
	       print "EST $testest\n" if $testgene =~ $try;
	       print "Genestart:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n" if $testgene =~ $try;
	     }
	   }    
	   next ESTEXON;
	 }
	 elsif ($CDS_exon_start < $est_exon_start) {
	   map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$est_exon_end);
	   if ($verbose) {
	     foreach my $try (@trygene) { 
	       print "EST $testest\n" if $testgene =~ $try;
	       print "Genestart:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n" if $testgene =~ $try;
	     }
	   }    
	   next ESTEXON;
	 }   
	 else {
	   if ($verbose) {
	     foreach my $tryest (@tryest) {
	       if ($testest =~ $tryest) {
		 print "Missed and overlap for genestart:\tgene $testgene\tEST $testest\n";
		 print "Coordinates:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n\n";
	       }
	     }         
	   }
	 }
	 
       }
       
       # $w = 0 AND $w = $geneblocks (single exon gene)
       elsif ($w == $geneblocks) {
	 
	 # single exon gene completely covered by EST
	 if (($est_exon_start <= $CDS_exon_start ) && ($est_exon_end >= $CDS_exon_end)) {
	   map {$store_match{$testgene}->{$_}=1;} ($CDS_exon_start..$CDS_exon_end);
	   if ($verbose) {
	     foreach my $try (@trygene) { 
	       print "EST $testest\n" if $testgene =~ $try;
	       print "singleGene:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n" if $testgene =~ $try;
	     }
	   }    
	   next GENE; 
	 }   
	 
	 # single exon gene end covered
	 elsif (($est_exon_start >= $CDS_exon_start) && ($est_exon_end >= $CDS_exon_end) && ($v == 0)) {
	   map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$CDS_exon_end); 
	   if ($verbose) {
	     foreach my $try (@trygene) { 
	       print "EST $testest\n" if $testgene =~ $try;
	       print "singleGene:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n" if $testgene =~ $try;
	     }
	   }    
	   next EST; 
	 }
	 
	 # single exon gene start covered
	 elsif (($est_exon_start <= $CDS_exon_start) && ($est_exon_end <= $CDS_exon_end) && ($v == $estblocks)) {
	   map {$store_match{$testgene}->{$_}=1;} ($CDS_exon_start..$est_exon_end);  
	   if ($verbose) {
	     foreach my $try (@trygene) { 
	       print "EST $testest\n" if $testgene =~ $try;
	       print "singleGene:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n" if $testgene =~ $try;
	     }
	   }    
	   next EST; 
	 }
	 
	 # single exon EST inside single exon gene
	 elsif (($est_exon_start >= $CDS_exon_start) && ($est_exon_end <= $CDS_exon_end) && ($estblocks == 0)) {
	   map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$est_exon_end); 
	   if ($verbose) {
	     foreach my $try (@trygene) { 
	       print "EST $testest\n" if $testgene =~ $try;
	       print "singleGene:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n" if $testgene =~ $try;
	     }
	   }    
	   next EST; 
	 }
	 else {
	   if ($verbose) {
	     foreach my $tryest (@tryest) {
	       if ($testest =~ $tryest) {
		 print "Missed and overlap for single exon gene:\tgene $testgene\tEST $testest\n";
		 print "Coordinates:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n\n";
	       }
	     } 
	   }         
	 }
       }
     }
     # $w = $geneblocks (geneend) 
     elsif ($w == $geneblocks) {
       if (($CDS_exon_start == $est_exon_start) && ($CDS_exon_end > $est_exon_start)) {
	 if ($CDS_exon_end <= $est_exon_end) {
	   map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$CDS_exon_end);
	   if ($verbose) {
	     foreach my $try (@trygene) {    
	       print "EST $testest\n" if $testgene =~ $try;
	       print "Geneend:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n" if $testgene =~ $try;
	     }
	   }    
	   next EST;
	 }
	 elsif ($CDS_exon_end > $est_exon_end) {
	   map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$est_exon_end);
	   if ($verbose) {
	     foreach my $try (@trygene) {  
	       print "EST $testest\n" if $testgene =~ $try;
	       print "Geneend:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n" if $testgene =~ $try;
	     }
	   }    
	   next EST;
	 }
	 else {
	   if ($verbose) {
	     foreach my $tryest (@tryest) {
	       if ($testest =~ $tryest) {
		 print "Missed and overlap for geneend:\tgene $testgene\tEST $testest\n";
		 print "Coordinates:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n\n";
	       }
	     } 
	   }            
	 }
       }
     }  
     else {
       if ($verbose) {
	 foreach my $tryest (@tryest) {
	   if ($testest =~ $tryest) {
	     print "Missed an overlap:\tgene $testgene\tEST $testest\n";
	     print "Coordinates:\tgene $CDS_exon_start\tEST $est_exon_start\tgene $CDS_exon_end\tEST $est_exon_end\tEST $v\tgene $w\n\n";
	   }
	 }
       }
     }
   }  
   }
   }
 }
 }
  return \%store_match;
}    

################

#sub confirm  {
#
#    my $gene   = shift;
#    my $start  = shift;
#    my $end    = shift;
#    my $recref = shift;
#    
#    for (my $n = $start; $n <= $end; $n++) {
#	$recref->{$gene}->{$n} = 1;
#    }
#}  

################

sub printhelp {
    exec ('perldoc',$0);
}    

__END__

=pod

=head2   NAME - confirm_genes.pl

=head1 USAGE

=over 4 

=item  confirm_genes.pl -options

=back

confirm_genes.pl parses GFF files and looks for gene models that are completely supported by 
EST or cDNA matches, depending on your choice. These gene models get an appropriate tag.

mandatory arguments:

=over 4

=item -est for checking for EST matches

=item or

=item -mrna for checking for cDNA matches 

=back

optional arguments

=over 4

=item -verbose for being verbose

=item -help for help

=back

author: Kerstin Jekosch (kj2@sanger.ac.uk)

=cut












