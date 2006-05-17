#!/usr/local/bin/perl5.8.0 -w
#
# confirm_genes.pl
#
# by Dan Lawson
#
# Makes CDS status information by looking at transcript to exon mappings
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2006-05-17 08:58:19 $      


use strict;
use lib  $ENV{'CVS_DIR'};
use Wormbase;
use Carp;
use IO::Handle;
use Getopt::Long;
use Storable;

##############################
# command-line options       #
##############################

my ($help, $verbose, $test, $quicktest, $store, $debug); # the main options
my $load; # whether to load to autoace or not

GetOptions (
	    "help"         => \$help,
            "verbose"      => \$verbose,
	    "test"         => \$test,
	    "load"         => \$load,
	    "quicktest"    => \$quicktest,
	    "store:s"      => \$store,
	    "debug:s"      => \$debug
            );

##############################
# Other script variables     #
##############################
$test = $quicktest if $quicktest;

my $wormbase;
if ($store) {
  $wormbase = Storable::retrieve($store) or croak("cant restore wormbase from $store\n")
} else {
  $wormbase = Wormbase->new( -debug => $debug, -test => $test, );
}

my $log = Log_files->make_build_log($wormbase);
my $maintainers = "All";

my @chromosomes  = qw( I II III IV V X MtDNA);                    # all chromosomes
@chromosomes     = qw( III ) if $quicktest;

&printhelp() if ($help);

my %CDSstatus;			# status data for each CDS
my %unconf = ();

###############
# directories #
###############


my $gffdir = $wormbase->gff_splits;
my $outdir = $wormbase->acefiles;
my $output = "$outdir/gene_confirmation_status.ace"; # output acefile name & path


####################################################################################################################

#############
# get genes #
#############

# open output acefile 

open (ACE,">$output") || die "Cannot open $output\n";

foreach my $chromosome (@chromosomes) {

  $log->write_to("Processing chromosome $chromosome\n");
  
  ###################################################
  # get CDS names (@CDSlist) and coordinates (%CDS) #
  ###################################################

  print $wormbase->runtime, " : Getting coordinates of 'coding_exon' features for chromosome $chromosome\n" if ($verbose);
  my @CDS_GFF  = &read_gff('curated',$chromosome);
  my %CDS      = %{$CDS_GFF[0]}; 
  my @CDSlist  = @{$CDS_GFF[1]};
  print $wormbase->runtime, " : Finished getting coordinates\n" if ($verbose);    
    
  ########################
  # get EST and/or cDNAs #
  ########################
    
  print $wormbase->runtime, " : Getting BLAT_TRANSCRIPT_BEST feature coordinates for chromosome $chromosome\n"  if ($verbose);
  #create BLAT_TRANSCRIPT_BEST for ease of transition!
  $wormbase->run_command("cd $gffdir; cat CHROMOSOME_${chromosome}_BLAT_EST_BEST.gff CHROMOSOME_${chromosome}_BLAT_OST_BEST.gff CHROMOSOME_${chromosome}_BLAT_mRNA_BEST.gff >  CHROMOSOME_${chromosome}_BLAT_TRANSCRIPT_BEST.gff", $log);
  my @EST_GFF  = &read_gff('BLAT_TRANSCRIPT_BEST',$chromosome);
  my %EST     = %{$EST_GFF[0]};
  my @ESTlist = @{$EST_GFF[1]};
  print $wormbase->runtime, " : Finished getting coordinates\n\n"  if ($verbose);
  
  ###############
  # check exons #
  ###############
  
  print $wormbase->runtime, " : Checking exons\n" if ($verbose);
  my %confirm = %{&find_match(\%EST,\@ESTlist,\%CDS,\@CDSlist)}; # Time costly step...... should replace it with a GFF_overlap call
  print $wormbase->runtime, " : Calculated exon-EST hash\n" if ($verbose);
  
  my $EXONpass;			# vars for tracking the status of the exons
  my $EXONhalf;
  my $EXONfail;
  my $BASEfail;
  my $BASEpass;

  #check for exon confirmation
  foreach my $look (@CDSlist) {

    # calculate the number of exons in this CDS
    my $no_exons = scalar @{$CDS{$look}};

    # reset tracking vars for each exon
    $EXONpass = 0;
    $EXONfail = 0;
    $EXONhalf = 0;

    # loop over each exon....
    for (my $exon = 0; $exon < $no_exons; $exon++) { # loop through each exon

      my $first = $CDS{$look}->[$exon][0]; # Exon start coordinate
      my $last  = $CDS{$look}->[$exon][1]; # Exon stop coordinate

      print "$look [$first-$last] :" if ($verbose);

      # reset tracking vars for each base
      $BASEfail = 0;
      $BASEpass = 0;

      for (my $base = $first; $base <= $last; $base++) { # Loop through each base of the exon
		
	if ($confirm{$look}->{$base}) {
	  print "$confirm{$look}->{$base}" if ($verbose); 
	  $BASEpass = 1;
	} else {
	  print "0" if ($verbose);
	  $BASEfail = 1;
	}
      }				#_ next base to be checked
	    
      print "\n" if ($verbose);
	    
      # logic for exon status
      if ( ($BASEfail == 1) && ($BASEpass == 0) ) { # Fail if no bases are supported
	$EXONfail++;
      } elsif ( ($BASEfail == 1) && ($BASEpass == 1) ) { # Support if more than 1 base is supported and more than 1 base is unsupported
	$EXONhalf++;
      } else {			# Confirm if all bases are supported
	$EXONpass++;
      }

    }				#_ next exon in list

    # Assign status
	
    if ( ($EXONpass == $no_exons) && ($EXONfail == 0) ) {
      $CDSstatus{$look} = "Confirmed";
      $CDSstatus{$look} = "Confirmed -C \"exon status($no_exons:$EXONpass.$EXONhalf.$EXONfail)\"" if ($verbose);
    } elsif ( ($EXONpass > 0) || ($EXONhalf > 0) ) {
      $CDSstatus{$look} = "Partially_confirmed";
      $CDSstatus{$look} = "Partially_confirmed -C \"exon status($no_exons:$EXONpass.$EXONhalf.$EXONfail)\"" if ($verbose);
    } elsif ( ($EXONfail == $no_exons) && (($EXONpass == 0) && ($EXONhalf == 0))  ) {
      $CDSstatus{$look} = "Predicted";
      $CDSstatus{$look} = "Predicted -C \"exon status($no_exons:$EXONpass.$EXONhalf.$EXONfail)\"" if ($verbose);
    }
  }				#_ next CDS in @CDSlist

  print $wormbase->runtime, " : Finished checking\n\n" if ($verbose);

  # Print output 
  print $wormbase->runtime, " : Writing data to file for Chromosome_$chromosome\n" if ($verbose);
  foreach my $look (@CDSlist) {
    print ACE "CDS : \"$look\"\n";
    print ACE "Prediction_status $CDSstatus{$look}\n\n";
  }
  print $wormbase->runtime, " : Finished with Chromosome_$chromosome\n\n" if ($verbose);
}

close(ACE);

# load file to autoace?
$wormbase->load_to_database($wormbase->autoace, $output,"confirm_genes") if($load);

$log->mail();
exit(0);

####################################################################################################################

###############
# subroutines #
###############

sub read_gff {
    
  my $filename   = shift;
  my $chromosome = shift; 
  my $name     = "";
  my %hash     = ();
   
  # open the GFF file
  open (GFF, "$gffdir/CHROMOSOME_${chromosome}_$filename.gff") or die "Cannot open $gffdir/CHROMOSOME_${chromosome}_$filename.gff $!\n";
  while (<GFF>) {
	
    # discard header lines
    next if /\#/;

    my @f = split /\t/;
    undef $name;
    ($name) = ($f[8] =~ /\"(\S+)\"/)          if ($filename eq 'curated' and /coding_exon/); # CDS name
    ($name) = ($f[8] =~ /\"Sequence:(\S+)\"/) if ($filename eq 'BLAT_TRANSCRIPT_BEST'); # transcript name
    next unless $name;
    push @{$hash{$name}} , [$f[3],$f[4]]; 
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

sub find_match {    
  
  my %ESTs        = %{shift;};
  my @ESTlist     = @{shift;};
  my %genes       = %{shift;};
  my @genelist    = @{shift;};
  my $lastfail    = 0;
  my %store_match = ();

  print $wormbase->runtime, " : Start searching, modus exon\n" if ($verbose);
    
    
  ##################
  # loop over ESTs #
  ##################

 GENE: for (my $y = 0; $y < @genelist; $y++) {
      
    my $testgene   = $genelist[$y];
    next if $unconf{$testgene};
      
    my $geneblocks = (scalar @{$genes{$testgene}})-1; # -1 for array pos last exon
    my $genestart  = $genes{$testgene}->[0][0];
    my $geneend    = $genes{$testgene}->[$geneblocks][1];
      
  EST: for (my $x = $lastfail; $x < @ESTlist; $x++) {
      my $testest   = $ESTlist[$x];
      next unless $ESTs{$testest};
      my $estblocks = (scalar @{$ESTs{$testest}})-1; # -1 for array pos last exon
      my $eststart  = $ESTs{$testest}->[0][0];
      my $estend    = $ESTs{$testest}->[$estblocks][1];
	
      ##################
      # loop over ESTs #
      ##################

	
      #######################################
      # No overlap between HSP and CDS exon #
      #######################################
	
      # next gene if geneend is left of ESTstart
	
      # no match,            EXON  ----------------<========================--------========*=----------------
      #                      HSP                                                                  ##########
	
      if ($eststart > $geneend) {
	next GENE;
      }

	
      # next EST if genestart is right of ESTend
	
      # no match,            EXON  ----------------<========================--------========*=----------------
      #                      HSP   ###########

      elsif ($estend < $genestart) {
	$lastfail = $x;
	next EST;
      }
   


      # HSP overlaps with CDS exon somehow

      else {
	    
	# loop over gene exons/introns
      GENEEXON:  for (my $w = 0; $w <= $geneblocks; $w++) {
	      
	  my $CDS_exon_start = $genes{$testgene}->[$w][0];
	  my $CDS_exon_end   = $genes{$testgene}->[$w][1];
	      
	  # loop over EST exons/introns
	ESTEXON:  for (my $v = 0; $v <= $estblocks; $v++) {
		
	    my $est_exon_start = $ESTs{$testest}->[$v][0];
	    my $est_exon_end   = $ESTs{$testest}->[$v][1];
		
	    ################################################
	    # IDEAL OUTCOME, BLAT HSP matches exon exactly #
	    ################################################
		
	    #                                  EXON_start ¬             EXON_stop ¬
	    # complete match,      EXON  ------->===^^^^^^^========================^^^^^^^^^========*-------------
	    #                      HSP                  ---########################---
		
	    if ( ($est_exon_start == $CDS_exon_start) && ($est_exon_end == $CDS_exon_end) ) {
	      map {$store_match{$testgene}->{$_} = 1;} ($est_exon_start..$est_exon_end);   
	      next ESTEXON;                                                                
	    }


	    if ( ($est_exon_start < $CDS_exon_end) && ($est_exon_end == $CDS_exon_end)) { 
		    
	      #                                         HSP_start¬         EXON_stop¬
	      # partial match        EXON  ------->===^^^^^^^========================^^^^^^^^^=======*-------------
	      #                      HSP                          ###################

	      if ($est_exon_start >= $CDS_exon_start) {                                  
		map {$store_match{$testgene}->{$_} = 1;} ($est_exon_start..$est_exon_end);
		next ESTEXON;                                                             
	      }

	      #                                   EXON_start¬              EXON_stop¬
	      # complete match       EXON    ----->===^^^^^^^========================^^^^^^^^^========*-------------
	      #                      HSP                  ###########################

	      elsif ($est_exon_start < $CDS_exon_start) {                                      
		map {$store_match{$testgene}->{$_} = 1;} ($CDS_exon_start..$est_exon_end);   
		next ESTEXON;                                                                
	      }
	    }
		


	    if (($CDS_exon_start == $est_exon_start) && ($CDS_exon_end > $est_exon_start)) {   
		    
	      #                                                HSP_start¬               CDS_stop¬
	      # complete match    EXON  ----------->==============^^^^^^^========================^^^^^^^^^==============*----------
	      #                   HSP                                    ##########################

	      if ($CDS_exon_end <= $est_exon_end) {                                       
		map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$CDS_exon_end);
		next ESTEXON;                                                           
	      }

	      #                                                HSP_start¬           HSP_stop¬
	      # partial match     EXON  ----------->==============^^^^^^^========================^^^^^^^^^==============*----------
	      #                   HSP                                    ####################                 
		    
	      elsif ($CDS_exon_end > $est_exon_end) {                                                            
		map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$est_exon_end);                       
		next ESTEXON;                                                                                  
	      }
	    }
		
	    #################
	    # SPECIAL CASES #
	    #################
		
	    # First HSP of transcript sequence ($v = 0)
	    if ($v == 0) {
	      if (($est_exon_start < $CDS_exon_end) && ($est_exon_end == $CDS_exon_end)) { 

		#                                         HSP_start¬         EXON_stop¬
		# partial match        EXON  ------->===^^^^^^^========================^^^^^^^^^=======*-------------
		#                      HSP                          ###################

		if ($est_exon_start >= $CDS_exon_start) {                                        
		  map {$store_match{$testgene}->{$_} = 1;} ($est_exon_start..$est_exon_end);   
		  next GENEEXON;                                                               
		}

		#                                   EXON_start¬              EXON_stop¬
		# complete match       EXON    ----->===^^^^^^^========================^^^^^^^^^========*-------------
		#                      HSP                  ###########################

		elsif ($est_exon_start < $CDS_exon_start) {                                      
		  map {$store_match{$testgene}->{$_} = 1;} ($CDS_exon_start..$est_exon_end);   
		  next GENEEXON;                                                               
		}
	      }
	    }			#_ end of initial HSP processing
		    
	    # Single HSP for transcript  ($v = 0 AND $v = $estblocks)  i.e. single HSP EST
	    if ($v == $estblocks) {
			
	      # single exon EST within geneexon

	      #                                      HSP_start ¬        HSP_stop ¬
	      # complete match          EXON  ----->===^^^^^^^========================^^^^^^^^^========*-------------
	      #                         HSP                     ##################
		    
	      if (($est_exon_start >= $CDS_exon_start) && ($est_exon_end <= $CDS_exon_end)) {      
		map {$store_match{$testgene}->{$_} = 1;} ($est_exon_start..$est_exon_end); 
		next EST;                                                                  
	      }
		    
	      # single exon EST overlapping gene start

	      #                                   EXON_start ¬    HSP_stop ¬
	      # partial match           EXON  ---------------->=======================^^^^^^^^^========*-------------
	      #                         HSP             ####################

	      elsif (($est_exon_start <= $CDS_exon_start) && ($est_exon_end <= $CDS_exon_end) && ($w == 0)) {      
		map {$store_match{$testgene}->{$_} = 1;} ($CDS_exon_start..$est_exon_end);                       
		next EST;                                                                                        
	      }
		    
	      # single exon EST overlapping gene end

	      #                                           HSP_start ¬      EXON_stop ¬
	      # partial match           EXON  ---->====^^^^^^^=======================*-------------------------------
	      #                         HSP                          ####################

	      elsif (($est_exon_start >= $CDS_exon_start) && ($est_exon_end >= $CDS_exon_end) && ($w == $geneblocks)) { 
		map {$store_match{$testgene}->{$_} = 1;} ($est_exon_start..$CDS_exon_end);                            
		next EST;                                                                                             
	      }
			
	      # single exon EST covering whole single exon gene

	      #                                             EXON_start¬                  EXON_stop ¬
	      # complete match          EXON  ------------------------->===========================*-------------------------------
	      #                         HSP                        ####################################

	      elsif (($est_exon_start < $CDS_exon_start) && ($est_exon_end > $CDS_exon_end) && ($geneblocks == 0)) {    
		map {$store_match{$testgene}->{$_} = 1;} ($CDS_exon_start..$CDS_exon_end);                            
		next GENE;                                                                                            
	      }
	    }			#_ end of Single HSP processing

	    # Final HSP of transcribed sequence ($v = $estblocks)
	    if ($v == $estblocks) {
	      if ( ($est_exon_start == $CDS_exon_start) && ($est_exon_end > $CDS_exon_start) ) { 

		#                                                                                EXON_start¬  HSP_stop ¬
		# partial match     EXON  ----------->==============^^^^^^^========================^^^^^^^^^==============*----------
		#                   HSP                                                                     ############

		if ($est_exon_end <= $CDS_exon_end) {                                                       
		  map {$store_match{$testgene}->{$_} = 1;} ($est_exon_start..$est_exon_end);              
		  next EST;                                                                               
		}

		#                                                                                EXON_start¬     EXON_stop¬
		# complete match    EXON  ----------->==============^^^^^^^=======================^^^^^^^^^^==============*----------
		#                   HSP                                                                     #################
			
		elsif ($est_exon_end > $CDS_exon_end) {                                                     
		  map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$CDS_exon_end);                
		  next EST;                                                                               
		}  
	      }
	    }			#_ end of final HSP processing
		
	    # First CDS exon ($w = 0)
	    if ($w == 0) {
	      if ( ($CDS_exon_start < $est_exon_end) && ($CDS_exon_end == $est_exon_end) ) {
			
		#                         EXON_start¬     EXON_stop¬
		# partial match     EXON  ----------->==============^^^^^^^========================^^^^^^^^^==============*---------- 
		#                   HSP           ##################
			
		if ($CDS_exon_start >= $est_exon_start) {                                                   
		  map {$store_match{$testgene}->{$_}=1;} ($CDS_exon_start..$est_exon_end);                
		  next ESTEXON;                                                                           
		}

		#                              HSP_start¬ EXON_stop¬
		# partial match     EXON  ----------->==============^^^^^^^========================^^^^^^^^^==============*----------
		#                   HSP                  ###########

		elsif ($CDS_exon_start < $est_exon_start) {                                                 
		  map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$est_exon_end);                
		  next ESTEXON;                                                                           
		}   
	      }
	    }			#_ end of intial CDS exon processing
		
	    # Single exon CDS  ($w = 0 AND $w = $geneblocks)
	    if ($w == $geneblocks) {
			
	      # single exon gene completely covered by EST

	      #                                             EXON_start¬                   EXON_stop¬
	      # complete match    EXON  ------------------------->===========================*-------------------------------
	      #                   HSP                         ####################################

	      if (($est_exon_start <= $CDS_exon_start ) && ($est_exon_end >= $CDS_exon_end)) {                     
		map {$store_match{$testgene}->{$_}=1;} ($CDS_exon_start..$CDS_exon_end);                         
		next GENE;                                                                                       
	      }   
		    
	      # single exon gene end covered

	      #                                                  HSP_start¬               EXON_stop¬
	      # partial match     EXON  ------------------------->===========================*-------------------------------
	      #                   HSP                                ####################################

	      elsif (($est_exon_start >= $CDS_exon_start) && ($est_exon_end >= $CDS_exon_end) && ($v == 0)) {      
		map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$CDS_exon_end);                         
		next EST;                                                                                        
	      }
		    
	      # single exon gene start covered
		    
	      #                                             EXON_start¬                 HSP_stop¬
	      # partial match     EXON  ------------------------->===========================*-------------------------------
	      #                   HSP                  ####################################
		    
	      elsif (($est_exon_start <= $CDS_exon_start) && ($est_exon_end <= $CDS_exon_end) && ($v == $estblocks)) { 
		map {$store_match{$testgene}->{$_}=1;} ($CDS_exon_start..$est_exon_end);                             
		next EST;                                                                                            
	      }
		    
	      # single exon EST inside single exon gene

	      #                                                 HSP_start¬            HSP_stop¬
	      # complete match    EXON  ------------------------->===========================*-------------------------------
	      #                   HSP                               #####################
		    
	      elsif (($est_exon_start >= $CDS_exon_start) && ($est_exon_end <= $CDS_exon_end) && ($estblocks == 0)) {  
		map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$est_exon_end);                             
		next EST;                                                                                            
	      }
	    }
		 
	    if ($w == $geneblocks) {
	      if (($CDS_exon_start == $est_exon_start) && ($CDS_exon_end > $est_exon_start)) {   
		    
		#                                                HSP_start¬               CDS_stop¬
		# complete match    EXON  ----------->==============^^^^^^^========================^^^^^^^^^==============*----------
		#                   HSP                                    ##########################
			
		if ($CDS_exon_end <= $est_exon_end) {                                                                
		  map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$CDS_exon_end);                         
		  next EST;                                                                                        
		}
			
		#                                                HSP_start¬           HSP_stop¬
		# partial match     EXON  ----------->==============^^^^^^^========================^^^^^^^^^==============*----------
		#                   HSP                                    ####################                 

		elsif ($CDS_exon_end > $est_exon_end) {                                                              
		  map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$est_exon_end);     
		  next EST;                                                                    
		}
	      }
	    }
		
		
	  }			#_ ESTEXON
	}			#_ GENEEXON
      }
    }				#_ EST
  }				#_ GENE
    
  return \%store_match;
    
}


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
transcript matches. The output file contains a Prediction_status for each CDS along the 
following lines:

    Predicted : No transcript data overlaps the coding exons
    
    Partially_confirmed : Some overlap between transcript data and coding exons
 
    Confirmed : All coding exons are supported by transcript data


This is in line with the models for WS133 onward. Each CDS should have a status in the
resulting acefile.

MANDATORY arguments:

=over 4

=item none

=back

OPTIONAL arguments

=over 4

=item -verbose, Extensive output to STDOUT and makes an acefile which 

=item -         includes the exon status comments

=item -help,  this help

=back

author: Dan Lawson (dl1@sanger.ac.uk) after Kerstin Jekosch

=cut












