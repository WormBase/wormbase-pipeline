#!/usr/local/bin/perl5.8.0 -w
#
# confirm_genes.pl
#
# by Dan Lawson
#
# Makes CDS status information by looking at transcript to exon mappings
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2004-10-20 09:46:00 $      


use strict;
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
use Wormbase;
use Carp;
use IO::Handle;
use Getopt::Long;

##############################
# command-line options       #
##############################

my ($help, $verbose, $test, $quicktest); # the main options
my $load; # whether to load to autoace or not

GetOptions (
	    "help"         => \$help,
            "verbose"      => \$verbose,
	    "test"         => \$test,
	    "load"         => \$load,
	    "quicktest"    => \$quicktest
            );

##############################
# Other script variables     #
##############################

my $log = Log_files->make_build_log();
my $maintainers = "All";

my @chromosomes  = qw( I II III IV V X );                    # all chromosomes
@chromosomes     = qw( III ) if $quicktest;

&printhelp() if ($help);

my %CDSstatus;                # status data for each CDS
my %unconf = ();

###############
# directories #
###############

my $basedir = "/wormsrv2";
$basedir    = glob("~wormpub/TEST_BUILD") if ($test or $quicktest);

my $gffdir = "$basedir/autoace/GFF_SPLITS/GFF_SPLITS";
my $outdir = "$basedir/autoace/acefiles";
my $output = "$outdir/gene_confirmation_status.ace";                             # output acefile name & path


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

    print &runtime, " : Getting coordinates of 'coding_exon' features for chromosome $chromosome\n" if ($verbose);
    my @CDS_GFF  = &read_gff('coding_exon',$chromosome);
    my %CDS      = %{$CDS_GFF[0]}; 
    my @CDSlist  = @{$CDS_GFF[1]};
    print &runtime, " : Finished getting coordinates\n" if ($verbose);    
    
    ########################
    # get EST and/or cDNAs #
    ########################
    
    print &runtime, " : Getting BLAT_TRANSCRIPT_BEST feature coordinates for chromosome $chromosome\n"  if ($verbose);
    my @EST_GFF  = &read_gff('BLAT_TRANSCRIPT_BEST',$chromosome);
    my %EST     = %{$EST_GFF[0]};
    my @ESTlist = @{$EST_GFF[1]};
    print &runtime, " : Finished getting coordinates\n\n"  if ($verbose);

    ###############
    # check exons #
    ###############

    print &runtime, " : Checking exons\n" if ($verbose);
    my %confirm = %{&find_match(\%EST,\@ESTlist,\%CDS,\@CDSlist)};          # Time costly step...... should replace it with a GFF_overlap call
    print &runtime, " : Calculated exon-EST hash\n" if ($verbose);
    
    my $EXONpass;     # vars for tracking the status of the exons
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
	for (my $exon = 0; $exon < $no_exons; $exon++) {             # loop through each exon

	    my $first = $CDS{$look}->[$exon][0];                     # Exon start coordinate
	    my $last  = $CDS{$look}->[$exon][1];                     # Exon stop coordinate

	    print "$look [$first-$last] :" if ($verbose);

	    # reset tracking vars for each base
	    $BASEfail = 0;
	    $BASEpass = 0;

	    for (my $base = $first; $base <= $last; $base++) {       # Loop through each base of the exon
		
		if ($confirm{$look}->{$base}) {
		    print "$confirm{$look}->{$base}" if ($verbose); 
		    $BASEpass = 1;
		}
		else {
		    print "0" if ($verbose);
		    $BASEfail = 1;
		}
	    } #_ next base to be checked
	    
	    print "\n" if ($verbose);
	    
	    # logic for exon status
	    if ( ($BASEfail == 1) && ($BASEpass == 0) ) {           # Fail if no bases are supported
		$EXONfail++;
	    }
	    elsif  ( ($BASEfail == 1) && ($BASEpass == 1) ) {       # Support if more than 1 base is supported and more than 1 base is unsupported
		$EXONhalf++;
	    }
	    else {                                                  # Confirm if all bases are supported
		$EXONpass++;
	    }

	} #_ next exon in list

	# Assign status
	
	if ( ($EXONpass == $no_exons) && ($EXONfail == 0) ) {
	    $CDSstatus{$look} = "Confirmed";
	    $CDSstatus{$look} = "Confirmed -C \"exon status($no_exons:$EXONpass.$EXONhalf.$EXONfail)\"" if ($verbose);
	}
	elsif ( ($EXONpass > 0) || ($EXONhalf > 0) ) {
	    $CDSstatus{$look} = "Partially_confirmed";
	    $CDSstatus{$look} = "Partially_confirmed -C \"exon status($no_exons:$EXONpass.$EXONhalf.$EXONfail)\"" if ($verbose);
	}
	elsif ( ($EXONfail == $no_exons) && (($EXONpass == 0) && ($EXONhalf == 0))  ) {
	    $CDSstatus{$look} = "Predicted";
	    $CDSstatus{$look} = "Predicted -C \"exon status($no_exons:$EXONpass.$EXONhalf.$EXONfail)\"" if ($verbose);
	}
    } #_ next CDS in @CDSlist
    
    print &runtime, " : Finished checking\n\n" if ($verbose);
    
    # Print output 
    print &runtime, " : Writing data to file for Chromosome_$chromosome\n" if ($verbose);
    foreach my $look (@CDSlist) {
	print ACE "CDS : \"$look\"\n";
	print ACE "Prediction_status $CDSstatus{$look}\n\n";
    }
    print &runtime, " : Finished with Chromosome_$chromosome\n\n" if ($verbose);
    
}

close(ACE);

# load file to autoace?
if($load){
  $log->write_to("Loading file to autoace\n");
  my $command = "autoace_minder.pl -load /wormsrv2/autoace/acefiles/gene_confirmation_status.ace -tsuser interpro_motifs";

  my $status = system($command);
  if(($status >>8) != 0){
    $log->write_to("ERROR: Loading interpro_motifs.ace file failed \$\? = $status\n");
  }
}


$log->mail("$maintainers", "BUILD REPORT: $0");

# Ate logo

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
    open (GFF, "$gffdir/CHROMOSOME_$chromosome.$filename.gff") or die "Cannot open $gffdir/CHROMOSOME_$chromosome.$filename.gff $!\n";
    while (<GFF>) {
	
	# discard header lines
	next if /\#/;

	my @f = split /\t/;

	($name) = ($f[8] =~ /\"(\S+)\"/)          if ($filename eq 'coding_exon');    # CDS name
	($name) = ($f[8] =~ /\"Sequence:(\S+)\"/) if ($filename ne 'coding_exon');    # transcript name
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

    print &runtime, " : Start searching, modus exon\n" if ($verbose);
    
    
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

	
	# No overlap between HSP and CDS exon
	
	# next gene if geneend is left of ESTstart 
	if ($eststart > $geneend) {                                                                                                # no match,            EXON  ----------------<========================--------========*=----------------
	    next GENE;                                                                                                             #                      HSP                                                                  ##########
	}
	
	# next EST if genestart is right of ESTend
	elsif ($estend < $genestart) {                                                                                             # no match,            EXON  ----------------<========================--------========*=----------------
	    $lastfail = $x;                                                                                                        #                      HSP   ###########
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
		
		if ( ($est_exon_start == $CDS_exon_start) && ($est_exon_end == $CDS_exon_end) ) {                                  #                                EXON_start ¬             EXON_stop ¬
		    map {$store_match{$testgene}->{$_} = 1;} ($est_exon_start..$est_exon_end);                                     # complete match,      EXON  ------->===^^^^^^^========================^^^^^^^^^========*-------------
		    next ESTEXON;                                                                                                  #                      HSP                  ---########################---
		}
		
		if ( ($est_exon_start < $CDS_exon_end) && ($est_exon_end == $CDS_exon_end)) { 
		    if ($est_exon_start >= $CDS_exon_start) {                                                                      #                                         HSP_start¬         EXON_stop¬
			map {$store_match{$testgene}->{$_} = 1;} ($est_exon_start..$est_exon_end);                                 # partial match        EXON  ------->===^^^^^^^========================^^^^^^^^^=======*-------------
			next ESTEXON;                                                                                              #                      HSP                          ###################
		    }
		    elsif ($est_exon_start < $CDS_exon_start) {                                                                    #                                   EXON_start¬              EXON_stop¬
			map {$store_match{$testgene}->{$_} = 1;} ($CDS_exon_start..$est_exon_end);                                 # complete match       EXON    ----->===^^^^^^^========================^^^^^^^^^========*-------------
			next ESTEXON;                                                                                              #                      HSP                  ###########################
		    }
		}
		
		if (($CDS_exon_start == $est_exon_start) && ($CDS_exon_end > $est_exon_start)) {   
		    
		    if ($CDS_exon_end <= $est_exon_end) {                                                                         #                                                HSP_start¬               CDS_stop¬
			map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$CDS_exon_end);                                  # complete match    EXON  ----------->==============^^^^^^^========================^^^^^^^^^==============*----------
			next ESTEXON;                                                                                             #                   HSP                                    ##########################
		    }
		    elsif ($CDS_exon_end > $est_exon_end) {                                                                       #                                                HSP_start¬           HSP_stop¬
			map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$est_exon_end);                                  # partial match     EXON  ----------->==============^^^^^^^========================^^^^^^^^^==============*----------
			next ESTEXON;                                                                                             #                   HSP                                    ####################                 
		    }
		}

		#################
		# SPECIAL CASES #
		#################
		
		# First HSP of transcript sequence ($v = 0)
		if ($v == 0) {
		    if (($est_exon_start < $CDS_exon_end) && ($est_exon_end == $CDS_exon_end)) { 
			if ($est_exon_start >= $CDS_exon_start) {                                                                  #                                         HSP_start¬         EXON_stop¬
			    map {$store_match{$testgene}->{$_} = 1;} ($est_exon_start..$est_exon_end);                             # partial match        EXON  ------->===^^^^^^^========================^^^^^^^^^=======*-------------
			    next GENEEXON;                                                                                         #                      HSP                          ###################
			}
			elsif ($est_exon_start < $CDS_exon_start) {                                                                #                                   EXON_start¬              EXON_stop¬
			    map {$store_match{$testgene}->{$_} = 1;} ($CDS_exon_start..$est_exon_end);                             # complete match       EXON    ----->===^^^^^^^========================^^^^^^^^^========*-------------
			    next GENEEXON;                                                                                         #                      HSP                  ###########################
			}
		    }
		} #_ end of initial HSP processing
		    
		# Single HSP for transcript  ($v = 0 AND $v = $estblocks)  i.e. single HSP EST
		if ($v == $estblocks) {
			
		    # single exon EST within geneexon
		    if (($est_exon_start >= $CDS_exon_start) && ($est_exon_end <= $CDS_exon_end)) {                               #                                      HSP_start ¬        HSP_stop ¬
			map {$store_match{$testgene}->{$_} = 1;} ($est_exon_start..$est_exon_end);                                # complete match          EXON  ----->===^^^^^^^========================^^^^^^^^^========*-------------
			next EST;                                                                                                 #                         HSP                     ##################
		    }
			
		    # single exon EST overlapping gene start
		    elsif (($est_exon_start <= $CDS_exon_start) && ($est_exon_end <= $CDS_exon_end) && ($w == 0)) {               #                                   EXON_start ¬    HSP_stop ¬
			map {$store_match{$testgene}->{$_} = 1;} ($CDS_exon_start..$est_exon_end);                                # partial match           EXON  ---------------->=======================^^^^^^^^^========*-------------
			next EST;                                                                                                 #                         HSP             ####################
		    }
		    
		    # single exon EST overlapping gene end
		    elsif (($est_exon_start >= $CDS_exon_start) && ($est_exon_end >= $CDS_exon_end) && ($w == $geneblocks)) {     #                                           HSP_start ¬      EXON_stop ¬
			map {$store_match{$testgene}->{$_} = 1;} ($est_exon_start..$CDS_exon_end);                                # partial match           EXON  ---->====^^^^^^^=======================*-------------------------------
			next EST;                                                                                                 #                         HSP                          ####################
		    }
			
		    # single exon EST covering whole single exon gene
		    elsif (($est_exon_start < $CDS_exon_start) && ($est_exon_end > $CDS_exon_end) && ($geneblocks == 0)) {        #                                             EXON_start¬                  EXON_stop ¬
			map {$store_match{$testgene}->{$_} = 1;} ($CDS_exon_start..$CDS_exon_end);                                # complete match          EXON  ------------------------->===========================*-------------------------------
			next GENE;                                                                                                #                         HSP                        ####################################
		    }
		} #_ end of Single HSP processing

		# Final HSP of transcribed sequence ($v = $estblocks)
		if ($v == $estblocks) {
		    if ( ($est_exon_start == $CDS_exon_start) && ($est_exon_end > $CDS_exon_start) ) { 
			if ($est_exon_end <= $CDS_exon_end) {                                                                     #                                                                                EXON_start¬  HSP_stop ¬
			    map {$store_match{$testgene}->{$_} = 1;} ($est_exon_start..$est_exon_end);                            # partial match     EXON  ----------->==============^^^^^^^========================^^^^^^^^^==============*----------
			    next EST;                                                                                             #                   HSP                                                                     ############
			}
			elsif ($est_exon_end > $CDS_exon_end) {                                                                   #                                                                                EXON_start¬     EXON_stop¬
			    map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$CDS_exon_end);                              # complete match    EXON  ----------->==============^^^^^^^=======================^^^^^^^^^^==============*----------
			    next EST;                                                                                             #                   HSP                                                                     #################
			}  
		    }
		} #_ end of final HSP processing
		
		# First CDS exon ($w = 0)
		if ($w == 0) {
		    if ( ($CDS_exon_start < $est_exon_end) && ($CDS_exon_end == $est_exon_end) ) {
			if ($CDS_exon_start >= $est_exon_start) {                                                                  #                         EXON_start¬     EXON_stop¬
			    map {$store_match{$testgene}->{$_}=1;} ($CDS_exon_start..$est_exon_end);                              # partial match     EXON  ----------->==============^^^^^^^========================^^^^^^^^^==============*---------- 
			    next ESTEXON;                                                                                         #                   HSP           ##################
			}
			elsif ($CDS_exon_start < $est_exon_start) {                                                               #                              HSP_start¬ EXON_stop¬
			    map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$est_exon_end);                              # partial match     EXON  ----------->==============^^^^^^^========================^^^^^^^^^==============*----------
			    next ESTEXON;                                                                                         #                   HSP                  ###########
			}   
		    }
		} #_ end of intial CDS exon processing
		    
		# Single exon CDS  ($w = 0 AND $w = $geneblocks)
		if ($w == $geneblocks) {
			
		    # single exon gene completely covered by EST
		    if (($est_exon_start <= $CDS_exon_start ) && ($est_exon_end >= $CDS_exon_end)) {                              #                                             EXON_start¬                   EXON_stop¬
			map {$store_match{$testgene}->{$_}=1;} ($CDS_exon_start..$CDS_exon_end);                                  # complete match    EXON  ------------------------->===========================*-------------------------------
			next GENE;                                                                                                #                   HSP                         ####################################
		    }   
		    
		    # single exon gene end covered
		    elsif (($est_exon_start >= $CDS_exon_start) && ($est_exon_end >= $CDS_exon_end) && ($v == 0)) {               #                                                  HSP_start¬               EXON_stop¬
			map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$CDS_exon_end);                                  # partial match     EXON  ------------------------->===========================*-------------------------------
			next EST;                                                                                                 #                   HSP                                ####################################
		    }
		    
		    # single exon gene start covered
		    elsif (($est_exon_start <= $CDS_exon_start) && ($est_exon_end <= $CDS_exon_end) && ($v == $estblocks)) {      #                                             EXON_start¬                 HSP_stop¬
			map {$store_match{$testgene}->{$_}=1;} ($CDS_exon_start..$est_exon_end);                                  # partial match     EXON  ------------------------->===========================*-------------------------------
			next EST;                                                                                                 #                   HSP                  ####################################
		    }
		    
		    # single exon EST inside single exon gene
		    elsif (($est_exon_start >= $CDS_exon_start) && ($est_exon_end <= $CDS_exon_end) && ($estblocks == 0)) {       #                                                 HSP_start¬            HSP_stop¬
			map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$est_exon_end);                                  # complete match    EXON  ------------------------->===========================*-------------------------------
			next EST;                                                                                                 #                   HSP                               #####################
		    }
		}
		 
		if ($w == $geneblocks) {
		    if (($CDS_exon_start == $est_exon_start) && ($CDS_exon_end > $est_exon_start)) {   
		    
			if ($CDS_exon_end <= $est_exon_end) {                                                                     #                                                HSP_start¬               CDS_stop¬
			    map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$CDS_exon_end);                              # complete match    EXON  ----------->==============^^^^^^^========================^^^^^^^^^==============*----------
			    next EST;                                                                                             #                   HSP                                    ##########################
			}
			elsif ($CDS_exon_end > $est_exon_end) {                                                                   #                                                HSP_start¬           HSP_stop¬
			    map {$store_match{$testgene}->{$_}=1;} ($est_exon_start..$est_exon_end);                              # partial match     EXON  ----------->==============^^^^^^^========================^^^^^^^^^==============*----------
			    next EST;                                                                                             #                   HSP                                    ####################                 
			}
		    }
		}
		

	    } #_ ESTEXON
	  }   #_ GENEEXON
	}
    }         #_ EST
  }           #_ GENE
    
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












