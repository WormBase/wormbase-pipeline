#!/usr/local/bin/perl5.8.0 -w
#
# map_RNAi.pl
#
# Add information to RNAi objects based on overlaps in GFF files 
#
# by Kerstin Jekosch
#
# Last updated by: $Author: krb $                      
# Last updated on: $Date: 2004-10-19 11:14:40 $        

use strict;
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Cwd;
use Ace;


######################################
# variables and command-line options #
######################################

my $tace        = &tace;                                  # tace executable path
my $dbdir       = "/wormsrv2/autoace";                    # Database path
my $gffdir      = "/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS";        # GFF_SPLITS directory
my @chromosomes = qw( I II III IV V X );                  # chromosomes
my $db_version  = &get_wormbase_version_name;             # WS version name

my %output      = (); # for RNAi
my %output2     = (); # for Expression profiles

my %finaloutput  = (); # for RNAi
my %finaloutput2 = (); # for Expression profiles

my $maintainers = "All";

my $rundate = &rundate;
my $runtime = &runtime;

my $name;
my $help;       # Help perdoc
my $test;       # Test mode
my $debug;      # Debug mode, verbose output to user running script
my $verbose;

my $log = Log_files->make_build_log();

GetOptions ("debug=s"   => \$debug,
	    "verbose"   => \$verbose,
	    "test"      => \$test,
            "help"      => \$help);


# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

if( $test ) {
  $dbdir       = glob("~wormpub/TEST_BUILD/autoace");                    # Database path
  $gffdir      = glob("~wormpub/TEST_BUILD/autoace/GFF_SPLITS/GFF_SPLITS");        # GFF_SPLITS directory
  @chromosomes = qw( III );
}


##########################
# MAIN BODY OF SCRIPT
##########################



###########################################################
# get exons, RNAis and Expr_profiles out of the gff files #
########################################################### 
   
my @line;

our %genetype         = ();


foreach my $chromosome (@chromosomes) {

  print "Processing chromosome $chromosome\n";
  $log->write_to("Processing chromosome $chromosome\n");

  my %RNAi             = ();
  my %RNAicount        = ();
  my %genes            = ();
  my %exon             = ();
  my %exoncount        = ();
  my %expression       = ();
  my %expression_count = ();
  
  # loop through the split GFF RNAi file  
  print "Loop through RNAi GFF file CHROMOSOME_${chromosome}\n" if ($verbose);
  open (GFF, "<$gffdir/CHROMOSOME_${chromosome}.RNAi.gff") || die "Failed to open RNAi gff file\n\n";
  while (<GFF>) {
    chomp;
    s/^\#.*//;
    next unless /\S/;
    @line = split /\t/;
    
    my ($name) = ($line[8] =~ /\"(\S+.+)\"$/);
    $RNAicount{$name}++;
    my $RNAiname = $name.".".$RNAicount{$name};
    $RNAi{$RNAiname} = [$line[3],$line[4]];
    print "RNAi : '$name'\n" if ($verbose);
    
  }
  close(GFF);


  # loop through the split GFF Expr_profile file  
  print "Loop through Expr_profile GFF file CHROMOSOME_${chromosome}\n" if ($verbose);
  open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.Expr_profile.gff") || die "Failed to open Expr_profile gff file\n\n";
  while (<GFF_in>) {
    chomp;
    s/^\#.*//;
    next unless /\S/;
    @line = split /\t/;
    
    my ($name) = ($line[8] =~ /\"(.*)\"$/);
    $expression_count{$name}++;
    my $expression_name = $name.".".$expression_count{$name};
    $expression{$expression_name} = [$line[3],$line[4]];
    print "Expr_profile : '$name'\n" if ($verbose);    
  }
  close(GFF_in);


  # loop through the split GFF exon file  
  print "Loop through exon GFF file CHROMOSOME_${chromosome}\n" if ($verbose);
  open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.exon.gff") || die "Failed to open exon gff file\n\n";
  while (<GFF_in>) {
    chomp;
    s/^\#.*//;
    next unless /\S/;
    @line = split /\t/;

    my ($name) = ($line[8] =~ /\"(\S+.+)\"/);
    $exoncount{$name}++;
    my $exonname = $name.".".$exoncount{$name};
    $exon{$exonname} = [$line[3],$line[4]];
    $genetype{$name} = "CDS"        if ($line[1] eq "curated");
    print "Gene : '$name' [$genetype{$name}] exon $exoncount{$name}\n" if ($verbose);
  }
  close(GFF_in);
  

  # loop through the split GFF exon_pseudogene file  
  print "Loop through exon_pseudogene GFF file CHROMOSOME_${chromosome}\n" if ($verbose);
  open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.exon_pseudogene.gff") || die "Failed to open exon pseudogene gff file\n\n";
  while (<GFF_in>) {
    chomp;
    s/^\#.*//;
    next unless /\S/;
    @line = split /\t/;

    my ($name) = ($line[8] =~ /\"(\S+.+)\"/);
    $exoncount{$name}++;
    my $exonname = $name.".".$exoncount{$name};
    $exon{$exonname} = [$line[3],$line[4]];
    $genetype{$name} = "Pseudogene" if ($line[1] eq "Pseudogene");
    print "Gene : '$name' [$genetype{$name}] exon $exoncount{$name}\n" if ($verbose);
  }
  close(GFF_in);

  
  # loop through the split GFF exon_noncoding file  
  print "Loop through exon_noncoding GFF file CHROMOSOME_${chromosome}\n" if ($verbose);
  open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.exon_noncoding.gff") || die "Failed to open exon_noncoding gff file\n\n";
  while (<GFF_in>) {
    chomp;
    s/^\#.*//;
    next unless /\S/;
    @line = split /\t/;

    my ($name) = ($line[8] =~ /\"(\S+.+)\"/);
    $exoncount{$name}++;
    my $exonname = $name.".".$exoncount{$name};
    $exon{$exonname} = [$line[3],$line[4]];
    $genetype{$name} = "Transcript" if ($line[1] eq "Non_coding_transcript");
    print "Gene : '$name' [$genetype{$name}] exon $exoncount{$name}\n" if ($verbose);
  }
  close(GFF_in);
  


  print "Finished GFF loop\n" if ($verbose);
  
  #########################   
  # make exons into genes #
  #########################
  
  print "Turn exons into genes\n" if ($verbose);
  
  foreach my $name (sort keys %exoncount) {
    my $v = $exoncount{$name};
    my $w = $name.".".$v;
    $genes{$name} = [$exon{$name.".1"}->[0],$exon{$w}->[1]];
  }
  
  ###################
  # make indexlists #
  ###################
  
  print "Index exons,genes,RNAi and Expr_profiles\n" if ($verbose);
  
  my @exonlist       = sort { $exon{$a}->[0]  <=> $exon{$b}->[0]  } keys %exon;
  my @genelist       = sort { $genes{$a}->[0] <=> $genes{$b}->[0] } keys %genes;
  my @RNAilist       = sort { $RNAi{$a}->[0]  <=> $RNAi{$b}->[0]             || $a cmp $b } keys %RNAi;     
  my @expressionlist = sort { $expression{$a}->[0]  <=> $expression{$b}->[0] || $a cmp $b } keys %expression;
  
  #############
  # map RNAis #
  #############
  
  print "Find overlaps for RNAi\n" if ($verbose);
  
  my $lastfail = 0;
  
  for (my $x = 0; $x < @RNAilist; $x++) {
    my $testRNAi   = $RNAilist[$x];
    my $RNAistart  = $RNAi{$testRNAi}->[0];
    my $RNAistop   = $RNAi{$testRNAi}->[1];
    
    for (my $y = 0; $y < @genelist; $y++) {
      my $testgene = $genelist[$y];
      my $genestart= $genes{$testgene}->[0];
      my $genestop = $genes{$testgene}->[1];
      
      if ($RNAistart > $genestop) {
	$lastfail = $y;
	next; 
      }
      elsif ($RNAistop < $genestart) {
	last; 
      }
      
      else {
	for (my $z = 1; $z <= $exoncount{$testgene}; $z++) {
	  my $exon_start = $exon{"$testgene.$z"}->[0];
	  my $exon_stop  = $exon{"$testgene.$z"}->[1];
	  
	  if ( not (($RNAistart > $exon_stop) || ($RNAistop < $exon_start))) {
	    my ($RNAi) = ($testRNAi =~ /(\S+)\.\d+$/);
	    push @{$output{$RNAi}}, $testgene;
	  }
	}
      }                
    }
  }
  
  ###########################
  # map Expression profiles #
  ###########################
  
  print "Find overlaps for Expr_profile\n" if ($verbose);
  
  $lastfail = 0;
    
  for (my $x = 0; $x < @RNAilist; $x++) {
    my $testRNAi   = $RNAilist[$x];
    my $RNAistart  = $RNAi{$testRNAi}->[0];
    my $RNAistop   = $RNAi{$testRNAi}->[1];

    
    for (my $y = 0; $y < @expressionlist; $y++) {
      my $testexpression = $expressionlist[$y];

      my $expressionstart= $expression{$testexpression}->[0];
      my $expressionstop = $expression{$testexpression}->[1];

      
      if ($RNAistart > $expressionstop) {
	$lastfail = $y;
	next; 
      }
      elsif ($RNAistop < $expressionstart) {
	last; 
      }
      
      else {
	my ($RNAi) = ($testRNAi =~ /(\S+)\.\d+$/);
	push @{$output2{$RNAi}}, $testexpression;
	print "$RNAi mapped to $testexpression\n" if ($verbose);
      }                
    }
  }
}

############################
# sort the output for RNAi #
############################
 
foreach my $mess (sort keys %output) {
  @{$output{$mess}} = sort @{$output{$mess}};
  my $count = 0; 
  push @{$finaloutput{$mess}}, $output{$mess}->[0];
  for (my $m = 1; $m < (scalar @{$output{$mess}}); $m++) {
    if ($output{$mess}->[$count] ne $output{$mess}->[$m]) {
      push @{$finaloutput{$mess}}, $output{$mess}->[$m];
      $count = $m;    
    }    
  }
}

####################################
# sort the output for Expr_profile #
####################################
 
foreach my $mess (sort keys %output2) {
  @{$output2{$mess}} = sort @{$output2{$mess}};
  my $count = 0; 
  push @{$finaloutput2{$mess}}, $output2{$mess}->[0];
  for (my $m = 1; $m < (scalar @{$output2{$mess}}); $m++) {
    if ($output2{$mess}->[$count] ne $output2{$mess}->[$m]) {
      push @{$finaloutput2{$mess}}, $output2{$mess}->[$m];
      $count = $m;    
    }    
  }
}



########################
# produce output files #
########################

open (OUTACE, ">$dbdir/acefiles/RNAi_mappings.ace");

# Produce connections for RNAi->Genes

# remove existing connections
foreach my $mapped (sort keys %finaloutput) {   
  
  for (my $n = 0; $n < (scalar @{$finaloutput{$mapped}}); $n++) {
    print OUTACE "RNAi : $mapped\n";
    print OUTACE "-D Inhibits\n\n";
  }
} 

# open autoace connection with aceperl
my $db = Ace->connect(-path=>$dbdir, -program=>$tace) || die "Couldn't connect to $dbdir\n", Ace->error;

foreach my $mapped (sort keys %finaloutput) {
    
  # $mapped is the name of the RNAi object
  # %finaloutput is a hash of arrays containing the gene o/l data
 
  
  my $seq;
  my $locus;
  my $gene;
  
  for (my $n = 0; $n < (scalar @{$finaloutput{$mapped}}); $n++) {
    
    $gene = $finaloutput{$mapped}->[$n];
    print "'$mapped' is mapped to $gene\t" if ($verbose);
    
    print OUTACE "\nRNAi : \"$mapped\"\n";
    
    my $seq;
    my $gene;       # e.g. WBGene00001231
    my $worm_gene;  # CDS, Transcript, or Pseudogene name

    for (my $n = 0; $n < (scalar @{$finaloutput{$mapped}}); $n++) {
	
      $worm_gene = $finaloutput{$mapped}->[$n];
      print "'$mapped' is mapped to $worm_gene\t" if ($verbose);

      print OUTACE "\nRNAi : \"$mapped\"\n";
      
      # Does this CDS have a Gene object?
      if ($genetype{$worm_gene} eq "CDS") {
	print OUTACE "Predicted_gene \"$worm_gene\"\n";
	$seq = $db->fetch(-class=>'CDS',-name=>" $finaloutput{$mapped}->[$n]");
	if (defined($seq->at('Visible.Gene'))) {
	  ($gene) = ($seq->get('Gene'));
	  print OUTACE "Gene \"$gene\"\n";
	}
	print " which is a CDS\n" if ($verbose);
      }
      # Does this Pseudogene have a Gene object?
      elsif ($genetype{$worm_gene} eq "Pseudogene") {
	print OUTACE "Pseudogene \"$worm_gene\"\n";
	$seq = $db->fetch(-class=>'Pseudogene',-name=>" $finaloutput{$mapped}->[$n]");
	if (defined($seq->at('Visible.Gene'))){
	  ($gene) = ($seq->get('Gene'));
	  print OUTACE "Gene \"$gene\"\n";
	}
	print " which is a pseudogene\n" if ($verbose);
      }
      # Does this transcript have a Gene object?
      elsif ($genetype{$worm_gene} eq "Transcript") {
	print OUTACE "Transcript \"$worm_gene\"\n";
	$seq = $db->fetch(-class=>'Transcript',-name=>" $finaloutput{$mapped}->[$n]");
	if (defined($seq->at('Visible.Gene'))){
	  ($gene) = ($seq->get('Gene'));
	  print OUTACE "Gene \"$gene\"\n";
	}
	print " which is a transcript\n" if ($verbose);
      }
    }
  }
} 

$db->close;

print OUTACE "\n\n//Expression profiles\n";

# Produce connections for RNAi->Expr_profile
foreach my $mapped (sort keys %finaloutput2) {
  for (my $n = 0; $n < (scalar @{$finaloutput2{$mapped}}); $n++) {
    my ($expr_profile) = (@{$finaloutput2{$mapped}}->[$n] =~ /(\S+)\.\d+$/);
    print OUTACE "RNAi : $mapped\n";
    print OUTACE "Expr_profile $expr_profile\n\n";
  }
} 

close(OUTACE);


#########################################################
# read acefiles into autoace (unless running test mode) #
#########################################################

unless ($test) {

  my $command = "pparse /wormsrv2/autoace/acefiles/RNAi_mappings.ace\nsave\nquit\n";

  open (TACE,"| $tace -tsuser map_RNAi $dbdir") || die "Couldn't open tace connection to $dbdir\n";
  print TACE $command;
  close (TACE);
}

exit(0);

###############################
# Prints help and disappears  #
###############################

##############################################################
#
# Subroutines
#
##############################################################


sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

############################################

__END__

=pod

=head2 NAME - map_RNAi.pl

=head1 USAGE

=over 4

=item map_RNAi.pl [-options]

=back

map_RNAi.pl calculates the overlap between the genomic regions used in RNAi
experiments and the CDS, transcript and pseudogene coordinates in the WS
database release. It will generate an acefile which will remove any existing
connections and make new ones. It wil check the current database and make Gene
connections where valid and attach expression_profiles as needed. This acefile
is then loaded into /wormsrv2/autoace

map_RNAi mandatory arguments:


=over 4

=item none

=back

map_RNAi optional arguments:

=over 4

=item -debug, Verbose/Debug mode

=item -test, Test mode, generate the acefile but do not upload them 

=item -help, Help pages

=back

=cut
