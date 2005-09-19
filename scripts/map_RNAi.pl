#!/usr/local/bin/perl5.8.0 -w
#
# map_RNAi.pl
#
# Add information to RNAi objects based on overlaps in GFF files 
#
# by Kerstin Jekosch
#
# Last updated by: $Author: gw3 $                      
# Last updated on: $Date: 2005-09-19 12:38:51 $        

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

my $tace        = &tace;                                            # tace executable path
my $dbdir       = "/wormsrv2/autoace";                              # Database path
my $gffdir      = "/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS";        # GFF_SPLITS directory
my @chromosomes = qw( I II III IV V X );                            # chromosomes
my $db_version  = &get_wormbase_version_name;                       # WS version name

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
my $load;       # use to automatically load file to autoace

my $log = Log_files->make_build_log();

GetOptions ("debug=s"   => \$debug,
	    "verbose"   => \$verbose,
	    "test"      => \$test,
            "help"      => \$help,
	    "load"      => \$load);


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
  #
  # New RNAi lines : CHROMOSOME_I    RNAi_primary    RNAi_reagent    1681680 1683527 .       .       .       Target "RNAi:WBRNAi00004820" 1 1848


  print "Loop through primary RNAi GFF file CHROMOSOME_${chromosome}\n" if ($verbose);
  open (GFF, "<$gffdir/CHROMOSOME_${chromosome}.RNAi_primary.gff") || die "Failed to open RNAi gff file\n\n";
  while (<GFF>) {
    chomp;
    s/^\#.*//;
    next unless /\S/;
    @line = split /\t/;
    
    my ($name) = ($line[8] =~ /\"RNAi:(\S+.+)\"\s+\d+\s+\d+$/);
    $RNAicount{$name}++;
    my $RNAiname = $name.".".$RNAicount{$name};
    # NB store type of RNAi (primary) with the start/end positions 
    $RNAi{$RNAiname} = [$line[3],$line[4], "primary"];
    print "RNAi : '$name'\n" if ($verbose);
    
  }
  close(GFF);

  # add the seondary RNAi hits to the same data structure
  # note which is secondary by adding "secondary" to the gene mapped to
  print "Loop through secondary RNAi GFF file CHROMOSOME_${chromosome}\n" if ($verbose);
  open (GFF, "<$gffdir/CHROMOSOME_${chromosome}.RNAi_secondary.gff") || die "Failed to open RNAi gff file\n\n";
  while (<GFF>) {
    chomp;
    s/^\#.*//;
    next unless /\S/;
    @line = split /\t/;
    
    my ($name) = ($line[8] =~ /\"RNAi:(\S+.+)\"\s+\d+\s+\d+$/);
    $RNAicount{$name}++;
    my $RNAiname = $name.".".$RNAicount{$name};
    # NB store type of RNAi (secondary) with the start/end positions 
    $RNAi{$RNAiname} = [$line[3],$line[4], "secondary"];
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
	      # strip the $RNAicount value from the end of the RNAi name
	    my ($RNAi) = ($testRNAi =~ /(\S+)\.\d+$/);
	    # store the gene the RNAi maps to and its type (primary/secondary)
	    push @{$output{$RNAi}}, [$testgene, $RNAi{$testRNAi}->[2]];
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

# store some statistics
# count number of primary RNAis (after removing duplicate mapping to same gene)
my $no_of_rnais_primary = 0; 
# count number of secondary RNAis (after removing duplicate mappings to same gene)
my $no_of_rnais_secondary = 0; 
# count number of RNAis which map to the same gene more than once
my $no_of_multiple_mappings = 0;
# count secondary RNAis which map to a gene which is already mapped by that RNAi as a primary
my $no_of_duplicate_secondaries = 0;

############################
# sort the output for RNAi #
############################

print "Remove duplicates for RNAi\n" if ($verbose);

foreach my $RNAiID (sort keys %output) {
  # sort by the gene names
  # with secondary sort key term to retain primary RNAis in preference to secondary ones
  # this was done by adding a term to also sort by "primary" before "secondary": 
  # 'or  $a->[1] cmp $b->[1]'
  # so that primary RNAi mapping is sorted to before any secondary mapping

  # the %output values are a list of lists - each hash value contains a list of:
  # [the names of genes, primary / secondary RNAi mapped to this gene]
  # we now wish to sort this list of pairs of values by gene name
  my @array = @{$output{$RNAiID}};
  @{$output{$RNAiID}} = sort { $a->[0] cmp $b->[0]  or  $a->[1] cmp $b->[1] } @array;

  my $count = 0; 
  push @{$finaloutput{$RNAiID}}, $output{$RNAiID}->[0];
  for (my $m = 1; $m < (scalar @{$output{$RNAiID}}); $m++) {
    # compare gene names to get a unique, sorted set of genes
    if ($output{$RNAiID}->[$count][0] ne $output{$RNAiID}->[$m][0]) {
      push @{$finaloutput{$RNAiID}}, $output{$RNAiID}->[$m];
      $count = $m;
    } else {
      # get some statistics on duplicate mappings
      $no_of_multiple_mappings++;
      if ($output{$RNAiID}->[$count][1] ne $output{$RNAiID}->[$m][1]) {
	$no_of_duplicate_secondaries++;
      }
    }
  }
}

####################################
# sort the output for Expr_profile #
####################################

print "Remove duplicates for Expr_profile\n" if ($verbose);

foreach my $RNAiID (sort keys %output2) {
  @{$output2{$RNAiID}} = sort @{$output2{$RNAiID}};
  my $count = 0; 
  push @{$finaloutput2{$RNAiID}}, $output2{$RNAiID}->[0];
  for (my $m = 1; $m < (scalar @{$output2{$RNAiID}}); $m++) {
    if ($output2{$RNAiID}->[$count] ne $output2{$RNAiID}->[$m]) {
      push @{$finaloutput2{$RNAiID}}, $output2{$RNAiID}->[$m];
      $count = $m;    
    }    
  }
}



########################
# produce output files #
########################

print "Produce output file\n" if ($verbose);

# store the gene names and the RNAis that hit it, together with the RNAi's primary/secondary type
# so that we can later populate the ?Gene model RNAi data with primary/secondary information
my %inverse = ();


open (OUTACE, ">$dbdir/acefiles/RNAi_mappings.ace") || die "Failed to open RNAi_mappings.ace file\n";

# Produce connections for RNAi->Genes

# remove existing connections
foreach my $mapped (sort keys %finaloutput) {     
  print OUTACE "RNAi : $mapped\n";
  print OUTACE "-D Inhibits\n\n";
}

# open autoace connection with aceperl
my $db = Ace->connect(-path=>$dbdir, -program=>$tace) || die "Couldn't connect to $dbdir\n", Ace->error;

foreach my $mapped (sort keys %finaloutput) {
    
  # $mapped is the name of the RNAi object
  # %finaloutput is a hash of arrays containing the gene overlap data

    
  my $seq;
  my $gene;       # e.g. WBGene00001231
  my $worm_gene;  # CDS, Transcript, or Pseudogene name

  # $type holds whether the RNAi mapping to this gene is primary or secondary
  my $type;

  for (my $n = 0; $n < (scalar @{$finaloutput{$mapped}}); $n++) {

    $worm_gene = $finaloutput{$mapped}->[$n][0];
    $type      = $finaloutput{$mapped}->[$n][1];

    print "$type RNAi '$mapped' is mapped to $worm_gene\t" if ($verbose);

    if ($type eq "primary") {
      $no_of_rnais_primary++;
    } else {
      $no_of_rnais_secondary++;
    }

    print OUTACE "\nRNAi : \"$mapped\"\n";
      
    # Does this CDS have a Gene object?
    if ($genetype{$worm_gene} eq "CDS") {
      $seq = $db->fetch(-class=>'CDS',-name=>" $worm_gene");
      if (defined $seq) { 
	print OUTACE "Predicted_gene \"$worm_gene\" Inferred_automatically \"RNAi_$type\"\n";
	if (defined($seq->at('Visible.Gene'))) {
	  ($gene) = ($seq->get('Gene'));
	  print OUTACE "Gene \"$gene\" Inferred_automatically \"RNAi_$type\"\n";
	  push @{$inverse{$gene}}, [$mapped, $type];
	}
	print " which is a CDS\n" if ($verbose);
      } else {
	print "*** WARNING - skipping missing gene $worm_gene\n";
      }
    }
    # Does this Pseudogene have a Gene object?
    elsif ($genetype{$worm_gene} eq "Pseudogene") {
      $seq = $db->fetch(-class=>'Pseudogene',-name=>" $worm_gene");
      if (defined $seq) {
	print OUTACE "Pseudogene \"$worm_gene\" Inferred_automatically \"RNAi_$type\"\n";
	if (defined($seq->at('Visible.Gene'))){
	  ($gene) = ($seq->get('Gene'));
	  print OUTACE "Gene \"$gene\" Inferred_automatically \"RNAi_$type\"\n";
	  push @{$inverse{$gene}}, [$mapped, $type];
	}
	print " which is a pseudogene\n" if ($verbose);
      } else {
	print "*** WARNING - skipping missing gene $worm_gene\n";  
      }
    }
    # Does this transcript have a Gene object?
    elsif ($genetype{$worm_gene} eq "Transcript") {
      $seq = $db->fetch(-class=>'Transcript',-name=>" $worm_gene");
      if (defined $seq) {
	print OUTACE "Transcript \"$worm_gene\" Inferred_automatically \"RNAi_$type\"\n";
	if (defined($seq->at('Visible.Gene'))){
	  ($gene) = ($seq->get('Gene'));
	  print OUTACE "Gene \"$gene\" Inferred_automatically \"RNAi_$type\"\n";
	  push @{$inverse{$gene}}, [$mapped, $type];
	}
	print " which is a transcript\n" if ($verbose);
      } else {
	print "*** WARNING - skipping missing gene $worm_gene\n";
      }
    }
  }
} 

$db->close;

print OUTACE "\n\n//Expression profiles\n";

# Produce connections for RNAi->Expr_profile
foreach my $mapped (sort keys %finaloutput2) {
  for (my $n = 0; $n < (scalar @{$finaloutput2{$mapped}}); $n++) {
    my ($expr_profile) = ($finaloutput2{$mapped}->[$n] =~ /(\S+)\.\d+$/);
    print OUTACE "\nRNAi : \"$mapped\"\n";
    print OUTACE "Expr_profile $expr_profile\n";
  }
} 


# Write primary/secondary type evidence for RNAi tag in ?Gene model
# (This is done explicitly because you can't (easily) make a XREF between #Evidence tags in the models)
foreach my $gene (sort keys %inverse) {
  for (my $n = 0; $n < (scalar @{$inverse{$gene}}); $n++) {
    my $RNAi = $inverse{$gene}->[$n][0];
    my $type = $inverse{$gene}->[$n][1];
    print "Gene: $gene RNAi: $RNAi type: $type\n" if ($verbose);
    print OUTACE "\nGene : \"$gene\"\n";
    print OUTACE "Experimental_info RNAi_result  \"$RNAi\" Inferred_automatically \"RNAi_$type\"\n";
  }
}

close(OUTACE);


#########################################################
# read acefiles into autoace (unless running test mode) #
#########################################################
if($load){
  $log->write_to("Loading file to autoace\n");
  my $command = "autoace_minder.pl -load $dbdir/acefiles/RNAi_mappings.ace -tsuser RNAi_mappings";
  
  my $status = system($command);
  if(($status >>8) != 0){
    $log->write_to("ERROR: Loading RNAi_mappings.ace file failed \$\? = $status\n");
  }
}

####################################
# print some statistics to the log #
####################################

# count secondary RNAis which map to a gene which is already mapped by that RNAi as a primary
$log->write_to("\n\nStatistics\n");
$log->write_to("----------\n\n");
$log->write_to("No. of primary   RNAi links to genes written to database: $no_of_rnais_primary\n");
$log->write_to("No. of secondary RNAi links to genes written to database: $no_of_rnais_secondary\n\n");
#$log->write_to("The following duplicates have been detected and not written to the database...\n");
#$log->write_to("No. of RNAis found mapping more than once to the same gene: $no_of_multiple_mappings\n");
#$log->write_to("No. of secondary RNAis found mapping to a gene already mapped by that RNAi as a primary: $no_of_duplicate_secondaries\n");

$log->mail("$maintainers","BUILD REPORT: $0");

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
connections and make new ones. It will check the current database and make Gene
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

=item -load, loads file to autoace

=back

=cut
