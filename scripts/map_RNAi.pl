#!/usr/local/bin/perl5.6.1 -w
#
# map_RNAi.pl
#
# Add information to RNAi objects based on overlaps in GFF files 
#
# by Kerstin Jekosch
#
# Last updated by: $Author: krb $                      
# Last updated on: $Date: 2003-03-14 09:29:03 $        


$|=1;
use strict;
use lib "/wormsrv2/scripts/";
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
my $gffdir      = "/wormsrv2/autoace/CHROMOSOMES";        # CHROMOSOME dir
my @chromosomes = qw( I II III IV V X );                  # chromosomes
my $db_version  = &get_wormbase_version_name;             # WS version name

my %output      = (); # for RNAi
my %output2     = (); # for Expression profiles

my %finaloutput  = (); # for RNAi
my %finaloutput2 = (); # for Expression profiles

my $maintainers = "All";
my $rundate = `date +%y%m%d`; chomp $rundate;
my $runtime = `date +%H:%M:%S`; chomp $runtime;
my $name;
my $help;       # Help perdoc
my $debug;      # Debug mode, verbose output to user running script
our $log;

GetOptions ("debug=s"   => \$debug,
            "help"      => \$help);


# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

&create_log_files;


##########################
# MAIN BODY OF SCRIPT
##########################



###########################################################
# get exons, RNAis and Expr_profiles out of the gff files #
########################################################### 
   
my @line;

foreach my $chromosome (@chromosomes) {
  my %RNAi             = ();
  my %RNAicount        = ();
  my %genes            = ();
  my %exon             = ();
  my %exoncount        = ();
  my %expression       = ();
  my %expression_count = ();
  open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.gff") || die "Failed to open gff file\n\n";
  while (<GFF_in>) {
    chomp;
    s/^\#.*//;
    next unless /\S/;
    @line = split /\t/;
    
    if ($line[1] eq "RNAi") {
      my ($name) = ($line[8] =~ /\"(.*)\"$/);
      $RNAicount{$name}++;
      my $RNAiname = $name.".".$RNAicount{$name};
      $RNAi{$RNAiname} = [$line[3],$line[4]];
    }

    elsif ($line[1] eq "Expr_profile") {
      my ($name) = ($line[8] =~ /\"(.*)\"$/);
      $expression_count{$name}++;
      my $expression_name = $name.".".$expression_count{$name};
      $expression{$expression_name} = [$line[3],$line[4]];
    }

    elsif (($line[2] eq "exon") && (($line[1] eq "curated") || 
	   ($line[1] eq "Pseudogene") || ($line[1] eq "provisional"))) {
      my ($name) = ($line[8] =~ /\"(\S+)\"/);
      $exoncount{$name}++;
      my $exonname = $name.".".$exoncount{$name};
      $exon{$exonname} = [$line[3],$line[4]];
    }

  }
  close(GFF_in);
  
  #########################   
  # make exons into genes #
  #########################
  
  foreach my $name (sort keys %exoncount) {
    my $v = $exoncount{$name};
    my $w = $name.".".$v;
    $genes{$name} = [$exon{$name.".1"}->[0],$exon{$w}->[1]];
  }
  
  ###################
  # make indexlists #
  ###################
  
  my @exonlist       = sort { $exon{$a}->[0]  <=> $exon{$b}->[0]  } keys %exon;
  my @genelist       = sort { $genes{$a}->[0] <=> $genes{$b}->[0] } keys %genes;
  my @RNAilist       = sort { $RNAi{$a}->[0]  <=> $RNAi{$b}->[0]             || $a cmp $b } keys %RNAi;     
  my @expressionlist = sort { $expression{$a}->[0]  <=> $expression{$b}->[0] || $a cmp $b } keys %expression;


  
  #############
  # map RNAis #
  #############
  
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
	    print LOG "$RNAi mapped to $testgene\n";
	  }
	}
      }                
    }
  }

  ###########################
  # map Expression profiles #
  ###########################
  
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
	print LOG "$RNAi mapped to $testexpression\n";
#	print "\t\t$RNAi mapped to $testexpression\n";
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

open (OUT,    ">/wormsrv2/autoace/MAPPINGS/RNAi_mappings.$db_version");
open (OUTACE, ">/wormsrv2/autoace/MAPPINGS/RNAi_mappings.$db_version.ace");

# Produce connections for RNAi->Genes
foreach my $mapped (sort keys %finaloutput) {
    print OUT "$mapped\t@{$finaloutput{$mapped}}\n";
    for (my $n = 0; $n < (scalar @{$finaloutput{$mapped}}); $n++) {
        print OUTACE "RNAi : $mapped\n";
        print OUTACE "-D Predicted_gene\n\n";
    }
} 



#open autoace connection with aceperl
my $db = Ace->connect(-path=>$dbdir, -program=>$tace) || die "Couldn't connect to $dbdir\n", Ace->error;

foreach my $mapped (sort keys %finaloutput) {
  print OUT "$mapped\t@{$finaloutput{$mapped}}\n";
  for (my $n = 0; $n < (scalar @{$finaloutput{$mapped}}); $n++) {
    # Does this sequence have a locus?
    my $seq = $db->fetch(-class=>'Sequence',-name=>" $finaloutput{$mapped}->[$n]");
    my $locus;
    if(defined($seq->at('Visible.Locus_genomic_seq'))){
      ($locus) = ($seq->get('Locus_genomic_seq'));
    }

    print OUTACE "RNAi : \"$mapped\"\n";
    print OUTACE "Predicted_gene \"$finaloutput{$mapped}->[$n]\"\n";
    print OUTACE "Locus \"$locus\"\n" if ($locus);
    print OUTACE "\n";
  }
} 

$db->close;


# Produce connections for RNAi->Expr_profile
foreach my $mapped (sort keys %finaloutput2) {
    print OUT "$mapped\t@{$finaloutput2{$mapped}}\n";
    for (my $n = 0; $n < (scalar @{$finaloutput2{$mapped}}); $n++) {
      my ($expr_profile) = (@{$finaloutput2{$mapped}}->[$n] =~ /(\S+)\.\d+$/);
      print OUTACE "RNAi : $mapped\n";
      print OUTACE "Expr_profile $expr_profile\n\n";
    }
} 


print LOG "Script ended at ",`date`,"\n";
close(LOG);
close(OUT);
close(OUTACE);


##############################
# read acefiles into autoace #
##############################

my $command =<<END;
pparse /wormsrv2/autoace/MAPPINGS/RNAi_mappings.$db_version.ace
save
quit
END

open (TACE,"| $tace -tsuser map_RNAi $dbdir") || die "Couldn't open tace connection to $dbdir\n";
print TACE $command;
close (TACE);

exit(0);

###############################
# Prints help and disappears  #
###############################




##############################################################
#
# Subroutines
#
##############################################################

sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate     = `date +%y%m%d`; chomp $rundate;
  $log        = "/wormsrv2/logs/$script_name.$rundate.$$";

  open (LOG, ">$log") or die "cant open $log";
  print LOG "$script_name\n";
  print LOG "started at ",`date`,"\n";
  print LOG "=============================================\n";
  print LOG "\n";

}

##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

############################################

