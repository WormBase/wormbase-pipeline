#!/usr/local/bin/perl5.6.1 -w
#
# map_RNAi.pl
#
# Add information to RNAi objects based on overlaps in GFF files 
#
# by Kerstin Jekosch
#
# Last updated by: $Author: dl1 $                      
# Last updated on: $Date: 2003-08-21 13:44:05 $        


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
my $test;       # Test mode
my $debug;      # Debug mode, verbose output to user running script
our $log;

GetOptions ("debug=s"   => \$debug,
	    "test"      => \$test,
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

our %genetype         = ();


foreach my $chromosome (@chromosomes) {
    my %RNAi             = ();
    my %RNAicount        = ();
    my %genes            = ();
    my %exon             = ();
    my %exoncount        = ();
    my %expression       = ();
    my %expression_count = ();

    # loop through the GFF file

    print "Loop through GFF file CHROMOSOME_${chromosome}\n" if ($debug);
    open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.gff") || die "Failed to open gff file\n\n";
    while (<GFF_in>) {
	chomp;
	s/^\#.*//;
	next unless /\S/;
	@line = split /\t/;
	
	if ($line[1] eq "RNAi") {
	    my ($name) = ($line[8] =~ /\"(\S+.+)\"$/);
	    $RNAicount{$name}++;
	    my $RNAiname = $name.".".$RNAicount{$name};
	    $RNAi{$RNAiname} = [$line[3],$line[4]];
	    print "RNAi : '$name'\n" if ($debug);
	}
	elsif ($line[1] eq "Expr_profile") {
	    my ($name) = ($line[8] =~ /\"(.*)\"$/);
	    $expression_count{$name}++;
	    my $expression_name = $name.".".$expression_count{$name};
	    $expression{$expression_name} = [$line[3],$line[4]];
	    print "Expr_profile : '$name'\n" if ($debug);
	}
	elsif (($line[2] eq "exon") && (($line[1] eq "curated")       || 
					($line[1] eq "Pseudogene")    ||
					($line[1] eq "transcript")    ||
					($line[1] eq "provisional")))  {
	    my ($name) = ($line[8] =~ /\"(\S+.+)\"/);
	    $exoncount{$name}++;
	    my $exonname = $name.".".$exoncount{$name};
	    $exon{$exonname} = [$line[3],$line[4]];
	    $genetype{$name} = "CDS"        if ($line[1] eq "curated");
	    $genetype{$name} = "Transcript" if ($line[1] eq "transcript");
	    $genetype{$name} = "Pseudogene" if ($line[1] eq "Pseudogene");
	    print "Gene : '$name' [$genetype{$name}] exon $exoncount{$name}\n" if ($debug);
	}
    }
    close GFF_in;
    
    print "Finished GFF loop\n" if ($debug);
    
    #########################   
    # make exons into genes #
    #########################
    
    print "Turn exons into genes\n" if ($debug);
    
    foreach my $name (sort keys %exoncount) {
    my $v = $exoncount{$name};
    my $w = $name.".".$v;
    $genes{$name} = [$exon{$name.".1"}->[0],$exon{$w}->[1]];
  }
  
  ###################
  # make indexlists #
  ###################
  
    print "Index exons,genes,RNAi and Expr_profiles\n" if ($debug);
    
    my @exonlist       = sort { $exon{$a}->[0]  <=> $exon{$b}->[0]  } keys %exon;
    my @genelist       = sort { $genes{$a}->[0] <=> $genes{$b}->[0] } keys %genes;
    my @RNAilist       = sort { $RNAi{$a}->[0]  <=> $RNAi{$b}->[0]             || $a cmp $b } keys %RNAi;     
    my @expressionlist = sort { $expression{$a}->[0]  <=> $expression{$b}->[0] || $a cmp $b } keys %expression;
    
  #############
  # map RNAis #
  #############
  
    print "Find overlaps for RNAi\n" if ($debug);
    
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
	    print "$RNAi mapped to $testgene\n" if ($debug);
	  }
	}
      }                
    }
  }

  ###########################
  # map Expression profiles #
  ###########################
  
    print "Find overlaps for Expr_profile\n" if ($debug);
   
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
	print "$RNAi mapped to $testexpression\n" if ($debug);
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

# remove existing connections
foreach my $mapped (sort keys %finaloutput) {
    
    print OUT "$mapped\t@{$finaloutput{$mapped}}\n";
    
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
    
    print OUT "$mapped\t@{$finaloutput{$mapped}}\n";

    my $seq;
    my $locus;
    my $gene;

    for (my $n = 0; $n < (scalar @{$finaloutput{$mapped}}); $n++) {
	
	$gene = $finaloutput{$mapped}->[$n];
	print "'$mapped' is mapped to $gene\t" if ($debug);

	print OUTACE "\nRNAi : \"$mapped\"\n";

	# Does this sequence have a locus?
	if ($genetype{$gene} eq "CDS") {
	    print OUTACE "Predicted_gene \"$gene\"\n";
	    $seq = $db->fetch(-class=>'Sequence',-name=>" $finaloutput{$mapped}->[$n]");
	    if (defined($seq->at('Visible.Locus_genomic_seq'))) {
		($locus) = ($seq->get('Locus_genomic_seq'));
		print OUTACE "Locus \"$locus\"\n";
	    }
	    print " which is a CDS\n" if ($debug);
	}
	# Does this Pseudogene have a locus?
	elsif ($genetype{$gene} eq "Pseudogene") {
	    print OUTACE "Pseudogene \"$gene\"\n";
	    $seq = $db->fetch(-class=>'Pseudogene',-name=>" $finaloutput{$mapped}->[$n]");
	    if (defined($seq->at('Genetics.Locus'))){
		($locus) = ($seq->get('Locus'));
		print OUTACE "Locus \"$locus\"\n";
	    }
	    print " which is a pseudogene\n" if ($debug);
	}
	# Does this transcript have a locus?
	elsif ($genetype{$gene} eq "Transcript") {
	    print OUTACE "Transcript \"$gene\"\n";
	    $seq = $db->fetch(-class=>'Transcript',-name=>" $finaloutput{$mapped}->[$n]");
	    if (defined($seq->at('Visible.Locus'))){
		($locus) = ($seq->get('Locus'));
		print OUTACE "Locus \"$locus\"\n";
	    }
	    print " which is a transcript\n" if ($debug);
	}
    }
} 

$db->close;

print OUTACE "\n\n//Expression profiles\n";

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


#########################################################
# read acefiles into autoace (unless running test mode) #
#########################################################

unless ($test) {

    my $command = "pparse /wormsrv2/autoace/MAPPINGS/RNAi_mappings.$db_version.ace\nsave\nquit\n";

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
connections and make new ones. It wil check the current database and make Locus
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
