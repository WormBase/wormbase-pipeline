#!/usr/local/bin/perl5.8.0 -w
#
# map_Oligo_set
#
# Cronjob integrity check controls for generic ACEDB database.
#
# Usage: map_Oligo_set.pl [-options]
#
# 010927 :  dl : modified output print line to include speech marks. this prevents any acefile
#                parsing problems
# 010717 :  kj : modified to (hopefully) run faster and produce the final output files in one go
#
#############################################################################################


#############
# variables #
#############

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Cwd;

##########################
# Script variables (run) #
##########################

my $maintainers = "All";
my $rundate     = &rundate;
my $runtime     = &runtime;
my %output      = ();
my %finaloutput = ();

########################
# command-line options #
########################

my $help;       # Help perdoc
my $test;       # Test mode
my $debug;      # Debug mode, output only goes to one user
my $verbose;    # verbose mode, more command line outout

GetOptions ("debug=s"   => \$debug,
	    "verbose"   => \$verbose,
	    "test"      => \$test,
            "help"      => \$help
	    );

# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}


#############
# Paths etc #
#############

my $tace        = &tace;                                     # tace executable path
my $dbdir       = "/wormsrv2/autoace";                       # Database path
my $gffdir      = "/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS"; # GFF splits directory
my @chromosomes = qw( I II III IV V X );                     # chromosomes to parse 
my $db_version  = &get_wormbase_version_name;                # WS version number
our %genetype   = ();                                        # gene type hash

################
# Open logfile #
################
my $log = Log_files->make_build_log();


#####################################################
# get exons and Oligo_set info out of the gff files #
#####################################################        
   
foreach my $chromosome (@chromosomes) {
  $log->write_to("Processing chromosome $chromosome\n");
  print "\nProcessing chromosome $chromosome\n" if ($verbose);
  my %exoncount  = ();
  my %Oligocount = ();
  my %genes      = ();
  my %Oligo      = ();
  my %exon       = ();
  my @f          = ();
  
  # Get Oligo set info from split GFF file
  open (GFF, "<$gffdir/CHROMOSOME_${chromosome}.Oligo_set.gff") || die "Failed to open Oligo_set gff file\n\n";
  while (<GFF>) {
    chomp;	     
    s/\#.*//;
    next unless /\S/;
    @f = split /\t/;
    
    my ($name) = ($f[8] =~ /Oligo_set \"(.*)\"$/);
    unless ($name) {
      $log->write_to("ERROR: Cant get name from $f[8]\n");
      next;
    }
    $Oligocount{$name}++;
    my $Oligoname = $name.".".$Oligocount{$name};
    $Oligo{$Oligoname} = [$f[3],$f[4]];
    print "Oligo_set : '$name'\n" if ($verbose);
  }
  close(GFF);


  # Get exon info from split exon GFF files
  open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.exon.gff") || die "Failed to open exon gff file\n\n";
  while (<GFF_in>) {
    chomp;	     
    s/\#.*//;
    next unless /\S/;
    @f = split /\t/;

    my ($name) = ($f[8] =~ /\"(\S+)\"/);
    $exoncount{$name}++;
    my $exonname = $name.".".$exoncount{$name};
    $exon{$exonname} = [$f[3],$f[4]];
    $genetype{$name} = "CDS"        if ($f[1] eq "curated");    
    print "Gene : '$name' [$genetype{$name}] exon $exoncount{$name}\n" if ($verbose);
  }
  close(GFF_in);


  # Get exon info from split pseudogene exon GFF files
  open (GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.exon_pseudogene.gff") || die "Failed to open exon_pseudogene gff file\n\n";
  while (<GFF_in>) {
    chomp;	     
    s/\#.*//;
    next unless /\S/;
    @f = split /\t/;

    my ($name) = ($f[8] =~ /\"(\S+)\"/);
    $exoncount{$name}++;
    my $exonname = $name.".".$exoncount{$name};
    $exon{$exonname} = [$f[3],$f[4]];
    $genetype{$name} = "Pseudogene" if ($f[1] eq "Pseudogene");
    print "Gene : '$name' [$genetype{$name}] exon $exoncount{$name}\n" if ($verbose);    
  }
  close(GFF_in);


  # Get exon info from split transcript exon GFF file
  open (GFF, "<$gffdir/CHROMOSOME_${chromosome}.exon_noncoding.gff") || die "Failed to open exon_noncoding gff file\n\n";
  while (<GFF>) {
    chomp;	     
    s/\#.*//;
    next unless /\S/;
    @f = split /\t/;

    my ($name) = ($f[8] =~ /\"(\S+)\"/);
    $exoncount{$name}++;
    my $exonname = $name.".".$exoncount{$name};
    $exon{$exonname} = [$f[3],$f[4]];
    $genetype{$name} = "Transcript" if ($f[1] eq "Non_coding_transcript");
    print "Gene : '$name' [$genetype{$name}] exon $exoncount{$name}\n" if ($verbose);
  }
  close(GFF);
  

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
  
  print "Index exons, genes, and oligos\n" if ($verbose);
  
  my @exonlist  = sort { $exon{$a}->[0]  <=> $exon{$b}->[0]  } keys %exon;
  my @genelist  = sort { $genes{$a}->[0] <=> $genes{$b}->[0] } keys %genes;
  my @Oligolist = sort { $Oligo{$a}->[0] <=> $Oligo{$b}->[0] || $a cmp $b } keys %Oligo;
  
  ##########
  # map it #
  ##########
  
  print "Find overlaps for Oligo_set\n" if ($verbose);
  
  my $lastfail = 0;
  
  for (my $x = 0; $x < @Oligolist; $x++) {
    my $testOligo   = $Oligolist[$x];
    my $Oligostart  = $Oligo{$testOligo}->[0];
    my $Oligostop   = $Oligo{$testOligo}->[1];
    
    for (my $y = 0; $y < @genelist; $y++) {
      my $testgene = $genelist[$y];
      my $genestart= $genes{$testgene}->[0];
      my $genestop = $genes{$testgene}->[1];
      
      if ($Oligostart > $genestop) {
	$lastfail = $y;
	next; 
      }
      
      elsif ($Oligostop < $genestart) {
	last; 
      }
      
      else {
	for (my $z = 1; $z <= $exoncount{$testgene}; $z++) {
	  my $exon_start = $exon{"$testgene.$z"}->[0];
	  my $exon_stop  = $exon{"$testgene.$z"}->[1];
	  
	  if ( not (($Oligostart > $exon_stop) || ($Oligostop < $exon_start))) {
	    my ($Oligo) = ($testOligo =~ /(\S+)\.\d+$/);
	    push @{$output{$Oligo}}, $testgene;
	  }
	}
      }                
    }
  }
}


###################
# sort the output #
###################
 
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

########################
# produce output files #
########################

open (OUTACE, ">/wormsrv2/autoace/acefiles/Oligo_mappings.ace");

foreach my $mapped (sort keys %finaloutput) {
    
  print "$mapped\t@{$finaloutput{$mapped}}\n" if ($verbose);
  
  for (my $n = 0; $n < (scalar @{$finaloutput{$mapped}}); $n++) {
    my $gene = $finaloutput{$mapped}->[$n];
    
    if ($genetype{$gene} eq "CDS") {
      print OUTACE "CDS : \"$gene\"\n";
      print OUTACE "Corresponding_oligo_set \"$mapped\"\n\n";
    }
    elsif ($genetype{$gene} eq "Pseudogene") {
      print OUTACE "Pseudogene : \"$gene\"\n";
      print OUTACE "Corresponding_oligo_set \"$mapped\"\n\n";
    }
    elsif ($genetype{$gene} eq "Transcript") {
      print OUTACE "Transcript : \"$gene\"\n";
      print OUTACE "Corresponding_oligo_set \"$mapped\"\n\n";
    }
  }
} 
close(OUTACE);

##############################
# read acefiles into autoace #
##############################

unless ($test) {
  
  my $command = "pparse /wormsrv2/autoace/acefiles/Oligo_mappings.ace\nsave\nquit\n";
  
  open (TACE,"| $tace -tsuser map_Oligo_products $dbdir") || die "Couldn't open tace connection to $dbdir\n";
  print TACE $command;
  close (TACE);
}

###############
# hasta luego #
###############

$log->mail("$maintainers","BUILD REPORT: $0");

exit(0);

##############################################################
##### Subroutines #####
##############################################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}
__END__

=pod

=head2 NAME - map_Oligo_products.pl

=head1 USAGE

=over 4

=item map_Oligo_products.pl [-options]

=back

map_Oligo_products.pl calculates the overlap between the mapped Oligo_products
and the CDS, transcript and pseudogene coordinates in the WS database release. 
It will generate an acefile for this data and upload into /wormsrv2/autoace

map_Oligo_products mandatory arguments:

=over 4

=item none

=back

map_Oligo_products optional arguments:

=over 4

=item -debug, debug mode email goes to user specified by -debug

=item -verbose, toggles extra output

=item -test, Test mode, generate the acefile but do not upload 

=item -help, Help pages

=back

=cut
