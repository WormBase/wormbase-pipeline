#!/usr/local/bin/perl5.8.0 -w
#
# find_short_genes.pl  
# 
# by Keith Bradnam, aged 12 and half
#
# A script to find (and classify) potentially short, spurious genes (default = <50 aa)
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2004-03-24 09:57:04 $     


use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Getopt::Long;


##############################
# misc variables             #
##############################

my $maintainers  = "All";           # who will receive log file
my $ws_version   = &get_wormbase_version_name;
my $dir          = "/wormsrv2/tmp"; # default output location
my $file         = "short_spurious_genes.$ws_version.csv"; # initial output file


##############################
# command-line options       #
##############################

my $help;                # Help/Usage page
my $build;               # flag to say whether this is a build script or run outside of build
my $verbose;             # turn on extra output
my $debug;               # For sending output to just one person
my $database;            # which database to use?
my $cutoff;              # what length of gene should you use as the cutoff?
my $test;                # use test environment
my $html;                # Also produce html files for website?

GetOptions ("database=s" => \$database,
            "verbose"    => \$verbose,
	    "cutoff=i"   => \$cutoff,
            "debug=s"    => \$debug,
            "help"       => \$help,
	    "build"      => \$build,
	    "html"       => \$html,
	    "test"       => \$test
            );


# Help pod if needed
exec ('perldoc',$0) if ($help);


# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}


# database/file paths and locations
my $basedir     = "/wormsrv2";
$basedir        = glob("~wormpub")."/TEST_BUILD" if ($test);
 


# set default cutoff to 50 amino acids if not specified on command line
if(!defined($cutoff)){
  $cutoff = 50;        
}


# specify different path if using build mode
$dir = "$basedir/autoace/CHECKS" if ($build);



#####################################################################
#
# Open log file, output file, and database connection
#
#####################################################################

my $log = Log_files->make_build_log();  
$log->write_to("Looking for potentially spurious genes shorter or equal to $cutoff amino acids in length\n");
print "Looking for potentially spurious genes shorter or equal to $cutoff amino acids in length\n" if ($verbose);

$log->write_to("Output file(s) will be written to $dir\n");


# open initial output file
open(OUT,">$dir/$file") || die "Can't open output file\n";


# open a database connection and grab list of genes
my $db = Ace->connect(-path  =>  "$database")  || die "Cannot open $database",Ace->error;
my @genes = $db->fetch(-class => "elegans_CDS");


# used for testing purposes
my $counter = 0;


################################################
#
# M A I N   L O O P  -  loop through each gene
#
################################################

foreach my $gene (@genes){

  # when debugging can reduce this if you want less output to speed up script
  last if ($counter > 50000);  $counter++;
  
  # get protein length information from translation of gene
  my $peptide = $gene->asPeptide();

  # trim FASTA header and remove new lines, and grab length
  $peptide =~ /^>\S+\s+([\w\*].*)/s;
  my $peptide_seq = $1 ; 
  $peptide_seq =~ s/\s//g;
  my $length = length($peptide);

  # ignore proteins longer than cutoff value
  next if ($length > $cutoff);


  # Get confirmed/partially_confirmed etc. status
  my $status;
  if   ($gene->at("Properties.Coding.Confirmed_by")){$status = "Confirmed";}	
  elsif($gene->at("Visible.Matching_cDNA"))         {$status = "Partially_confirmed";}
  else                                              {$status = "Predicted";}			
  

  # get lab
  my $lab = $gene->From_laboratory;


  # get RNAi info, status defaults to 'N/A' when there are no RNAi results
  # otherwise status is set to WT or non-WT
  my $rnai_result = "N/A";
  my @RNAi = $gene->RNAi_result;
  foreach my $item (@RNAi){
    $rnai_result = "WT";
    my $RNAi = $db->fetch(RNAi => "$item");
    my $phenotype = $RNAi->Phenotype;
    if ($phenotype ne "WT"){
      $rnai_result = "non-WT";
      $RNAi->DESTROY();
      last;
    }
    $RNAi->DESTROY();
  }

  
  # look for associated PCR products that do or do not amplify
  # status is 'N/A' if there are no associated PCR products
  my $pcr_result = "N/A";
  my @PCRs = $gene->Corresponding_PCR_product;
  foreach my $item (@PCRs){
    $pcr_result = "Does not amplify";
    my $pcr = $db->fetch(PCR_product => "$item");
    my $amplify_status = $pcr->Amplified;
    if ((defined($amplify_status)) && ($amplify_status == 1)){
      $pcr_result = "Amplifies";
      $pcr->DESTROY();
      last;
    }
    $pcr->DESTROY();
  }



  # flag output status depending on how bad we think the prediction is...
  # 1 - no cDNA evidence AND two types of evidence AGAINST gene prediction (RNAi AND PCR)
  # 2 - no cDNA evidence AND one type of evidence AGAINST gene prediction (RNAi OR PCR) and
  #     other type of evidence not available
  # 3 - no cDNA evidence AND no available RNAi AND PCR information
  # 4 - no cDNA evidence AND contradictory RNAi AND PCR information (wild type AND amplifies, OR 
  #     non-wild type AND doesn't amplify)
  # 5 - no cDNA evidence BUT both RNAi and PCR information confirm it is a real gene
  # 6 - cDNA evidence

  my $evidence;

  if($status eq "Predicted"){
    if(($rnai_result eq "WT") && ($pcr_result eq "Does not amplify")){
      $evidence = "1";
    }
    elsif((($rnai_result eq "WT")  && ($pcr_result eq "N/A")) ||
	  (($rnai_result eq "N/A") && ($pcr_result eq "Does not amplify"))){
      $evidence = "2";
    }
    elsif(($rnai_result eq "N/A") && ($pcr_result eq "N/A")){
      $evidence = "3";
    }
    elsif((($rnai_result eq "WT")     && ($pcr_result eq "Amplifies")) || 
	  (($rnai_result eq "non-WT") && ($pcr_result eq "Does not amplify"))){
      $evidence = "4";
    }
    else{ 
      $evidence = "5";
    }
  }
  else{
    $evidence = "6";
  }
  
  print OUT "EVIDENCE $evidence,$lab,$gene,$length,$status,$rnai_result,$pcr_result\n";
  print     "EVIDENCE $evidence,$lab,$gene,$length,$status,$rnai_result,$pcr_result\n" if ($verbose);



  # kill AcePerl objects
  $gene->DESTROY();
}

close(OUT);
$db->close;





############################################################
#
# Process main output file into files for web (if required)
#
############################################################

if($html){
  
  $log->write_to("Producing html file\n");
  print "Producing html file\n" if ($verbose);

  # sort output file by various factors using sort command
  my $status = system("sort -t, -k 2,2 -k 1,1 -k 5,5r -k 4n,4 $dir/$file -o $dir/$file.tmp");
  $log->write_to("ERROR: Couldn't run sort command") if ($status != 0);
  
  open(IN,"<$dir/$file.tmp")      || die "Couldn't open $dir/$file.tmp for reading";
  open(OUT, ">$dir/$file.html")   || die "Couldn't write to $dir/$file.html";
  
  # set once to be able to flag when you reach RW part of output file
  my $lab_flag = "HX";

  print OUT "<TR bgcolor=lightgrey>\n";
  print OUT "<TD colspan=6 align=center><B> Sanger genes </B></TD>\n";
  print OUT "</TR>\n";

  while(<IN>){
    my @line = split(/,/);
    my $evidence = $line[0];
    my $lab      = $line[1];
    my $gene     = $line[2];
    my $length   = $line[3];
    my $status   = $line[4];
    my $rnai     = $line[5];
    my $pcr      = $line[6];

    # drop evidence prefix for web files, just show number
    $evidence =~ s/EVIDENCE //;


    # print web output
    if(($lab_flag eq "HX") && ($lab eq "RW")){
      # now set flag to RW so this only get executed once
      $lab_flag = "RW";
      print OUT "<TR bgcolor=lightgrey>\n";
      print OUT "<TD colspan=6 align=center><B> St. Louis genes </B></TD>\n";
      print OUT "</TR>\n";
    }

    print OUT "<TR bgcolor=lightblue>\n";
    print OUT "<TD align=center><B>$evidence</B></TD>\n";
    print OUT "<TD align=center><a href=\"http://dev.wormbase.org/db/gene/gene?name=${gene};class=Sequence\">$gene</a></TD>\n";
    print OUT "<TD align=center>$length</TD>\n";
    print OUT "<TD align=center>$status</TD>\n";
    print OUT "<TD align=center>$rnai</TD>\n";    
    print OUT "<TD align=center>$pcr</TD>\n";
    print OUT "</TR>\n";

  }
  close(IN);
  close(OUT);

  # remove temp file
  $status = system("rm $dir/$file.tmp");
  $log->write_to("ERROR: Couldn't remove temp file") if ($status != 0);

}




################################
# Tidy up and exit             #
################################

# only mail if running as part of build or if debugging
$log->mail("$maintainers") if ($build || $debug);
exit;





__END__

=pod

=head2   NAME - find_short_genes.pl

=head1 USAGE

=over 4

=item find_short_genes.pl -[options]

=back

=head1 DESCRIPTION

This script will check the 'elegans_CDS' subclass in any target database
and find all genes less than a certain size (defaults to 50 amino acids).
It will then check each of those genes for the presence/absence of cDNA
information and also whether there are associated RNAi experiments that give
a wild-type phenotype.  Finally, it also looks to see if there are associated
PCR products that don't amplify.  The various combinations of these different
types of evidence allows the script to categorise and prioritise the genes,
such that curators can inspect these short genes to see if they are real.

Can be run as part of build or not.  If so it will write the output file to
/wormsrv2/autoace/CHECKS/ else it will write to /wormsrv2/tmp

=back

=head1 MANDATORY arguments: -database

=over 4

=item -database <path to valid acedb database>

Database will be expected to have a valid 'elegans_CDS' subclass

=back

=head1 OPTIONAL arguments: -cutoff, -build, -debug, -verbose, -help, -html

=over 4

=item -cutoff <integer>

Specify the length (in amino acids) of genes to consider by this script.  Will
look at every gene equal to or less than the specified length


=item -build

If specified will assume that this script is being run as part of the build,
the only differences being that the output file will be placed in /wormsrv2/autoace/CHECKS
rather than in /wormsrv2/tmp and the log file will be emailed to everyone


=item -html

Additionally produce a html file that can be added to the website as part of the 
data consistency checks page

=item -debug <user>

Send log report to specified user only

=item -help

This help.

=back

=item -verbose

Toggles a little extra output on the command line when running the script

=back


=item -test

Use test environment (~wormpub/TEST_BUILD) rather than /wormsrv2

=back

=head1 AUTHOR Keith Bradnam (krb@sanger.ac.uk) 

=back

=cut

