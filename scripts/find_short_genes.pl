#!/usr/local/bin/perl5.8.0 -w
#
# find_short_genes.pl  
# 
# by Keith Bradnam, aged 12 and half
#
# A script to find (and classify) potentially short, spurious genes (<100 aa)
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2003-10-29 13:20:03 $     


use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Ace;
use Getopt::Long;


##############################
# misc variables             #
##############################

my $maintainers = "All";  # who will receive log file
my $ws_version   = &get_wormbase_version_name;



##############################
# command-line options       #
##############################

my $help;                # Help/Usage page
my $build;               # flag to say whether this is a build script or run outside of build
my $verbose;             # turn on extra output
my $debug;               # For sending output to just one person
my $database;            # which database to use?
my $cutoff;              # what length of gene should you use as the cutoff?

GetOptions ("database=s" => \$database,
            "verbose"    => \$verbose,
	    "cutoff=i"   => \$cutoff,
            "debug=s"    => \$debug,
            "help"       => \$help,
	    "build"      => \$build
            );


# Help pod if needed
exec ('perldoc',$0) if ($help);


# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

# set default cutoff to 100 amino acids if not specified on command line
if(!defined($cutoff)){
  $cutoff = 100;        
}


#####################################################################
#
# Open log file, output file, and database connection
#
#####################################################################

my $log = Log_files->make_build_log();  
$log->write_to(&runtime,": starting script\n===============\n\n");
$log->write_to("Looking for spurious genes shorter or equal to $cutoff amino acids in length\n");
print "Looking for spurious genes shorter or equal to $cutoff amino acids in length\n" if ($verbose);


# choice of output file location depends on if you are running this as part of build
if($build){
  open(OUT,">/wormsrv2/autoace/CHECKS/short_spurious_genes.$ws_version.csv") || die "Can't open output file\n";
}
else{
  open(OUT,">/wormsrv2/tmp/short_spurious_genes.$ws_version.csv") || die "Can't open output file\n";
}


# open a database connection and grab list of genes
my $db = Ace->connect(-path  =>  "$database")  || die "Cannot open $database",Ace->error;
my @genes = $db->fetch(-class => "Predicted_gene");




################################################
#
# M A I N   L O O P  -  loop through each gene
#
################################################

foreach my $gene (@genes){
  
  # get protein information, and make sure that there is a protein object!
  my ($protein) = $gene->at("Visible.Corresponding_protein");
  next if !($protein);  
  $protein = $db->fetch(Protein => "$protein") || die "Cannot fetch protein\n";
  my $length = $protein->at("Peptide")->right(2);

  # ignore proteins longer than cutoff value
  if ($length > $cutoff){
    $protein->DESTROY();
    next;
  }


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
  $protein->DESTROY();
}


################################
# Tidy up and exit             #
################################

close(OUT);
$db->close;
$log->write_to(&runtime,": finishing script\n");
# only mail if running as part of build
$log->mail("$maintainers") if ($build);
exit;





__END__

=pod

=head2   NAME - find_short_genes.pl

=head1 USAGE

=over 4

=item find_short_genes.pl -[options]

=back

=head1 DESCRIPTION

This script will check the 'Predicted_gene' subclass in any target database
and find all genes less than a certain size (defaults to 100 amino acids).
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

Database will be expected to have a valid 'Predicted_gene' subclass

=back

=head1 OPTIONAL arguments: -cutoff, -build, -debug, -verbose, -help

=over 4

=item -cutoff <integer>

Specify the length (in amino acids) of genes to consider by this script.  Will
look at every gene equal to or less than the specified length


=item -build

If specified will assume that this script is being run as part of the build,
the only differences being that the output file will be placed in /wormsrv2/autoace/CHECKS
rather than in /wormsrv2/tmp and the log file will be emailed to everyone

=item -debug <user>

Send log report to specified user only

=item -help

This help.

=back

=item -verbose

Toggles a little extra output on the command line when running the script

=back


=head1 AUTHOR Keith Bradnam (krb@sanger.ac.uk) 

=back

=cut

