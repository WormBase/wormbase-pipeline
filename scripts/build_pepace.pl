#!/usr/local/bin/perl5.6.0 -w                   
# build_pepace.pl                           
# 
# by Anthony Rogers
#Ge
# This creates an acefile that can be loaded in to an empty 
# database to completely recreate what was pepace. This is based
# solely in the wormpep.history file.
# 
#
# Last updated by: $Author: ar2 $                     
# Last updated on: $Date: 2002-07-24 10:11:59 $     


use strict;                                     
use lib "/wormsrv2/scripts/";                  
use Wormbase;

##############
# variables  #                                                                   
##############

# Produce a log file that is a) emailed to us all 
# and b) copied to /wormsrv2/logs

my $maintainers = "All";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
our $log        = "/wormsrv2/logs/build_pepace.$rundate";

my $ver = 80;#&get_wormbase_version();
my $wormpepdir = "/wormsrv2/WORMPEP/wormpep$ver";

open( LOG, ">$log") || die "cant open $log";
print LOG "build_pepace log file $rundate $runtime using wormpep$ver and perl version $]\n---------------------------------------------------\n\n";


open (HISTORY, "$wormpepdir/wormpep.history$ver") || die "wormpep.history$ver";

#read file in 
our ($gene, $CE, $in, $out);
my %CE_history;
my %CE_gene;
my %gene_CE;
my %CE_live;

my $stem;
my $isoform;
my $existingCE;
my $existingGene;

my $count;
while(<HISTORY>)
      {
	my @data = split(/\s+/,$_);
	($gene, $CE, $in, $out) = @data;


	if( defined( $CE_gene{$CE} ) )
	  {
	    $existingGene = $CE_gene{$CE};
	    if( "$gene" eq "$existingGene" )
	      {
		if( $CE_live{$CE} == 1 ) {
		    print LOG "$CE is being replaced by $gene when not dead\n";
		  }
		else {
		  #reappeared protein
		  &reappearedPeptide;
		}

	      }
	    #is this an isoform of a pre-exisiting gene?
	    elsif( $gene =~ m/((\w+\.\d+)\w*)/)
	      {
		$stem = $2;     #eg FK177.8
		$isoform = $1;  #eg FK177.8a
		
		if( $existingGene =~ m/^($stem)\w*/ )
		  {
		    #$gene is isoform
		    if( $CE_live{$CE} == 1 )
		      {
			#Became isoform
			# ZK177.8  CE02097 8 
			# ZK177.8a CE02097 11
			&becameIsoform; 
		      }
		    else
		      {
			#Reappeared as isoform to 
			# ZK177.8  CE02097 8 11
			# ZK177.8a CE02097 12
			&reappearedAsIsoform;
		      }		    
		  }
		else
		  {
		    #peptide coded by multiple genes
		    print LOG "$CE coded by multiple genes\n"
		  }		
	      }
	  }
	#processing entry where CE not known
	else
	  {
	    if( defined( $gene_CE{$gene}) )
	      {
		# if peptide coded by >1 genes
		# do something clever
		# else
		$existingCE = $gene_CE{$gene};
		&addNewPeptide;
		&replacePeptide;
	      }
	  
	    elsif( $gene =~ m/((\w+\.\d+)\w*)/ ) 
	      {
		$stem = $2;     #eg FK177.8
		$isoform = $1;  #eg FK177.8a
		
		if( defined( $gene_CE{$stem}) )
		  {
		    #$gene is isoform of already entered gene
		    &addNewPeptide;
		    $existingCE = $gene_CE{$stem};
		    #add to history of $CE that it is isoform og $existingCE
		  }
		else{
		  &addNewPeptide;
		}
	      }
	  }
	$count++;
	
#	if( $count == 3000 ){
#	  last;
#      }
      }
close HISTORY;

foreach my $key(sort keys %CE_history){
  print "$key\t$CE_history{$key}\n";
}






close LOG;
#### use Wormbase.pl to mail Log ###########
my $name = "$0";
$maintainers = "ar2\@sanger.ac.uk";
#&mail_maintainer ($name,$maintainers,$log);
#########################################
exit(0);


sub addNewPeptide
  {
    # 1st occurance of peptide
    $CE_history{$CE} .= "Created $in\t" ;
    $CE_gene{$CE} = $gene;
    $CE_live{$CE} = 1;   #assume live when put in unless explicitly killed
    $gene_CE{$gene} = $CE;
    if (defined ($out) ){
      $CE_history{$CE} .= "Removed $out\t" ;
      $CE_live{$CE} = 0;
    } 
  }

sub replacePeptide
  {
    $CE_live{$existingCE} = 0;
    $CE_history{$existingCE} .= "$in Replaced_by $CE\t";
    $CE_history{$CE} .= "$in Replaces $existingCE\t";
  }

sub reappearedPeptide
  {
    $CE_live{$CE} = 1;
    $CE_history{$CE} .= "$in Reappeared\t";
    if( $out )
      {
	$CE_live{$CE} = 0;
	$CE_history{$CE} .= "$in Removed\t";
      }
    #gene is same as was previously if this routine called
  }

sub reappearedAsIsoform
  {
    $CE_live{$CE} = 1;
    $CE_history{$CE} .= "$in Reappeared as isoform to $stem \t"; 
    $CE_gene{$CE} = $gene;
    $gene_CE{$CE} = $CE;
    if( $out )
      {
	$CE_live{$CE} = 0;
	$CE_history{$CE} .= "$in Removed\t";
      }
  }

sub becameIsoform
  {
    $CE_live{$CE} = 1;
    $CE_history{$CE} .= "$in became isoform to $stem \t"; 
    $CE_gene{$CE} = $gene;
    $gene_CE{$CE} = $CE;
    if( $out )
      {
	$CE_live{$CE} = 0;
	$CE_history{$CE} .= "$in Removed\t";
      }   
  }




# Add perl documentation in POD format
# This should expand on your brief description above and add details of any options
# that can be used with the program.  Such documentation can be viewed using the perldoc
# command.


__END__

=pod

=head2 NAME - script_template.pl

=head1 USAGE

=over 4

=item script_template.pl  [-options]

=back

This script:

 1)
 2)
 3)
 4)
 5)

script_template.pl MANDATORY arguments:

=over 4

=item none

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
