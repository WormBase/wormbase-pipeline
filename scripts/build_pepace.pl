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
# Last updated on: $Date: 2002-07-25 13:10:53 $     


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
my %multicodedPeps;

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
		    if( $CE_live{$CE} == 1 )
		      {
			#peptide coded by multiple genes
			if( defined( $out ) ) {
			  print LOG "$CE temporarily mulitcoded $in - $out\n";
			}
			else {
			  &addMultiCoded;
			}
		      }
		    else
		      {
			#peptide was previously coded by a different gene
			&changePepGene;
		      }	
	
		  }			
	      }
	    else
	      {
		if ( &oldStyleName($gene) )
		  {
		    if( $CE_live{$CE} == 1 )
		    {
		      #peptide coded by multiple genes
		      print LOG "$CE strange oldstyle name stuff $gene\n"
		    }
		    else
		      {
			#peptide was previously coded by a different gene
			&changePepGene;
		      } 
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
	  
	    else
	      {
		if( $gene =~ m/((\w+\.\d+)\w*)/ ) 
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
		else 
		  {
		    if( $gene =~ m/\w+\.\w+/ )#one of the pain in the arse old named genes eg Y53F4B.AA orY53F4B.A - both exist! 
		      {
			&addNewPeptide; #just add it
		      }
		  }			
	      }
	  }
	$count++;
	
#	if( $count == 3000 ){
#	  last;
#      }
      }
close HISTORY;
my $release;
foreach my $key(sort keys %CE_history)
{
  print "$key";
  foreach $release(sort byRelease keys %{ $CE_history{$key} })
    {
      print "\t$release $CE_history{$key}{$release}\n";
    } 
  if( $CE_live{$key} == 1 ){
    print "\tLive \n\n";
  }
  else{
    print "\n";
  }
}


###check multicodedPeps is true# # 
print "About to check multis\n\n\n\n\n\n\n\n\n\n\n\n";
my $db =Ace->connect('/wormsrv2/current_DB') || die "cant connect to current_DB\n\n";
print "connected\n";
foreach my $multigene( keys %multicodedPeps ){
  my $testgene = "WP:"."$multigene";
  print $testgene;
  my $DBpep = $db->fetch(Protein => "$testgene");
  if( defined($DBpep) )
    {
      my @corresponding_DNA = $DBpep->at('Visible.Corresponding_DNA');
      my $corresponding_DNA_COUNT = @corresponding_DNA;
      if ($corresponding_DNA_COUNT > 1){
	print LOG "$multigene correctly id'd as multicoded ($corresponding_DNA_COUNT X)\n";
      }
      else{
	print LOG "$multigene has $corresponding_DNA_COUNT coresponding_DNA's\n";
      }
    }
}


$db->close;
close LOG;
#### use Wormbase.pl to mail Log ###########
my $name = "$0";
$maintainers = "ar2\@sanger.ac.uk";
#&mail_maintainer ($name,$maintainers,$log);
#########################################
exit(0);

sub byRelease
  {
    #used by sort 
    $a <=> $b
  }

sub addNewPeptide
  {
    # 1st occurance of peptide
    $CE_history{$CE}{$in} = "created"; #.= "Created $in\t" ;
    $CE_gene{$CE} = $gene;
    $CE_live{$CE} = 1;   #assume live when put in unless explicitly killed
    $gene_CE{$gene} = $CE;
    if (defined ($out) ){
      $CE_history{$CE}{$out} = "removed";#.= "Removed $out\t" ;
      $CE_live{$CE} = 0;
    } 
  }

sub replacePeptide
  {
    $CE_live{$existingCE} = 0;
    $CE_history{$existingCE}{$in} = "replaced_by $CE";# .= "$in Replaced_by $CE\t";
    $CE_history{$CE}{$in} = "Created to replace $existingCE";#.= "$in Replaces $existingCE\t";
  }

sub reappearedPeptide
  {
    $CE_live{$CE} = 1;
    $CE_history{$CE}{$in} = "reappeared" ;#.= "$in Reappeared\t";
    if( $out )
      {
	$CE_live{$CE} = 0;
	$CE_history{$CE}{$out} = "removed";# .= "$out Removed\t";
      }
    #gene is same as was previously if this routine called
  }

sub reappearedAsIsoform
  {
    $CE_live{$CE} = 1;
    #check if becoming isoform is same release as removal - if so modify history to show conversion rather than reappearance
    if( (defined("$CE_history{$CE}{$in}") ) && ("$CE_history{$CE}{$in}" eq "removed") ) { 
      $CE_history{$CE}{$in} = "converted to isoforom of $stem";
    }
    else{
      $CE_history{$CE}{$in} = "reappeared as isoform to $stem";# .= "$in Reappeared as isoform to $stem \t"; 
    }
    
    $CE_gene{$CE} = $gene;
    $gene_CE{$CE} = $CE;
    if( $out )
      {
	$CE_live{$CE} = 0;
	$CE_history{$CE}{$out} = "removed";# .= "$out Removed\t";
      }
  }

sub becameIsoform
  {
    $CE_live{$CE} = 1;
    $CE_history{$CE}{$in} = "became isoform to $stem";#  .= "$in became isoform to $stem \t"; 
    $CE_gene{$CE} = $gene;
    $gene_CE{$CE} = $CE;
    if( $out )
      {
	$CE_live{$CE} = 0;
	$CE_history{$CE}{$out} = "removed";# .= "$out Removed\t";
      }   
  }

sub changePepGene
  {
    $CE_live{$CE} = 1;
    $CE_gene{$CE} = $gene;
    $CE_history{$CE}{$in} = "reappeared coded by another gene";# .= "$in reappeared coded by another gene\t";
    if( $out )
      {
	$CE_live{$CE} = 0;
	$CE_history{$CE}{$out} = "removed";# .= "$out Removed\t";
      } 
  }
 sub oldStyleName
  {
    if( $gene =~ m/\w+\.\p{IsAlpha}+/ ) {
      return 1;
    }
    else {
      return 0;
    }
  }

sub addMultiCoded
  {
    $multicodedPeps{$CE}++;
    $CE_gene{$CE} .= " $gene";
    $gene_CE{$gene} = $CE;
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
