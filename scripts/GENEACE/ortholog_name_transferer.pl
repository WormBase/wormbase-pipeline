#!/software/bin/perl -w


use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
use NameDB_handler;


######################################
# variables and command-line options #
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $user);
my $database;
my $sanitycheck;
my $explain_gene;
my $explain_gene_output='';

my $acefile = "/nfs/wormpub/DATABASES/geneace/CGC/orthos.ace";
my $batchfile  = "/nfs/wormpub/DATABASES/geneace/CGC/batch_load";
my $deletebatchfile  = "/nfs/wormpub/DATABASES/geneace/CGC/delete_batch_load";

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "test"       => \$test,
            "verbose"    => \$verbose,
            "store:s"    => \$store,
	    "database:s" => \$database,
	    "ace:s"      => \$acefile, # ace file to be loaded into geneace
	    "deletebatch:s"    => \$deletebatchfile, # list of genes with names to be deleted using batch_delName.pl
	    "batch:s"    => \$batchfile, # list of genes with names to be added using batch_pname_update.pl
	    "user:s"     => \$user, #manditory option as required for .ace output.
	    "sanitycheck" => \$sanitycheck, # check everything is OK
	    "explain_gene:s" => \$explain_gene, # provide a detailed explanation of what has been done to a GeneID
            );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} 
else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			   );
}

my $log = Log_files->make_build_log($wormbase);


my %CGC_species = ('briggsae' => 'Cbr',
		   'remanei'  => 'Cre',
		   'japonica' => 'Cjp',
		   'brenneri' => 'Cbn',
		   'pacificus'=> 'Ppa',
		   'malayi'   => 'Bma', # NB not 'brugia'
                   'volvulus' => 'Ovo', # NB not 'ovolvulus'
                   'ratti'    => 'Sra', # NB not 'tmuris'
		   'muris'    => 'Tmu', # NB not 'sratti'
		  );

$database = $database || $wormbase->database('geneace');
$log->write_to("Database : $database\n\n");

my $ace;
my $namedb;
my $deletenamedb;
open($ace, ">$acefile")      || $log->log_and_die("cant write acefile - $acefile : $!\n");
open($namedb, ">$batchfile") || $log->log_and_die("cant write batch load file - $batchfile: $!\n");
open($deletenamedb, ">$deletebatchfile") || $log->log_and_die("cant write batch load file - $deletebatchfile: $!\n");


my %person = (
	      'pad' => 'WBPerson1983',
	      'klh' => 'WBPerson3111',
	      'gw3' => 'WBPerson4025',
	      'mh6' => 'WBPerson4055',
	     );

unless (defined $user) {
  print "\nERROR: please specify a user ID - e.g. pad\n\n[INPUT]:";
  my $tmp = <STDIN>;
  chomp $tmp;
  if (($tmp ne 'pad') && ($tmp ne 'gw3') && ($tmp ne 'mh6')){
    $log->log_and_die("UNKNOWN USER.....TERMINATE\n\n");
    print "UNKNOWN USER.....TERMINATE\n\n" if $debug;
  }
  else {
    $user = "$tmp";
  }
}


##################################################################
# Main Body
##################################################################

#Get maximum number of evidence databases. 
#Count number of matches
#If have one matching ortholog and this ortholog also has one match that is reciprocal then score as ‘1 to 1’. 
#If either side has more than one match and the matches are all reciprocal then score as ‘1 to many’ or 'many to 1'
#1 to 1 - check if already assigned ok. If so and it is different then report it as needing work. If not assigned then assign it. 
#1 to many - then construct a split name or a merged name. check if already assigned. If so and it is different then report it as needing work. If not assigned then assign it. 

if (defined $explain_gene) {$explain_gene_output = "EXPLANATION OF GENE $explain_gene\n"}

my $acedb = Ace->connect('-path' => $database) || Ace->error; # geneace

# get elegans genes with a CGC_name
$log->write_to("Fetch elegans genes\n");
my @elegans_genes = $acedb->fetch (-query => 'FIND Gene WBGene* WHERE Live AND Species = "Caenorhabditis elegans" AND Ortholog');
if (defined $explain_gene && grep /^$explain_gene$/, @elegans_genes) {$explain_gene_output .= "It is a Live elegans gene with Orthologs\n"}

# get other species genes
# Currently the policy is not to transfer names to T. muris, S. ratti and O. volvulus because MB doesn't wish it.
$log->write_to("Fetch non-elegans non-ovolvulus non-tmuris non-sratti genes\n");
my @other_genes = $acedb->fetch (-query => 'FIND Gene WBGene* WHERE Live AND Species != "Caenorhabditis elegans" AND Species != "Onchocerca volvulus" AND Species != "Trichuris muris" AND Species != "Strongyloides ratti"');
if (defined $explain_gene && grep /^$explain_gene$/, @other_genes) {$explain_gene_output .= "It is a Live non-elegans gene\n"}

# get other species DEAD genes
# We wish to be able to display a warning if we assign a CGC_name and there is another gene that is already DEAD but which still has the CGC_name - it needs to be manually removed from the DEAD gene.
$log->write_to("Fetch DEAD non-elegans non-ovolvulus non-tmuris non-sratti genes\n");
my @dead_other_genes = $acedb->fetch (-query => 'FIND Gene WBGene* WHERE Dead AND Species != "Caenorhabditis elegans" AND Species != "Onchocerca volvulus" AND Species != "Trichuris muris" AND Species != "Strongyloides ratti"');
if (defined $explain_gene && grep /^$explain_gene$/, @dead_other_genes) {$explain_gene_output .= "It is a Dead non-elegans gene\n"}

if (defined $explain_gene && $explain_gene_output eq 'EXPLANATION OF GENE $explain_gene\n') {$explain_gene_output .= "It is not a gene that has been searched for - possibly a species we are not considering\n"}



$log->write_to("Ignore genes with paper_evidence for the name\n");
my %paper_bork_list = get_paper_bork_list(@other_genes);
if (defined $explain_gene && exists $paper_bork_list{$explain_gene}) {$explain_gene_output .= "It has a paper giving evidence for its manually-assigned gene name\n"}

$log->write_to("Ignore genes with a ncRNA Biotype\n");
my %RNA_bork_list = get_RNA_bork_list(@other_genes);
if (defined $explain_gene && exists $RNA_bork_list{$explain_gene}) {$explain_gene_output .= "It is an ncRNA gene, so has a manually-assigned gene name\n"}

# store existing non-elegans gene CGC names 
# list of existing names to check that all have been set up correctly
$log->write_to("Store the existing genes with names\n");
my %existing = store_existing(@other_genes);
my %assigned_names;
my %correct_already;
if (defined $explain_gene && exists $existing{$explain_gene}) {$explain_gene_output .= "It already has the existing gene name $existing{$explain_gene}\n"}


# store existing DEAD non-elegans gene CGC names 
# list of existing DEAD names to check that all have been set up correctly
$log->write_to("Store the existing DEAD genes with names\n");
my %dead_existing = store_existing(@dead_other_genes);
if (defined $explain_gene && exists $dead_existing{$explain_gene}) {$explain_gene_output .= "It has the existing gene name $dead_existing{$explain_gene}\n"}

if ($sanitycheck) {
  $log->write_to("Sanity check\n");
  sanity_check(@other_genes);
  
} else {
  
  # go through elegans genes getting $elegans_ortholog{$elegans_gene}{$species}{$other_gene} = $number_of_evidence_databases; where this is the maximum number
  # go through other genes getting $other_ortholog{$other_gene}{$species}{$elegans_gene} = $number_of_evidence_databases; where this is the maximum number
  $log->write_to("Get elegans orthologs\n");
  my %elegans_ortholog = get_orthologs(@elegans_genes);
  $log->write_to("Get non-elegans orthologs\n");
  my %other_ortholog = get_orthologs(@other_genes);
  

  # delete unused data
  @elegans_genes = ();
  @other_genes = ();
  
  # classify as '1-to-1', 1-to-many', 'many-to-1'
  $log->write_to("Classify ortholog groups\n");
  my ($one2one, $many2one, $one2many) = classify(\%elegans_ortholog, \%other_ortholog);
  
  # delete unused data
  %elegans_ortholog = ();
  %other_ortholog = ();
  
  # get CGC names of the non-elegans gene
  #   - report ones with an existing incorrect name
  #   - output data files for names to transfer
  
  $log->write_to("Update 1:1 orthologs\n");
  update_one2one($one2one);
  $log->write_to("Update 1:many orthologs\n");
  update_one2many($one2many);
  $log->write_to("Update many:1 orthologs\n");
  update_many2one($many2one);
  
  # check to see if there are any CGC names that are set but should not be.
  $log->write_to("Check existing CGC names\n");
  check_existing(\%dead_existing, \%existing, \%assigned_names, \%correct_already);
 

}

close $ace;
close $namedb;
close $deletenamedb;

if (defined $explain_gene) {$log->write_to("\n==================================\n$explain_gene_output\n==================================\n\n")}


$log->write_to("\nCreated orthos.ace and batch_load and delete_batch_load under /nfs/wormpub/DATABASES/geneace/CGC/\n\n
Dont forget to load the ace file in to Geneace.(default is orthos.ace)\n
\nDelete names from the Nameserver before creating new ones.\n
Use batch_delName.pl -type CGC -user \$USER -pass pass -file delete_batch_load\n
Use batch_pname_update.pl -newonly -cgc -domain Gene -user \$USER -pass pass -file batch_load\n
\nFinished\n");
$log->mail;
print "\nCreated orthos.ace and batch_load and delete_batch_load under /nfs/wormpub/DATABASES/geneace/CGC/\n\n
Dont forget to load the ace file in to Geneace.(default is orthos.ace)\n
\nDelete names from the Nameserver before creating new ones.\n
Use batch_delName.pl -type CGC -user \$USER -pass pass -file delete_batch_load\n
Use batch_pname_update.pl -newonly -cgc -domain Gene -user \$USER -pass pass -file batch_load\n
\nFinished\n" if $debug;

exit;


#######################################################################################
# takes a list of gene objects and returns a hash of the 'best' orthologs in different species
# we define the "best" ortholog as the one with the most number of ortholog databases used in evidence
sub get_orthologs {
  my (@genes) = @_;

  my %store;
  foreach my $gene (@genes) {
    my $gene_name = $gene->name;
    if ($gene_name =~ /^WBGene/) {
      my @orthos = $gene->Ortholog;
      #    print $orthos[0]->name."\n";
      foreach my $orth (@orthos){
	my $ortholog_gene_name = $orth->name;
	if (!defined $orth->right) {$log->write_to( "ERROR - no species name in Gene $gene_name Ortholog $ortholog_gene_name\n"); next}
	my ($species) = $orth->right->name =~ /(\w+)$/;
	my $number = scalar $orth->right->right->col;
	# ensure that we only store the highest-scoring orthologs in any one species (where the 'score' is the number of evidence databases)
	if (exists $store{$gene_name}{$species}) {
	  my @keys = keys %{$store{$gene_name}{$species}};
	  if ($store{$gene_name}{$species}{$keys[0]} <= $number) {
	    if ($store{$gene_name}{$species}{$keys[0]} < $number) {
	      delete $store{$gene_name}{$species};
	    }
	    $store{$gene_name}{$species}{$ortholog_gene_name} = $number;
	  }
	} else {
	  $store{$gene_name}{$species}{$ortholog_gene_name} = $number;
	}
      }
      #if (! exists $store{$gene_name}) {$log->write_to( "ERROR - no data stored for $gene_name\n");}
    }
  }

  return %store;
}
#######################################################################################
# classify as '1-to-1', 1-to-many', 'many-to-1', (and we ignore any 'many-to-many')
sub classify {
  my ($elegans_ortholog, $other_ortholog) = @_;

  my (%one2one, %many2one, %one2many);


  # explain a gene stuff
  if (defined $explain_gene) {
    if (exists $elegans_ortholog->{$explain_gene}) {
      $explain_gene_output .= "It has the following orthologs:\n";
      foreach my $species (keys %{$elegans_ortholog->{$explain_gene}}) {
	$explain_gene_output .= "\t$species\n";
	foreach my $explain_ortholog (keys %{$elegans_ortholog->{$explain_gene}{$species}}) {
	  $explain_gene_output .= "\t\t$explain_ortholog\twith $elegans_ortholog->{$explain_gene}{$species}{$explain_ortholog} supporting databases\n";
	}
      }      
    } elsif (exists $other_ortholog->{$explain_gene}) {
      $explain_gene_output .= "It has the following elegans orthologs:\n";
      foreach my $species (keys %{$other_ortholog->{$explain_gene}}) {
	if ($species eq 'elegans') {
	  foreach my $explain_ortholog (keys %{$other_ortholog->{$explain_gene}{$species}}) {
	    $explain_gene_output .= "\t\t$explain_ortholog\twith $other_ortholog->{$explain_gene}{$species}{$explain_ortholog} supporting databases\n";
	    # take a look at the elegans orthologs
	    foreach my $explain_species (keys %{$elegans_ortholog->{$explain_ortholog}}) {
	      foreach my $explain_elegans_ortholog (keys %{$elegans_ortholog->{$explain_ortholog}{$explain_species}}) {
		$explain_gene_output .= "\twith orthologs back:\t$explain_species\t$explain_elegans_ortholog\twith $elegans_ortholog->{$explain_ortholog}{$explain_species}{$explain_elegans_ortholog} supporting databases\n";
		
	      }

	    } 

	  }
	}
      }      
      
    } else {
      $explain_gene_output .= "It has no orthologs, so gene names cannot be transferred using this gene.\n";
    }
  }
  

  # classify stuff
  foreach my $elegans_gene_name (keys %{$elegans_ortholog}) {
    foreach my $species (keys %{$elegans_ortholog->{$elegans_gene_name}}) {
      my $number = scalar keys %{$elegans_ortholog->{$elegans_gene_name}{$species}}; # number of orthologs to this species
      my @ortholog_gene_names  = keys %{$elegans_ortholog->{$elegans_gene_name}{$species}};
      if (defined $explain_gene && $elegans_gene_name eq $explain_gene) {$explain_gene_output .= "There are $number $species orthologs for this elegans gene.\n"}

      if ($number == 1) {
	my $ortholog_gene_name = $ortholog_gene_names[0];
	my $ortholog_number = scalar keys %{$other_ortholog->{$ortholog_gene_name}{elegans}};
	my @elegans_ids = keys %{$other_ortholog->{$ortholog_gene_name}{elegans}};
	
	if ($ortholog_number == 1) {
	  # 1-to-1
	  if ($elegans_ids[0] eq $elegans_gene_name) {
	    $one2one{$elegans_gene_name}{$species}{$ortholog_gene_name} = 1; # NB %one2one has the elegans gene ID as the first key
	    if (defined $explain_gene && ($elegans_gene_name eq $explain_gene  || $ortholog_gene_name eq $explain_gene )) {$explain_gene_output .= "There is a 1:1 orthology between elegans $elegans_gene_name and $species $ortholog_gene_name\n"}
	  }
	  
	} else {
	  # many-to-1
	  # check all the elegans orthologs of the other species gene are only orthologs of the one other species gene
	  my $all_match = 1;
	  foreach my $reverse_elegans (@elegans_ids) {
	    if (scalar keys %{$elegans_ortholog->{$elegans_gene_name}{$species}} != 1) {$all_match = 0}
	  }
	  if ($all_match) {
	    foreach my $reverse_elegans (@elegans_ids) {
	      $many2one{$ortholog_gene_name}{$species}{$reverse_elegans} = 1; # NB %many2one has the non-elegans gene ID as the first key
	      if (defined $explain_gene && ($elegans_gene_name eq $explain_gene  || $ortholog_gene_name eq $explain_gene )) {$explain_gene_output .= "There is a many:1 orthology between elegans $elegans_gene_name and $species $ortholog_gene_name\n"}
	    }
	  }
	  
	  
	}
      } else {
	# 1-to-many
	# check to see that the other species genes have a orthology back to this elegans gene
	# and check to see that the other species genes have only one orthology back to any elegans gene
	my $all_match = 1;
	foreach my $ortholog_gene_name (@ortholog_gene_names) {
	  # if no homology back or have homologies to other elegans genes then this is not a 1-to-many, it is some sort of many-to-many and we ignore those
	  if (!exists $other_ortholog->{$ortholog_gene_name}{elegans}{$elegans_gene_name}) {$all_match = 0}
	  if (defined $explain_gene && ($elegans_gene_name eq $explain_gene || (grep /^$explain_gene$/, @ortholog_gene_names) )) {
	    $explain_gene_output .= "Testing for a 1:many orthology between elegans $elegans_gene_name and $species $ortholog_gene_name - the $species  $ortholog_gene_name gene ";
	    if ($all_match == 1) {$explain_gene_output .= "has some orthology to elegans $elegans_gene_name - so it could be OK\n"}
	    if ($all_match != 1) {$explain_gene_output .= "has no orthology to elegans $elegans_gene_name - so it is not a 1:many - we ignore complex matches like this\n"}
	  }
	  my $ortholog_number = scalar keys %{$other_ortholog->{$ortholog_gene_name}{elegans}};
	  if ($ortholog_number != 1) {$all_match = 0}
	  if (defined $explain_gene && ($elegans_gene_name eq $explain_gene || (grep /^$explain_gene$/, @ortholog_gene_names) )) {
	    $explain_gene_output .= "Testing for a 1:many orthology between elegans $elegans_gene_name and $species $ortholog_gene_name - the $species  $ortholog_gene_name gene has $ortholog_number matches to elegans $elegans_gene_name ";
	    if ($ortholog_number == 1) {$explain_gene_output .= "so it could be OK\n"}
	    if ($ortholog_number != 1) {$explain_gene_output .= "so it is not a 1:many - we ignore complex matches like this\n"}
	  }
	}	
	if ($all_match) {
	  foreach my $ortholog_gene_name (@ortholog_gene_names) {
	    $one2many{$elegans_gene_name}{$species}{$ortholog_gene_name} = 1; # NB %one2many has the elegans gene ID as the first key
	    if (defined $explain_gene && ($elegans_gene_name eq $explain_gene  || $ortholog_gene_name eq $explain_gene )) {$explain_gene_output .= "There is a 1:many orthology between elegans $elegans_gene_name and $species $ortholog_gene_name\n"}
	  }
	} else {
	    if (defined $explain_gene && $elegans_gene_name eq $explain_gene ) {$explain_gene_output .= "There is a many:many orthology between elegans $elegans_gene_name and $species genes - we ignore complex matches like this\n"}
	}
      }
    }
  }


  return (\%one2one, \%many2one, \%one2many);
}
#######################################################################################
# NB %one2one has the elegans gene ID as the first key
sub update_one2one {
  my ($one2one) = @_;

  foreach my $elegans_gene_name (keys %{$one2one}) {
    my $gene_obj = $acedb->fetch('Gene', $elegans_gene_name);
    if ($gene_obj->CGC_name) {
      my $cgc = $gene_obj->CGC_name->name;
          if ($cgc =~ /tag/) {next;}
      foreach my $species (keys %{$one2one->{$elegans_gene_name}}) {
	if (! exists $CGC_species{$species}) {$log->log_and_die("Species '$species' is not found in the hash \%CGC_species\n")}
	my $new_name = $CGC_species{$species}."-$cgc";

	foreach my $ortholog_gene_name (keys %{$one2one->{$elegans_gene_name}{$species}}) {
	  my $ortholog_obj = $acedb->fetch('Gene', $ortholog_gene_name);
	  # write out the new CGC name if there is not already a name in the ortholog or if the ortholog name is different to the new name
	  delete $existing{$ortholog_gene_name};
	  if ($ortholog_obj->CGC_name) {
	    my $ortholog_cgc = $ortholog_obj->CGC_name->name;
	    if ($ortholog_cgc ne $new_name) {
	      if (exists $paper_bork_list{$ortholog_gene_name}) {
		$log->write_to("WARNING - $ortholog_gene_name should be assigned a name $new_name but it already has a name $ortholog_cgc that has Paper_evidence - no action taken, but this should be looked at\n");
	      } elsif (exists $RNA_bork_list{$ortholog_gene_name}) {
		$log->write_to("WARNING - $ortholog_gene_name should be assigned a name $new_name but it already has a name $ortholog_cgc and is a ncRNA - no action taken, but this should be looked at\n");
	      } else {
		&write_new_orthology($elegans_gene_name, $ortholog_gene_name, $new_name); #gene ortho cgc
	      }
	    } else {
	      $correct_already{$ortholog_gene_name} = $ortholog_cgc;
	    }
	  } else {
	    &write_new_orthology($elegans_gene_name, $ortholog_gene_name, $new_name); #gene ortho cgc
	  }
	}
      }
    }
  }

}
#######################################################################################
# NB %one2many has the elegans gene ID as the first key
sub update_one2many {
  my ($one2many) = @_;

  foreach my $elegans_gene_name (keys %{$one2many}) {
    my $gene_obj = $acedb->fetch('Gene', $elegans_gene_name);
    if ($gene_obj->CGC_name) {
      my $cgc = $gene_obj->CGC_name->name;

      foreach my $species (keys %{$one2many->{$elegans_gene_name}}) {
	if (! exists $CGC_species{$species}) {$log->log_and_die("Species '$species' is not found in the hash \%CGC_species\n")}
	my $new_name = $CGC_species{$species}."-$cgc";
	
	my %used_names;
	# get prefix of new CGC name and see if any ortholog genes already use CGC names based on this
	foreach my $ortholog_gene_name (keys %{$one2many->{$elegans_gene_name}{$species}}) {
	  my $ortholog_obj = $acedb->fetch('Gene', $ortholog_gene_name);
	  if ($ortholog_obj->CGC_name) {
	    my $ortholog_cgc = $ortholog_obj->CGC_name->name;
	    if ($ortholog_cgc =~ /^$new_name\./) {
	      $used_names{$ortholog_gene_name} = $ortholog_cgc;
	    }
	  }
	}
	
	# write out the new CGC name if there is not already a name in the ortholog or if the ortholog name is different to the new name
	my $count = 1;
	my $new_name_paralog = $new_name . "." . $count;
	foreach my $ortholog_gene_name (keys %{$one2many->{$elegans_gene_name}{$species}}) {
	  delete $existing{$ortholog_gene_name};
	  my $ortholog_obj = $acedb->fetch('Gene', $ortholog_gene_name);
	  if ($ortholog_obj->CGC_name) {
	    my $ortholog_cgc = $ortholog_obj->CGC_name->name;
	    if (grep /^$ortholog_cgc$/, values %used_names) {
	      # this ortholog is already named correctly - no changes required
	      $correct_already{$ortholog_gene_name} = $ortholog_cgc;
	      next;
	    }
	    while (grep /^$new_name_paralog$/, values %used_names) {
	      $count++;
	      $new_name_paralog = $new_name . "." . $count;
	    }
	    $used_names{$ortholog_gene_name} = $new_name_paralog;
	    if (exists $paper_bork_list{$ortholog_gene_name}) {
	      $log->write_to("WARNING - $ortholog_gene_name should be assigned a name $new_name_paralog but it already has a name $ortholog_cgc that has Paper_evidence - no action taken\n");
	    } elsif (exists $RNA_bork_list{$ortholog_gene_name}) {
	      $log->write_to("WARNING - $ortholog_gene_name should be assigned a name $new_name_paralog but it already has a name $ortholog_cgc and is a ncRNA - no action taken\n");
	    } else {
	      &write_new_orthology($elegans_gene_name, $ortholog_gene_name, $new_name_paralog); #gene ortho cgc
	    }
	    
	  } else {
	    while (grep /^$new_name_paralog$/, values %used_names) {
	      $count++;
	      $new_name_paralog = $new_name . "." . $count;
	    }
	    $used_names{$ortholog_gene_name} = $new_name_paralog;
	    &write_new_orthology($elegans_gene_name, $ortholog_gene_name, $new_name_paralog); #gene ortho cgc
	  }
	}
      }
    }
  }
  
  
  
}
#######################################################################################
# NB %many2one has the non-elegans gene ID as the first key
sub update_many2one {
  my ($many2one) = @_;
  
  # Look at all of the elegans genes that are orthologs of this other species gene.
  # If the elegans genes all have CGC names that are formed from a basic name, 
  # then we accept them as a suitable origin for the new CGC name.
  #
  # It might be thought better to use the elegans Paralog data, but that results
  # in groups of paralogs that share names that often have no common gene class.
  # This results in no easy way to transfer the many CGC names.
  # At least this method will track changes in the nomenclature of the elegans genes.

  foreach my $ortholog_gene_name (keys %{$many2one}) {
    my $ortholog_gene_obj = $acedb->fetch('Gene', $ortholog_gene_name);
    my $ortholog_cgc;
    if ($ortholog_gene_obj->CGC_name) {
      $ortholog_cgc = $ortholog_gene_obj->CGC_name->name;
    }

    foreach my $species (keys %{$many2one->{$ortholog_gene_name}}) { # there should only be one species for this $ortholog_gene_name
 
      # check through the CGC names of the elegans genes that are orthologs of the non-elegans gene: $ortholog_gene_name
      # we want them all to have a common root CGC name with paralog .digit added
      my $all_match = 1;
      my $cgc_root = '';
      my $elegans_gene_name; # define this outside the foreach loop so that we have an example name to pass to &write_new_orthology()
      foreach $elegans_gene_name (keys %{$many2one->{$ortholog_gene_name}{$species}}) {
	my $elegans_gene_obj = $acedb->fetch('Gene', $elegans_gene_name);
	if ($elegans_gene_obj->CGC_name) {
	  my $elegans_cgc = $elegans_gene_obj->CGC_name->name;
	  my ($root) = $elegans_cgc =~ /^(\S+)\.d+/;
	  if (!defined $root) {$all_match = 0; $root = ''}
	  if ($cgc_root eq '') {$cgc_root = $root}
	  if ($cgc_root ne $root) {$all_match = 0}
	} else {
	  $all_match = 0; # if any of the elegans genes has no CGC name, then don't transfer a name - better safe than sorry
	}
      }
     
      if ($all_match) {
	if (! exists $CGC_species{$species}) {$log->log_and_die("Species '$species' is not found in the hash \%CGC_species\n")}
	my $new_name = $CGC_species{$species}."-$cgc_root";
	if (!defined $ortholog_cgc || $ortholog_cgc ne $new_name) {
	  delete $existing{$ortholog_gene_name};
	  if (defined $ortholog_cgc) {
	    if (exists $paper_bork_list{$ortholog_gene_name}) {
	      $log->write_to("WARNING - $ortholog_gene_name should be assigned a name $new_name but it already has a name $ortholog_cgc that has Paper_evidence - no action taken\n");
	    } elsif (exists $RNA_bork_list{$ortholog_gene_name}) {
	      $log->write_to("WARNING - $ortholog_gene_name should be assigned a name $new_name but it already has a name $ortholog_cgc and is a ncRNA - no action taken\n");
	    } else {
	      &write_new_orthology($elegans_gene_name, $ortholog_gene_name, $new_name); #gene ortho cgc
	    }
	  } else {
	    &write_new_orthology($elegans_gene_name, $ortholog_gene_name, $new_name); #gene ortho cgc
	  }
	} elsif (defined $ortholog_cgc) {
	  $correct_already{$ortholog_gene_name} = $ortholog_cgc;
	}
      }
    }
  }
  
  

}
#######################################################################################
sub write_new_orthology {
    my $gene = shift;
    my $ortholog = shift;
    my $new_name = shift;
    my $geneObj = $acedb->fetch('Gene',$ortholog);
    my $version = $geneObj->Version->name;
    $version++;

    $assigned_names{$new_name} = $ortholog;
    print $ace "\nGene : $ortholog\nVersion $version\n";
    print $ace "CGC_name $new_name From_analysis Inferred_from_orthology\n";
    print $ace "Public_name $new_name\n";
    print $ace "Version_change $version now $person{$user} Name_change CGC_name $new_name\n";
    my ($class) = $new_name =~ /-(\w+)-/;
    print $ace "Gene_class $class\n";
    $log->write_to("Transfering CGC_name: $new_name from $gene to $ortholog");

    if($geneObj->CGC_name) {
	#get existing evidence
	my $evidence = $geneObj->CGC_name->right(2);
	my $old_name = $geneObj->CGC_name->name;
	if ($new_name eq $old_name) {$log->log_and_die("ERROR - Attempt to assign the same existing CGC_name $new_name to $ortholog\n")}  # 
	my $old_class = $geneObj->Gene_class->name;
	print $ace "Version_change ",$version," now $person{$user} Name_change Other_name $old_name\n";
	print $ace "Other_name $old_name\n";
	#print old name as Other_name with evidence transferred.
	foreach ($geneObj->CGC_name(2)){
	    print $ace "Other_name $old_name ", $_->name."\t\"".$_->right->name."\"\n";
	}

	#gene class update
	print $ace "\nGene_class : $old_class\nOld_member $old_name\n";
	$log->write_to(" and replacing $old_name");

	# If there is an Other_name that is the same as the new CGC name, then remove it.
	if ($geneObj->Other_name) {
	  my @other_name = $geneObj->Other_name;
	  foreach my $other_name (@other_name) {
	    if ($other_name->name eq $new_name) {
	      print $ace "\nGene : $ortholog\n";
	      print $ace "-D Other_name $other_name\n";
	    }
	  }
	}

    }

    print $namedb "$ortholog\t$new_name\n";
    $log->write_to("\n");
}

#######################################################################################
sub delete_incorrect_orthology {
    my ($ortholog, $existing_CGC_name, $remark) = @_;

    my $geneObj = $acedb->fetch('Gene', $ortholog);
    my $seqname = $geneObj->Sequence_name->name;
    my $version = $geneObj->Version->name;
    $version++;

#    $log->write_to("Deleting CGC_name: $existing_CGC_name in $ortholog\n");

    my $date = date();

    print $ace "\nGene : $ortholog\n";
    print $ace "-D CGC_name\n";
    print $ace "-D Public_name $existing_CGC_name\n";
    print $ace "-D Gene_class\n";

    print $ace "\nGene : $ortholog\n";
    print $ace "Version $version\n";
    print $ace "Version_change $version now $person{$user} Name_change CGC_name \"Removed\"\n";
    print $ace "Public_name \"$seqname\"\n";
    print $ace "Remark \"[$date $user] $remark\" Curator_confirmed $person{$user}\n" if ($remark);

    print $deletenamedb "$ortholog\t$existing_CGC_name\n";
}


#######################################################################################
# store existing non-elegans genes to check that all have been correctly set
sub store_existing {
  my (@genes) = @_;

  my %existing;
  foreach my $gene (@genes) {
    my $gene_name = $gene->name;
    if ($gene_name =~ /^WBGene/) {
      if ($gene->CGC_name) {
	my $cgc = $gene->CGC_name->name;
	$existing{$gene_name} = $cgc;
      }
    }
  } 
  return %existing;
}
#######################################################################################
# report any genes which are DEAD but which have a CGC name that we have just assigned to another gene
# report any non-elegans genes that have a CGC name that we did not set
sub check_existing {
  my ($dead_existing, $existing, $assigned_names, $correct_already) = @_;

  foreach my $dead_gene (keys %{$dead_existing}) {
    my $name = $existing->{$dead_gene};
    if (exists $assigned_names->{$name}) {
      $log->write_to("The Dead gene $dead_gene has a name $name but this is bogus: this name has now correctly been assigned to $assigned_names->{$name} - no action has been taken, but it needs attention\n");
    }
  }


  foreach my $gene (keys %{$existing}) {
    my $name = $existing->{$gene};
    if (exists $assigned_names->{$name}) {
      if (exists $paper_bork_list{$gene}) {
	$log->write_to("The gene $gene has a name $name but the only evidence for this is Paper_evidence: this name has correctly been assigned to $assigned_names->{$name} - no action has been taken, but it needs attention\n");
      } elsif (exists $RNA_bork_list{$gene}) {
	$log->write_to("The gene $gene has a name $name and is a ncRNA: this name has correctly been assigned to $assigned_names->{$name} - no action has been taken, but it needs attention\n");
      } else {
	$log->write_to("The gene $gene has a name $name but this is bogus: it will be deleted: this name has correctly been assigned to $assigned_names->{$name}\n");
	delete_incorrect_orthology($gene, $name, "The CGC name $name was deleted as there was no evidence to support it. This name has correctly been assigned to $assigned_names->{$name}");
      }      
    } else {
      if (exists $paper_bork_list{$gene}) {
	$log->write_to("The gene $gene has a name $name but the only evidence for this is Paper_evidence - no action has been taken\n");
      } elsif (exists $RNA_bork_list{$gene}) {
	$log->write_to("The gene $gene has a name $name and is a ncRNA - no action has been taken\n");
      } else {
	$log->write_to("The gene $gene has a name $name but this is bogus: it will be deleted\n");
	delete_incorrect_orthology($gene, $name, "The CGC name $name was deleted as there was no evidence to support it.");
      }
    }
  }

  if ($verbose) {
    foreach my $gene (keys %{$correct_already}) {
      my $name = $correct_already->{$gene};
      $log->write_to("The gene $gene already has a correct CGC name $name\n");
    }
  } else {
    my $count = scalar keys %{$correct_already};
    $log->write_to("\nThere are $count genes that already have a correct CGC name\n");
  }
}

#######################################################################################
# we do not wish to change the names of Tier II genes with paper_evidence
# for example, there are some Brugia-specific gene classes that have been hand-curated by Michael Paulini - these all have Paper_evidence
sub get_paper_bork_list {

  my @paper_genes = $acedb->fetch (-query => 'FIND Gene WBGene* WHERE CGC_name AND NEXT AND NEXT = "Paper_evidence"');
  my %paper_bork;

  foreach my $paper_gene (@paper_genes) {
    $paper_bork{$paper_gene->name} = 1;
  }

  return %paper_bork;
}

#######################################################################################
# we do not wish to change the names of Tier II genes which are ncRNA genes - these will have had their names explicitly set.
sub get_RNA_bork_list {

  #                                                                              miRNA                     ncRNA                     ncRNA_gene               miRNA                      siRNA                     snoRNA                   snRNA                     tRNA_gene                 rRNA_gene                 piRNA_gene               lincRNA_gene
  my @RNA_genes = $acedb->fetch (-query => 'FIND Gene WBGene* WHERE Biotype = "SO:0001265" OR Biotype = "SO:0000655" OR Biotype = "SO:0001263" OR Biotype = "SO:0001265" OR Biotype = "SO:0001266" OR Biotype = "SO:0001267" OR Biotype = "SO:0001268" OR Biotype = "SO:0001272" OR Biotype = "SO:0001637" OR Biotype = "SO:0001638" OR Biotype = "SO:0001641"');
  my %RNA_bork;

  foreach my $RNA_gene (@RNA_genes) {
    $RNA_bork{$RNA_gene->name} = 1;
  }

  return %RNA_bork;
}

#######################################################################################

sub remove_genes_with_paper_evidence {
  my ($paper_bork_list, $genes) = @_;
  my @OK;

  foreach my $gene (@{$genes}) {
    if (exists $paper_bork_list->{$gene->name}) {
      my $gene_name = $gene->name;
      $log->write_to("Ignore gene $gene_name with Paper_evidence\n"); 
      next;
    }
    push @OK, $gene;
  }


  return @OK;
}
#######################################################################################
sub date {
  my ($day, $mon, $yr)  = (localtime)[3,4,5];
  my $date = sprintf("%02d%02d%02d",$yr-100, $mon+1, $day);
}
#######################################################################################
# do a sanity check of the non-elegans data
sub sanity_check {
  my (@genes) = @_;
  
  my $PASS = $user;
  my $USER = $user;
  my $DOMAIN = 'Gene';
  
  my $DB = 'nameserver_live;web-wwwdb-core-02:3449';
  my $db = NameDB_handler->new($DB,$USER,$PASS);
  $db->setDomain($DOMAIN);



  foreach my $gene (@genes) {
    my $gene_name = $gene->name;

    # fix missing Public_name - commented out as not required at present
    if (0) {    
      if (!defined $gene->Sequence_name) {
	$log->write_to("WARNING - $gene_name has no Sequence_name - probably an uncloned non-elegans Gene?\n");
      } else {
	my $seqname = $gene->Sequence_name->name;
	if (!defined $gene->Public_name) {
	  if ($gene->CGC_name) {
	    my $cgc = $gene->CGC_name->name;
	    print $ace "\nGene : $gene_name\n";
	    print $ace "Public_name \"$cgc\"\n";
	  } else {
	    print $ace "\nGene : $gene_name\n";
	    print $ace "Public_name \"$seqname\"\n";
	  }
	}
      }
    }
    
    # check against Nameserver
    my $public = $gene->Public_name;
    if (defined $public) {
      my $ID = $db->idGetByTypedName('Public_name'=>$public)->[0];
      if (defined $ID) {
	if ($gene_name ne $ID) {
	  $log->write_to("ERROR - $gene_name has a Public_name $public but this is the name of $ID in the Nameserver\n");
	}
      } else {
	$log->write_to("ERROR - $gene_name has a Public_name $public but this is not found in the Nameserver\n");
      }
      
    }
  }
}

