#!/usr/local/bin/perl5.8.0 -w
#
# inherit_GO_terms.pl
#
# map GO_terms to ?Sequence objects from ?Motif and ?Phenotype
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2013-01-10 10:00:43 $      

use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Ace;

##############################
# Script variables (run)     #
##############################

my ($help, $debug, $motif, $phenotype, $tmhmm, $store);
my $verbose;             # for toggling extra output
my $maintainers = "All"; # who receives emails from script
my $noload;              # generate results but do not load to autoace
my $database;
my $species;
my $test;

##############################
# command-line options       #
##############################

GetOptions ("help"      => \$help,
            "debug=s"   => \$debug,
            "phenotype" => \$phenotype,
	    "tmhmm"	=> \$tmhmm,
	    "motif"     => \$motif,
	    "noload"    => \$noload,
    	    "store:s"   => \$store,
    	    "database:s" => \$database,
    	    "species:s"  => \$species,
	    "test"       => \$test,
    	);

# Display help if required
&usage("Help") if ($help);

# recreate configuration 
my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("cant restore wormbase from $store\n"); 
} else { 
  $wormbase = Wormbase->new( -debug => $debug, 
                             -test => $test, 
                             -organism => $species );
}

# Variables Part II (depending on $wormbase) 
$debug = $wormbase->debug if $wormbase->debug;    # Debug mode, output only goes to one user

my $log=Log_files->make_build_log($wormbase);

##############################
# Paths etc                  #
##############################

my $tace      = $wormbase->tace;      # tace executable path
# Database path

my $dbpath    = $wormbase->orgdb;   
$dbpath = $database if (defined $database);

my %cds2gene;
$wormbase->FetchData('cds2wbgene_id',\%cds2gene);                                

my $out=$wormbase->acefiles."/inherited_GO_terms.ace";

open (OUT,">$out") or $log->log_and_die("cant open $out :$!\n");

########################################
# Connect with acedb server            #
########################################


&motif($dbpath)     if ($motif);
&phenotype($dbpath) if ($phenotype);
&tmhmmGO($dbpath)   if ($tmhmm);

close OUT;

##############################
# read acefiles into autoace #
##############################

$wormbase->load_to_database($dbpath,$out,'inherit_GO_terms', $log) unless ($noload || $debug) ;

##############################
# mail $maintainer report    #
##############################
$log->mail();

##############################
# hasta luego                #
##############################

exit(0);

########################################################################################
####################################   Subroutines   ###################################
########################################################################################

########################################################################################
# motif to sequence mappings                                                           #
########################################################################################

sub motif {
  my ($dbpath) = @_;
  my $def = "$dbpath/wquery/SCRIPT:inherit_GO_terms.def";
  
  # these GO terms should not be attached to the Gene or CDS
  my @stopterms = (
		   'sporulation',
		   'forespore',
		   'photosynthesis',
		   'photosynthetic',
		   'chlorophyll',
		  );
  my @stopGO = (
		'GO:0009772', # photosynthetic electron transport in photosystem II
		'GO:0045282', # plasma membrane succinate dehydrogenase complex (only_in_taxon Bacteria)
		'GO:0009288', # bacterial-type flagellum
		'GO:0007391', # dorsal closure (only_in_taxon Insecta)
		'GO:0009103', # lipopolysaccharide biosynthetic process (only in prokaryotes)
      );
  
  # get the GO terms
  my $term_def = &write_GO_def;
  my $term_query = $wormbase->table_maker_query($dbpath,$term_def);
  my %terms;
  while(<$term_query>) {
    chomp;
    s/\"//g;  #remove "
    next if (/acedb/ or /\/\//);
    my @data = split("\t",$_);
    my ( $GO, $term) = @data;
    $terms{$GO} = $term;
  }
  
  my $query = $wormbase->	table_maker_query($dbpath, $def);
  while(<$query>) {
    s/\"//g;#"
    next if (/acedb/ or /\/\//);
    my($motif,$GO,$protein,$cds,$gene) = split;
    next if (! defined $gene || ! defined $cds || ! defined $GO || ! defined $motif);
    next if (&matching(\@stopterms, \@stopGO, $terms{$GO}, $GO, $motif, $protein, $gene));
    print OUT "\nGene : $gene\nGO_term \"$GO\" IEA inferred_automatically \"$motif\"\n";
    print OUT "\nCDS  : \"$cds\"\nGO_term \"$GO\" IEA inferred_automatically \"$motif\"\n";
  }
}

########################################################################################
# phenotype to sequence mappings                                                       #
########################################################################################

sub matching {
  my ($stopterms_aref, $stopGO_aref, $GO_term, $GO, $motif, $protein, $gene) = @_;
  if (! defined $GO_term) {return 0;}
  foreach my $term (@{$stopterms_aref}) {
    if ($GO_term =~ /\b$term\b/) {
      $log->write_to("The invalid term '$term' was found in the description '$GO_term' of GO-term $GO ($motif) from protein $protein and will not be attached to $gene\n");
      return 1;
    }
  }
  foreach my $term (@{$stopGO_aref}) {
    if ($GO eq $term) {
      $log->write_to("The invalid GO ID '$term' ('$GO_term' $motif) from protein $protein was found and will not be attached to $gene\n");
      return 1;
    }
  }
  return 0;
}



sub tmhmmGO {
  my ($dbpath) = @_;

  my $mydb = Ace->connect(-path=>$dbpath,
                          -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};
  

  my $query = "find protein where species = \"".$wormbase->full_name."\" where Feature AND NEXT = \"Tmhmm\"; follow Corresponding_CDS; follow Gene";
  my $genes = $mydb->fetch_many(-query => $query);
  while(my $gene = $genes->next){
    print OUT "\nGene : ".$gene->name."\nGO_term \"GO:0016021\"\tIEA\tInferred_automatically\t\"CBS:TMHMM\"\n";
  }

  $mydb->close;
}

sub phenotype {
  my ($dbpath) = @_;
   
  my %phenotype_names;

  my $mydb = Ace->connect(-path=>$dbpath,
                          -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};
  
  my $query = "find protein where species = \"".$wormbase->full_name."\" where Feature AND NEXT = \"Tmhmm\"; follow Corresponding_CDS; follow Gene";
  my $phenos = $mydb->fetch_many(-class => 'Phenotype');
  while(my $pheno = $phenos->next){
    if ($pheno->Primary_name) {
      $phenotype_names{$pheno->name} = $pheno->Primary_name;
    }
  }
  $mydb->close;

 
  my $def = &write_phenotype_def;
  my %include_list; #only papers in this list should be included.
  READARRAY: while (<DATA>) {
    chomp;
    $include_list{$_} =1;
  }

  my $tm_query = $wormbase->table_maker_query($dbpath,$def);
  while(<$tm_query>) {
    s/\"//g;  #remove "
    next if (/acedb/ or /\/\//);
    my @data = split("\t",$_);
    my ( $rnai, $paper, $cds, $gene, $phenotype_id,$go,$not) = @data;
    chomp $not;
    next if (($paper !~ /WBPaper/) or (! defined $phenotype_id) or ((defined $not)and ($not eq 'Not')));
    next unless ($include_list{$paper});
    my $phenotype; 
    if($phenotype_id =~ /WBPheno/) {
      $phenotype = exists $phenotype_names{$phenotype_id} ? $phenotype_names{$phenotype_id} : $phenotype_id;
    }
    else {
      next;
    }
    unless($gene and $phenotype and $go) {
      $log->write_to("bad data (causes a phenotype, but doesn't affect a gene - inform CalTech) $_");
      next;
    }
    print OUT "\nCDS : \"$cds\"\nGO_term \"$go\" IMP Inferred_automatically \"$phenotype ($phenotype_id|$rnai)\"\n" if $cds ;
    print OUT "\nGene : \"$gene\"\nGO_term \"$go\" IMP Inferred_automatically \"$phenotype ($phenotype_id|$rnai)\"\n" if $gene;
  }
  #tidy up
  $wormbase->run_command("rm -f $def", $log);
}


######################################################################


sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}


# this will write out an acedb tablemaker defn to a temp file
sub write_phenotype_def {
	my $def = '/tmp/inherit_GO.def';
	open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
	my $txt = <<END;
Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class RNAi 
From 1 
  
Colonne 2 
Width 12 
Mandatory
Visible 
Class 
Class Paper 
From 1 
Tag Reference 

Colonne 3 
Width 12 
Optional 
Visible 
Class 
Class CDS 
From 1 
Tag Predicted_gene  
 
Colonne 4 
Width 12 
Optional 
Visible 
Class 
Class Gene 
From 1
Tag Gene  
 
Colonne 5
Width 12 
Optional 
Visible 
Class 
Class Phenotype 
From 1 
Tag Phenotype  
 
Colonne 6 
Width 12 
Mandatory 
Visible 
Class 
Class GO_term 
From 5
Tag GO_term 
 
Colonne 7
Width 12 
Optional 
Visible 
Next_Tag 
Right_of 5 

END

	print TMP $txt;
	close TMP;
	return $def;
}


# this will write out an acedb tablemaker defn to a temp file
sub write_GO_def {
	my $def = '/tmp/inherit_GO_term.def';
	open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
	my $txt = <<END;
Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class GO_term 
From 1 
  
Colonne 2 
Width 120
Mandatory
Visible 
Class 
Class Text 
From 1 
Tag Term


END

	print TMP $txt;
	close TMP;
	return $def;
}


__DATA__
WBPaper00004402
WBPaper00004403
WBPaper00004540
WBPaper00004651
WBPaper00004769
WBPaper00005599
WBPaper00005654
WBPaper00006395
WBPaper00024497
WBPaper00024925
WBPaper00025054
WBPaper00026763
WBPaper00005736
WBPaper00026593
WBPaper00028783
WBPaper00029254
WBPaper00030951




__END__
=pod

=head2   NAME - inherit_GO_terms.pl


=head1 USAGE

=over 4

=item inherit_GO_terms.pl [-options]

=back

inherit_GO_terms.pl assigns GO terms to sequences based on Interpro motifs
and RNAi phenotypes. The resulting acefile will be loaded into the database
as part of the script run (default database is autoace).

inherit_GO_terms.pl mandatory arguments:

=over 4

=item none, (but it won\'t do anything)

=back

inherit_GO_terms.pl OPTIONAL arguments:

=over 4

=item -motif, parse Interpro motif data

=item -phenotype, parse phenotype data

=item -noload, do not upload results to autoace

=item -debug, debug (results not loaded into autoace)

=item -help, help

=item -verbose, toggle extra output to screen

=item -store <storable_file>, specifiy stored commandline options

=back

=cut
