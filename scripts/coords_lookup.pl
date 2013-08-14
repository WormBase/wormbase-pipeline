#!/software/bin/perl -w
#
# ilk.pl                           
# 
# by Gary Williams                         
#
# This finds information about the protein coded for by the specified gene/CDS
# and finds information about the homologous wormpep proteins
#
# Last updated by: $Author: pad $     
# Last updated on: $Date: 2013-08-14 16:12:51 $      

# To do:
# accept wormpep ID as input
# align the genes
# put all images in one canvas
# scroll the zoom of the canvas
# scroll up and down
# display mass-spec peptide positions
# display the descriptions of the protein domains
# 
# 



use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
#use Sequence_extract;
use Coords_converter;
use Tk;
use Tk::Dialog; 
use Tk::FileDialog;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($input, $textonly, $group);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "input:s"    => \$input,
	    "group:s"    => \$group,
	    "textonly"   => \$textonly,	# don't show the interface
	    );

# always in test and debug mode
$test = 1;
$debug = "gw3";


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

##########################
# MAIN BODY OF SCRIPT
##########################


# TK global variables
my $status;			# label to display status and other information
my $VERSION = '1.1.0';
my $MW;				# Main Window handle
my $quick_entry;		# entry field for quick frame
my $compare_cds;		# entry field for comparing a CDS and the target to the homologs

my %data;

# open an ACE connection to parse details
my $tace            = $wormbase->tace;        # TACE PATH

# get the type of grouping
if (!defined $group || $group eq "") {$group = "same"}
if ($group ne "all" && 
    $group ne "same" && 
    $group ne "other" &&
    $group ne "elegans" &&
    $group ne "briggsae" &&
    $group ne "remanei") {
  die "-group should be set to one of:\n'same' (only same species), 'all' (all species), 'other' (only other species), or a species name.\nThe default is 'same'\n";
}

print "Connecting to Ace\n";
my $db = Ace->connect (-path => $wormbase->database('current'),
                       -program => $tace) || die "cannot connect to database at $wormbase->database('current')\n";


# get CDS name
if (! defined $input || $input eq "") {
  if ($textonly) {
    die "please specify the cds/gene name\n";
  }
} else {

# get the data on the cds
  %data = &get_data($input, $group);
}

# show the interface
if (! $textonly) {
  &interface(%data);
}


# close the ACE connection
$db->close;

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);






##############################################################
#
# Subroutines
#
##############################################################

##########################################
# start the TK interface
# this never returns
sub interface {

  my (%data) = @_;

    
  $MW = MainWindow->new;
  $MW->optionAdd('*BorderWidth' => 1); # make the style better with a thinner border

  # if there is no data, get it
  if (! %data) {		
    my $cds_dialog =   $MW->Dialog(-title   => 'Enter CDS name', 
                                 -text    => "Enter CDS name:", 
                                 -buttons => ['OK'], 
                                 );
    # +++ %data = get_data($cds, $group);
  }

  &init(%data); 

  &MainLoop; 
}

##########################################
# get the data on the CDS

sub get_data {
  my ($cds, $group) = @_;
  my %data;


  $cds = uc($cds);

# get the protein product
  my $protein = &get_protein($cds);
  print "Protein: $protein\n";
  $data{protein} = $protein;

# get CDS exon details
  my ($gene, $cgc_name, $clone, $lab, $rnai, $cds_start, $cds_end, $exons_start_ref, $exons_end_ref) = &get_cds_details($cds);
  print "Position: $clone:$cds_start..$cds_end\n";
  print "Gene name: $gene\n";
  print "CGC name: $cgc_name\n" if (defined $cgc_name);
  print "Lab: $lab\n";
  $data{cds_name} = $cds;	# NB only the target has this key/value pair
  $data{gene} = $gene;
  $data{cgc_name} = $cgc_name;
  $data{clone} = $clone;
  $data{lab} = $lab;
  $data{cds_start} = $cds_start;
  $data{cds_end} = $cds_end;
  $data{exons_start_ref} = $exons_start_ref;
  $data{exons_end_ref} = $exons_end_ref;

  for (my $i = 0; $i < @{$exons_start_ref}; $i++) {
    print "Exon ",$i+1,": $exons_start_ref->[$i]..$exons_end_ref->[$i] length: ",$exons_end_ref->[$i] - $exons_start_ref->[$i] + 1,"\n";
  }

# get the RNAi phenotype
  $data{phenotype} = &get_phenotype($rnai);
  print "Phenotype: $data{phenotype}\n";

# get domains and homologs in the protein product
  my ($domains_href, $homologs_href, $mass_spec_href, $signalp) = &get_domains($protein, "target", $group);
  $data{domains_href} = $domains_href;
  $data{mass_spec_href} = $mass_spec_href;
  $data{homologs_href} = $homologs_href;
  $data{signalp} = $signalp;
  
# print the domains
  print "\nDomains of $protein\n";
  foreach my $domain (keys %{$domains_href}) {
    print "$domain $domains_href->{$domain}->{'start'} $domains_href->{$domain}->{'end'}\n";
  }

# print the mass spec
  print "\nMass-spec of $protein\n";
  foreach my $mass_spec (keys %{$mass_spec_href}) {
    print "$mass_spec @{$mass_spec_href->{$mass_spec}->{'start'}} @{$mass_spec_href->{$mass_spec}->{'end'}}\n";
  }

# get the wormpep homologs
  print "\nHomologs of $protein\n";
  foreach my $homolog (keys %{$homologs_href}) {

    print "\nHomolog: $homolog\n";

# get the CDS names of the homologs
    my $cds_name_of_homolog = &get_cds_name($homolog);
    if (! defined $cds_name_of_homolog) {
      delete $homologs_href->{$homolog};
      next;
    }
    print "\tCDS name: $cds_name_of_homolog\n";
    $data{homologs}{$homolog}{cds_name_of_homolog} = $cds_name_of_homolog;

# get the homolog CDS exon details
    my ($hom_gene, $hom_cgc_name, $hom_clone, $hom_lab, $hom_rnai, $hom_cds_start, $hom_cds_end, $hom_exons_start_ref, $hom_exons_end_ref) = &get_cds_details($cds_name_of_homolog);
    print "\tPosition: $hom_clone:$hom_cds_start..$hom_cds_end\n";
    print "\tGene name: $hom_gene\n";
    print "\tLab: $hom_lab\n";
    $data{homologs}{$homolog}{gene} = $hom_gene;
    $data{homologs}{$homolog}{cgc_name} = $hom_cgc_name;
    $data{homologs}{$homolog}{clone} = $hom_clone;
    $data{homologs}{$homolog}{lab} = $hom_lab;
    $data{homologs}{$homolog}{cds_start} = $hom_cds_start;
    $data{homologs}{$homolog}{cds_end} = $hom_cds_end;
    $data{homologs}{$homolog}{exons_start_ref} = $hom_exons_start_ref;
    $data{homologs}{$homolog}{exons_end_ref} = $hom_exons_end_ref;

# get the RNAi phenotype
    $data{homologs}{$homolog}{phenotype} = &get_phenotype($hom_rnai);
    print "Phenotype: $data{homologs}{$homolog}{phenotype}\n";

    my @starts       = @{ $homologs_href->{$homolog}->{'start'} };
    my @ends         = @{ $homologs_href->{$homolog}->{'end'} };
    my @other_starts = @{ $homologs_href->{$homolog}->{'other_start'} };
    my @other_ends   = @{ $homologs_href->{$homolog}->{'other_end'} };
    my $score        = $homologs_href->{$homolog}->{'score'};

    # store the homologous regions of the target to this homolog
    $data{homologs}{$homolog}{homology_starts} = [ @starts ];
    $data{homologs}{$homolog}{homology_ends} = [ @ends ];
    $data{homologs}{$homolog}{homology_other_starts} = [ @other_starts ];
    $data{homologs}{$homolog}{homology_other_ends} = [ @other_ends ];
    $data{homologs}{$homolog}{score} = $score;

    foreach my $start (@starts) {
      my $end = shift @ends;
      my $other_start = shift @other_starts;
      my $other_end = shift @other_ends;
      print "\tHomology: $start..$end homolog_positions: $other_start..$other_end Score: $score\n";
    }

# get domains in the homologs
    my ($homolog_domains_href, $homologs_homologs_href, $homolog_mass_spec_href, $signalp) = &get_domains($homolog, "homolog", $group);
    $data{homologs}{$homolog}{domains_href} = $homolog_domains_href;
    $data{homologs}{$homolog}{mass_spec_href} = $homolog_mass_spec_href;
    $data{homologs}{$homolog}{homologs_href} = $homologs_homologs_href;	# this is the set of homologous matches of this homolog to other proteins - probably not useful
    $data{homologs}{$homolog}{signalp} = $signalp;

# print the domains in the homologs
    foreach my $domain (keys %{$homolog_domains_href}) {
      print "\tDomain: $domain $homolog_domains_href->{$domain}->{'start'} $homolog_domains_href->{$domain}->{'end'}\n";
    }

# print the mass-spec in the homologs
    foreach my $mass_spec (keys %{$homolog_mass_spec_href}) {
      print "\tMass-spec: $mass_spec @{$homolog_mass_spec_href->{$mass_spec}->{'start'}} @{$homolog_mass_spec_href->{$mass_spec}->{'end'}}\n";
    }

    for (my $i = 0; $i < @{$hom_exons_start_ref}; $i++) {
      print "\tExon ",$i+1,": $hom_exons_start_ref->[$i]..$hom_exons_end_ref->[$i] length: ",$hom_exons_end_ref->[$i] - $hom_exons_start_ref->[$i] + 1,"\n";
    }

  }

  # get the description of the domains - most domains will be duplicated
  # across the proteins, so we only have to get their descriptions
  # once here and not for each homolog.
  my %domains_desc;
  foreach my $homolog (keys %{$data{homologs}}) {
    foreach my $domain (keys %{$data{homologs}{$homolog}{domains_href}}) {
      if (! exists $domains_desc{$domain}) {
	my $motif_obj = $db->fetch("Motif" => $domain);
	if (defined $motif_obj) {
	  my $title = $motif_obj->Title; # get the description
	  if (defined $title) {
	    print "$domain:\t $title\n";
	    $domains_desc{$domain} = $title;
	  }
	}
      }
    }
  }
  $data{domain_desc_href} = \%domains_desc;

  return %data;
}

##########################################
# my ($gene, $clone, $lab, $rnai, $cds_start, $cds_end, $exons_start_ref, $exons_end_ref) = &get_cds_details($cds_name);
# get the clone, clone start position and exons list from a CDS
sub get_cds_details {
  my ($cds_name) = @_;

  my $cds_obj = $db->fetch("CDS" => $cds_name);
  # debug
  if (! defined $cds_obj) {
    print "Can't fetch CDS object for $cds_name\n";
    return (undef, undef, undef, undef, undef, undef, undef);
  }

  my $gene = $cds_obj->Gene;
  if (! defined $gene) {
    print "Can't fetch gene name for $cds_name\n";
    return (undef, undef, undef, undef, undef, undef, undef);
  }

  # get the CGC name, if any
  my $cgc_name = $gene->CGC_name;
  

  my $clone = $cds_obj->Sequence;
  if (! defined $clone) {
    print "Can't fetch clone for $cds_name\n";
    return (undef, undef, undef, undef, undef, undef, undef);
  }

  my $lab = $cds_obj->From_laboratory;
  if (! defined $lab) {
    print "Can't fetch laboratory for $cds_name\n";
    $lab = '';
  }

  # if there is no RNAi phenotype, then $rnai is undef
  my $rnai = $cds_obj->RNAi_result;

  # get the named object
  # then get the data following a tag in that object
  # and the NEXT data following the tag
  my @exons_start = $cds_obj->Source_exons;
  my @exons_end = $cds_obj->Source_exons(2);
  $cds_obj->DESTROY();
  if (! @exons_start || ! @exons_end) {
    print "Can't fetch exons for $cds_name\n";
    return (undef, undef, undef, undef, undef, undef, undef);
  }

  # get the start position of the CDS in the clone's SMap
  #print "Clone = $clone\n";
  my $clone_obj = $db->fetch("Sequence" => $clone);
  if (! defined $clone_obj) {
    print "Can't fetch clone object for $clone\n";
    return (undef, undef, undef, undef, undef, undef, undef);
  }

  my @cds_names = $clone_obj->CDS_child;
  my @cds_starts = $clone_obj->CDS_child(2);
  my @cds_ends = $clone_obj->CDS_child(3);
  $clone_obj->DESTROY();
  if (! @cds_names || ! @cds_starts || ! @cds_ends) {
    print "Can't fetch CDS starts/ends for $clone\n";
    return (undef, undef, undef, undef, undef, undef, undef);
  }

  # find the start of the CDS we want
  my $count = 0;
  foreach my $next_cds_name (@cds_names) {
    if ($next_cds_name eq $cds_name) {last;}
    $count++;
  }
  my $cds_start = $cds_starts[$count];
  my $cds_end = $cds_ends[$count];

  # if the CDS is on the reverse sense, then $cds_end > $cds_start
  return ($gene, $cgc_name, $clone, $lab, $rnai, $cds_start, $cds_end, \@exons_start, \@exons_end);

}



##########################################
# get the protein product of a CDS

sub get_protein {
  my ($cds_name) = @_;

  my $cds_obj = $db->fetch("CDS" => $cds_name);
  # debug
  if (! defined $cds_obj) {
    print "Can't fetch CDS object for $cds_name\n";
    return (undef);
  }

  my $protein = $cds_obj->Corresponding_protein;
  $cds_obj->DESTROY();

  return $protein;
}

##########################################
# get the domain and homolog and mass-spec information for this protein
# my ($domains_href, $homologs_aref, $mass_spec_aref) = &get_domains($protein, "target");

sub get_domains {
  my ($protein, $type, $group) = @_;

  my (%domains, %homologs, %mass_spec, $signalp);

  # table of first two letters of protein names of various species
  my %species_letters = ('elegans'  => 'WP',
			 'remanei'  => 'RP',
			 'briggsae' => 'BP',
			 );

  # get the first two letters of the protein name - other members of
  # this group will share this name.
  my $group_start = substr($protein, 0, 2);

  my $protein_obj = $db->fetch("Protein" => $protein);
  if (! defined $protein_obj) {
    print "Can't fetch Protein object for $protein\n";
    return (undef, undef);
  }

  # get the signalp status
  $signalp = defined $protein_obj->at('Feature.Signalp')?1:0;

# this skips over the repeated domains
  my @domain_names = $protein_obj->Motif_homol;
  my @domain_starts = $protein_obj->Motif_homol(4);
  my @domain_ends = $protein_obj->Motif_homol(5);

# doesn't return anything
#  my @domain_names = $protein_obj->at('Motif_homol[1]');
#  my @domain_starts = $protein_obj->at('Motif_homol[4]');
#  my @domain_ends = $protein_obj->at('Motif_homol[5]');

# this skips over the repeated domains
#  my @domain_names = $protein_obj->get('Motif_homol');
#  my @domain_starts = $protein_obj->get('Motif_homol', 4);
#  my @domain_ends = $protein_obj->get('Motif_homol', 5);

  # make a hash of hashes, keyed on the name of the interpro/pfam domain
  foreach my $name (@domain_names) {
    # check that there is a start and end position - some INTERPRO entries have this missing
    if (defined $domain_starts[0] && defined $domain_ends[0]) {
      $domains{$name}->{'start'} = shift @domain_starts;
      $domains{$name}->{'end'} = shift @domain_ends;
    }
  }


#
# now get the homologies
#

  # just using constructs like 'my @homol_starts =
  # $protein_obj->Pep_homol(4)' doen't work because you get a list of
  # all of the start positions with no indication as to which belong
  # to which homology.

  # You have to creep up on the data you want by using the $obj->at()
  # structure, in which you can specify exactly the path of tags and
  # data through the object that you wish to traverse and then you can
  # get the list of resulting positions under the sub-path that you
  # want. This means that you must break up the query into a loop
  # looking at every score value of every wormbase homology and then
  # you can get the positions under that score.

  my @homol_names = $protein_obj->Pep_homol;
  foreach my $homol_name (@homol_names) {
    # we only want the homologs if this is the main target - ignore the homologs of the homologs
    # and we only want the homologs from the same species if $group eq 'same'
    if ($type eq "target") { 
      if ($homol_name !~ /^MSP:/) {
	if (($group eq "same" && $homol_name =~ /^$group_start/) ||
	    ($group eq "other" && $homol_name !~ /^$group_start/) ||
	    (exists $species_letters{$group} && $homol_name =~ /$species_letters{$group}/) ||
	    ($group eq "all")
	    ) { 
	  my $max_score = 0;
	  my @scores = $protein_obj->at("Homol.Pep_homol.$homol_name"."[2]");
	  foreach my $score (@scores) {
	    if ($score > $max_score) {
	      $max_score = $score;
	    }
	    $score =~ s/\./\\./;	# quote the decimal points
	    #print "want = Homol.Pep_homol.$homol_name.wublastp_worm.$score"."[1]\n";
	    my @starts = $protein_obj->at("Homol.Pep_homol.$homol_name.wublastp_worm.$score"."[1]");
	    my @ends   = $protein_obj->at("Homol.Pep_homol.$homol_name.wublastp_worm.$score"."[2]");
	    my @other_starts = $protein_obj->at("Homol.Pep_homol.$homol_name.wublastp_worm.$score"."[3]");
	    my @other_ends   = $protein_obj->at("Homol.Pep_homol.$homol_name.wublastp_worm.$score"."[4]");
	    push @{$homologs{$homol_name}->{'start'}}, @starts;
	    push @{$homologs{$homol_name}->{'end'}},   @ends;
	    push @{$homologs{$homol_name}->{'other_start'}}, @other_starts;
	    push @{$homologs{$homol_name}->{'other_end'}},   @other_ends;
	  }
	  $homologs{$homol_name}->{'score'} = $max_score;
	}
      }
    } elsif ($homol_name =~ /^MSP:/) { # mass-spec peptide mapping
      #print "want = Homol.Pep_homol.$homol_name.wublastp_worm.1"."[1]\n";
      my @starts = $protein_obj->at("Homol.Pep_homol.$homol_name.mass-spec.1"."[1]");
      my @ends   = $protein_obj->at("Homol.Pep_homol.$homol_name.mass-spec.1"."[2]");
      my @other_starts = $protein_obj->at("Homol.Pep_homol.$homol_name.mass-spec.1"."[3]");
      my @other_ends   = $protein_obj->at("Homol.Pep_homol.$homol_name.mass-spec.1"."[4]");
      push @{$mass_spec{$homol_name}->{'start'}}, @starts;
      push @{$mass_spec{$homol_name}->{'end'}},   @ends;
      push @{$mass_spec{$homol_name}->{'other_start'}}, @other_starts;
      push @{$mass_spec{$homol_name}->{'other_end'}}, @other_ends;
    }
  }

  $protein_obj->DESTROY();



  return (\%domains, \%homologs, \%mass_spec, $signalp);
}

##########################################
# get the rnai phenotype
#  $data{phenotype} = &get_phenotype($rnai);

sub get_phenotype {
  my ($rnai) = @_;

  my $result = "";
  # $rnai is undef if no phenotype has been scored
  if (! defined $rnai) {return $result;}

  my $rnai_obj = $db->fetch("RNAi" => $rnai); 

  if (! defined $rnai_obj) {die "Can't find an object for $rnai";}

  my @pheno = $rnai_obj->Phenotype(2);

  foreach my $phenotype (@pheno) {
    if ($result eq "" && $phenotype eq "Not") {
      $result = "Not"
    }
    if ($phenotype ne "Not") {
      $result = "OK";
    }
  }

  return $result;
}


##########################################
# get the CDS name of a protein

sub get_cds_name {
  my ($protein) = @_;

  my $protein_obj = $db->fetch("Protein" => $protein);
  my @corr_cds = $protein_obj->Corresponding_CDS;
  foreach my $cds (@corr_cds) {
    if ($cds !~ /\:wp/) {return $cds;}
  }
  print "no valid CDS found for protein $protein\n";
  return undef;
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

#####################################################
# build the contents of the menu bar
 sub build_menubar { 
 
 
     # Create the menubar and File and Quit menubuttons.  Note 
     my $menubar = $MW->Menu; 
     $MW->configure(-menu => $menubar); 
     my $file = $menubar->cascade(-label => '~File'); 
     my $help = $menubar->cascade(-label => '~Help', -tearoff => 0); 
     # Create the menuitems for each menu.  First, the File menu item.
     $file->command(-label => "New", -command => \&new_cds);
     $file->command(-label => "Save As", -command => \&export);
     $file->command(-label => "Quit", -command => \&quit); 

     # Finally, the Help menuitems. 
     $help->command(-label => 'Version'); 
     $help->separator; 
     $help->command(-label => 'About'); 
     my $ver_dialog =   $MW->Dialog(-title   => 'Version', 
                                 -text    => "ilk\n\nVersion $VERSION", 
                                 -buttons => ['OK'], 
                                 );
     my $about_dialog = $MW->Dialog(-title   => 'About ilk', 
                                 -text    => 'ilk is a script to extract information about the homologs of a protein product of a CDS and to display graphical information about those homologs.

more text', 
                                 -buttons => ['OK']); 
     my $menu = $help->cget('-menu'); 
     $menu->entryconfigure('Version', -command => [$ver_dialog   => 'Show']); 
     $menu->entryconfigure('About',   -command => [$about_dialog => 'Show']); 
     return $menubar;
 } # end build_menubar 

#####################################################
# menu command for Export
# write the protein sequence to a file
 sub export { 

     # get file to export to
   my $file = entry_dialog("Save As", "Enter file name");
   if (defined $file && length($file)) {
     #print "filename = $file\n";
     # get protein name of target
     #print "export target protein $data{protein}\n";
     system("pfetch $data{protein} > $file");
     # run mfetch to write it to the file
     # get the protein sequence of each homolog
     foreach my $homolog (keys %{ $data{homologs_href} }) {
       #print "export protein $homolog\n";
       #   run pfetch to append it to the file
       system("pfetch $homolog >> $file");
     }
   }


 } # end export 

#####################################################
# menu command for New
# get a new CDS name and display its details
 sub new_cds { 

   # get CDS name
   my $cds = entry_dialog("New", "Enter CDS name");
   if (defined $cds && length($cds)) {
     #print "CDS = $cds\n";
     %data = &get_data($cds, $group);
     &interface(%data);
   }

 } # end export 
#####################################################
# menu command for Quit
# Exit from the program
 
sub quit { 
  $db->close;			# close the Ace connection
  exit 0; 
} # end quit

#####################################################
# set the interface up
 
sub init { 
  my (%data) = @_;

  $MW->title("ilk $VERSION"); 
  my $menubar = build_menubar; 

  # create the quick-entry bar
  &add_quick_entry(\%data);

  # create the status line
  &add_status_line;
  # display the initial $status line
  &status_line("Status Display");

  &make_canvases(%data);

  # add a scrollbar
  #my $srl_y = $MW -> Scrollbar(-orient=>'v',-command=>[yview => $MW]);
  #$MW->configure(-yscrollcommand=>['set', $srl_y])


} # end init 

#####################################################
# work on the data to produce a set of canvases displaying 
# the target and the homologs

sub make_canvases {
  my (%data) = @_;


  # find the longest homolog/target sequence
  my $max_cds_length = 0;
  foreach my $homolog (keys %{ $data{homologs} }) {
    my $cds_length = 0;
    #print "homolog = $homolog\n";
    my $hom_exons_start_ref = $data{homologs}{$homolog}{exons_start_ref};
    my $hom_exons_end_ref = $data{homologs}{$homolog}{exons_end_ref};
    for (my $i = 0; $i < @{$hom_exons_start_ref}; $i++) {
      my $len = $hom_exons_end_ref->[$i] - $hom_exons_start_ref->[$i] + 1;
      $cds_length += $len;
      #print "\tExon ",$i+1,": $hom_exons_start_ref->[$i]..$hom_exons_end_ref->[$i] length: $len\n";
    }
    if ($cds_length > $max_cds_length) {
      $max_cds_length = $cds_length;
    }
  }

  #print "max length of homologs in bases = $max_cds_length\n";

  # get the length of the target
  my $exons_start_ref = $data{exons_start_ref};
  my $exons_end_ref = $data{exons_end_ref};

  my $cds_length = 0;
  for (my $i = 0; $i < @{$exons_start_ref}; $i++) {
    my $len = $exons_end_ref->[$i] - $exons_start_ref->[$i] + 1;
    $cds_length += $len;
    #print "Exon ",$i+1,": $exons_start_ref->[$i]..$exons_end_ref->[$i] length: $len\n";
  }
  if ($cds_length > $max_cds_length) {
    $max_cds_length = $cds_length;
  }
  #print "max length of homologs & target in bases = $max_cds_length\n";

  # get the scaling factor
  my $scale = 1800/$max_cds_length;
  #print "scale = $scale\n";

  my %colours;			# holds randomly generated colours keyed on the names of the domains

  # set up the canvas for the target
  my ($target_canv, %homologous_regions) = &add_canvas($scale, "target", \%data, \%colours, undef);

  # set up the Frame for the homologs
  # this frame wraps up all of the homologs and is scrollable
    use Tk::Pane;
    my $scrolled = $MW->Scrolled("Frame", 
				 -scrollbars => "osoe",
				 -height => 500,
				 -width  => 750,
				 )->pack(-side => 'top',
					 -fill => 'x',
					 );



  # set up the canvases for the homologs
  # sort them by descending Blast scores
  foreach my $homolog (sort { $data{homologs}{$b}{score} <=>  $data{homologs}{$a}{score} } keys %{ $data{homologs} }) {
      &add_canvas($scale, $homolog, $data{homologs}{$homolog}, \%colours, $scrolled, $target_canv, $homologous_regions{$homolog}, $data{domain_desc_href});
  }

}


#####################################################
# make the canvas for a target or homolog

sub add_canvas {
  # parameters:
  # $target_canv = canv object created by the target, here so that homologs can raise/lower the homologous regions
  # $target_homologous_regions_aref = aref of this homolog's homologous region rectangle ojects in target's canvas
  my ($scale, $id, $info, $colours, $scrolled, $target_canv, $target_homologous_regions_aref, $target_domain_desc_href) = @_;

  my %info = %{$info};		# the bit of the %data global variable dealing specifically with this homolog
  my $have_target = 0;		# flag set true if we have the target CDS


  my $protein;
  my %domain_desc;

  if ($id eq 'target') {
    $have_target = 1;
    $id = $info{cds_name};	# NB only the target has this key/value pair
    $protein = $info{protein};
    %domain_desc = %{$info{domain_desc_href}};
  } else {
    $protein = $id;
    $id = $info{cds_name_of_homolog};
    %domain_desc = %{$target_domain_desc_href};
  }
  my $gene = $info{gene};
  my $cgc_name = $info{cgc_name};
  my $clone = $info{clone};
  my $lab = $info{lab};
  my $phenotype = $info{phenotype};
  my $cds_start = $info{cds_start};
  my $cds_end = $info{cds_end};
  my $exons_start_ref = $info{exons_start_ref};
  my $exons_end_ref = $info{exons_end_ref};
  my $domains_href = $info{domains_href};
  my $mass_spec_href = $info{mass_spec_href};
  my $homologs_href = $info{homologs_href};
  my $signalp = $info{signalp};
  my $score = $info{score};
  if (! defined $score) {$score = 100.0;} # the target doesn't have a score to itself, so make a dummy one

  # top holds all of the stuff created here
  my $top;
  if (defined $scrolled) {	# NB only the target has this undefined
    $top = $scrolled->Frame(-background => 'gray')->pack(-side => 'top',
							    -fill => 'x'); 
  } else {
    $top = $MW->Frame(-background => 'gray')->pack(-side => 'top',
							    -fill => 'x'); 
  }

  # count the number of domains so we know how much space is required
  my $domain_count = keys %{$domains_href};

  # count one other line for the mass-spec
  if (keys %{$mass_spec_href}) {$domain_count++;}

  # canv is the canvas that we will draw on
  my $canv = $top->Canvas(-background => 'gray',
			  -width  => 850,
			  -height => 50 + ($domain_count * 5),
			  )->pack(-side => 'left',
				  -expand => 1,
				  -fill => "x",
				  );

  # display some information about the gene
  my $lab1 = $canv->createText(5, 10, -text => $id, -anchor => 'w');

  $canv->bind($lab1,		# display some information about the gene on the status line
	      			# and raise the homologous region in the target's canvas
	      '<Enter>',
              sub {
		&status_line("$id Lab: $lab Gene: $gene Protein: $protein Homology blast score: $score");
		foreach my $h (@{ $target_homologous_regions_aref }) {
		  $target_canv->raise($h);
		}
	      }
	      );
  $canv->bind($lab1,		# lower the homologous region in the target's canvas
	      '<Leave>',
              sub {
		foreach my $h (@{ $target_homologous_regions_aref }) {
		  $target_canv->lower($h);
		}
	      }
	      );
  $canv->bind($lab1,		# make the acedb FMAP display the gene
	      '<Button-1>',
              sub {
		&goto_location($clone, $cds_start, $cds_end, '+', 10);
	      }
	      );

  # display the lab 
  my $lab_colour = $lab eq "HX" ? "darkblue" : "darkgreen";
  my $lab2 = $canv->createText(5, 25, -text => $lab, -anchor => 'w', -fill => $lab_colour);

  # display the phenotype
  my $pheno_colour = $phenotype eq "Not" ? "darkred" : "blue";
  if ($phenotype ne "") {
    my $lab3 = $canv->createText(5, 40, -text => "Phenotype $phenotype", -anchor => 'w', -fill => $pheno_colour);
  }
  
  if (defined $info{cgc_name}) {
    my $lab4 = $canv->createText(80, 10, -text => "$info{cgc_name}", -anchor => 'w', -fill => 'red');
  }

  # display the signalp region
  if ($signalp) {
    my $sig = $canv->createRectangle(95, 20, 
				     105, 35,
				     -fill => 'white',
				     );
    $canv->bind($sig,
		'<Enter>',
		sub {
		  &status_line("$id Lab: $lab Gene: $gene Protein: $protein has a SignalP region");
		});
  }

  # display the exons
  my $x = 100;
  my $y = 25;
  my $height = 5;
  for (my $i = 0; $i < @{$exons_start_ref}; $i++) {
    my $len = $exons_end_ref->[$i] - $exons_start_ref->[$i] + 1;
    #print "\tExon ",$i+1,": $hom_exons_start_ref->[$i]..$hom_exons_end_ref->[$i] length: $len\n";

    my $aa_len = ($len/3) * $scale; # get length in residues and scale to fit

    my $exon = $canv->createRectangle($x, $y, 
				      $x+$aa_len, $y+$height,
				      -fill => 'blue',
				      -tags => "$id.$i"); # tag made from ID and exon number

    my $exon_count = $i + 1;
    $canv->bind("$id.$i",
		'<Enter>',
		sub {
		  &status_line("$id Lab: $lab Gene: $gene Protein: $protein Exon: $exon_count Length: $len bp");
		});

    $canv->bind("$id.$i",		
	      '<Button-1>',
              sub {
		&status_line("Going to: $id");
		$input = $id;
		&quick_next_cds;		
	      }
	      );


    # mark the exon boundaries
    if ($i != 0) {
      my $exon_mark = $canv->createLine($x, $y-2,
					$x, $y+$height+3,
					-fill => 'green',
					-width => 2);
    }

    $x += $aa_len;		# get the start of the next exon
  }

  my %homologous_regions;	# hash of arrays of objects of homologous regions in target

  if (! $have_target) {		# i.e if this is a homolog
    # display the regions of homology - and gaps in homology
    # get the homologous regions of the target to this homolog
    my @starts       = @{ $info{homology_starts} };
    my @ends         = @{ $info{homology_ends} };
    my @other_starts = @{ $info{homology_other_starts} };
    my @other_ends   = @{ $info{homology_other_ends} };

    foreach my $start (@starts) {
      my $end = shift @ends;
      my $other_start = shift @other_starts;
      my $other_end = shift @other_ends;
      #print "\tHomology: $start..$end homolog_positions: $other_start..$other_end\n";

      my $height = 3;
      my $y = 26;
      my $x = 100 + ($other_start * $scale);
      my $len = $other_end - $other_start + 1;
      my $aa_len = $len * $scale; # get length in residues and scale to fit
      my $homologous_region = $canv->createRectangle($x, $y, 
						     $x+$aa_len, $y+$height,
						     -fill => 'red');
      $canv->bind($homologous_region,		
		  '<Button-1>',
		  sub {
		    &status_line("Going to: $id");
		    $input = $id;
		    &quick_next_cds;		
		  }
		  );

    }
  } else {
    # this is the target
    # display the regions of homology of this target to the homologs
    # and then lower the display so that it is not visible until it is raise
    # by a mouse event
    foreach my $homolog (keys %{ $data{homologs} }) {
      my @target_starts = @{ $data{homologs}{$homolog}{homology_starts} };
      my @target_ends = @{ $data{homologs}{$homolog}{homology_ends} };
      foreach my $target_start (@target_starts) {
	my $target_end = shift @target_ends;
	my $height = 3;
	my $y = 26;
	my $x = 100 + ($target_start * $scale);
	my $len = $target_end - $target_start + 1;
	my $aa_len = $len * $scale; # get length in residues and scale to fit
	my $homologous_region = $canv->createRectangle($x, $y, 
						     $x+$aa_len, $y+$height,
						     -fill => 'red');
	$canv->lower($homologous_region);
	push @{ $homologous_regions{$homolog} }, $homologous_region;

      }
    }
  }


  # display the mass-spec
  my $count = 0;
  foreach my $mass_spec (keys %{$mass_spec_href}) {
    #print "$mass_spec $mass_spec_href->{$mass_spec}{'start'} $mass_spec_href->{$mass_spec}{'end'}\n";
    my $height = 4;
    my $y = 32 + ($count * 5);
    my $x = 100 + (@{$mass_spec_href->{$mass_spec}{'start'}}[0] * $scale);
    my $len = @{$mass_spec_href->{$mass_spec}{'end'}}[0] - @{$mass_spec_href->{$mass_spec}{'start'}}[0] + 1;
    my $aa_len = $len * $scale; # get length in residues and scale to fit

    my $colour = "#FFCC22";

    my $mass_spec_obj = $canv->createRectangle($x, $y, 
					$x+$aa_len, $y+$height,
					-fill => "$colour");
    $canv->bind($mass_spec_obj,
		'<Enter>',
		sub {
		  &status_line("Mass_Spec: $mass_spec");
		});
  }

  if (keys %{$mass_spec_href}) {$count++;} # use only one line to display the mass-spec

  # display the domains
  foreach my $domain (keys %{$domains_href}) {
    #print "$domain $domains_href->{$domain}{'start'} $domains_href->{$domain}{'end'}\n";
    my $height = 4;
    my $y = 32 + ($count * 5);
    my $x = 100 + ($domains_href->{$domain}{'start'} * $scale);
    my $len = $domains_href->{$domain}{'end'} - $domains_href->{$domain}{'start'} + 1;
    my $aa_len = $len * $scale; # get length in residues and scale to fit

    my $colour;
    if (exists $colours->{$domain}) {
      $colour = $colours->{$domain};
    } else {
      # construct a colour from the domain ID name
      my @colour;
      $colour[0] = 0;
      $colour[1] = 0;
      $colour[2] = 0;
      foreach (my $i = 0; $i < length $domain; $i++) {
	my $chr = substr($domain, $i, 1);
	if ($chr =~ /\d/) {
	  $colour[$i % 3] .= $chr;
	}
      }
      foreach my $c (@colour) {
	if ($c < 12) {$c *= 10;}
      }
      $colour = sprintf( "#%02x%02x%02x", $colour[0]+127, $colour[1]+127, $colour[2]+127); # light random colour
      $colours->{$domain} = $colour;
    }

    my $domain_obj = $canv->createRectangle($x, $y, 
					$x+$aa_len, $y+$height,
					-fill => "$colour");
    my $domain_desc = $domain_desc{$domain};
    if (!defined $domain_desc) {$domain_desc = "";}
    $canv->bind($domain_obj,
		'<Enter>',
		sub {
		  &status_line("Domain: $domain\t$domain_desc");
		});

    $count++;
  }

  return ($canv, %homologous_regions);
}

#####################################################
# add the quick entry frame and label
sub add_quick_entry {
  my $data = shift;

  my $quick = $MW->Frame(-borderwidth => 3,
			-relief => 'groove',
			-background => 'lightyellow')->pack(-side => 'top',
						     -fill => 'x');
  $quick->Label( 
	    -text => 'Next CDS: ',
	    -background => 'lightyellow')->pack(-side => 'left');		

  # entry field
  $quick_entry = $quick->Entry( 
				   -width => 25,
				   -textvariable=> \$input,
				   -background => 'grey')->pack(-side => 'left');
  
  # make Return and Enter submit CDS 
  $quick_entry->bind("<Return>",[ \&quick_next_cds]);
  $quick_entry->bind("<KP_Enter>",[ \&quick_next_cds]);

  # OK button
  my $make = $quick->Button( -text => "OK",
			     -command => [\&quick_next_cds],
			     -background => 'lightyellow',
			     )->pack(-side => 'left',
				     -pady => '1',
				     -padx => '6',
				     -anchor => "w"
				     );

  # compare entry field - compare this CDS's protein to all the other
  # homologs and compare the target protein to all the homologs to see
  # which one is best
  
  $quick->Label( 
	    -text => ' Compare to CDS: ',
	    -background => 'lightyellow')->pack(-side => 'left');		

  # entry field
  my $compare_entry = $quick->Entry( 
				   -width => 25,
				   -textvariable=> \$compare_cds,
				   -background => 'grey')->pack(-side => 'left');
  
  # make Return and Enter submit CDS 
  $compare_entry->bind("<Return>",[ \&compare_cds, $data]);
  $compare_entry->bind("<KP_Enter>",[ \&compare_cds, $data]);

  # OK button
  $quick->Button( -text => "OK",
		     -command => [\&compare_cds, $data],
		     -background => 'lightyellow',
		     )->pack(-side => 'left',
			     -pady => '1',
			     -padx => '6',
			     -anchor => "w"
			     );


  # groups pop-up menu
  my $groups=$quick->Menubutton(-text=>'Groups',
  -background=>'lightyellow',
  -activebackground=>'grey')->pack(-side=>'right');
 

  $groups->command(-label=>'Same species',
		   -background=>'lightyellow',
		   -command=>[\&change_group, 'same']);
  $groups->command(-label=>'Other species',
		   -background=>'lightyellow',
		   -command=>[\&change_group, 'other']);
  $groups->command(-label=>'All species',
		   -background=>'lightyellow',
		   -command=>[\&change_group, 'all']);
  $groups->separator();
  $groups->command(-label=>'elegans only',
		   -background=>'lightyellow',
		   -command=>[\&change_group, 'elegans']);
  $groups->command(-label=>'briggsae only',
		   -background=>'lightyellow',
		   -command=>[\&change_group, 'briggsae']);
  $groups->command(-label=>'remanei only',
		   -background=>'lightyellow',
		   -command=>[\&change_group, 'remanei']);
 


}




##############################################################
sub quick_next_cds {
  my $cds = $input;
  if (defined $cds && length($cds)) {
    print "displaying CDS = $cds\n";
    %data = &get_data($cds, $group);
    &interface(%data);
  }
}
##############################################################
sub change_group {
  ($group) = @_; 
  print "Group = $group\n";
  &quick_next_cds;
}

##############################################################
#
# Compare the target gene to another CDS by taking their protein
# sequences and for both of them doing a global alignment against each
# of the homolog proteins. The scores of the alignments are displayed.

sub compare_cds {
  my $data_href = shift;
  my $data = %{$data_href};

  if (defined $compare_cds && length($compare_cds)) {

    # get the target protein data
    my $target_cds = $data{cds_name};	# NB only the target has this key/value pair
    my $target_protein = $data{protein};
    my $target_obj = $db->fetch(Protein => $target_protein);
    my $target_seq = $target_obj->asPeptide;
    my $target_file = "/tmp/target_${target_cds}";
    open (TARGET, ">$target_file") || die "Can't open $target_file\n";
    print TARGET "$target_seq";
    close(TARGET);

    # get the data for the CDS to compare the target to
    my $compare_obj = $db->fetch(CDS => $compare_cds);
    my $compare_seq = $compare_obj->asPeptide;
    my $compare_file = "/tmp/compare_${compare_cds}";
    open (COMPARE, ">$compare_file") || die "Can't open $compare_file\n";
    print COMPARE "$compare_seq";
    #print "$compare_seq";
    close(COMPARE);

    print "comparing target $target_cds to the CDS $compare_cds\n\n";
    
    my $result_file = "/tmp/${target_cds}_comparison";
    open (RES, ">$result_file") || die "Can't open $result_file\n";
    foreach my $homolog (sort { $data{homologs}{$b}{score} <=>  $data{homologs}{$a}{score} } keys %{ $data{homologs} }) {
      my $homolog_cds = $data{homologs}{$homolog}{cds_name_of_homolog};
      my $obj = $db->fetch(Protein => $homolog);
      my $hom_file = "/tmp/homolog_${homolog_cds}";
      open (HOM, ">$hom_file") || die "Can't open $hom_file\n";
      print HOM $obj->asPeptide . "\n";
      close(HOM);
      my $target_score;
      my $compare_score;
      open (NEEDLE, "needle $target_file $hom_file -auto -stdout|");
      while (my $line = <NEEDLE>) {
	if ($line =~ /Score:\s+([\d\.]+)/) {$target_score = $1;}
	print RES $line;
      }
      open (NEEDLE, "needle $compare_file $hom_file -auto -stdout|");
      while (my $line = <NEEDLE>) {
	if ($line =~ /Score:\s+([\d\.]+)/) {$compare_score = $1;}
	print RES $line;
      }
      unlink $hom_file;

      # report the results
      my $difference = $compare_score - $target_score;
      print "Homolog: $homolog_cds $homolog\ttarget: $target_score\tcompare: $compare_score\tdiff: $difference\n";

    }
    close(RES);
    print "\nAlignments are in $result_file\n";

    # tidy up
    unlink $target_file;
    unlink $compare_file;

  }
}
#####################################################
# add the status frame and label
sub add_status_line {
  my $line = $MW->Frame(-borderwidth => 3,
			-relief => 'groove',
			-background => 'cyan')->pack(-side => 'top',
						     -fill => 'x');

  $status = $line->Entry(-background => 'cyan',
			 -width => 120)->pack(-side => 'left',
					       -fill => 'x'); # $status is a global variable
			
}

#####################################################
# set the status line to display the argument
sub status_line {
  #$status->configure(-text => $_[0]);
  $status->delete(0, 'end');
  $status->insert(0, $_[0]);
}

############################################################################
# make the acedb FMAP display the specified location
sub goto_location
  {
    my $seqobj  = shift;	# clone or superlink name
    my $x       = shift;	# start position
    my $y       = shift;	# end position
    my $sense   = shift;	# sense '+' or '-'
    my $zoomout = shift;	# about of sequence to add around the edges

    if (!defined $zoomout) {
      $zoomout = 500;
    }
    $x -= $zoomout;
    $y += $zoomout;

    my $reversed_command = "";
    if (defined $sense && $sense eq '-') {
      $reversed_command = "; seqactions -rev_comp";
    }
    my $command = "xremote -remote 'gif seqget $seqobj -coords $x $y; seqdisplay $reversed_command'";

    print "command=$command\n";
    my $return_status = system($command);
    if ( ( $return_status >> 8 ) != 0 ) {
      &error_warning("WARNING", "X11 connection appears to be lost");
    }

  }

############################################################################
# display a dialog box to inform the user of an error
sub error_warning
  {
    my $title = shift;
    my $text = shift;
    my $D = $MW->Dialog(
			      -title => "$title",
			      -text  => "$text",
			      -default_button => 'OK',
			      -buttons        => ['OK']
			     );
    $D->configure(-background => 'red',
		  -foreground => 'black');

    my $choice = $D->Show;
  }
############################################################################
# display a dialog box to get some text from the user
sub entry_dialog
  {
    my ($title, $label) = @_;

    my $D = $MW->Dialog(
			-title => "$title",
			-default_button => 'OK',
			-buttons        => ['OK', 'Cancel']
			);
    $D->configure(-background => 'lightyellow',
		  -foreground => 'black');
    $D->add("Label", 
	    -text => $label,
	    -background => 'lightyellow')->pack();
    my $entry = $D->add("Entry", 
			-width => 25,
			-background => 'lightyellow')->pack();
    my $choice = $D->Show;
    if ($choice eq "OK") {
      return $entry->get;
    } else {
      return undef;
    }
  }
