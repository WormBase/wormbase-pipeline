#!/software/bin/perl -w
#
# Make enhancer region Features.
#
# Specify the WBPaper ID
# Specify a gene that is affected (by WBGene, sequence name of CGC name).
# Get the regions by either an offset from the START codon of the specified gene  or primer sequences or both (they check each other).
# Get any overlapping pre-existing Features and report them in case this region has been worked on already.
# Get the type of region (promoter, regulatory_region, TF_binding_site_region, TF_binding_site, silencer, enhancer, nc_conserved_region)
# make the Feature by populating the tags:
#
#(WBsf ID) ?Feature
#(clone)                          Sequence UNIQUE ?Sequence XREF Feature_object
#(clone)                          Mapping_target UNIQUE ?Sequence 
#(flanking sequences) Flanking_sequences UNIQUE Text UNIQUE Text
#(name it is called in the paper) Name Public_name UNIQUE ?Text
#(C. elegans)                     Species UNIQUE ?Species 
#(short description)              Description ?Text
#(WBPaper ID)                     Defined_by_paper    ?Paper XREF Feature #Evidence
#(Gene)                           Associated_with_gene       ?Gene       XREF Associated_feature #Evidence 
#(TF)                             Bound_by_product_of ?Gene XREF Gene_product_binds #Evidence 
#(TF)                             Associated_with_transcription_factor UNIQUE ?Transcription_factor XREF Binding_site 
#(Caltech)                        Associated_with_expression_pattern         ?Expr_pattern XREF Associated_feature #Evidence
#(date tagged remark)             Remark ?Text #Evidence
#(enhancer, silencer, promoter, TF_binding_site, binding_site) Method UNIQUE ?Method
#
#Feature : ""
#Sequence 
#Mapping_target 
#Flanking_sequences 
#DNA_text
#Description "."
#Remark ". [2013-07-23 gw3]"
#Associated_with_gene 
# Bound_by_product_of 
# Associated_with_transcription_factor 
#Method  enhancer
#        TF_binding_site
#        regulatory_region
#        promoter 
#SO_term "SO:0000165" // enhancer
#        "SO:0000235" // TF_binding_site
#        "SO:0005836" // regulatory_region
#        "SO:0000167" // promoter 
#        "SO:0000334" // nc_conserved_region
#        "SO:0000409" // binding_site
#Defined_by_paper 
#Public_name 


# 
# by Gary Williams                     
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2014-09-02 14:09:14 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Carp;
use Coords_converter;
use Modules::Remap_Sequence_Change;
use Sequence_extract;
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use Ace;
use Modules::Overlap;
use Feature_mapper;



my ($help, $debug, $test, $verbose, $store, $wormbase, $output, $gene_name, $paper, $feature_id);


GetOptions ("help"         => \$help,
            "debug=s"      => \$debug,
            "test"         => \$test,
            "verbose"      => \$verbose,
            "store:s"      => \$store,
	    "paper:s"      => \$paper,
	    "gene:s"       => \$gene_name,
	    "feature_id:s" => \$feature_id,
	    "output:s"     => \$output,
            );

$debug = $ENV{USER};
$test = 0;

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


if (! defined $output) {$log->log_and_die ("No -output acefile specified\n");}
if (! defined $paper) {$log->log_and_die ("No -paper WBPaperID specified\n");}
if (! defined $gene_name) {$log->log_and_die ("No -gene genename specified\n");}
#if (! defined $feature_id) {$log->log_and_die ("No -feature_id WBsfID specified\n");}


my $currentdb = $wormbase->database('current');

print "connect to Ace\n";
my $db = Ace->connect(-path => $currentdb) || die "cannot connect to database at $currentdb\n";

print "init coords converter\n";
my $coords = Coords_converter->invoke($currentdb, 0, $wormbase);

print "init sequence extract\n";
my $seq_obj = Sequence_extract->invoke($currentdb, undef, $wormbase);

print "init feature mapper\n";
my $mapper = Feature_mapper->new($currentdb, undef, $wormbase);

print "load common data\n";
my %seq2gene = $wormbase->FetchData("worm_gene2geneID_name", undef, "$currentdb/COMMON_DATA"); # 'AC3.3' => 'WBGene00000024'
my %cgc2gene = $wormbase->FetchData("cgc_name2gene", undef, "$currentdb/COMMON_DATA");         # 'abu-1' => 'WBGene00000024',

# open an ACE connection to parse details for mapping to genome



# set up defaults
my %so_term = (
	       TF_binding_site     => "SO:0000235", # A region of a nucleotide molecule that binds a Transcription Factor or Transcription Factor complex 
	       enhancer            => "SO:0000165", # A cis-acting sequence that increases the utilization of (some) eukaryotic promoters, and can function in either orientation and in any location (upstream or downstream) relative to the promoter.
	       silencer            => "SO:0000625", # A regulatory region which upon binding of transcription factors, suppress the transcription of the gene or genes they control.
	       regulatory_region   => "SO:0005836", # A region of sequence that is involved in the control of a biological process.
	       promoter            => "SO:0000167", # A regulatory_region composed of the TSS(s) and binding sites for TF_complexes of the basal transcription machinery.
	       nc_conserved_region => "SO:0000334", # Non-coding region of sequence similarity by descent from a common ancestor.
	       binding_site        => "SO:0000409", # A biological_region of sequence that, in the molecule, interacts selectively and non-covalently with other molecules.
	      );

print "get transcription factors\n";
my ($TFs, $TF_genes) = get_TFs();
my %TFs = %{$TFs};
my %TF_genes = %{$TF_genes};

print "get gene ID\n";
my $gene_id = get_gene_id($gene_name);

print "get gene location\n";
my ($chromosome, $gene_start, $gene_sense) = get_gene_start($gene_id, $gene_name);



open (ACE, "> $output") || die "Can't open $output\n";

while (1) {
  print "\n\nNext Feature.\n\n";

  # get the user to indicate the location
  my ($quit, $start, $end, $sense);
  ($quit, $chromosome, $start, $end, $sense) = get_location($gene_id, $chromosome, $gene_start, $gene_sense);
  if ($quit) {last}

  # get flanking sequence
  my ($target_sequence, $sequence_start, $sequence_end, $left_flank, $right_flank) = get_flanking_sequences($chromosome, $start, $end, $sense);
  
  # get sequence if less than 100 bp long
  my $sequence_text = get_sequence_text($chromosome, $start, $end, $sense);

  # get details
  my ($method, $public_name, $desc, $TFID, $TF_gene, $TF_name, $remark) = get_details($gene_name, $TFs, $TF_genes);

  # check for overlapping features
  check_overlapping_features($method, $target_sequence, $sequence_start, $sequence_end, $sense);
  
  # get feature ID
  $feature_id = get_feature_id($feature_id);

  # write feature
  write_feature($target_sequence, $sequence_start, $sequence_end, $left_flank, $right_flank, $method, $public_name, $desc, $TFID, $TF_gene, $TF_name, $remark, $sequence_text);
}

close (ACE);

# close the ACE connection
$db->close;


$log->mail();
print "Finished.\n";
exit(0);



##############################################################
#
# Subroutines
#
##############################################################



##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################################################

sub write_feature {
  my ($target_sequence, $sequence_start, $sequence_end, $left_flank, $right_flank, $method, $public_name, $desc, $TFID, $TF_gene, $TF_name, $remark, $sequence_text) = @_;

  if (defined $public_name) {$public_name =~ s/,//;}
  if (defined $desc) {$desc =~ s/,//;}
  if (defined $remark) {$remark =~ s/,//;}

  print ACE "\n";
  print ACE "Feature $feature_id\n";
  print ACE "Sequence $target_sequence\n";
  print ACE "Mapping_target $target_sequence\n";
  print ACE "Flanking_sequences $left_flank $right_flank\n";
  print ACE "Public_name \"$public_name\"\n";
  print ACE "DNA_text $sequence_text\n" if ($sequence_text ne '');
  print ACE "Species \"Caenorhabditis elegans\"\n";
  print ACE "Description \"$desc\"\n" if (defined $desc && $desc ne '');
  print ACE "Method $method\n";
  print ACE "SO_term \"",$so_term{$method},"\"\n";
  print ACE "Defined_by_paper \"$paper\"\n";
  print ACE "Associated_with_gene $gene_id\n";
  print ACE "Associated_with_transcription_factor $TFID\n" if (defined $TFID && $TFID ne '');
  foreach my $gene (@{$TF_gene}) {
    print ACE "Bound_by_product_of $gene\n" if (defined $gene);
  }
  print ACE "Remark \"$remark\"\n" if (defined $remark && $remark ne '');
  print ACE "\n";
  print ACE "Sequence : \"$target_sequence\"\n";
  print ACE "Feature_object \"$feature_id\" $sequence_start $sequence_end\n";
  print ACE "\n";
  $feature_id++;
    
}

##########################################################################

sub get_flanking_sequences {
  my ($chromosome, $start, $end, $sense) = @_;

# get the chromosome location of the Feature
# map the Feature onto the chromosome
# get the flanking sequence start & end and find the minimal sequence object that will allow these positions to map to it (LocateSpan) and return this

  # returned data
  my ($left_flank, $right_flank);


  # get the left and right flanking sequences for this position
  # fourth argument says this is not a zero-length feature
  # seventh argument allows non-unique short sequences at the contig sequence ends
  my $is_zero_length = 0;
  my $min_flank_length = 30;
  my $no_unique_check = 0;
  my $short_end_flanks_allowed = 0; # don't allow truncated flanks at the end of the clone

  my @clone_coords = $coords->LocateSpan($chromosome, $start-$min_flank_length, $end+$min_flank_length ); 
  my ($target_clone, $clone_start, $clone_end) = @clone_coords;
  $clone_start += $min_flank_length;
  $clone_end -= $min_flank_length;

  # use reverse sense coords
  if ($sense eq '-') {
    ($clone_end, $clone_start) = ($clone_start, $clone_end);
    $clone_end--;
    $clone_start--;
  }

  ($left_flank, $right_flank) = $mapper->get_flanking_sequence_for_feature($target_clone, $clone_start, $clone_end, $is_zero_length, $min_flank_length, $no_unique_check, $short_end_flanks_allowed);

  if (defined $left_flank && defined $right_flank) {
    return ($target_clone, $clone_start, $clone_end, $left_flank, $right_flank)
  }
   


  # if the search for unique flanks in the clone fails, search in the whole chromosome
  # this can be quite slow, which is why we look for possible unique flanks in the clone first

  # get the left and right flanking sequences for this position
  # fourth argument says this is not a zero-length feature
  # seventh argument allows non-unique short sequences at the contig sequence ends
  $is_zero_length = 0;
  $min_flank_length = 30;
  $no_unique_check = 0;
  $short_end_flanks_allowed = 1; # allow truncated flanks at the end of the chromosoem/contig
  ($left_flank, $right_flank) = $mapper->get_flanking_sequence_for_feature($chromosome, $start, $end, $is_zero_length, $min_flank_length, $no_unique_check, $short_end_flanks_allowed); 
  
  if (! defined $left_flank) {
    $log->write_to( "ERROR when trying to find flanking regions at $chromosome, $start, $end\n");
    return (undef, undef, undef);
  }

  return ($chromosome, $start, $end, $left_flank, $right_flank);
    
}
##########################################################################
# get the stuff to describe this region
sub  get_details {
  my ($gene_name, $TFs, $TF_genes) = @_;

  my $method                     = get_method(); # Method
  my $public_name                = get_public_name(); # name it is called in the paper
  my ($TFID, $TF_gene, $TF_name) = get_TF($method, $TFs, $TF_genes); # TF ID, gene ID, name
  my $desc                       = get_desc($method, $TF_name, $public_name, $gene_name); # short description
  my $remark                     = get_remark($method, $TF_name, $public_name); # date tagged remark

  return ($method, $public_name, $desc, $TFID, $TF_gene, $TF_name, $remark);
}
##########################################################################
# gets the Method

sub get_method {

  my $candidate = "enhancer";
  my $method;
  my @methods = ('enhancer', 'silencer', 'TF_binding_site', 'regulatory_region', 'nc_conserved_region', 'binding_site', 'promoter');

  do {
    print "Method [$candidate] (",join ', ',@methods,") > ";
    my $input =  <STDIN>;
    chomp ($input);
    $input =~ s/,//;
    if ($input eq '') {
      $method = $candidate;
    } else {
      if (grep /$input/, @methods) {
	$method = $input;
      } else {
	print "WARNING: '$input' is an unknown Method\n";
      }
    }
  } until (defined $method);

  return $method;
}
##########################################################################
# gets the public name - name it is called in the paper
sub get_public_name {
  my $public_name;

  my $candidate = "";

  do {
    print "Public name of region in paper (e.g. 'enh1', '350bp enhancer', 'CES-1 site') > ";
    my $input =  <STDIN>;
    chomp ($input);
    if ($input eq '') {
	print "WARNING: no public name entered.\n";
    } else {
	$public_name = $input;
    }
  } until (defined $public_name);

  return $public_name;
}
##########################################################################
# get TF ID, gene ID, name
sub get_TF {
  my ($method) = @_;

  if ($method ne 'TF_binding_site') {return undef}

  my %TFs = %{$TFs};
  my %TF_genes = %{$TF_genes};

  my ($TFID, $TF_name, $TF_genes);

  my $candidate = '';

  do {
    if ($candidate eq '') {
      print "TF name > ";
    } else {
      print "TF name [$candidate] > ";      
    }
    my $input =  <STDIN>;
    chomp ($input);
    if ($input eq '') {
      if ($candidate ne '') {
	$TF_name = $candidate;
	$TFID = $TFs{$candidate};
	$TF_genes = $TF_genes{$candidate};
      }
    } else {

      if (exists $TFs{$input}) {
	$TFID = $TFs{$input};
	$TF_name = $input;
	$TF_genes = $TF_genes{$input}; # ref to array of gene IDs
      } else {
	print "WARNING: '$input' is an unknown TF\n";
	# make a guess based on what the user input
	print "Did you mean any of the following?\n";
	
	# see if it is a partial TF name
	foreach my $TF_name (sort keys %TFs) { 
	  if ($TF_name !~ /$input/) {next}
	  $candidate = $TF_name;
	  print "\t$TF_name\n";
	}
	
	# see if it is a gene CGC name that matches a TF name
	if (exists $cgc2gene{$input}) {
	  my $gene_id = $cgc2gene{$input};
	  foreach my $TF_name (sort keys %TF_genes) {
	    if (!grep /$gene_id/, @{$TF_genes{$TF_name}}) {next}
	    $candidate = $TF_name;
	    print "\t$TF_name (@{$TF_genes{$TF_name}})\n";
	  } 
	}
	
	# see if it is a gene sequence name that matches a TF name
	if (exists $seq2gene{$input}) {
	  my $gene_id = $cgc2gene{$input};
	  foreach my $TF_name (sort keys %TF_genes) {
	    if (!grep /$gene_id/, @{$TF_genes{$TF_name}}) {next}
	    $candidate = $TF_name;
	    print "\t$TF_name (@{$TF_genes{$TF_name}})\n";
	  } 
	}
	
	# see if it is a WBGene ID
	foreach my $TF_name (sort keys %TF_genes) { 
	  if (!grep /$input/, @{$TF_genes{$TF_name}}) {next}
	  $candidate = $TF_name;
	  print "\t$TF_name (@{$TF_genes{$TF_name}})\n";
	}
      }
    }
  } until (defined $TFID);

  return ($TFID, $TF_genes, $TF_name);

}
##########################################################################
# short description
sub get_desc {
  my $desc;
  my ($method, $TF_name, $public_name, $gene_name) = @_;

  my $candidate='';
  if ($method eq 'TF_binding_site') {
    $candidate = "This is a $TF_name transcription factor binding site for $gene_name.";
  } elsif ($method eq 'binding_site') {
    $candidate = "This is the '$public_name' binding site for $gene_name.";
  } elsif ($method eq 'enhancer') {
    $candidate = "This is the '$public_name' enhancer region for $gene_name.";
  } elsif ($method eq 'silencer') {
    $candidate = "This is the '$public_name' silencer region for $gene_name.";
  } elsif ($method eq 'regulatory_region') {
    $candidate = "This is the '$public_name' regulatory region for $gene_name.";
  } elsif ($method eq 'nc_conserved_region') {
    $candidate = "This is the conserved '$public_name' region for $gene_name.";
  } elsif ($method eq 'promoter') {
    $candidate = "This is the '$public_name' promoter region for $gene_name.";
  }

  do {
    print "Description [$candidate] > ";
    my $input = <STDIN>;
    chomp ($input);
    if ($input eq '') {
      $desc = $candidate;
    } else {
      $desc = $input;
    }
  } until (defined $desc);

  return $desc;
}
##########################################################################
# next Feature ID
sub get_feature_id {
  my ($feature_id) = @_;
  my $candidate = $feature_id;
  my $valid = 0;

  do {
    if (defined $candidate && $candidate ne '') {
      print "Feature ID [$candidate] > ";
    } else {
      print "Feature ID > ";
    }
    my $input = <STDIN>;
    chomp ($input);
    if ($input eq '') {
      if (defined $candidate && $candidate ne '') {
	$feature_id = $candidate;
	$valid = 1;
      }
    } else {
      $feature_id = $input;
      $valid = 1;
    }
    if (defined $feature_id) {$feature_id =~ s/\s//g};
  } until ($valid && $feature_id =~ /^WBsf\d+$/);


  return $feature_id;
}
##########################################################################
# get the date-tagged remark
sub get_remark {
  my ($method, $TF_name, $public_name) = @_;
  my $remark;

  my ($day, $mon, $yr)  = (localtime)[3,4,5];
  my $date = sprintf("%04d-%02d-%02d",$yr+1900, $mon+1, $day);


  print "Remark (anything else?) > ";
  my $input =  <STDIN>;
  chomp ($input);
  if ($input ne '') {
    $remark = $input . " [$date $ENV{USER}]";
  }

  return $remark;
}
##########################################################################
# details of the transcription factor objects
sub get_TFs {

  my %TFs;
  my %TF_genes;


  #load TF details from table maker query
  my $table = $wormbase->table_maker_query($currentdb, &write_def_file);
  while(<$table>) {
    s/\"//g; #"
    next if (/acedb/ or /\/\//);
    chomp;
    my ($TFID, $TF_name, $TF_gene) = split(/\t/,$_);
    if (defined $TF_gene) {
      $TF_name =~ s/\\//g;
      $TFs{$TF_name} = $TFID;
      push @{$TF_genes{$TF_name}}, $TF_gene;
    }
  }

  return (\%TFs, \%TF_genes);
}

##########################################################################


sub write_def_file {
        my $def = '/tmp/make_enhancer_features.def';
        open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
        my $txt = <<END;
// Spread sheet definition for the ACeDB software 
// User: gw3
// Date: 2014-06-30_16:55:23

// %n (%%n in the graphic) are parameter to be given on the command line in tace
// or by default by the Parameters given in this file
// \%n (%n in the graphic) are substituted by the value of column n at run time
// Line starting with // are ignored, starting with # are comments

Sortcolumn 1

Colonne 1 
Subtitle TFID 
Width 35 
Optional 
Visible 
Class 
Class Transcription_factor 
From 1 
 
Colonne 2 
Subtitle Name 
Width 60 
Optional 
Visible 
Text 
From 1 
Tag Name 
 
Colonne 3 
Subtitle Gene 
Width 20 
Optional 
Visible 
Class 
Class Gene 
From 1 
Tag Product_of 
 
 

// End of these definitions
END

        print TMP $txt;
        close TMP;
        return $def;
}

##########################################################################
# get the user to indicate the location
sub get_location {
  my ($gene_id, $chromosome, $gene_start, $gene_sense) = @_;
  my ($start, $end, $sense);
  my $valid = 0;
  my $quit = 0;

  #print "gene location: $gene_id, $chromosome, $gene_start, $gene_sense\n";

  do {
    print "\n";
    print "New Feature location\n";
    print "offsets: a pair of numbers, -ve for upsteam, +ve for downstream of the START\n";
    print "flanking sequences: acgtacgt agctagct\n";
    print "sequence of region: aaggcctt\n";
    print "quit: quit\n";
    print "\n";
    print "Location (offsets, sequence, flanking sequences, or quit) > ";
    my $input = <STDIN>;
    chomp ($input);
    if ($input eq 'q' || $input eq 'quit') {$quit = 1; last}
    if ($input eq '') {
      print "WARNING: no location given\n";
    } else {
      # a pair of numbers, -ve for upsteam, +ve for downstream of the START
      if ($input =~ /([+-]*\d+)[\s.,]+([+-]*\d+)/) {
	my $offset1 = $1;
	my $offset2 = $2;
	if ($gene_sense eq '+') {
	  $start = $gene_start + $offset1;
	  $end = $gene_start + $offset2;
	} else {
	  $start = $gene_start - $offset1;
	  $end = $gene_start - $offset2;
	}
	$sense = $gene_sense;
	$valid = 1;


      # a pair of flanking sequences
      } elsif ($input =~ /([ACGTacgt]+)\s+([ACGTacgt]+)/) {
	my $flank1 = lc $1;
	my $flank2 = lc $2;
      # search within 30Kb of the START of the gene in both senses
	my $area = 30000;
	my $seq = $seq_obj->Sub_sequence($chromosome, $gene_start-$area, $area * 2);
	if ($gene_sense eq '-') {
	  $seq = $seq_obj->DNA_revcomp($seq);
	}

	my $start_looking = 0;
	my @matches1;
	while ((my $loc = index($seq, $flank1, $start_looking)) >= 0) {
	  push @matches1, $loc + length($flank1);
	  $start_looking = $loc + 1;
	}
	if ($#matches1 > 0) {
	  print "WARNING: there is more than 1 match to the left flank, increase the length of this.\n";
	} elsif ($#matches1 == -1) {
	  print "WARNING: the left flank was not found.\n";
	} else {
	  $start_looking = 0;
	  my @matches2;
	  while ((my $loc = index($seq, $flank2, $start_looking)) >= 0) {
	    push @matches2, $loc;
	    $start_looking = $loc + 1;
	  }
	  if ($#matches2 > 0) {
	    print "WARNING: there is more than 1 match to the right flank, increase the length of this.\n";
	  } elsif ($#matches2 == -1) {
	    print "WARNING: the right flank was not found.\n";
	  } else {
	    if ($gene_sense eq '+') {
	      $start = $gene_start - $area + $matches1[0] + 1;
	      $end = $gene_start - $area + $matches2[0];
	    } else {
	      $start = $gene_start - $area - $matches2[0] + 1;
	      $end = $gene_start - $area + $matches1[0];
	    }
	    $sense = $gene_sense;
	    $valid = 1;
	  }
	}


      # single sequence
      } elsif ($input =~ /([ACGTacgt]+)/) {
	my $flank1 = lc $1;
      # search within 10Kb of the START of the gene in both senses
	my $area = 10000;
	my $seq = $seq_obj->Sub_sequence($chromosome, $gene_start-$area, $area * 2);
	if ($gene_sense eq '-') {
	  $seq = $seq_obj->DNA_revcomp($seq);
	}

	my $start_looking = 0;
	my @matches1;
	while ((my $loc = index($seq, $flank1, $start_looking)) >= 0) {
	  push @matches1, $loc;
	  $start_looking = $loc + 1;
	}
	if ($#matches1 > 0) {
	  print "WARNING: there is more than 1 match to this sequence, use flanking sequences instead.\n";
	} elsif ($#matches1 == -1) {
	  print "WARNING: the sequence was not found, use flanking sequences instead.\n";
	} else {
	  if ($gene_sense eq '+') {
	    $start = $gene_start - $area + $matches1[0] + 1;
	    $end = $start + length($flank1) - 1;
	  } else {
	    $end = $gene_start - $area + $matches1[0];
	    $start = $end - length($flank1) + 1;
	  }
	  $sense = $gene_sense;
	  $valid = 1;
	}
      } else {
	print "WARNING: invalid input\n";
	print "Enter one of:\n";
	print "\tstart and end offset from START codon first base (upstream is -ve)\n";
	print "\tflanking sequences\n";
	print "\tthe sequence\n";
      }
    }
  } until ($valid);

  #print "result: $chromosome, $start, $end, $sense\n";

  return($quit, $chromosome, $start, $end, $sense);
}
##########################################################################
sub get_gene_start {
  my ($gene_id, $gene_name) = @_;

  my ($chromosome, $gene_start, $gene_end, $gene_sense);
  my $cds;

  # get the CDS from the gene
  # if more than one, go for the 'a' isoform as this will have been in existence longer than the others
  my $gene_obj = $db->fetch(Gene => "$gene_id");
  if (defined $gene_obj) {
    my @cds = $gene_obj->Corresponding_CDS;
    if (grep /^$gene_name$/, @cds) {
      $cds = $gene_name; # the name on the command-line is an exact match to an isoform, so use this
    } else {
      foreach my $this_cds (@cds) { # see if we can find the 'a' isoform, else use whatever is there
	$cds = $this_cds;
	$cds =~ /(\w)$/;
	if ($1 eq 'a') {last} # get the 'a' isoform if it exists
      }
      if ($#cds > 0) {print "WARNING: there is more than one isoform in this gene.\nUsing isoform '$cds' for offsets.\n";}
    }
  } else {
    die "ERROR: $gene_id is not a coding gene";
  }

  my $cds_obj = $db->fetch(CDS => "$cds");
  if (defined $cds_obj) {
    my $clone_start;
    my $clone_end;
    my $clone = $cds_obj->Sequence;
    my @clone_CDSs = $clone->CDS_child;
    foreach my $CDS ( @clone_CDSs ) {
      next unless ($CDS->name eq $cds);
      $clone_start = $CDS->right->name;
      $clone_end = $CDS->right->right->name;
      last;
    }
    ($chromosome, $gene_start) = $coords->Coords_2chrom_coords($clone, $clone_start);
    ($chromosome, $gene_end) = $coords->Coords_2chrom_coords($clone, $clone_end);
    $gene_sense = '+';
    if ($gene_start > $gene_end) {
      $gene_sense = '-';
    }
  }
  return ($chromosome, $gene_start, $gene_sense)
}
##########################################################################
# get the gene id from whatever was put on the command-line
sub get_gene_id {

  my ($gene_name) = @_;

  my $gene_id;

# %seq2gene  # 'AC3.3' => 'WBGene00000024'
# %cgc2gene  # 'abu-1' => 'WBGene00000024',


  if ($gene_name =~ /^WBGene\d+$/) {
    $gene_id = $gene_name;
  } elsif (exists $seq2gene{$gene_name}) {
    $gene_id = $seq2gene{$gene_name};
  } elsif (exists $cgc2gene{$gene_name}) {
    $gene_id = $cgc2gene{$gene_name};
  }

  if (!defined $gene_id) {die "ERROR: gene name '$gene_name' is not recognised.\n"}
  return $gene_id;
}
##########################################################################
# get sequence if less than 100 bp long
sub get_sequence_text {
  my ($chromosome, $start, $end, $sense) = @_;

  my $seq = $seq_obj->Sub_sequence($chromosome, $start-1, $end-$start+1);
  if ($gene_sense eq '-') {
    $seq = $seq_obj->DNA_revcomp($seq);
  }
  if (length($seq) > 100) {$seq = ''}

  return $seq;
}
##########################################################################
# this looks in the current DB for overlapping Feature objects and
# reports them in case the user wishes to use them as the Feature ID

sub check_overlapping_features {

  my ($method, $sequence, $start, $end, $sense) = @_;

  my @exact_matches;
  my @overlapping;

  my $seq_obj = $db->fetch(Sequence => "$sequence");
  if (defined $seq_obj) {
    my @features = $seq_obj->Feature_object;
    foreach my $ft ( @features ) {
      my $ft_name = $ft->name;
      my $clone_start = $ft->right->name;
      my $clone_end = $ft->right->right->name;
      my $ft_desc = $ft->Description;
      if (!defined $ft_desc) {$ft_desc = ''}
      my $ft_method = $ft->Method;
      my $stars = '   ';
      if ($method =~ /bind/ && $ft_method =~ /bind/) {$stars = '** '} # two stars if they are both some sort of binding method
      if ($method eq $ft_method) {$stars = '***'} # three stars if the method is the same
      my $len = abs($end-$start);
      my $TF = $ft->Associated_with_transcription_factor;
      my $TransFac = '';
      if (defined $TF) {
	foreach my $TF_name (keys $TFs) {
	  my $TF_id = $TFs{$TF_name};
	  if ($TF_id eq $TF) {
	    $TransFac = "TF: $TF_name\t";
	  }
	}
      }
      my $line = "$stars $ft_name\tMethod: $ft_method\t${TransFac}Length: $len\tDescription: $ft_desc";
      if ($start == $clone_start && $end == $clone_end) {
	push @exact_matches, $line;
      } elsif ($sense eq '+' && ($start >= $clone_start && $start <= $clone_end || $end >= $clone_start && $end <= $clone_end) || 
	       $sense eq '-' && ($start >= $clone_end && $start <= $clone_start || $end >= $clone_end && $end <= $clone_start)) {
	push @overlapping, $line;	
      }
      

    }
  } else {
    print "WARNING: Can't find Sequence $sequence in CurrentDB\n";
  }
  
  if ($#exact_matches > -1) {
    print "\nEXACT Feature matches:\n";
    foreach my $line (@exact_matches) {
      print "$line\n";
    }
    print "\n";
  }
  
  if ($#overlapping > -1) {
    print "\nOverlapping Features:\n";
    foreach my $line (@overlapping) {
      print "$line\n";
    }
    print "\n";
  }


}
