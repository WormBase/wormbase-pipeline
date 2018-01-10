#!/software/bin/perl -w
#
# changegene.pl
#
# by Keith Bradnam
#
# simple script for changing class of gene objects (e.g. CDS->Pseudogene)
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2013-09-27 10:51:30 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;


###################################################
# misc variables and command line options         # 
###################################################

my $input;                   # when loading from input file
my $seq;                     # sequence name of gene to change
my $id;                      # gene ID of gene to change, ignore trailing zeros
my $class;                   # two or more letters to indicated old and new class (e.g. CDS Transcript)
my $who;                     # Person ID for new genes being created (defaults to mt3 = WBPerson2970)
my $person = "WBPerson2970"; # default
my $load;                    # load results to geneace (default is to just write an ace file)
my $verbose;                 # toggle extra (helpful?) output to screen
my $store;
my $species='elegans';

# hash for class lookups
my %change = ("CP" => ["CDS", "Pseudogene"],
	      "CR" => ["CDS", "Transcript"],
	      "CT" => ["CDS", "Transposon"],
	      "PR" => ["Pseudogene", "Transcript"],
	      "PT" => ["Pseudogene", "Transposon"],
	      "PC" => ["Pseudogene", "CDS"],
	      "RC" => ["Transcript", "CDS"],
	      "RP" => ["Transcript", "Pseudogene"],
	      "RT" => ["Transcript", "Transposon"],
	      "NULL" => ["NULL"]);
    
GetOptions ("input=s"   => \$input,
            "seq=s"     => \$seq,
	    "id=i"      => \$id,
	    "class=s"   => \$class,
	    "who=i"     => \$who,
	    "load"      => \$load,
	    "verbose"   => \$verbose,
	    "store:s"   => \$store,
	    "species:s" => \$species
);

my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new('-organism' => $species);
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

&check_command_line_options;

############################################################
# set database path, open connection and open output file
############################################################

my $tace = $wormbase->tace;
my $database = $wormbase->database('geneace');

my $db = Ace->connect(-path  => $database,
		      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

my $outdir = $database."/NAMEDB_Files/";
my $backupsdir = $outdir."BACKUPS/";
my $outname;
if (!$id) {$outname = "changegene_".&seq2gene($seq).".ace";}
elsif ($id) {$outname = "changegene_".$id.".ace";}
else {$outname = "changegene_".$seq.".ace";}

my $outfile = "$outdir"."$outname";
if (-e $outfile) {print "Warning this gene has probably already been processed.\n";}


open(OUT, ">$outfile") || die "Can't write to output file $outfile \n";

$log->write_to("changing gene ");
# get gene ID if -seq was specified
# else create a valid Gene object name based on numerical id from -id option
if ($seq){
  $log->write_to("$seq " );
  $id = &seq2gene($seq);
}
elsif($id){
  $id = "WBGene" . sprintf("%08d",$id);
}

$log->write_to("$id\n") if $id;



#######################################################################################
# Process list of genes if -input is specified, else just process command line options
#######################################################################################

if ($input){
  open(IN, "<$input") || die "Could not open $input\n";

  # process each gene in file, warning for errors
  while(<IN>){
    my($seq, $class) = split(/\s+/, $_);

    # set CGC to NULL if not specified
    $class = "NULL" if (!$class);

    # skip bad looking class fields
    if (!exists ($change{$class})){
      print "ERROR Not valid arguments for -class option\n";
      next;
    }

    # skip bad looking sequence names
    if ($seq !~ m/^\w+\.(\d{1,2}|\d{1,2}[a-z])$/){
      next;
    }

    # look up gene ID
    ($id = &seq2gene($seq)) if ($seq);
    next if (!defined($id));
    &process_gene($id,$class);
  }
  close(IN);
}
else{
  &process_gene($id,$class);
}


###################
# tidy up and exit
###################

$db->close;
close(OUT);

# load information to geneace if -load is specified
$wormbase->load_to_database($database, "$outfile", 'changegene', $log, undef, 1 ) if $load;
$wormbase->run_command("mv $outfile $backupsdir/$outname\n") if $load;
print "Output file has been cleaned away like a good little fellow\n" if $load;
$log->mail;

exit(0);

###############################################
#
# The main subroutine
#
###############################################


sub process_gene{
  my $id = shift;
  my $class = shift;

  # Look up gene based on ID
  my ($gene) = $db->fetch(-query=>"Find Gene $id");
  if(!$gene){
    print "ERROR: No gene corresponding to $id\n";
    last;
  }

  # get version of gene
  my ($version) = $gene->Version;
  print "ERROR: No gene version for $id\n" if (!$version);
  my $new_version = $version+1;


  # new version number
  print OUT "Gene $gene\n";
  print OUT "Version $new_version\n";


  # get class details from change hash
  my $old = $change{$class}[0];
  my $new = $change{$class}[1];

  # need to handle transposons differently
  if ($new eq "Transposon"){
    print OUT "History Version_change $new_version now $person Event Transposon_in_origin\nTransposon_in_origin\nSuppressed\nRemark \"This gene was determined to be of Transposon in origin so has been suppressed from the protein set. Detailed information about the origin of this gene can be found in the corresponding Transposon object associated with this gene.\"\n\n";
    
    print OUT "Gene $gene\n-D Gene_info\n\n";
    
  }
  else{
    print OUT "History Version_change $new_version now $person Event Changed_class $old $new\n\n";
  }

  print "Gene exists:  $gene (version $version), making version $new_version ($old->$new) \n" if ($verbose);

}



#####################################################
# warn about incorrect usage of command line options
#####################################################

sub check_command_line_options{
  die "-seq option not valid if -input is specified\n"     if ($input && $seq);
  die "Please specify either -input <file> or -seq <sequence> or -id <gene ID>\n" if (!$seq && !$id && !$input);
  die "-who option must be an integer\n"                   if ($who && ($who !~ m/^\d+$/));
  die "can't use -id option if processing input file\n"    if ($id && $input);
  die "Use -id or -seq option, not both\n"                 if ($id && $seq);
  die "can't use -class option if processing input file\n" if ($class && $input);
  my $cds_regex = $wormbase->cds_regex;
  die "-seq option is not a valid type of sequence name\n" if ($seq && ($seq !~ m/$cds_regex/));

  # look up -class option to see if it exists in hash
  if (!defined $class) {die "No arguments supplied for -class option....ending\n"}
  if ($class && (!exists ($change{$class}))){
    die "Not valid arguments for -class option\n";
  }

  # set person ID for curator if -who was specified
  $person = "WBPerson$who" if ($who);

}

############################################
#
# Look up Gene ID based on sequence name
#
############################################

sub seq2gene{
  my $seq = shift;

  my ($gene_name) = $db->fetch(-query=>"Find Gene_name $seq");
  if(defined($gene_name) && $gene_name->Sequence_name_for){
    $id = $gene_name->Sequence_name_for->name;
    print "\n$seq = $id\n" if ($verbose);
    return($id);
  }
  else{
    print "ERROR: no gene ID exists for $seq\n";
    return(0);
  }

}



=pod
                                                                                           
=head2   NAME - changegene.pl
                                                                                           
=head1 USAGE
                                                                                           
=over 4
                                                                                           
=item changegene.pl -[options]
  
=back
  
=head1 DESCRIPTION
  
A script designed to change the status of genes in geneace.  This is specifically to
cater for cases where a gene 'changes class', e.g. when a CDS becomes a Pseudogene or 
a CDS is made into a Transposon.  The result of this is to add one line of history
information to the gene object and to increment the version number.  For transposons,
additional changes are needed to 'kill' the gene object, so the following tags are removed:

Live
Map_info
Method
Sequence_name

It is possible that the Gene object has other information attached which may need to be
removed.

Just like the related script (newgene.pl), the script can process lists of genes if 
stored in an input file.  The -class option allows for two letters to denote the old and
new class:

C = CDS
P = Pseudogene
R = RNA gene (Transcript)
T = Transposon
 
Example
changegene.pl -seq AH6.4 -who 2970 -class CP -load
 
 
This would produce the following acefile
/nfs/disk100/wormpub/DATABASES/geneace/changegene_WBGene00005027.ace
and attempt to load it into geneace:
 
Gene WBGene00023428
Live
Version 2
History Version_change 2 now WBPerson2970 Event Changed_class CDS Pseudogene




=head2 MANDATORY arguments:

=over 4

=item -seq <name> or -id <number>

The -seq option must specify a valid CDS/Pseudogene/Transcript name.  This must correspond to an 
existing gene, else script will warn you.  Alternatively, specify the gene ID directly if you know 
it using the -id option.  In this case, just specify the significant digits and not any leading zeros,
e.g. for WBGene00001323 just specify -id 1323

=back

=head2 OPTIONAL arguments:
                                                                                           
=over 4


=item -who <number>

Where number should correspond to a person ID...if this number doesn't match anyone then 
the script will assume that it is mt3
                                                                                           

=item -verbose

writes extra output to screen
                                                                                           
=item -class

Must specify two letters which each correspond to one of four classes (C = CDS, R = RNA
gene (Transcript), P = Pseudogene, T = Transposon).  E.g. -class CP indicates a change
from a CDS to a Pseudogene

=item -input <file>

if input file has tab separated fields of sequence_name and class change (one pair 
per line) then script will process file in a batch style

=item -load

will attempt to load the acefile into geneace (need to have write access!)

=item -species <species>

species working with eg japonica
                                                                                           
                                                                                           
=head1 AUTHOR Keith Bradnam (mt3@sanger.ac.uk)
                                                                                           
=back
                                                                                           
=cut
