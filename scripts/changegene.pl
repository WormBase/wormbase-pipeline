#!/usr/local/bin/perl5.8.0 -w
#
# changegene.pl
#
# by Keith Bradnam
#
# simple script for changing class of gene objects (e.g. CDS->Pseudogene)
#
# Last edited by: $Author: krb $
# Last edited on: $Date: 2004-09-13 10:49:47 $

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;

###################################################
# misc variables and command line options         # 
###################################################

my $input;                   # when loading from input file
my $seq;                     # sequence name of gene to change
my $id;                      # gene ID of gene to change
my $class;                   # two or more letters to indicated old and new class (e.g. CDS Transcript)
my $who;                     # Person ID for new genes being created (defaults to krb = WBPerson1971)
my $person = "WBPerson1971"; # default
my $load;                    # load results to geneace (default is to just write an ace file)
my $verbose;                 # toggle extra (helpful?) output to screen

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
	    "verbose"   => \$verbose);

&check_command_line_options;


############################################################
# set database path, open connection and open output file
############################################################

my $tace = &tace;
my $database = "/wormsrv1/geneace";

my $db = Ace->connect(-path  => $database,
		      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

open(OUT, ">/wormsrv1/geneace/fix.ace") || die "Can't write to output file\n";


# get gene ID if -seq was specified
# else add preceding part of name if specifying an ID
if ($seq){
  $id = &seq2gene($seq);
}
elsif($id){
  $id = "WBGene$id";
}




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

    print "\n\n$seq - $class\n" if ($verbose);

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
if ($load){
  my $command = "pparse /wormsrv1/geneace/fix.ace\nsave\nquit\n";
  open (GENEACE,"| $tace -tsuser \"krb\" /wormsrv1/geneace") || die "Failed to open pipe to /wormsrv1/geneace\n";
  print GENEACE $command;
  close GENEACE;
}



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
    print OUT "History Version_change $new_version now $person Event Made_into_transposon\n";
    print OUT "-D Live\n";
    print OUT "WARNING: Extra gene information needs to be removed from $gene\n";
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
  die "You must specify either -input <file> or -seq <sequence> or -id <gene ID>\n" if (!$seq && !$input);
  die "-who option must be an integer\n"                   if ($who && ($who !~ m/^\d+$/));
  die "can't use -id option if processing input file\n"    if ($id && $input);
  die "Use -id or -seq option, not both\n"                 if ($id && $seq);
  die "can't use -class option if processing input file\n" if ($class && $input);
  die "-seq option is not a valid type of sequence name\n" if ($seq && ($seq !~ m/^\w+\.\d{1,2}$/));

  # look up -class option to see if it exists in hash
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
    print "$seq = $id\n" if ($verbose);
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
information to the gene object and to increment the versio number.  For transposon,
the Live tag is removed (other data should also be removed from genes which are actually
transposon encoded CDSs but the script doesn't do this yet).

Just like the related script (newgene.pl), the script can process lists of genes if 
stored in an input file.  The -class option allows for two letters to denote the old and
new class:

C = CDS
P = Pseudogene
R = RNA gene (Transcript)
T = Transposon
 
Example
changegene.pl -seq AH6.24 -who 1971 -class CP -load
 
 
This would produce the following acefile at /wormsrv1/geneace/fix.ace and attempt to
load it into geneace:
 
Gene WBGene00023428
Live
Version 2
History Version_change 2 now WBPerson1971 Event Changed_class CDS Pseudogene




=head2 MANDATORY arguments:

=over 4

=item -seq <name> or -id <number>

The -seq option must specify a valid CDS/Pseudogene/Transcript name.  This must correspond to an 
existing gene, else script will warn you.  Alternatively, specify the gene ID directly if you know 
it using the -id option.

=back

=head2 OPTIONAL arguments:
                                                                                           
=over 4


=item -who <number>

Where number should correspond to a person ID...if this number doesn't match anyone then 
the script will assume that it is krb
                                                                                           

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
                                                                                           
                                                                                           
=head1 AUTHOR Keith Bradnam (krb@sanger.ac.uk)
                                                                                           
=back
                                                                                           
=cut
