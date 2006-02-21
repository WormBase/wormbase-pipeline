#!/usr/local/bin/perl5.8.0 -w
#
# newgene.pl
#
# by Keith Bradnam
#
# simple script for creating new (sequence based) Gene objects 
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2006-02-21 14:30:40 $

use strict;
use lib $ENV{'CVS_DIR'};
use lib "/nfs/WWWdev/SANGER_docs/cgi-bin/Projects/C_elegans/lib";

use Wormbase;
use Getopt::Long;
use NameDB_handler;

###################################################
# command line options                            # 
###################################################

my $input;       # when loading from input file
my $seq;         # sequence name for new/existing gene
my $cgc;         # cgc name for new/existing gene
my $who;         # Person ID for new genes being created (defaults to mt3 = WBPerson2970)
my $p_clone;     # positive clone name for new/existing gene
my $id;          # force creation of gene using set ID
my $gene_id;     # stores highest gene ID
my $email;       # email new Gene IDs back to users to person who requested it
my $load;        # load results to geneace (default is to just write an ace file)
my $verbose;     # toggle extra (helpful?) output to screen
my $test;
my $update_nameDB;

GetOptions ("input=s"   => \$input,
            "seq=s"     => \$seq,
	    "cgc=s"     => \$cgc,
	    "who=i"     => \$who,
	    "id=i"      => \$id,
	    "email"     => \$email,
	    "load"      => \$load,
	    "verbose"   => \$verbose,
	    "test"      => \$test,
	    "namedb"    => \$update_nameDB

	    );

my $wormbase = Wormbase->new();

#####################################################
# warn about incorrect usage of command line options
#####################################################

die "-seq option not valid if -input is specified\n"     if ($input && $seq);
die "-cgc option not valid if -input is specified\n"     if ($input && $cgc);
die "-cgc option not valid if -seq is not specified\n"   if ($cgc && !$seq);
die "You must specify either -input <file> or -seq <sequence> -cgc <cgc name>\n" if (!$seq && !$input);
die "-cgc option is not a valid type of CGC name\n"      if ($cgc && ($cgc !~ m/^[a-z]{3,4}\-\d{1,3}$/));
die "-who option must be an integer\n"                   if ($who && ($who !~ m/^\d+$/));
die "can't use -id option if processing input file\n"    if ($id && $input);
die "-seq option is not a valid type of sequence name\n" if ($seq && ($seq !~ m/^\w+\.\d{1,2}$/));

# set CGC field to null string if not specified
$cgc = "NULL" if (!$cgc);



######################################
# set person ID for curator
######################################
my $person;

if($who){
  $person = "WBPerson$who";
}
else{
  # defaults to mt3
  $person = "WBPerson2970";
}

my $namedb;
if( $update_nameDB ) {
######################################
# setup NameDB connection
######################################

  my $DB   = 'wbgene_id;mcs2a';
  my $USER = "$person";
  my $PASS = 'wormpub';

#verify if valid name

  $namedb = NameDB_handler->new($DB,$USER,$PASS);
  $namedb->setDomain('Gene');
}
############################################################
# set database path, open connection and open output file
############################################################

my $tace = $wormbase->tace;
my $database = $wormbase->database('geneace');
$database = glob("~wormpub/DATABASES/geneace") if $test;

my $db = Ace->connect(-path  => $database,
		      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

open(OUT, ">$database/newgene.ace") || die "Can't write to output file\n";

# find out highest gene number in case new genes need to be created
my $gene_max = $db->fetch(-query=>"Find Gene");




#######################################################################################
# Process list of genes if -input is specified, else just process command line options
#######################################################################################

if ($input){
  open(IN, "<$input") || die "Could not open $input\n";

  # process each gene in file, warning for errors
  while(<IN>){
    my($seq, $cgc) = split(/\s+/, $_);

    # set CGC to NULL if not specified
    $cgc = "NULL" if (!$cgc);

    print "\n\n$seq - $cgc\n" if ($verbose);

    # skip bad looking sequence names
    if ($seq !~ m/^\w+\.(\d{1,2}|\d{1,2}[a-z])$/){
      print "ERROR: Bad sequence name, skipping\n";
      next;
    }
     
    &process_gene($seq,$cgc);
  }
  close(IN);
}
else{
  &process_gene($seq,$cgc);
}


###################
# tidy up and exit
###################

$db->close;
close(OUT);

# load information to geneace if -load is specified
if ($load){
  my $command = "pparse $database/newgene.ace\nsave\nquit\n";
  open (GENEACE,"| $tace -tsuser \"mt3\" $database") || die "Failed to open pipe to $database\n";
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
  my $seq = shift;
  my $cgc = shift;

  # flag to check whether gene already exists
  my $exists = 0;

  # Look up gene based on sequence name
  my $gene;
  my ($gene_name) = $db->fetch(-query=>"Find Gene_name $seq");
  
  # create positive clone name from sequence name
  my $p_clone = $seq;
  $p_clone =~ s/\.\S+$//;   

  # get gene object if sequence name is valid, else need to make new gene
  if(defined($gene_name) && $gene_name->Sequence_name_for){
    $gene = $gene_name->Sequence_name_for;
    my ($version) = $gene->Version;
    my $new_version = $version+1;
    print "Gene exists:  $gene (version $version)\n" if ($verbose);
    $exists = 1;

    # If gene exists but -cgc not specified, then we can't do anything else!
    if($cgc eq "NULL"){
      print "ERROR: $seq already exists as $gene\n";
    }
    
    # check that CGC name doesn't already exist if -cgc has been specified
    elsif($cgc && $gene->CGC_name){
      my $cgc_name = $gene->CGC_name;
      print "ERROR: $seq($gene) already has a CGC name ($cgc_name)\n";
    }

    # can now process CGC name
    else{

      # new version number
      print "Creating version $new_version: CGC_name = $cgc\n" if ($verbose);
      print OUT "Gene $gene\n";
      print OUT "Version $new_version\n";
      print OUT "History Version_change $new_version now $person Name_change CGC_name $cgc\n";
      print OUT "CGC_name $cgc\n";

      # need to also Gene_class link unless it already exists
      my $gene_class;
      if($gene->Gene_class){
	$gene_class = $gene->Gene_class;
	print "WARNING: $gene already connects to $gene_class\n";	
      }
      else{
	($gene_class) = ($cgc =~ m/^(\w+)\-\.*/);
	print OUT "Gene_class $gene_class\n";
	print "Creating new link to Gene_class $gene_class\n" if ($verbose);
      }

      # finally print out new public name field
      print OUT "Public_name $cgc\n\n";

    }
  }

  # gene object doesn't exist need to make it!
  else{
   
    # get new gene ID, unless specified by -id
    if($id){
      $gene_id = "WBGene" . sprintf("%08d",$id);
    }
    else{
      $gene_max++;
      $gene_id = "WBGene" . sprintf("%08d",$gene_max);
    }

    print "$seq does not exist, creating new Gene object $gene_id\n";
    
    print OUT "Gene : $gene_id\n";
    print OUT "Live\n";
    print OUT "Version 1\n";
    print OUT "Sequence_name $seq\n";
    print OUT "Species \"Caenorhabditis elegans\"\n";
    print OUT "History Version_change 1 now $person Event Created\n";
    print OUT "Method Gene\n";
    print OUT "Positive_clone $p_clone Inferred_automatically \"From sequence, transcript, pseudogene data\"\n";

    # set CGC name if it exists and set public name based on CGC name or sequence name
    if($cgc && ($cgc ne "NULL")){
      print OUT "CGC_name $cgc\n";
      print OUT "Public_name $cgc\n\n";
    }
    else{
      print OUT "Public_name $seq\n\n";      
    }
  }



    ######################################
    # email user to notify of new gene ID
    ######################################
  
  if($email){
    # set default address to mt3 in case wrong user ID used
    my $address = "mt3\@sanger.ac.uk";
    
    $address = "ar2\@sanger.ac.uk"          if ($person eq "WBPerson1847");
    $address = "gw3\@sanger.ac.uk"          if ($person eq "WBPerson4025");
    $address = "pad\@sanger.ac.uk"          if ($person eq "WBPerson1983");
    $address = "dblasiar\@watson.wustl.edu" if ($person eq "WBPerson1848");
    $address = "tbieri\@watson.wustl.edu"   if ($person eq "WBPerson1849");
    $address = "pozersky\@watson.wustl.edu" if ($person eq "WBPerson1867");
    
    my $email;
    if($exists){
      $email = "\n\nYou requested a new gene ID for $seq, but this gene already exists as $gene\n\n";
    }
    else{
      $email = "\n\nYou requested a new gene ID for $seq, this Gene ID is $gene_id\n\n";
    }
    $email .= "This email was generated automatically, please reply to mt3\@sanger.ac.uk\n";
    $email .= "if there are any problems\n";

    my $subject;
    if($exists){
      $subject = "WormBase Gene ID request for $seq:  FAILED";
    }
    else{
      $subject = "WormBase Gene ID request for $seq:  SUCCESSFUL";
    }
    open (MAIL,  "|/bin/mailx -r \"mt3\@sanger.ac.uk\" -s \"$subject\" $address");
    print MAIL "$email";
    close (MAIL);

    print "$address was emailed regarding gene ID for $seq\n";
  }

  # NameDB query / update

  if ( $update_nameDB ) {
    undef $gene_id;
    if ( $seq and $cgc ) {
      #make CGC name based gene
      my $name = $cgc;
      my $type = 'CGC';
      $namedb->validate_name($name, $type);
      $namedb->check_pre_exists($name, $type);
      $namedb->make_new_obj($name, $type);
      $$gene_id = $namedb->idGetByTypedName('CGC',$cgc);

      #add the sequence name
      $name = $seq;
      $type = 'CDS';
      $namedb->validate_name($name, $type);
      $namedb->check_pre_exists($name, $type);
      $namedb->isoform_exists($name, $type);
      $namedb->addName($id,"CDS",$seq);
    } elsif ( $seq ) {
      my $name = $seq;
      my $type = 'CDS';
      $namedb->validate_name($name, $type);
      $namedb->check_pre_exists($name, $type);
      $namedb->isoform_exists($name, $type);
      $namedb->make_new_obj($name, $type);
    } elsif ( $cgc ) {
      #make CGC name based gene
      my $name = $cgc;
      my $type = 'CGC';
      $namedb->validate_name($name, $type);
      $namedb->check_pre_exists($name, $type);
      $namedb->make_new_obj($name, $type);
      $$gene_id = $namedb->idGetByTypedName('CGC',$cgc);
    }

    $gene_name = $cgc if $cgc;
    $gene_name .= " $seq" if $seq;
    open (MAIL,  "|/bin/mailx -r \"ar2\@sanger.ac.uk\" -s \"$gene_name $gene_id\" \"ar2\@sanger.ac.uk");
    print MAIL "$email";
    close (MAIL);
  }
}

=pod
                                                                                           
=head2   NAME - newgene.pl
                                                                                           
=head1 USAGE
                                                                                           
=over 4
                                                                                           
=item newgene.pl -[options]
  
=back
  
=head1 DESCRIPTION
  
A script designed to create new gene objects to load into geneace.  Mainly written to
save time from adding all the mandatory tags that each new object needs.  Just supply
a sequence name, person ID of curator providing the information and a new Gene object
ID.  Resulting acefile will be made in /nfs/disk100/wormpub/DATABASES/geneace/newgene.ace

More powerfully the script can additionally assign CGC names to genes as it creates
them, or just assign CGC names to pre-existing genes.  Finally, the script can process
lists of genes if stored in an input file.
 
Example 1 
newgene.pl -seq AH6.24 -who 2970 -id 23428 -load
 
 
This would produce the following acefile at /nfs/disk100/wormpub/DATABASES/geneace/newgene.ace and attempt to
load it into geneace:
 
Gene WBGene00023428
Live
Version 1
Sequence_name AH6.24
Public_name AH6.24
Species "Caenorhabditis elegans"
History Version_change 1 now WBPerson2970 Event Created
Method Gene


Example 2
newgene.pl -seq AH6.24 -load

This would achieve the same effect (assuming that 23428 in the previous example is the
next available gene ID).  Here the script automatically looks up the highest gene ID
and adds 1 to get the new gene ID and assumed mt3 to the be the default option for -who



=head2 MANDATORY arguments:

=over 4

=item -seq

must specify a valid CDS/Pseudogene/Transcript name.  Script will tell you if it corresponds
to an existing gene, else will assume it will be a new gene

=back

=head2 OPTIONAL arguments:
                                                                                           
=over 4

=item -id <number>

Where the number is the new gene ID (ignore leading zeros).  If -id is not specified then the
script will look to see what the next available gene ID is

=item -who <number>

Where number should correspond to a person ID...if this number doesn't match anyone then 
the script will assume that it is mt3
                                                                                           
=item -email

person corresponding to -who option will be emailed notification, email goes to
mt3@sanger.ac.uk if -who option doesn't correspond to a curator

=item -verbose

writes extra output to screen
                                                                                           
=item -cgc

will also add CGC name details to gene as it is being created, you can also use this
to add CGC names to existing genes (in which case it will increment the version number of the gene).
                                                                                           
=item -input <file>

if input file has tab separated fields of sequence_name and cgc_name (one pair 
per line) then script will process file in a batch style

=item -load

will attempt to load the acefile into geneace (need to have write access!)
                                                                                           
                                                                                           
=head1 AUTHOR Keith Bradnam (krb@sanger.ac.uk)
                                                                                           
=back
                                                                                           
=cut
