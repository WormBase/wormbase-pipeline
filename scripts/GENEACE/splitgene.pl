#!/usr/local/bin/perl5.8.0 -w
#
# splitgene.pl
#
# by Keith Bradnam
#
# simple script for creating new (sequence based) Gene objects when splitting 
# existing gene 
#
# Last edited by: $Author: krb $
# Last edited on: $Date: 2004-09-20 13:23:20 $

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;

###################################################
# command line options                            # 
###################################################

my $old;         # sequence name for existing gene
my $new;         # sequence name for new gene
my $who;         # Person ID for new genes being created (defaults to krb = WBPerson1971)
my $id;          # force creation of gene using set ID
my $gene_id;     # stores highest gene ID
my $email;       # email new Gene IDs back to users to person who requested it
my $load;        # load results to geneace (default is to just write an ace file)
my $verbose;     # toggle extra (helpful?) output to screen

GetOptions ("old=s"     => \$old,
            "new=s"     => \$new,
	    "who=i"     => \$who,
	    "id=i"      => \$id,
	    "email"     => \$email,
	    "load"      => \$load,
	    "verbose"   => \$verbose);


#####################################################
# warn about incorrect usage of command line options
#####################################################

die "must specify -old and -new options\n"               if (!$old || !$new);
die "-who option must be an integer\n"                   if ($who && ($who !~ m/^\d+$/));
die "-old option is not a valid type of sequence name\n" if ($old !~ m/^\w+\.\d{1,2}$/);
die "-new option is not a valid type of sequence name\n" if ($new !~ m/^\w+\.\d{1,2}$/);


######################################
# set person ID for curator
######################################
my $person;

if($who){
  $person = "WBPerson$who";
}
else{
  # defaults to krb
  $person = "WBPerson1971";
}


############################################################
# set database path, open connection and open output file
############################################################

my $tace = &tace;
my $database = "/wormsrv1/geneace";

my $db = Ace->connect(-path  => $database,
		      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

open(OUT, ">/wormsrv1/geneace/fix.ace") || die "Can't write to output file\n";

# find out highest gene number in case new genes need to be created
my $gene_max = $db->fetch(-query=>"Find Gene");


# process gene information
&process_gene;


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

  # checks for old and new (old gene should exists, new gene shouldn't)
  my $old_exists = 1;
  my $new_exists = 0;
  my $split_gene;
  my $split_gene_name;

  # Look up existing gene based on sequence name

  my ($old_gene_name) = $db->fetch(-query=>"Find Gene_name $old");

  # proceed if this is valid
  if(defined($old_gene_name) && $old_gene_name->Sequence_name_for){
    my $old_gene = $old_gene_name->Sequence_name_for;
    my ($old_version) = $old_gene->Version;

    # now check that new split gene doesn't already exist
    ($split_gene_name) = $db->fetch(-query=>"Find Gene_name $new");

    if(defined($split_gene_name) && $split_gene_name->Sequence_name_for){
      $split_gene = $split_gene_name->Sequence_name_for;
      print "ERROR: $new already exists as $split_gene\n";
      $new_exists = 1;
    }
    else{

      # new version number
      my $new_version = $old_version+1;

      # get new gene ID, unless specified by -id
      if($id){
	$gene_id = "WBGene" . sprintf("%08d",$id);
      }
      else{
	$gene_max++;
	$gene_id = "WBGene" . sprintf("%08d",$gene_max);
      }

      print "Splitting $old into $new\n" if ($verbose);
      print "$old = $old_gene, version $old_version becomes version $new_version\n" if ($verbose);
      print "$new = $gene_id, version 1\n" if ($verbose);


      # update existing gene info
      print OUT "Gene : $old_gene\n";
      print OUT "Version $new_version\n";
      print OUT "History Version_change $new_version now $person Event Split_into $gene_id\n";
      print OUT "Split_into $gene_id\n\n";
     
      # write split gene info
      print OUT "Gene : $gene_id\n";
      print OUT "Live\n";
      print OUT "Version 1\n";
      print OUT "Sequence_name $new\n";
      print OUT "Public_name $new\n";
      print OUT "Species \"Caenorhabditis elegans\"\n";
      print OUT "History Version_change 1 now $person Event Split_from $old_gene\n";
      print OUT "Split_from $old_gene\n";
      print OUT "Method Gene\n\n";
    }

  }
  else{
    $old_exists = 0;
    print "ERROR: $old does not exist in the database\n";
  }
    
  
  ######################################
  # email user to notify of new gene ID
  ######################################
  
  if($email){
    # set default address to krb in case wrong user ID used
    my $address = "krb\@sanger.ac.uk";
    
    $address = "ar2\@sanger.ac.uk"          if ($person eq "WBPerson1847");
    $address = "dl1\@sanger.ac.uk"          if ($person eq "WBPerson1846");
    $address = "pad\@sanger.ac.uk"          if ($person eq "WBPerson1983");
    $address = "dblasiar\@watson.wustl.edu" if ($person eq "WBPerson1848");
    $address = "tbieri\@watson.wustl.edu"   if ($person eq "WBPerson1849");
    $address = "pozersky\@watson.wustl.edu" if ($person eq "WBPerson1867");
    
    # write email
    my $text;
    my $subject;

    if($new_exists==1){
      $text = "\n\nYou requested a new gene ID for $new (split from $old), but $new already exists as $split_gene\n\n";
      $subject = "WormBase Gene ID request for split gene $new:  FAILED";
    }
    elsif($old_exists==0){
      $text = "\n\nYou requested a new gene ID for $new (split from $old), but $old doesn't exist\n\n";
      $subject = "WormBase Gene ID request for split gene $new:  FAILED";
    } 
    else{
      $text = "\n\nYou requested a new gene ID for $new (split from $old), this Gene ID is $gene_id\n\n";
      $subject = "WormBase Gene ID request for split gene $new:  SUCCESSFUL";
    }
    $text .= "This email was generated automatically, please reply to krb\@sanger.ac.uk\n";
    $text .= "if there are any problems\n";
    
    open (MAIL,  "|/bin/mailx -r \"krb\@sanger.ac.uk\" -s \"$subject\" $address");
    print MAIL "$text";
    close (MAIL);
    
    print "$address was emailed regarding gene ID for $new\n";
  } 
} 



=pod
                                                                                           
=head2   NAME - splitgene.pl
                                                                                           
=head1 USAGE
                                                                                           
=over 4
                                                                                           
=item splitgene.pl -[options]
  
=back
  
=head1 DESCRIPTION
  
A script designed to create new gene objects to load into geneace by splitting an
existing gene.  Just supply an old (existing) sequence name, a new sequence name 
for the result of the split, a person ID of curator providing the information and 
optionally a new Gene object ID.  Resulting acefile will be made in /wormsrv1/geneace/fix.ace

Example 1 
splitgene.pl -old AH6.3 -new AH6.11 -who 1971 -id 2342 -load
 
 
=head2 MANDATORY arguments:

=over 4

=item -old

must specify a valid CDS/Pseudogene/Transcript name.  Script will tell you if it corresponds
to an existing gene, else will warn

=back

=item -new

must specify a new sequence name for split gene.  Script will warn if this already exists

=back

=head2 OPTIONAL arguments:
                                                                                           
=over 4

=item -id <number>

Where the number is the new gene ID (ignore leading zeros).  If -id is not specified then the
script will look to see what the next available gene ID is

=item -who <number>

Where number should correspond to a person ID...if this number doesn't match anyone then 
the script will assume that it is krb
                                                                                           
=item -email

person corresponding to -who option will be emailed notification, email goes to
krb@sanger.ac.uk if -who option doesn't correspond to a curator

=item -verbose

writes extra output to screen
                                                                                           
=item -load

will attempt to load the acefile into geneace (need to have write access!)
                                                                                           
                                                                                           
=head1 AUTHOR Keith Bradnam (krb@sanger.ac.uk)
                                                                                           
=back
                                                                                           
=cut
