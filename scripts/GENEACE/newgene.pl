#!/software/bin/perl -w
#
# newgene.pl
#
# by Keith Bradnam
#
# simple script for creating new (sequence based) Gene objects 
#
# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2015-04-27 13:30:28 $

use strict;
use lib $ENV{'CVS_DIR'};
use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";

use Log_files;
use Wormbase;
use Getopt::Long;
use NameDB_handler;
use Storable;

###################################################
# command line options                            # 
###################################################

my $input;       # when loading from input file
my $seq;         # sequence name for new/existing gene
my $cgc;         # cgc name for new/existing gene
my $who;         # Person ID for new genes being created (should specify your ID eg. 1983 but will request one.)
my $p_clone;     # positive clone name for new/existing gene
my $id;          # force creation of gene using set ID
my $gene_id;     # stores highest gene ID
my $email;       # email new Gene IDs back to users to person who requested it
my $load;        # load results to geneace (default is to just write an ace file)
my $verbose;     # toggle extra (helpful?) output to screen
my $test;        # this switches between the live and test nameserver.
my $update_nameDB; # allows for update od various IDs as well as new ID requests from commandline.
my $species = 'elegans';  #default to elegans if not specified
my $debug;
my $store;
my ($USER,$PASS,); # provide your mysql username and password to request new WBGeneIDs
my $sneak;       # option to store the ID::Gene data in a pseudo .ace format.
my $outdir;      # specify your own output directory for the files to load.
my $bio;         # Mandatory biotype

GetOptions ("input=s"     => \$input,
            "seq=s"       => \$seq,
	    "cgc=s"       => \$cgc,
	    "who=i"       => \$who,
	    "id=s"        => \$id,
	    "email"       => \$email,
	    "load"        => \$load,
	    "verbose"     => \$verbose,
	    "test"        => \$test,
	    "namedb"      => \$update_nameDB,
	    "species=s"   => \$species,
	    "debug=s"     => \$debug,
	    "store:s"     => \$store,   #
	    "user:s"      => \$USER,
	    "password:s"  => \$PASS,
	    "pseudoace:s" => \$sneak,
	    "out:s"       => \$outdir,
	    "bio:s"         => \$bio,
	    );

my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
			     -organism => $species,
                           );
}

my $log = Log_files->make_build_log($wormbase);

#####################################################
# warn about incorrect usage of command line options
#####################################################

my %species_data;
$species_data{'elegans'}->{'regex'}  = '^\w+\.\d{1,2}$';
$species_data{'briggsae'}->{'regex'} = '^CBG\d{5}$';
$species_data{'remanei'} = '1';
$species_data{'brenneri'} = '1';
$species_data{'japonica'} = '1';
$species_data{'brugia'} = '1';
$species_data{'ovolvulus'} = '1';
$species_data{'sratti'} = '1';


unless( $species_data{"$species"} ) {
	my @list = keys %species_data;
	die "$species not valid\nTry one of these: @list\n";
}
die "-seq option not valid if -input is specified\n"     if ($input && $seq);
die "-cgc option not valid if -input is specified\n"     if ($input && $cgc);

die "-cgc option not valid if -seq is not specified\n"   if ($cgc && !$seq);
die "You must specify either -input <file> or -seq <sequence> -cgc <cgc name>\n" if (!$seq && !$input);

die "-cgc option is not a valid type of CGC name\n"      if ($cgc && ($cgc !~ m/[a-z]{3,4}\-\d{1,3}(\.\d+)?$/));

die "-who option must be an integer\n"                   if ($who && ($who !~ m/^\d+$/));
die "can't use -id option if processing input file\n"    if ($id && $input);

die "-seq option is not a valid type of sequence name\n" if ($seq && ($seq !~ $wormbase->cds_regex));
die "-bio is a mandatory option\n" unless ($bio);

# set CGC field to null string if not specified
$cgc = "NULL" if (!$cgc);


# SO:0001217 CDS
# SO:0001263 Transcript
# SO:0000111 Transposon
# SO:0000336 Pseudogene

unless (($bio eq "CDS") || ($bio eq "Transcript") || ($bio eq "Transposon") || ($bio eq "Pseudogene")) {
  die "-bio option $bio is not valid, please use CDS/Transcript/Transposon/Pseudogene\n";
}

my $SO;

if ($bio eq "CDS") {$SO = "0001217";}
if ($bio eq "Transcript") {$SO = "0001263";}
if ($bio eq "Transposon") {$SO = "0000111";}
if ($bio eq "Pseudogene") {$SO = "0000336";}

######################################
# set person ID for curator
######################################
my $person;

if($who){
  $person = "WBPerson$who";
}
else {
  print "\nERROR: please specify a user ID - e.g. 1983\n\n[INPUT]:";
  my $tmp = <STDIN>;
  chomp $tmp;
  if (($tmp ne '1983') && ($tmp ne'2970') && ($tmp ne '4025') && ($tmp ne '4055') && ($tmp ne '3111') && ($tmp ne '615')){
    $log->log_and_die("UNKNOWN USER.....TERMINATE\n\n");
    print "UNKNOWN USER.....TERMINATE\n\n";
  }
  $person = "WBPerson$tmp";
}

my $namedb;
if( $update_nameDB ) {
  ######################################
  # setup NameDB connection
  ######################################
  my $DB;
  # Live or Test
  if ($test) {
    $DB = 'test_wbgene_id;utlt-db:3307';
    print "Using the TEST server $DB\n";
  }
  else {
    $DB   = 'wbgene_id;shap;3303';
    print "Using the LIVE server $DB\n";
  }
  unless (defined$USER && defined$PASS) {
    die "You must specify both a username and password when using the nameDB option!\n\n";
  }
  #verify if valid name
  $namedb = NameDB_handler->new($DB,$USER,$PASS);
  $namedb->setDomain('Gene');
}

############################################################
# set database path, open connection and open output file
############################################################

my $tace = $wormbase->tace;
my $database = $wormbase->database('geneace');
$database = glob("~wormpub/DATABASES/TEST_DATABASES/geneace") if $test;

$log->log_and_die("BAILING: someone has a writelock on the database and you specified to load it in\n") if $load && -e "$database/database/lock.wrm";


my $db = Ace->connect(-path  => $database,
		      -program =>$tace) || do { $log->write_to("tace Connection failure\n");
						print "Connection failure: ",Ace->error; die();};
unless ($outdir) {
  $outdir = $database."/NAMEDB_Files/";
}
my $backupsdir = $outdir."BACKUPS/";
my $outname;
if (defined$id){
  $outname = "newgene_".$id.".ace";
}
else {
  $outname = "newgene_".$seq.".ace";
}
my $outfile = "$outdir"."$outname";

if (defined $sneak){
  system ("/bin/touch $sneak");
  print "Creating $sneak file\n";
}

print "Creating $outfile\n";
if (-e $outfile) {print "Warning this split has probably already been processed.\n";}

open(OUT, ">$outfile") || die "Can't write to output file\n";
if ($sneak){
  open(SN,  ">>$sneak") || die "Can't write to $sneak file\n";
}

# find out highest gene number in case new genes need to be created
my $gene_max = $db->fetch(-query=>"Find Gene");
my $override;



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
#    if ($seq !~ m/^\w+\.(\d{1,2}|\d{1,2}[a-z])$/){
#      print "ERROR: Bad sequence name, skipping\n";
#      $log->write_to("ERROR: Bad sequence name, skipping\n");
#      next;
#    }
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
  my $command = "pparse $outfile\nsave\nquit\n";
  open (GENEACE,"| $tace -tsuser newgene $database") || die "Failed to open pipe to $database\n";
  print GENEACE $command;
  close GENEACE;
  $wormbase->run_command("mv $outfile $backupsdir"."$outname\n");
  print "Output file has been cleaned away like a good little fellow\n";
  $log->write_to("Output file has been cleaned away like a good little fellow\n");
  $log->write_to("Finished!!!!\n");
}
$log->mail();
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
    $log->write_to("Gene exists:  $gene (version $version)\n");
    $exists = 1;

    # If gene exists but -cgc not specified, then we can't do anything else!
    if($cgc eq "NULL"){
      print "ERROR: $seq already exists as $gene\n";
      $log->log_and_die("ERROR: $seq already exists as $gene\n");
    }
    
    # check that CGC name doesn't already exist if -cgc has been specified
    elsif($cgc && $gene->CGC_name){
      my $cgc_name = $gene->CGC_name;
      print "ERROR: $seq($gene) already has a CGC name ($cgc_name)\n";
      $log->write_to("ERROR: $seq($gene) already has a CGC name ($cgc_name)\n");
    }

    # can now process CGC name
    else{

      # new version number
      print "Creating version $new_version: CGC_name = $cgc\n" if ($verbose);
      print OUT "\nGene $gene\n";
      print OUT "Version $new_version\n";
      print OUT "History Version_change $new_version now $person Name_change CGC_name \"$cgc\"\n";
      print OUT "CGC_name \"$cgc\"\n";

      # need to also Gene_class link unless it already exists
      my $gene_class;
      if($gene->Gene_class){
	$gene_class = $gene->Gene_class;
	print "WARNING: $gene already connects to $gene_class\n";
	$log->write_to("WARNING: $gene already connects to $gene_class\n");
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
      $gene_id = $id =~ /WBGene/ ? $id : "WBGene" . sprintf("%08d",$id);
    }
    else{
      $gene_id = $namedb->new_gene($seq, 'Sequence', $species);
      $override = "1";
      print SN "Object : $seq\nGene $gene_id\n\n" if ($sneak);
    }
    
    print "$seq does not exist, creating new Gene object $gene_id\n";
    $log->write_to("$seq does not exist, creating new Gene object $gene_id\n");
    
    print OUT "\nGene : $gene_id\n";
    print OUT "Live\n";
    print OUT "Version 1\n";
    print OUT "Biotype \"SO:$SO\"\n";
    print OUT "Sequence_name \"$seq\"\n";
    if ($species eq "elegans"){
      print OUT "Other_name \"CELE_${seq}\"\n";
    }
    print OUT "Species \"${\$wormbase->full_name}\"\n";
    print OUT "History Version_change 1 now $person Event Created\n";
    print OUT "Method Gene\n";
    print OUT "Positive_clone $p_clone Inferred_automatically \"From sequence, transcript, pseudogene data\"\n" if ($species eq 'elegans');
    # set CGC name if it exists and set public name based on CGC name or sequence name
    if($cgc && ($cgc ne "NULL")){
      print OUT "CGC_name \"$cgc\"\n";
      print OUT "Public_name \"$cgc\"\n\n";
    }
    else{
      print OUT "Public_name \"$seq\"\n\n";
    }
  }



    ######################################
    # email user to notify of new gene ID
    ######################################
  
  if($email){
    # set default address to wormpub in case wrong user ID used
    my $address = "wormpub\@sanger.ac.uk";
    
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
    $email .= "This email was generated automatically, please reply to wormpub\@sanger.ac.uk\n";
    $email .= "if there are any problems\n";

    my $subject;
    if($exists){
      $subject = "WormBase Gene ID request for $seq:  FAILED";
    }
    else{
      $subject = "WormBase Gene ID request for $seq:  SUCCESSFUL";
    }
    open (MAIL,  "|/bin/mailx -r \"wormpub\@sanger.ac.uk\" -s \"$subject\" $address");
    print MAIL "$email";
    close (MAIL);

    print "$address was emailed regarding gene ID for $seq\n";
    $log->write_to("$address was emailed regarding gene ID for $seq\n");
  }

  # NameDB query / update

  if ($update_nameDB && !defined$override)  {
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
ID.  Resulting acefile will be made in /nfs/wormpub/DATABASES/geneace/newgene_WBGeneid.ace

More powerfully the script can additionally assign CGC names to genes as it creates
them, or just assign CGC names to pre-existing genes.  Finally, the script can process
lists of genes if stored in an input file.
 
Example 1 
newgene.pl -seq AH6.4 -who 2970 -id 5027 -load
 
 
This would produce the following acefile at /nfs/wormpub/DATABASES/geneace/newgene_WBGene00005027.ace 
and attempt to load it into geneace:
 
Gene WBGene00005027
Live
Version 1
Sequence_name AH6.4
Public_name AH6.4
Species "Caenorhabditis elegans"
History Version_change 1 now WBPerson2970 Event Created
Method Gene


Example 2
newgene.pl -seq AH6.4 -load

This would achieve the same effect (assuming that 5027 in the previous example is the
next available gene ID).  Here the script automatically looks up the highest gene ID
and adds 1 to get the new gene ID and requests a valid who option

Example 3
newgene.pl -seq Y116F11B.27 -who 1983 -species elegans -namedb -test -user <who> -password <pass> -pseudoace /tmp/data.txt

This would request a new WBGeneID from the test nameserver for Y116F11B.27 and output the data as other example but would also create the additional file data.txt which would be pseudo.ace in the format

object : "Y116F11B.27"
Gene WBGene00000123


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
the script will ask for one to be input to STDIN
                                                                                           
=item -email <string>

person corresponding to -who option will be emailed notification, email goes to
wormpub@sanger.ac.uk if -who option doesn't correspond to a curator

=item -verbose

writes extra output to screen
                                                                                           
=item -cgc

will also add CGC name details to gene as it is being created, you can also use this
to add CGC names to existing genes (in which case it will increment the version number of the gene).

=item -pseudoace <file>

will create a pseudo .ace file that will contain the object::Gene connection, useful for easy manipulation
prior to loading the data into a sequence database.

=item -user <string>

this is your mysql username

=item -password <string>

this is your mysql password

=item -test

This switches between the live nameserver and the test nameserver by switching the server info.
                                                                                           
=item -input <file>

if input file has tab separated fields of sequence_name and cgc_name (one pair 
per line) then script will process file in a batch style

=item -load

will attempt to load the acefile into geneace (need to have write access!)
                                                                                           
                                                                                           
=head1 AUTHOR Keith Bradnam (krb@sanger.ac.uk)
                                                                                           
=back
                                                                                           
=cut
