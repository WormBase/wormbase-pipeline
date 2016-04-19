#!/software/bin/perl -w
#
# splitgene.pl
#
# by Keith Bradnam
#
# simple script for creating new (sequence based) Gene objects when splitting 
# existing gene 
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2014-03-13 13:53:05 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

###################################################
# command line options                            # 
###################################################

my $old;         # sequence name for existing gene
my $new;         # sequence name for new gene
my $who;         # Person ID for new genes being created
my $id;          # force creation of gene using set ID
my $gene_id;     # stores highest gene ID
my $email;       # email new Gene IDs back to users to person who requested it
my $load;        # load results to geneace (default is to just write an ace file)
my $verbose;     # toggle extra (helpful?) output to screen
my $p_clone;     # positive clone name for new gene
my $species;
my $debug;
my $test;

GetOptions ("old=s"     => \$old,
            "new=s"     => \$new,
	    "who=i"     => \$who,
	    "id=s"      => \$id,
	    "test"      => \$test,
	    "email"     => \$email,
	    "load"      => \$load,
	    "verbose"   => \$verbose,
            "species=s" => \$species,
	    "debug=s"   => \$debug,);


my $wormbase = Wormbase->new( -debug   => $debug,
			      -organism => $species,
			   );

my $log = Log_files->make_build_log($wormbase);
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
  if (($tmp ne '1983') && ($tmp ne '2970') && ($tmp ne '4025') && ($tmp ne '4055') && ($tmp ne '3111')){
    $log->log_and_die("UNKNOWN USER.....TERMINATE\n\n");
    print "UNKNOWN USER.....TERMINATE\n\n";
  }
  $person = "WBPerson$tmp";
}

############################################################
# set database path, open connection and open output file
############################################################
my $tace = $wormbase->tace;
my $database = $wormbase->database('geneace');
$database = glob("~wormpub/DATABASES/TEST_DATABASES/geneace") if $test;

#####################################################
# warn about incorrect usage of command line options
#####################################################

die "must specify -old and -new options\n"               if (!$old || !$new);
die "-who option must be an integer\n"                   if ($who && ($who !~ m/^\d+$/));
die "-old option is not a valid type of sequence name\n" unless( ($old =~ $wormbase->cds_regex) or ($old =~ /WBGene\d{8}/) );
die "-new option is not a valid type of sequence name\n" if ($new !~ $wormbase->cds_regex);

##############################
# warn/notify on use of -load.
##############################

if (!defined$load) {print "\n\nWarning: Update will not be loaded into $database as -load option has not been specified.\n";}
elsif (defined$load) { print "Update will be loaded into $database\n";}

######################################
# open connection and open output file
######################################
my $db = Ace->connect(-path  => $database,
		      -program =>$tace) || do { $log->write_to("ERROR:tace Connection failure\n");print "Connection failure: ",Ace->error; die();};

my $outdir = $database."/NAMEDB_Files/";
my $backupsdir = $outdir."BACKUPS/";
my $outname = "splitgene_".$id.".ace";
my $outfile = "$outdir"."$outname";


if (-e $outfile) {$log->write_to("Warning this split has probably already been processed.\n");
		  print "Warning this split has probably already been processed.\n";}

open(OUT, ">$outfile") || die "Can't write to output file\n";
$log->write_to("Output file: $outfile\n");


# find out highest gene number in case new genes need to be created
my $gene_max = $db->fetch(-query=>"Find Gene");


# process gene information
&process_gene;


###################
# tidy up and exit
###################

$db->close;
close(OUT);
$log->write_to("Output file: $outfile\n");
# load information to geneace if -load is specified
$wormbase->load_to_database($database, "$outfile", 'split_gene') if $load;
$wormbase->run_command("mv $outfile $backupsdir/$outname\n") if $load;
print "Output file has been cleaned away like a good little fellow\n" if $load;
$log->write_to("Output file has been cleaned away like a good little fellow\n") if $load;
print "Finished!!!!\n";
$log->write_to("Finished!!!!\n");
$log->mail;
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

    # Look up existing gene based on sequence name or gene id
    my $old_gene;

    #if gene id passed
    if( $old =~ /WBGene\d{8}/ ) {
	$old_gene = $db->fetch(Gene => $old );
	die "$old doesn't exist\n" unless $old;
    }
    # if seq name passed
    else {
	my ($old_gene_name) = $db->fetch(-query=>"Find Gene_name $old");

  	# proceed if this is valid
	if(defined($old_gene_name) && $old_gene_name->Sequence_name_for){
	    $old_gene = $old_gene_name->Sequence_name_for;
	} 
	else{
	    $old_exists = 0;
	    print "ERROR: $old does not exist in the database\n";
	}
    }
    #carry on with $old_gene regardless of how we got it
    my ($old_version) = $old_gene->Version;
    
    # now check that new split gene doesn't already exist
    ($split_gene_name) = $db->fetch(-query=>"Find Gene_name $new");

    if(defined($split_gene_name) && $split_gene_name->Sequence_name_for){
	$split_gene = $split_gene_name->Sequence_name_for;
	$log->log_and_die("ERROR: $new already exists as $split_gene\n");
	$new_exists = 1;
    }
    else{
	
	# create new positive clone name from sequence name
	my $p_clone = $new;
	$p_clone =~ s/\.\S+$//;

	# new version number
	my $new_version = $old_version+1;

	# get new gene ID, unless specified by -id
	if($id){
	    $gene_id = $id =~ /WBGene/ ? $id : "WBGene" . sprintf("%08d",$id); #if the passed id includes WBGene - dont add it !
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
	if ($species eq "elegans"){
	  print OUT "Other_name \"CELE_$new\"\n";
	}
	print OUT "Species \"".$wormbase->full_name."\"\n";
	  print OUT "Positive_clone $p_clone Inferred_automatically \"From sequence, transcript, pseudogene data\"\n" if ($species eq 'elegans');
	print OUT "History Version_change 1 now $person Event Split_from $old_gene\n";
	print OUT "Split_from $old_gene\n";
	print OUT "Method Gene\n\n";
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
	$text .= "This email was generated automatically, please reply to hinxton\@wormbase.org\n";
	$text .= "if there are any problems\n";
	
	open (MAIL,  "|/bin/mailx -r \"wormpub\@sanger.ac.uk\" -s \"$subject\" $address");
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
    optionally a new Gene object ID.  Resulting acefile will be made in /nfs/wormpub/DATABASES/geneace/
    the file will be called splitgene_<newgeneid>.ace

    Example 1 
    splitgene.pl -old AH6.3 -new AH6.11 -who 2970 -id 2342 -load
    
    
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
    the script will ask for one to be input to STDIN
    
    =item -email

    person corresponding to -who option will be emailed notification, email goes to
    wormpub@sanger.ac.uk if -who option doesn't correspond to a curator

    =item -verbose

    writes extra output to screen
    
    =item -load

    will attempt to load the acefile into geneace (need to have write access!)
    
    
    =head1 AUTHOR Keith Bradnam (krb@sanger.ac.uk)
    
    =back
    
    =cut
