#status!/usr/local/bin/perl5.8.0 -w
# 
# cgc_strain_parser.pl
#
# by Keith Bradnam
#
# Script to convert cgc strain file into ace file for geneace
# Page download and update upload to geneace has been automated [ck1]

# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2015-05-05 11:38:21 $

use strict;
use lib $ENV{'CVS_DIR'};
use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
#use lib '/nfs/users/nfs_g/gw3/Nameserver-API';

use lib '/software/worm/lib/perl';

use Wormbase;
use GENEACE::Geneace;
use Log_files;
use NameDB_handler;
use Ace;

use Storable;
use Getopt::Long;

#######################
# check user is wormpub
#######################

my ($test, $help, $debug, $verbose, $store, $wormbase, $species, $logfile);
my ($verbose, $load,$ndbUser,$ndbPass, $path, $input_file,$pg);

GetOptions (
	    "test"              => \$test,
	    'help'              => \$help,
            'debug=s'           => \$debug,
            'verbose'           => \$verbose,
            'logfile:s'         => \$logfile,
            'load'              => \$load,
            'ndbuser=s'         => \$ndbUser,
            'ndbpass=s'         => \$ndbPass,
            'path=s'            => \$path,
            'cgcfile=s'         => \$input_file,
            'pgpassword=s'      => \$pg,
       )||die($@);


&usage if ($help);
my $user = `whoami`; chomp $user;
if ($user ne "wormpub"){
  die("You have to be wormpub to run this script!\n");
}

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -organism => $species,
			     -test => $test,
                             );
}


# establish log file.
my $log = $logfile ? Log_files->make_log($logfile, $debug) : Log_files->make_build_associated_log($wormbase);

$log->log_and_die("input file $input_file doesn't exist\n") unless -e $input_file;

#######################
# misc variables
#######################
my $geneace_dir = $wormbase->database('geneace');
my $tace        = $wormbase->tace;
my $rundate     = $wormbase->rundate;

$path = $geneace_dir.'/STRAIN_INFO' if not defined $path;
$input_file = "$path/cgc_strain_list_$rundate" if not defined $input_file;

############################################
# get hash to convert CGC name to Gene ID
############################################
my $ga = init Geneace($wormbase);
my %Transgene_ids = %{$ga->transgene_ids($pg)};
my %Gene_info = %{$ga -> gene_info()};
my $last_gene_id_number = $ga ->get_last_gene_id();

########################################################################
# Set up various output files
######################################################################## 
my $current_strain_ace      = "$path/cgc_strain_info_$rundate.ace";
my $delete_strain_ace       = "$path/cgc_strain_info_$rundate.delete.ace";
my $missingAuthors          = "$path/missingStrainAuthors_$rundate.txt";
my $gene_allele_connections = "$path/gene_allele_connections.$rundate.ace";
my $potential_new_genes     = "$path/potential_new_genes.$rundate.ace";
my $backup_file             = "$path/strain_class_backup.$rundate.ace";
my $allelefluff             = "$path/allele_public_name.$rundate.ace";
my $transgene_report        = "$path/transgene_report.$rundate.txt";

open(STRAIN,       ">$current_strain_ace") || die "cant create output file $current_strain_ace\n";
open(DELETE_STRAIN,">$delete_strain_ace") || die "can't create $delete_strain_ace\n";
open(GENE2ALLELE,  ">$gene_allele_connections") || die "\nCan't open $gene_allele_connections\n";
open(NEWGENES,     ">$potential_new_genes") || die "\nCan't open $potential_new_genes\n";
open(ALLELEFLUFF,  ">$allelefluff") || die "\nCan't open $allelefluff\n";
open(MISSINGPERSON,">$missingAuthors") || die "\nCan't open $missingAuthors\n";
open(TRANSGENEREPORT, ">$transgene_report") or die "\nCould not open $transgene_report for writing\n";

print NEWGENES "// This file should *ONLY* be loaded to geneace when it has been fully checked\n";
print NEWGENES "// by hand.  If these Gene objects are ok, then they will need Gene IDs added.\n";
print NEWGENES "// Some of these gene names might already exist if they are from non-elegans species.\n\n";

############################################
# loop through strain data making ace files
############################################

# use following record separator
$/ = "--------------------";

my $big_counter=0;

# Count how many strains to loop through
my $strain_count = `grep Strain: $input_file | wc -l`;

open(INPUT, $input_file) || die "Can't open inputfile!"; 

# setup the nameserver
my $db = NameDB_handler->new($wormbase, $test);
my $geneAceDB = Ace->connect(-path => $geneace_dir) or die Ace->error;

while(<INPUT>){
  # drop out of loop before you reach last line of file
  last if $big_counter == $strain_count;
  $big_counter++;

  #remove carriage returns and new lines
  s/\015//g;
  s/\012//g;

  m/\s+Strain: (.*?)Species:/;
  my $strain = $1;
  print $_ if ! defined($strain) && $verbose;
  print "$big_counter : $strain\n" if $verbose;
  # skip object if no strain name
  next unless ($strain =~ /\w/);

  $strain =~ s/\s+$//g;
  if ($strain){
    print STRAIN "\n\nStrain : \"$strain\"\n";
    print DELETE_STRAIN "\n\nStrain : \"$strain\"\n";
  }
  m/\s+Species: (.+)Genotype:/;
  my $species = $1;
  $species =~ s/\s+$//g;
  if ($strain){
    print STRAIN "Species \"$species\"\n";
    print DELETE_STRAIN  "Species \"$species\"\n";
  }

  m/Genotype: (.*?)Description:/;
  my $genotype = $1;
  my $original_genotype = $genotype; # will be used later
  $genotype =~ s/\s{2,}/ /g; # condense whitespace to single gap between words
  $genotype =~ s/\s+$//g; # remove trailing whitespace
  print STRAIN "Genotype \"$genotype\"\n" unless ($genotype eq "");
  print DELETE_STRAIN  "-D Genotype \n" unless ($genotype eq "");

  m/Description: (.*?)Mutagen:/;
  my $description = $1;
  $description =~ s/\s{2,}/ /g;
  $description =~ s/\s+$//g;
  $description =~ s/\"//g; # get rid of any quotation marks
  # change any URLs present else the double back slash will be treated as a comment
  $description =~ s/http:\/\//URL: /g;

  if ($description =~ /Attribution: (.*) (Paper_evidence WBPaper\d+)/){
      my $attribution="$1";
      my $paper ="$2";
      $attribution =~ s/\s+/ /g;
      chomp $attribution;
      print STRAIN "Remark \"$attribution\" $paper\n";
      $description =~ s/\s*Attribution: .* Paper_evidence WBPaper\d+//;
  }
  print STRAIN "Remark \"$description\" Inferred_automatically \"From CGC strain data\"\n" unless ($description eq '');
  print DELETE_STRAIN  "-D Remark \n" unless ($description eq ""); 

  # find simple locus allele combinations e.g. spt-3(hc184)
  my $reg_exp=qr/^([Ca-z\-]{3,6}\-\d+\.{0,1}\d*)\(([a-z]{1,2}\d+)\)/;
  while($genotype =~ m/$reg_exp/){
    my $gene = $1;
    my $allele = $2;
    print STDERR "simple combination: $genotype\n" if $verbose;
    &check_details($gene,$allele,$strain,$species);
    $genotype =~ s/$reg_exp//;
  }

  # find transposon insertions
  $reg_exp=qr/([a-z]+(Si|Ti)\d+)/;
  while($genotype =~ m/$reg_exp/) {
    my $allele = $1;
    &check_details(undef, $allele, $strain, $species);
    $genotype =~ s/$reg_exp//;
  }

  # find chromosomal aberrations e.g. szT1
  $reg_exp=qr/\s*([a-z]{1,3}(Dp|Df|In|T|C)\d+)/;
  while($genotype =~ m/$reg_exp/){
    my $rearrangement = $1;
    print STRAIN "Rearrangement \"$rearrangement\"\n";
    print DELETE_STRAIN  "-D Rearrangement \"$rearrangement\"\n";
    $genotype =~ s/$reg_exp//;
  }

  # find transgenes e.g. zhEx11
  $reg_exp=qr/\s*([a-z]{1,3}(Ex|Is)\d+)/;
  while($genotype =~ m/$reg_exp/){
    my $transgene = $1;

    if (exists $Transgene_ids{$transgene}) {
      print STRAIN "Transgene \"$Transgene_ids{$transgene}\"\n";
      # delete all transgene references next time round, for safety
      print DELETE_STRAIN  "-D Transgene\n";
    } else {
      print TRANSGENEREPORT "\nTransgene : $transgene\n";
      print TRANSGENEREPORT "Strain : $strain\n";
      print TRANSGENEREPORT "Description : $description\n";
      print TRANSGENEREPORT "WBID : NONE\n"; 
    }
    $genotype =~ s/$reg_exp//;
  }

  # find double barrelled alleles (revertants) e.g. daf-12(rh61rh412) 
  $reg_exp=qr/([Ca-z\-]{3,6}\-\d+\.{0,1}\d*)\(([a-z]{1,2}\d+)([a-z]{1,2}\d+)\)/;
  while($genotype =~ m/$reg_exp/){
    my $gene = $1;
    # need to split up allele name into two fields
    my $allele1 = $2;
    my $allele2 = $3;
    
    print STDERR "double barreled allele ($allele1 / $allele2): $genotype\n" if $verbose;

    &check_details($gene,$allele1,$strain,$species);
    &check_details($gene,$allele2,$strain,$species);
    $genotype =~ s/$reg_exp//;
  }

  # alleles affecting 2 genes like arf-1.1&F45E4.7(ok1840)
  $reg_exp=qr/([\w\.\-]+)&([\w\.\-]+)\(([a-z]{1,2}\d+)\)/;
  while($genotype =~ m/$reg_exp/){
    my $gene1 = $1;
    my $gene2 = $2;
    my $allele = $3;
    print STDERR "allele affecting 2 genes: $genotype\n" if $verbose;
    &check_details($gene1,$allele,$strain,$species);   
    &check_details($gene2,$allele,$strain,$species);   
    $genotype =~ s/$reg_exp//;
  }



  # find alleles attached to non-approved, or unusual gene names e.g. let-?(h661)
  while($genotype =~ m/([\w\.\-\?]+)\(([a-z]{1,2}\d+)\)/){
    my $gene = $1;
    my $allele = $2;
    print STDERR "allele attached to non-approved gene: $genotype\n" if $verbose;
    &check_details($gene,$allele,$strain,$species);
    $genotype =~ s/\([a-z]{1,2}\d+\)//;
  }

  # find any skulking gene names missed by steps above, these are often where there is no allele name
  # or the allele name is wild-type, e.g. unc-24(+)
  while($genotype =~ m/([a-z]{3,4}\-\d+)/){
    my $gene = $1;
    # can't use check_details subroutine as there is no allele name to pass so add Strain->Gene connections here
    if (exists $Gene_info{$gene}{'Gene'}){
      print STRAIN "Gene \"$Gene_info{$gene}{'Gene'}\"\n";
      print DELETE_STRAIN  "-D Gene \"$Gene_info{$gene}{'Gene'}\"\n";
    }
    $genotype =~ s/[Ca-z\-]{3,6}\-\d+//;
  }

  m/Mutagen: (.*?)Outcrossed:/;
  my $mutagen = $1;
  $mutagen =~ s/\s{2,}/ /g;
  $mutagen =~ s/\s+$//g;
  print STRAIN "Mutagen \"$mutagen\"\n" unless ($mutagen eq '');
  print DELETE_STRAIN  "-D Mutagen \"$mutagen\"\n" unless ($mutagen eq '');


  m/Outcrossed: (.*?)Reference:/;
  my $outcrossed = $1;
  $outcrossed =~ s/\s{2,}/ /g;
  $outcrossed =~ s/\s+$//g;
  print STRAIN "Outcrossed\t\"$outcrossed\"\n" unless ($outcrossed eq "");
  print DELETE_STRAIN  "-D Outcrossed\n" unless ($outcrossed eq "");

  my $made_by;
  m/Made by: (.*?)Received:/;
  $made_by = $1;
  $made_by =~ s/\s{2,}/ /g;
  $made_by =~ s/\s+$//g;
  
  my $wperson = &find_author($made_by);

  if ($wperson =~ (/WBPerson\d{1,5}/)) {    
    print STRAIN "Made_by $wperson\n";
    print DELETE_STRAIN  "-D Made_by $wperson\n";
  }
  else {
    print STRAIN "Remark \"Made_by: $made_by\" CGC_data_submission\n" if $made_by;
    print MISSINGPERSON "\"$made_by\" $strain\n";
    $log->write_to("$wperson is not a valid WBPerson\n\n");
  }

  if (m/Received: (\d+\/\d+\/\d+)/){
    my @dates = split(/\//,$1);
    if($dates[2] > 60){
      print STRAIN "CGC_received \"19$dates[2]-$dates[0]-$dates[1]\"\n";
      print DELETE_STRAIN  "-D CGC_received \"19$dates[2]-$dates[0]-$dates[1]\"\n";
    }
    else{
      print STRAIN "CGC_received \"20$dates[2]-$dates[0]-$dates[1]\"\n";
      print DELETE_STRAIN  "-D CGC_received \"20$dates[2]-$dates[0]-$dates[1]\"\n";
    }
  }

  # always add CGC lab details
  if ($strain){
    print STRAIN "Location \"CGC\"\n";
    print DELETE_STRAIN  "-D Location \"CGC\"\n";
  }

}

$geneAceDB->close();

close(INPUT);
close(STRAIN);
close(DELETE_STRAIN);
close(GENE2ALLELE);
close(NEWGENES);
close(ALLELEFLUFF);
close(MISSINGPERSON);
close(TRANSGENEREPORT);

##################################################
# 1. backup strain class with timestamp
# 2. grep last delete.ace file and load to Geneace
# 3. load updated CGC strain genotype to Geneace
##################################################

my @deleteACE;

opendir(DIR, "$path/") || die "Can't read directory $path\n";
my @dir=readdir DIR;
closedir (DIR);
foreach (@dir){
  if ($_ =~ /^cgc_strain_info_\d+.delete.ace/){
    push(@deleteACE, $_);
  }
}
@deleteACE=sort @deleteACE;
my $last_delete_ace=$deleteACE[-2];
print "\n\nDelete_ace file from last update: $last_delete_ace\n\n";

my $command=<<END;
find strain "*"
show -a -T -f $backup_file
pparse $last_delete_ace
pparse $current_strain_ace
pparse $gene_allele_connections
pparse $allelefluff
save
quit
END

# only load this data if -load specified     
if($load){
  open (FH,"| $tace -tsuser \"CGC_strain_update\" $geneace_dir") || die "Failed to upload to test_Geneace";
  print FH $command;
  close FH;
}

print "\nThe script has run to completeion and is now going to end.  Goodnight\n\n";

$log->mail('maryann.tuli@wormbase.org');
$log->mail;
exit(0);

################
# subroutines
################


###########################
# function to get a variation from GeneAce by Other_name
#   - as that one is not stored in the NameDB
#   - will only return *one* variationid at max
#
sub _getVariation_byOtherName {
   my ($othername) = @_;
   my ($variationid) = $geneAceDB->fetch(-query => "Variation; Other_name = $othername");
   return "$variationid";
}

###########################
#
sub _get_variationId {
    my ($id)=@_;

    my $var;
    my $var_ref = $db->find_variations($id);
    foreach my $hash (@{$var}) {
      if ($hash->{'variation/name'} eq '$id') {$var = $hash->{'variation/id'}}
    }

    print STDERR "found: $id -> $var\n" if ($var && $verbose);

    return $var if $var;

    my ($new_ids, $batch_id) = $db->new_variations([$id]);
    my ($new_id) = keys %{$new_ids};
    print STDERR "found: $id -> $new_id\n" if $verbose;

    return $new_id;
}

######################################################################
# function to find a WBPerson by searching through a collection of tags

sub find_author {
    my ($searchterm)=@_;
    my ($wbperson) = $geneAceDB->aql("select all class Person where ->Standard_name like \"$searchterm\"".
                    " or ->Full_name like \"$searchterm\" or ->Also_known_as like \"$searchterm\"");
    $$wbperson[0]||=$searchterm;
    return "$$wbperson[0]";
}


###########################################################################################
# subroutine to do some basic checking of allele and gene details from strain genotype
# field to determine whether a new Gene object needs to be made or whether a Gene->Allele
# connection can be made
###########################################################################################

sub check_details {
  my ($gene,$_allele,$strain,$species)= @_;

  my $variationId=_get_variationId($_allele);

  # First thing is to make Strain->Allele connection (this assumes that the allele name
  # will link to a valid ?Allele object)
  print STRAIN "Variation \"$variationId\"\n";
  print DELETE_STRAIN  "-D Variation \"$variationId\"\n";  

  # if the gene name corresponds to a valid Gene object, add a Gene->Allele and Strain->Gene connections
  if(defined $gene and defined($Gene_info{$gene}{'Gene'})){
    print GENE2ALLELE "Gene : $Gene_info{$gene}{'Gene'}\n";
    print GENE2ALLELE "Allele $variationId Inferred_automatically \"From strain object: $strain\"\n\n";

    print STRAIN "Gene \"$Gene_info{$gene}{'Gene'}\"\n";
    print DELETE_STRAIN "-D Gene \"$Gene_info{$gene}{'Gene'}\"\n";  
  }

  # otherwise we might need to make a new Gene object (to be checked by curator)
  # Only want to consider gene names which look sensible (i.e. can't do much with 'let-?')
  elsif(defined $gene and $gene =~ m/^(\w|\-)+\d+$/){

    print NEWGENES "Gene : WBGene\n";
    print NEWGENES "Evidence Inferred_automatically \"Gene name parsed from strain object: $strain\"\n";
    print NEWGENES "Version  1\n";
    print NEWGENES "Version_change 1 now \"WBPerson2970\" Created\n";
    print NEWGENES "Live\n";
    print NEWGENES "Species  \"$species\"\n";
    print NEWGENES "Strain \"$strain\"\n";

    # many of these will probably be tag genes which we know will generally be ok to CGC approve
    if($gene =~ m/tag\-\d+/){
      print NEWGENES "CGC_name \"$gene\"\n";
      print NEWGENES "Public_name \"$gene\"\n";
      print NEWGENES "Gene_class \"tag\"\n";
    }
    # otherwise just set Other_name field for now
    else{
      print NEWGENES "Other_name \"$gene\"\n";
      print NEWGENES "Public_name \"$gene\"\n";
    }
    print NEWGENES "Allele \"$variationId\" Inferred_automatically \"From strain object: $strain\"\n\n";

  }
  print ALLELEFLUFF "Variation : $variationId\nPublic_name \"$_allele\"\nStatus Live\n\n" unless $geneAceDB->fetch(Variation => $variationId); 
}

##################################################################################################

sub usage {
  system ('perldoc',$0);
  exit;
}


__END__

=pod

=head1 NAME - cgc_strain_parser.pl

=back


=head1 USAGE

=over 4

=item cgc_strain_parser.pl

=back

This script will convert the CGC file of strain information into ace format.
It should be run against the file available at

http://www.cbs.umn.edu/CGC/strains/gophstrnt.txt

The script will write two ace files to your current directory, one to be loaded 
into geneace, and a second to be archived in /nfs/disk100/wormpub/DATABASES/geneace which will have 
delete instructions for removing all the data you have just added.  Some files 
will be loaded automatically (strain objects in one, and Gene->Allele connections 
in another).  

A third ace file (potential_new_genes.ace) needs to be checked by hand and modified 
to add Gene IDs (if data is ok).  This file may contain allele names (that appear in 
the strain data file) that we do not yet have in geneace.  These are probably KO 
consortium alleles that will be in our next allele update, but which have strain 
data already submitted to the CGC.

For these, you will probably find the name of the CDS that the allele deletes in the 
strain file.  This will then allow you to determine the relevant gene ID and then 
connect the allele name and strain name to that gene ID.  If you believe that we 
should have received details of these alleles already, chase this up with Mark Edgley.

Other items in this file to check by hand will concern alleles of C. briggsae genes
which will already have a gene ID.  

=over 4

=back

=over 4

=item OPTIONAL arguments:

-help (this help page)

-load (loads data to geneace, else just writes files)

=head1 AUTHOR - Keith Bradnam

Email krb@sanger.ac.uk

=cut
