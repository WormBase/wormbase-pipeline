#!/software/bin/perl -w
#
# Script to write out an Analysis XML file for submitting a BAM file to the ENA.
# 
# by Gary Williams
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2013-09-27 10:08:55 $      

# Submission guides:
# http://www.ncbi.nlm.nih.gov/books/NBK49167/
# - this describes BAM file submission to the NCBI at the bottom. The requirements seem to be quite exact and restrictive, especially 
#   Unmapped reads included 
#   Duplicate reads marked
#   Vendor quality filtered reads identified
#   Nonredundancy 
# 3.    One record per sequencing production unit starting with @RG, to be configured as follows:
# ID: an arbitrary ID used to link reads back to the read group header
# PL: the sequencing platform that generated the reads.
# PU: the "platform unit" - a unique identifier which tells you what run/experiment created the data.  For Illumina, please follow this convention: Illumina flowcell barcode suffixed with a period and the lane number (and further suffixed with period followed by sample member name for pooled runs). If referencing an existing already archived run, then please use the run alias in the SRA.
# LB: the unique identifier of the sequencing library that was sequenced. This should correspond to the SRA library name for already-archived runs.
# DT: the run start date of the instrument run. Please use ISO-8601 format.
# SM: the sample identifier. This should be the sample alias loaded in the SRA or in the metadata being submitted to the SRA.
# CN: the sequencing center that produced the data (This should be the INSDC short name for the Center.)



use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
use Time::Piece;
use Net::FTP;
use LWP::UserAgent;



######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $species, $srx, $validate, $modify, $justtesting, $justsubmit);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "srx:s"      => \$srx,
	    "validate"   => \$validate, # validate the submission, but don't submit it for real
	    "modify"     => \$modify, # To update objects, you must use the MODIFY action instead of ADD action. The MODIFY action should point to XML documents containing objects to be updated. When MODIFY action is given it is verified that the objects being updated must already have been accessioned. This is done by comparing the object alias or the previously assigned accession number against objects already submitted.
	    "justtesting" => \$justtesting, # use the test server (it is a copy of the production archive but it gets re-written everyday so it does not matter if you make a mistake)
	    "justsubmit" => \$justsubmit, # Only do a re-submission with the existing files in the tophat directory. This helps when testing.
	    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species,
			     );
}

if (! defined $species) {$species = $wormbase->species}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

if (!defined $srx) {
  $log->log_and_die("-srx not specified\n");
}


##########################
# MAIN BODY OF SCRIPT
##########################

# the analysis.xml validator wants a '.version' after the genome
# assembly accession but they don't appear to have version numbers, so
# I have just appended '.1' to all of them

my %assembly_accessions = (
			 'elegans'  => 'GCA_000002985.1',
                         'briggsae' => 'GCA_000004555.1',
			 'brenneri' => 'GCA_000143925.1',
			 'japonica' => 'GCA_000147155.1',
			 'remanei'  => 'GCA_000149515.1',
			 'pristionchus' => 'GCA_000180635.1',
			 'brugia'   => 'GCA_000002995.1',
);

my ($Study, $Sample, $Experiment, $Run, $Instrument, $Library_Name);
my $RG_bamfile;

my $db = Ace->connect(-path => $wormbase->autoace);

my $dir = "/nfs/nobackup/ensemblgenomes/wormbase/BUILD/RNASeq/$species/SRA/$srx/tophat_out";
if (!-e $dir) {$log->log_and_die("The tophat_out directory for $srx was not found\n")}
chdir $dir;


my $bamfile = "accepted_hits.bam";
my $analysis_xmlfile = "$srx.analysis.xml";
my $submission_xmlfile = "$srx.submission.xml";
my $accession_file = "Analysis_accession.dat";

# do we want to just do a submission with the existing files for testing purposes?
if (!$justsubmit) {

  # check to see if there is an analysis_accession file already in this directory
  # if there is then there will have been a previous submission of a BAM file for this experiment and we can MODIFY instead of SUBMIT
  if (-e $accession_file && ! $validate) {
    $log->log_and_die("File '$accession_file' exists, so it looks like you should be using '-modify'\n");
  }

  $log->write_to("Get SRA $srx details.\n");
  ($Study, $Sample, $Experiment, $Run, $Instrument, $Library_Name) = get_SRA_details($srx);

  $log->write_to("Add the \@RG tag to the BAM file header.\n");
  $RG_bamfile = add_RG_header_to_BAM_file($srx, $bamfile);

  $log->write_to("Make the Analysis XML file.\n");
  make_analysis_xml($analysis_xmlfile, $Study, $Sample, $Experiment, $Run, $RG_bamfile);

  $log->write_to("Make the Submission XML file.\n");
  make_submission_xml($validate, $modify, $submission_xmlfile, $analysis_xmlfile);

  $log->write_to("Upload the data.\n");
  data_upload($RG_bamfile);

}

$log->write_to("Submit the data.\n");
my ($result, $text, $analysis_accession) = data_submission($submission_xmlfile, $analysis_xmlfile);
if ($result) {
  $log->write_to("Submission worked OK\n");
  $log->write_to("\n$text\n");
  $log->write_to("Save the analysis accession to file.\n");
  save_analysis_accession($analysis_accession);

} else {
  $log->write_to("Problem with the $srx submission:\n$text\n");
}

$db->close;

system("rm -f $RG_bamfile"); # tidy up

$log->write_to("Finished.\n");

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);






##############################################################
#
# Subroutines
#
##############################################################




##########################################
# get some details from the SRA
sub get_SRA_details {

  my ($srx) = @_;

# get lines like:

#Study	Sample	Experiment	Run	Analysis	Organism	Instrument Platform	Instrument Model	Library Name	Library Layout	Library Strategy	Library Source	Library Selection	Run Read Count	Run Base Count	File Name	File Size	md5	Ftp
#SRP001005	SRS004375	SRX006984	SRR019721	null	Caenorhabditis elegans	ILLUMINA	Illumina Genome Analyzer	N2_emb_rep 1	SINGLE	OTHER	GENOMIC	cDNA	13202626	448889284	SRR019721.fastq.gz	558Mb	fc4efcf7b219dcc946e2c1e05abd0d70	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR019/SRR019721/SRR019721.fastq.gz

  my @cols;
  my @Run;

  open(ENTRY, "/sw/arch/bin/wget -q -O - 'http://www.ebi.ac.uk/ena/data/view/reports/sra/fastq_files/$srx' |") || $log->log_and_die("Can't get information on SRA entry $srx\n");
  while (my $line = <ENTRY>) {
    chomp $line;
    if ($line =~ /^Study/) {next}
    @cols = split /\t+/, $line;   
    push @Run, $cols[3];
  }
  close(ENTRY);
  my ($Study, $Sample, $Experiment, $Run, $Instrument, $Library_Name) = ($cols[0], $cols[1], $cols[2], $cols[3], $cols[6], $cols[8]);

# sometimes we get the details returned like:
#Study	Sample	Experiment	Run	Analysis	Organism	Instrument Platform	Instrument Model	Library Name	Library Layout	Library Strategy	Library Source	Library Selection	Run Read Count	Run Base Count	File Name	File Size	md5	Ftp
#PRJNA33023	SRS001788	SRX001872	SRR006511	null	Caenorhabditis elegans	ILLUMINA	Illumina Genome Analyzer	L2_ce0109_rw001	SINGLE	EST	TRANSCRIPTOMIC	cDNA	6246633	224878788	SRR006511.fastq.gz	273Mb	ea41e72037e64a1f63f6063fea19add4	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR006/SRR006511/SRR006511.fastq.gz
#PRJNA33023	SRS001788	SRX001872	SRR006512	null	Caenorhabditis elegans	ILLUMINA	Illumina Genome Analyzer	L2_ce0109_rw001	SINGLE	EST	TRANSCRIPTOMIC	cDNA	60460517	2176578612	SRR006512.fastq.gz	2Gb	f393846e4662ca95ae43d5f37098a615	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR006/SRR006512/SRR006512.fastq.gz
#PRJNA33023	SRS001788	SRX001872	SRR006513	null	Caenorhabditis elegans	ILLUMINA	Illumina Genome Analyzer	L2_ce0109_rw001	SINGLE	EST	TRANSCRIPTOMIC	cDNA	49311467	1775212812	SRR006513.fastq.gz	2Gb	cd256a6d6d8e6b22c3a77daaea5a5fa2	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR006/SRR006513/SRR006513.fastq.gz

# These Study IDs like 'PRJNA33023' don't give the expected '<STUDY alias' lines in get_study_details()
# so we see if we can get a Secondary Study ID like 'SRP000401'
# output should be like:

#run_accession	secondary_study_accession
#SRR006511	SRP000401	

  if ($Study !~ /^SRP/ && $Study !~ /^ERP/) {
    open(ENTRY, "/sw/arch/bin/wget -q -O - 'http://www.ebi.ac.uk/ena/data/warehouse/search?query=%22study_accession=%22${Study}%22%22&result=read_run&limit=1&length=1&offset=0&display=report&fields=secondary_study_accession' |") || $log->log_and_die("Can't get information on Study entry $Study\n");
    while (my $line = <ENTRY>) {
      chomp $line;
      if ($line =~ /^run_accession/) {next}
      @cols = split /\t+/, $line;  
    }
    $Study = $cols[1];
    if (!defined $Study) {$log->log_and_die("Can't get a Study ID like SRPxxxxx for $srx\n")}
  }

  return ($Study, $Sample, $Experiment, \@Run, $Instrument, $Library_Name);
}

##########################################
# add the @RG header to the BAM file
# returns the name of the file created with the @RG head in

sub add_RG_header_to_BAM_file {
  my ($srx, $bamfile) = @_;

  my $Software = "/nfs/panda/ensemblgenomes/wormbase/software/packages";
  my $samtools = "$Software/samtools/samtools";
  my $status;
  my $outfile = "RG_$srx.bam";


## want to write:
## @RG ID:$Run SM:$Sample LB:$Library_Name PG:$Instrument PL:$Instrument
#open (HEAD, ">$dir/header.file") || die "Can't open $dir/header\n";
#print HEAD "\@RG\tID:$Sample\tLB:$Library_Name\tPG:$Instrument\tPL:$Instrument\n";
#close (HEAD);
#
#$wormbase->run_command("$samtools merge -r -h $dir/header.file $dir/$srx.bam $dir/accepted_hits.bam", $log);

  # use Picard to add @RG to the header
  my $AddOrReplaceReadGroups = "$Software/picard/AddOrReplaceReadGroups.jar";
  # the default java on acedb points to jdk1.5, ao explicitly use jdk1.6
  $status = $wormbase->run_command("/sw/arch/pkg/jdk1.6/bin/java -Xmx2g -jar $AddOrReplaceReadGroups  INPUT=$dir/$bamfile OUTPUT=$outfile RGID=$srx RGLB='$Library_Name' RGPL=$Instrument RGPU=$Instrument RGSM=$Sample", $log);
  $log->log_and_die("Problem running AddOrReplaceReadGroups on $dir/$bamfile\n") if ($status != 0);
  return $outfile;
}

##########################################
# make the analysis XML file
# See: http://www.ebi.ac.uk/ena/about/sra_preparing_metadata#bam

sub make_analysis_xml {

  my ($analysis_xmlfile, $Study, $Sample, $Experiment, $Run, $RG_bamfile) = @_;

  my ($study_alias, $study_center) = get_study_details($Study);
  my ($experiment_alias, $center_name) = get_experiment_details($Experiment);

  # write out the start of the XML file
  open (OUT, ">$analysis_xmlfile") || $log->log_and_die("Can't open $analysis_xmlfile\n");
  
  my $ftp_dir  = "RNASeq";

  # want ISO date like: 2011-09-29T13:59:40.0Z
  my $t = localtime;  
  my $date = $t->datetime . '.0Z';
  
  my $g_species = $wormbase->full_name('-g_species' => 1);
  my $version_name = $wormbase->get_wormbase_version_name;
  
  my $UNIQUE_NAME_FOR_ANALYSIS = 'WormBase_'. $g_species . '_' . $srx;
  my $center_name_abbreviation = "EBI";
  my $a_descriptive_title_for_the_analysis_shown_in_search_results = "RNASeq reads mapped using Tophat to the $g_species reference genome $version_name";
  my $a_detailed_description_of_the_analysis = "RNASeq reads from the SRA experiment ID $srx mapped to the $g_species reference genome $version_name as part of the WormBase database Build process using TopHat version tophat-2.0.5.Linux_x86_64";
  my $STUDY_ALIAS_OF_RELEVANT_STUDY_OBJECT = $Study;
  my $INSDC_assembly_accession = $assembly_accessions{$species};
  my $FILENAME = "$RG_bamfile";
  
  
  # get the md5 checksum for the file
  my $md5 = `md5sum $FILENAME`;
  my $CHECKSUM = (split /\s+/, $md5)[0];
  
  
  print OUT qq[<?xml version="1.0" encoding="UTF-8"?>
<ANALYSIS_SET>
    <ANALYSIS alias="$UNIQUE_NAME_FOR_ANALYSIS" 
        center_name="$center_name_abbreviation"
        broker_name="$center_name_abbreviation"
        analysis_center="$center_name_abbreviation" analysis_date="$date">
        <TITLE>$a_descriptive_title_for_the_analysis_shown_in_search_results</TITLE>
        <DESCRIPTION>$a_detailed_description_of_the_analysis</DESCRIPTION>
        <STUDY_REF accession="$Study"/>
        <SAMPLE_REF accession="$Sample"/>
];

  foreach my $run (@{$Run}) {
    print OUT qq[        <RUN_REF accession="${run}"/>
];
  }

print OUT qq[
        <ANALYSIS_TYPE>
            <REFERENCE_ALIGNMENT>
                <ASSEMBLY>
                    <STANDARD accession="$INSDC_assembly_accession"/>
                </ASSEMBLY>
];


  # loop over all chromosome names to give INSDC vs WormBase chromosome IDs
  foreach my $chrom ($wormbase->get_chromosome_names(-mito => 1, -prefix=> 1)) {
    my $chrom_obj = $db->fetch(Sequence => $chrom);
    my $INSDC_sequence_accession = $chrom_obj->at("DB_info.Database.NDB.NDB_AC[1]");
    if (! defined $INSDC_sequence_accession) {$log->log_and_die("Can't find DB_info.Database.NDB.NDB_AC tag in Sequence $chrom\n")}
    my $version = ena_entry_version($INSDC_sequence_accession);
    my $reference_sequence_name_in_the_BAM_file = $chrom; # the chromosome name
    print OUT qq[                <SEQUENCE accession="${INSDC_sequence_accession}.${version}" label="$reference_sequence_name_in_the_BAM_file"/>
];
  }


  print OUT qq[
            </REFERENCE_ALIGNMENT>
        </ANALYSIS_TYPE>
        <FILES>
            <FILE filename="${ftp_dir}/${FILENAME}" filetype="bam" checksum_method="MD5"
                checksum="$CHECKSUM"/>
        </FILES>
    </ANALYSIS>
</ANALYSIS_SET>
];

}

##########################################
# get the EXPERIMENT details from the XML file at the ENA

sub get_experiment_details {

  my ($Experiment) = @_;

#<?xml version="1.0" encoding="UTF-8"?>
#<ROOT request="SRX001872&amp;display=xml">
#<EXPERIMENT alias="L2 polyA  RNAseq" accession="SRX001872" center_name="UWGS-RW" broker_name="NCBI">
#     <IDENTIFIERS>
#          <PRIMARY_ID>SRX001872</PRIMARY_ID>
#          <SUBMITTER_ID namespace="UWGS-RW">L2 polyA  RNAseq</SUBMITTER_ID>
#     </IDENTIFIERS>
#     <TITLE>Illumina Genome Analyzer sequencing; Illumina sequencing of Caenorhabditis elegans mid-L2 transcript fragment library</TITLE>
#     <STUDY_REF accession="SRP000401" refname="Caenorhabditis elegans transcriptome study RW0001" refcenter="UWGS-RW">
#          <IDENTIFIERS>
#               <PRIMARY_ID>SRP000401</PRIMARY_ID>
#               <SUBMITTER_ID namespace="UWGS-RW">Caenorhabditis elegans transcriptome study RW0001</SUBMITTER_ID>
#          </IDENTIFIERS>
#     </STUDY_REF>

  my ($experiment_alias, $center_name);

  open(ENTRY, "/sw/arch/bin/wget -q -O - 'http://www.ebi.ac.uk/ena/data/view/${Experiment}&display=xml' |") || $log->log_and_die("Can't get information on SRA entry $Experiment\n");
  while (my $line = <ENTRY>) {
    chomp $line;
    if ($line !~ /<EXPERIMENT alias/) {next}
    my ($experiment_alias) = ($line =~ /<EXPERIMENT alias=\"(.+?)\"/);
    my ($center_name) = ($line =~ /center_name=\"(\S+)\"/);
    
    return ($experiment_alias, $center_name);
  }
  close(ENTRY);
  return undef;

}

##########################################
# get the STUDY details from the XML file at the ENA

sub get_study_details {

  my ($study) = @_;

# get lines like:
#<?xml version="1.0" encoding="UTF-8"?>
#<ROOT request="SRP000401&amp;display=xml">
#<STUDY alias="Caenorhabditis elegans transcriptome study RW0001" center_name="UWGS-RW" accession="SRP000401" broker_name="NCBI" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
#     <IDENTIFIERS>
#          <PRIMARY_ID>SRP000401</PRIMARY_ID>
#          <SUBMITTER_ID namespace="UWGS-RW">Caenorhabditis elegans transcriptome study RW0001</SUBMITTER_ID>
#     </IDENTIFIERS>
#     <DESCRIPTOR>
#          <STUDY_TITLE>Deep sequencing of the Caenorhabditis elegans transcriptome using RNA isolated from various developmental stages under various experimental conditions RW0001</STUDY_TITLE>
#          <STUDY_TYPE existing_study_type="Transcriptome Analysis"/>
#          <STUDY_ABSTRACT/>


  my @cols;

  open(ENTRY, "/sw/arch/bin/wget -q -O - 'http://www.ebi.ac.uk/ena/data/view/$study&display=xml' |") || $log->log_and_die("Can't get information on SRA entry $study\n");
  while (my $line = <ENTRY>) {
    chomp $line;
    if ($line !~ /<STUDY alias/) {next}
    my ($alias) = ($line =~ /<STUDY alias=\"(.+?)\"/);
    my ($center) = ($line =~ /center_name=\"(\S+)\"/);
    
    return ($alias, $center);
  }
  close(ENTRY);
  return undef;

}

##########################################
# upload a file to the drop-box


sub make_submission_xml {
  my ($validate, $modify, $submission_xmlfile, $analysis_xmlfile) = @_;

  my $g_species = $wormbase->full_name('-g_species' => 1);
  my $version_name = $wormbase->get_wormbase_version_name;
  
  my $UNIQUE_NAME_FOR_ANALYSIS = 'WormBase_' . $g_species . '_' . $srx;
  my $center_name_abbreviation = "EBI";

  my $add_or_validate;
  if ($validate) { # 'VALIDATE' does a test submission
    $add_or_validate = "VALIDATE";
  } elsif ($modify) { # 'MODIFY' changes an existing submission
    $add_or_validate = "MODIFY";
  } else {
    $add_or_validate = "ADD";    # submits a new data file
  }


  # write out the start of the XML file
  open (OUT, ">$submission_xmlfile") || $log->log_and_die("Can't open $submission_xmlfile\n");
  
  print OUT qq[
<SUBMISSION alias="$UNIQUE_NAME_FOR_ANALYSIS" center_name="$center_name_abbreviation">
    <ACTIONS>
      <ACTION>
          <$add_or_validate source="$analysis_xmlfile" schema="analysis"/>
      </ACTION>
 
        <ACTION>
           <RELEASE/> <!-- immediate release -->
        </ACTION>

    </ACTIONS>
</SUBMISSION>
];

}

##########################################
# upload the data


sub data_upload {

  my ($RG_bamfile) = @_;

  my $ftp_user;
  my $ftp_pass;
  my $ftp_host = "ftp.era.ebi.ac.uk";
  my $ftp_dir  = "RNASeq";

  my $login_details_file = $wormbase->wormpub . "/ebi_resources/EBISRAFTP.s";
  open(my $infh, $login_details_file)
      or $log->log_and_die("Can't open secure account details file $login_details_file\n");
  while (<$infh>){
    /^USER_ID:(\S+)$/ and $ftp_user = $1;
    /^PASSWD:(\S+)$/ and $ftp_pass = $1;
  }
  close($infh);

  $log->log_and_die("Could not find both user name and password in login details file\n")
      if not defined $ftp_user or not defined $ftp_pass;

  ###################################################
  # Establish ftp connection                        #
  ###################################################
  my $ftp = Net::FTP->new($ftp_host, Debug => 0) 
    or $log->log_and_die("Cannot connect to $ftp_host: $@");
  $ftp->login($ftp_user,"$ftp_pass\@")
    or $log->log_and_die ("Cannot login to $ftp_host using WormBase credentials\n". $ftp->message);
  $ftp->cwd($ftp_dir) 
    or $log->log_and_die ("Cannot change into dir '$ftp_dir' for upload of files\n". $ftp->message);

  
  ##################################
  # Deposit the data on the FTP site
  ##################################
  $ftp->binary
    or $log->log_and_die ("FTP-binary failed: ".$ftp->message."\n");
  $ftp->put($RG_bamfile) 
    or $log->log_and_die ("FTP-put failed for $RG_bamfile: ".$ftp->message."\n");
  $ftp->quit;
  
  $log->write_to("\nFile: $RG_bamfile uploaded to ENA SRA ftp account\n");



}

##########################################
# submit the data

# All SRA REST submission functionality can be accessed programmatically using curl.
# 
# First, you need to create an authenticated URL. Please follow the steps given below.
# 
# Go to: https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/  (Please don't remove the last slash '/' in the URL)
# Authenticate using your drop box name and password. 
# After successful login you should see the message:  'Submission xml not specified'. 
# Copy the URL. This is your authenticated URL, for example:
# https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ERA%20era-drop...
# 
# This URL can be used in curl to access and autheticate against the SRA REST service.
# 
# For example, one or more studies and samples can be submitted using for following curl command:
# 
# curl -F "SUBMISSION=@submission.xml" -F "STUDY=@study.xml" -F"SAMPLE=@sample.xml"  "https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ERA%20era-drop-30%20blahblahblaho%3D"
# 
# The -F option instructs the SRA REST service to process one of more XML files. Allowed XML file types are: SUBMISSION (mandatory), SAMPLE, STUDY, EXPERIMENT, RUN, ANALYSIS, DAC (for EGA submissions only), POLICY (for EGA submissions only), DATASET (foe EGA submissions only) and PROJECT.
# 
# You will receive a receipt XML informing on the success of the submission with other information including assigned accession numbers.

sub data_submission {

  my ($submission_xmlfile, $analysis_xmlfile) = @_;

  my $result = 0;
  my $analysis_accession;
  my $text;
  my $http;
  if ($justtesting) {
    $http = "https://www-test.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ERA%20era-drop-196%20PKwv9wrv9S9Co5YFMdQhta8FWSU%3D";
  } else {
    $http = "https://www.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ERA%20era-drop-196%20PKwv9wrv9S9Co5YFMdQhta8FWSU%3D";
  }

  my $submit_cmd = "/usr/bin/curl --max-time 600 --retry 3 -F \"SUBMISSION=\@$submission_xmlfile\" -F \"ANALYSIS=\@$analysis_xmlfile\" $http ";
  print "submit command is:\n";
  print "$submit_cmd\n";
  open(SUBMISSION, "$submit_cmd |") || $log->log_and_die("Can't submit SRA entry $srx\n");
  while (my $line = <SUBMISSION>) {
    $text .= $line;
    chomp $line;
    if ($line =~ /success="true"/) {$result=1}
    if ($line =~ /ANALYSIS accession=\"(\S+)\"/) {$analysis_accession = $1}
  }
  close (SUBMISSION);
 
  return ($result, $text, $analysis_accession);
}

##########################################
# get the version number of an entry in ENA

sub ena_entry_version {
  my ($id) = @_;

  my $version;
  my $ua       = LWP::UserAgent->new;
  
  my $query      = "http://www.ebi.ac.uk/ena/data/view/$id&display=text";
  
  # Doing query

  my $qa1 = $ua->get($query);
  die("Can't get URL -- " . $qa1->status_line) unless $qa1->is_success;
  my $text = $qa1->content;

  ($version) = ($text =~ /^ID\s\s\s${id};\sSV\s(\d+);/);

  return $version


}

##########################################
# save the accession number in a file in case it is required in the future
sub save_analysis_accession {
  my ($analysis_accession) = @_;

  open (ACC, "> $accession_file") || $log->log_and_die("Can't open file $accession_file\n");
  print ACC "$analysis_accession\n";
  close(ACC);

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

##########################################




# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - script_template.pl

=head1 USAGE

=over 4

=item script_template.pl  [-options]

=back

This script does...blah blah blah

script_template.pl MANDATORY arguments:

=over 4

=item None at present.

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item xxx (xxx@sanger.ac.uk)

=back

=cut
