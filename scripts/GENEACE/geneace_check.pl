#!/software/bin/perl -w
# 
# geneace_check.pl
#
# Initiated by Keith Bradnam
#
# Script to run consistency checks on the geneace database
#

use strict;
use lib $ENV{"CVS_DIR"};
use Wormbase;
use Ace;
use Ace::Object;
use Carp;
use Getopt::Long;
use GENEACE::Geneace;
use File::Path;

###################################################
# command line options                            # 
###################################################
my ($help, $debug, $test, $class, @classes, $database, $ace, $verbose);
my (@skip_methods, $excludeprojects, $one_variation);

GetOptions ('help'          => \$help,
            'debug=s'       => \$debug,
	    'class=s'       => \@classes,
	    'database=s'    => \$database,
            'ace=s'           => \$ace,
	    'verbose'       => \$verbose,
            'test'          => \$test,
            'skipmethod=s@' =>  \@skip_methods,
	    'excludeprojects' => \$excludeprojects, # don't test the large Allele projects
	    'variation=s'   => \$one_variation, # specify one test variation
	   );

###################################################
# Miscellaneous important variables               # 
###################################################
my $wb = Wormbase->new(-test => $test, -debug => $debug );
# choose database to query: default is /nfs/wormpub/DATABASES/geneace
$database = '/nfs/production/panda/ensemblgenomes/wormbase/DATABASES/geneace' unless $database;
print "Using database $database.\n\n";

my $tace = $wb->tace;          # tace executable path
my $curr_db = '/nfs/production/panda/ensemblgenomes/wormbase/DATABASES/current_DB'; # Used for some cross checking with geneace
my $def_dir = "${database}/wquery";                          # where lots of table-maker definitions are kept

my $rundate = $wb->rundate;                                # Used by various parts of script for filename creation
my $maintainers = join (
    'hinxton\@wormbase.org',
    );

my $log_dir = "$database/logs";                            # some of the many output files are put here (ar2)
my $log;                                                   # main log file for most output
my (%L_name_F_WBP, %L_name_F_M_WBP);                       # hashes for checking Person and Author merging?


##################################################
# other variables                                #
##################################################

# Display help if required
if ($help){&usage('Help')}

# Use debug mode?
if($debug){
  print "\nDEBUG = \"$debug\"\n";
  ($maintainers = "$debug" . '\@ebi.ac.uk');
}

# Open file for acefile output?
my $acefile;
if ($ace){
    $acefile = $ace;
    open(ACE, ">>$acefile") || croak $!;
    system("chmod 777 $acefile");
}
else {
    mkdir "$database/CHECKS" unless ( -e "$database/CHECKS");
    $acefile = "$database/CHECKS/geneace_check.$rundate.$$.ace";  
    open(ACE, ">>$acefile") || croak $!;
    system("chmod 777 $acefile");
}

my $next_build_ver = $wb->get_wormbase_version() + 1 ; # next build number

# Set up other log files
&create_log_files;


######################################################################################################
######################################################################################################
##                                                                                                  ##
##                                       MAIN PART OF SCRIPT                                        ##
##                                                                                                  ##
######################################################################################################
######################################################################################################

# open a connection to database
my $db = Ace->connect(-path  => $database,
		      -program =>$tace) || do { print LOG "Connection failure: ",Ace->error; die();};


my $ga = init Geneace($wb);
# hash for converting locus/seq. names <-> gene_id
my $Gene_info = $ga -> gene_info($database, "seq2id");
my %Gene_info = %{$Gene_info};

#basic checks to always run
#bad WBPaper_IDs
foreach my $bpaper ($db->fetch(-query=>'Find Paper WHERE !WBPaper0*')) {
    unless ($bpaper =~ /^WBPaper\d{8}/){
	print LOG "ERROR: $bpaper is not valid\n";
	print "ERROR: $bpaper is not valid\n";
    }
}


#bad WBPerson_IDs
foreach my $bperson ($db->fetch(-query=>'Find Person !WBPerson*')){
    unless ($bperson =~ /WBPerson\d+/){
        print LOG "ERROR: $bperson is not valid\n";
	print "ERROR: $bperson is not valid\n";
    }
}



# Process separate classes if specified on the command line else process all classes, 
@classes = ("gene", "laboratory", "evidence", "allele", "strain", "rearrangement", "feature") if (!@classes);
foreach $class (@classes){
  if ($class =~ m/gene/i)          {&process_gene_class}
  elsif ($class =~ m/laboratory/i)    {&process_laboratory_class}
  elsif ($class =~ m/evidence/i)      {&check_evidence}
  elsif ($class =~ m/allele/i)        {&process_allele_class}
  elsif ($class =~ m/strain/i)        {&process_strain_class}
  elsif ($class =~ m/rearrangement/i) {&process_rearrangement}
  elsif ($class =~ m/paper/i)         {&process_paper_class}
  elsif ($class =~ m/feature/i)       {&process_feature_class}
#  if ($class =~ m/mapping/i)       {&check_genetics_coords_mapping}
#  if ($class =~ m/multipoint/i)    {&check_dubious_multipt_gene_connections}
  else {print LOG "\nNo valid class has been specified CLASS:$class\n";}
}


#######################################
# Tidy up and mail relevant log files #
#######################################

$db->close;
print LOG "\nEnded at ",`date`,"\n";
close(LOG)
or warn $! ? "Error closing log: $!"
                   : "Exit status $? log close";
print "Exit status $? log close\n\n" if ($debug);


close(ACE) if ($ace);
# email everyone specified by $maintainers
$wb->mail_maintainer('geneace_check:',$maintainers,$log);

exit(0);
#--------------------------------------------------------------------------------------------------------------------




                          ###########################################################
                          #         SUBROUTINES FOR -class gene option              #
                          ###########################################################

sub process_gene_class{

  print "Checking Gene class for errors:\n";
  print LOG "\nChecking Gene class for errors:\n--------------------------------\n";

  #check if any genes in operons have been killed
  &check_operons;

  # Can first check general errors by grabbing sets of genes for querying

  # check that there is a Version tag
  foreach my $gene ($db->fetch(-query=>'Find Gene WHERE NOT Version')){
    print LOG "ERROR: $gene ($Gene_info{$gene}{'Public_name'}) has no Version number\n";    
  }

  # check that there is a Status tag
  foreach my $gene ($db->fetch(-query=>'Find Gene WHERE NOT Status')){
    print LOG "ERROR: $gene ($Gene_info{$gene}{'Public_name'}) has no Status tag\n";    
  }

  # test for Other_name tag but no value
  foreach my $gene ($db->fetch(-query=>'Find Gene WHERE Other_name AND NOT NEXT')){
    print LOG "ERROR: $gene ($Gene_info{$gene}{'Public_name'}) has 'Other_name' tag without value\n";
  }

  # test for Public_name tag but no value          
  foreach my $gene ($db->fetch(-query=>'Find Gene WHERE Public_name AND NOT NEXT')){
      print LOG "ERROR: $gene ($gene) has 'Public_name' tag without value\n";
  }

  # test for Public_name tag in live WB genes 
  foreach my $gene ($db->fetch(-query=>'Find Gene "WBGene*" WHERE !Public_name AND Live')){
      print LOG "ERROR: $gene ($gene) is live but has no 'Public_name'\n";
  }


  # checks that when a Gene belongs to a Gene_class, it should have a CGC_name
  foreach my $gene ($db->fetch(-query=>'Find Gene WHERE Gene_class AND NOT CGC_name')){
    my $gc = $gene->Gene_class;
    if(!defined($gc)){
      print LOG "ERROR: $gene ($Gene_info{$gene}{'Public_name'}) has no CGC_name but has an unpopulated Gene_class tag\n";
    }
    else{
      print LOG "ERROR: $gene ($Gene_info{$gene}{'Public_name'}) has no CGC_name but links to Gene_class $gc\n";
    }
  }

  # checks existence of a CGC name but no gene_class
  foreach my $gene ($db->fetch(-query=>'Find Gene WHERE CGC_name AND NOT Gene_class')){
    my $cgc_name = $gene->CGC_name;
    print  LOG "ERROR: $gene has CGC name ($cgc_name) but is not linked to its Gene_class\n";
  }

  # checks Genes that do not have a map position nor an interpolated_map_position but has sequence info
  my $query = "Find Live_genes WHERE !(Map | Interpolated_map_position) & Sequence_name & Species=\"*elegans\" & !Positive_clone=\"MTCE\" & !Made_into_transposon";
  foreach my $gene ($db->fetch(-query=>"$query")){
    if( $Gene_info{$gene} ) {
      print LOG "ERROR: $gene ($Gene_info{$gene}{'Public_name'}) has neither Map nor Interpolated_map_position info but has Sequence_name\n";
    }
    else {
      print LOG "$gene not in Gene_info hash\n";
    }
  }

  # test for Map tag and !NEXT
  foreach my $gene ($db->fetch(-query=>"Find Gene WHERE Map AND NOT NEXT")){
    print LOG "ERROR: $gene ($Gene_info{$gene}{'Public_name'}) has a 'Map' tag without a value\n";
  }

  # test for Interpolated_map_position tag and !NEXT
  foreach my $gene ($db->fetch(-query=>"Find Gene WHERE Interpolated_map_position AND NOT NEXT")){
    print LOG "ERROR: $gene ($Gene_info{$gene}{'Public_name'}) has an 'Interpolated_map_position' tag without a value\n";
  }

  # test for Other_sequence tag and no value
  foreach my $gene ($db->fetch(-query=>"Find Other_sequence AND NOT NEXT")){
    print LOG "ERROR(a): $gene ($Gene_info{$gene}{'Public_name'}) has Other_sequence tag but no associated value\n";
    if ($ace){print ACE "\n\nGene : \"$gene\"\n-D Other_sequence\n";}
  }


  # A seq. should normally be linked to only one gene id: there are cases where .a is WBGenex and .b is WBGy (eg, ZC416.8a/b)
  # The query is to find sequences (CDS/Transcript/Pseudogene) that have more than one Sequence_name_for values
  # this tells you if a gene id is linked to > 1 sequences
  foreach my $gene_name ($db->fetch(-query=>"Find Gene_name WHERE COUNT Sequence_name_for > 1")){

    # skip hard-coded exceptions for eat-18/lev-10, cha-1/unc-17 & B0564.1/tin-9.2 
    next if ($gene_name eq "Y105E8A.7" || $gene_name eq "ZC416.8" || $gene_name eq "B0564.1");

    my @gene_ids = $gene_name->Sequence_name_for;
    print LOG "ERROR: $gene_name is connected to multiple gene IDs: @gene_ids\n";
  }

  # Look for missing Method tag for Live genes
  foreach my $gene ($db->fetch(-query=>"Find Gene WHERE Live AND NOT Method")){
    print LOG "ERROR: $gene is a Live gene but has no 'Method' tag\n";
  }

  # Look for Method tag but no Method field after it
  foreach my $gene ($db->fetch(-query=>"Find Gene WHERE Method AND NOT NEXT")){
    print LOG "ERROR: $gene has a 'Method' tag but no value\n";
  }


  # test for missing Species tag
  foreach my $gene ($db->fetch(-query=>"Find Gene WHERE NOT Species")){
    print LOG "ERROR(a): $gene ($Gene_info{$gene}{'CGC_name'}) has no 'Species' info\n";
    if ($ace){
      print ACE "\n\nGene : \"$gene\"\n";
      if ( $Gene_info{$gene}{'Public_name'} !~ /^Cb-|^Cr-|^Cv-/ ){
        print ACE "Species \"Caenorhabditis elegans\"\n";
      }
      if ( $Gene_info{$gene}{'Public_name'} =~/^Cb-.+/ ){
        print ACE "Species \"Caenorhabditis briggsae\"\n";
      }
      if ( $Gene_info{$gene}{'Public_name'} =~/^Cr-.+/ ){
        print ACE "Species \"Caenorhabditis remanei\"\n";
      }
      if ( $Gene_info{$gene}{'Public_name'} =~/^Cv-.+/ ){
        print ACE "Species \"Caenorhabditis vulgaris\"\n";
      }
    }
  }

  # Look for Species tag but no Species field after it
  foreach my $gene ($db->fetch(-query=>'Find Gene WHERE Species AND NOT NEXT')){
    print LOG "ERROR: $gene has a 'Species' tag but no value\n";
 }

  # checks that a gene with alleles are not dead (i.e. merged into something else)
  foreach my $gene ($db->fetch(-query=>'Find Gene WHERE Dead AND Allele')){
    print LOG "ERROR: Mama mia! $gene is dead but is still connected to an allele\n";
  }

  # checks that a gene with references are not dead (i.e. merged into something else)
  foreach my $gene ($db->fetch(-query=>'Find Gene WHERE Dead AND Reference')){
      print LOG "ERROR: Oh my sainted Aunt! $gene is dead but is still connected to a reference\n";
  }

  # checks that a gene with orthologs has not been merged into something else ie Dead. nb. It is OK for Transposon_CDSs to have Orthologs
  foreach my $gene ($db->fetch(-query=>'Find Gene WHERE NOT Live AND Merged_into AND Ortholog')){
      print LOG "ERROR: Zut alors! $gene has been merged (and is Dead) but is still connected to an ortholog\n";
  }
  
  # checks that a gene with mapping data is not dead.
  foreach my $gene ($db->fetch(-query=>'Find Gene WHERE Dead AND Map_info')){
    print LOG "ERROR: ERROR: Bleedin' Nora! $gene is dead but still has mapping data\n";
  }

  # checks that a Gene doesn't have both Map and Interpolated_map_position tags
  foreach my $gene ($db->fetch(-query=>'Find Gene WHERE Map AND Interpolated_map_position')){
    print LOG "ERROR: $gene has both Map and Interpolated_map_position tags, are you crazy?\n";
  }

  # checks for genes that have no Live tag but a split_from tag
  foreach my $gene ($db->fetch(-query=>'Find Gene WHERE Split_from AND Dead AND NOT Merged_into AND NOT Made_into_transposon')){
    print LOG "ERROR: $gene has Split_from tag but no Live tag\n" unless ($gene eq "WBGene00195119");
  }



  # get all genes
  my @gene_ids = $db->fetch(-class => 'Gene',
	                    -name  => 'WBGene*');

  # now loop through looking for specific errors that might be in any Gene object
  foreach my $gene_id (@gene_ids){
    # useful to see where you are in the script if running on command line

    my $species = $gene_id->Species;
    if ($verbose) {
      print "$gene_id:";
      if($species eq 'Caenorhabditis elegans' and exists $Gene_info{$gene_id}{'Public_name'}) {
        print $Gene_info{$gene_id}{'Public_name'};
      }
      print "\n";
    }
    &test_locus_for_errors($gene_id);
  }

  
  # check highest gene id number = gene obj. number
  my $last_gene_id = $gene_ids[-1];
  $last_gene_id =~ s/WBGene(0)+//;
  
  if ( scalar @gene_ids != $last_gene_id ){
    print LOG "ERROR: The highest gene id ($last_gene_id) is not equal to total number of Gene objects (", scalar @gene_ids, ")";
  }

}


sub test_locus_for_errors{
  my $gene_id = shift;
  if ($debug) {print "$gene_id\n";}
  my $warnings;


  # check the number of version and Version_change of a Gene is identical
  if ( defined $gene_id->Version && defined $gene_id->Version_change ){
    my @ver_changes = $gene_id->Version_change;
    my $ver = $gene_id->Version;
    if ( "$ver" ne "$ver_changes[-1]" ){
      my $dname = (exists $Gene_info{$gene_id}{'Public_name'}) ? $Gene_info{$gene_id}{'Public_name'} : "No public name";
      $warnings .= "ERROR: $gene_id ($dname) has Version problem: current ($ver) => history ($ver_changes[-1])\n";
    }
  }



  # test for Public_name is different from CGC_name
  if ( defined $gene_id->CGC_name && defined $gene_id->Public_name ){				
    my $cgc_name = $gene_id->CGC_name;
    my $pub_name = $gene_id->Public_name;
    if ( $cgc_name ne $pub_name ){
      $warnings .= "ERROR(a): $gene_id ($Gene_info{$gene_id}{'CGC_name'}) has a Public_name ($pub_name) different from its CGC_name ($cgc_name)\n";
      if ($ace){
	print ACE "\nGene : \"$gene_id\"\n";
	print ACE "Public_name \"$Gene_info{$gene_id}{'CGC_name'}\"\n";
      }
    }
  }

  # test for missing Public_name and assign one if so
  if( !defined $gene_id->Public_name && (defined $gene_id->CGC_name || defined $gene_id->Sequence_name || defined $gene_id->Other_name) ){
    if ($gene_id->at('Identity.Live')){
	$warnings .= "ERROR(a): $gene_id has no Public_name but has CGC/Sequence/Other_name\n";
	
	if ($ace){
	    print ACE "\nGene : \"$gene_id\"\n";
	    print ACE "Public_name \"$Gene_info{$gene_id}{'CGC_name'}\"\n" if exists $Gene_info{$gene_id}{'CGC_name'};
	    print  "Public_name \"$Gene_info{$gene_id}{'CGC_name'}\"=== E6\n" if exists $Gene_info{$gene_id}{'CGC_name'};
	    print ACE "Public_name \"$Gene_info{$gene_id}{'Sequence_name'}\"\n" if !exists $Gene_info{$gene_id}{'CGC_name'} && exists $Gene_info{$gene_id}{'Sequence_name'};
	    if ( !exists $Gene_info{$gene_id}{'CGC_name'} && !exists $Gene_info{$gene_id}{'Sequence_name'} && exists $Gene_info{$gene_id}{'Other_name'} ){
		print ACE "Public_name \"$Gene_info{$gene_id}{'Other_name'}\"\n"
	    }
	}
    }
    else {
	print "Warning: $gene_id has no Public_name but has CGC/Sequence/Other_name but is Dead so probably OK\n";
    }
  }

  if( !defined $gene_id->Public_name && !defined $gene_id->CGC_name && !defined $gene_id->Sequence_name && !defined $gene_id->Other_name ){
      $warnings .= "WARNING: $gene_id has no name information please check it was merged into a cgc_named gene (266 known issues)\n" unless  (($gene_id->Remark =~ /[MERGED INTO UNCLONED GENE]/) || ($gene_id->Merged_into));
  }

  # test for discrepancy betw. CGC_name and Gene_class name, ie, for aap-1, the gene_class should be aap
  # The idea is that when two Gene objects are merged, the latter one will overwrite the first one, 
  # as the gene_class field is UNIQUE

  if( defined $gene_id->Gene_class && defined $gene_id->CGC_name ){
    my $cgc_name = $gene_id->CGC_name;
    my $gc_name = $cgc_name;
    if ( $gc_name =~ /^\w+-\w+-/){
	$gc_name =~ s/^\w+-//;
       }
    $gc_name =~ s/-.+//;
    my $gc = $gene_id->Gene_class;

    if ($gc ne $gc_name){
      $warnings .= "ERROR: $gene_id ($cgc_name) has incorrect Gene_class $gc\n";
    }
  }


  # check that live gene id should not have wpxxx appended to its Sequence_name
  foreach my $tag ("Sequence_name", "Public_name"){
    if ( defined $gene_id->$tag && $gene_id->at('Identity.Live') ){
      my $history = $gene_id->$tag;
      if ( $history =~ /:wp\d+/ ){
	$warnings .= "ERROR: $gene_id has $tag with :wpxxx history name appended\n";
      }
    }
  }

  # Look for Genes with Live tag and no Positive_clone info but which can be derived from its sequence info
  if( !defined $gene_id->Positive_clone(1) && defined $gene_id->Sequence_name){
    # don't need to do this for C. briggsae genes
    my $species = $gene_id->Species;
    if ($species eq 'Caenorhabditis elegans'){
      my $seq = $gene_id->Sequence_name;
      
      # don't need to do this for Dead genes.
      my $identity = $gene_id->Status;
      if ($identity eq 'Live'){
	my $seq = $gene_id->Sequence_name;
	
	# need to chop off the ending to just get clone part
	$seq =~ s/\..*//;
	$warnings .= "ERROR(a): $gene_id ($Gene_info{$gene_id}{'Public_name'}) has no Positive_clone but info is available from Sequence_name $seq\n";
	
	# must use Inferred_automatically from Evidence hash for this type of info
	print ACE "\n\nGene : \"$gene_id\"\nPositive_clone \"$seq\" Inferred_automatically \"From sequence, transcript, pseudogene data\"\n" if ($ace);
      }
    }
  }



  # Look for genes with Live tag but also with Merged_into or Killed tags
  if ( defined $gene_id->at('Identity.Live')){

    my @merge =  $gene_id->Merged_into;
    if (@merge){
      $warnings .= "ERROR(a): $gene_id ($Gene_info{$gene_id}{'Public_name'}) has 'Live' tag but also has been merged into $merge[-1] ($Gene_info{$merge[-1]}{'Public_name'}) => Lose Live\n";
      if ($ace){
	print ACE "\nGene : \"$gene_id\"\n";
	print ACE "-D Live\n";
      }
    }
  }
    
  if ( $gene_id->Version_change ){
    my $tag = &get_event($gene_id);
    if (!defined($tag)){
      $warnings .= "ERROR: $gene_id has an Event tag but does not have any following tag\n";
    }
    elsif ($tag eq 'Killed'){
      if ( defined $gene_id->at('Identity.Live')){
	$warnings .= "ERROR(a): $gene_id ($Gene_info{$gene_id}{'Public_name'}) has 'Live' tag but also has 'Killed' tag => Lose Live\n";
	if ($ace){
	  print ACE "\nGene : \"$gene_id\"  \/\/20.2\n";
	  print ACE "-D Live\n";
	}
      }
    }
  }
  

  # checks that a Gene has identical CGC_name and Other_name
  if ( defined $gene_id->CGC_name && defined $gene_id->Other_name ){
    my @other_names =  $gene_id->Other_name;
    foreach my $o (@other_names) {
      if ( $o eq $gene_id->CGC_name ){
	$warnings .= "ERROR: $gene_id (" . $gene_id->CGC_name . ") has an identical Other_name\n";
      }
    }
  }

  print LOG "$warnings" if(defined($warnings));
}


sub get_event {
  my $gene = shift;

  my $version = $gene->Version;
  my $tag;

  my $vc = $gene->Version_change;
  while ( $vc->down ){
    $vc = $vc->down;
  }

  my $right_step = 0;
   do { 
     if( $vc->right ) {
       $vc = $vc->right;
     }
     else {
       print LOG $gene->name," has data missing from the Version_change info\n";
       return 
     }
     $right_step++;
   }while( $right_step < 4 );

  $tag = $vc->name;

  return $tag;
}






                          #############################################################
                          #         SUBROUTINES FOR -class evidence option            #
                          #############################################################

sub check_evidence {

  my $WBPerson_F_M_L_names="Table-maker -p \"$def_dir/WBperson_first_middle_last_names.def\"\nquit\n";
  my @WBPerson = &process_WBPerson_names($WBPerson_F_M_L_names, $curr_db);

  print "\n\nChecking misuse of Evidence and converting Author to Person / Non-Person to Author:\n\n";
  print LOG "\n\nChecking misuse of Evidence and converting Author to Person / Non-Person to Author:\n";
  print LOG "-----------------------------------------------------------------------------------\n";


# dump flat files with time stamps
my $dump_dir = "$database/CHECKS";

my $command=<<END;
Find Gene * 
show -a -T -f $dump_dir/gene_dump.ace

Find Variation *
show -a -T -f $dump_dir/variation_dump.ace

Find strain *
show -a -T -f $dump_dir/strain_dump.ace

Find gene_class *
show -a -T -f $dump_dir/geneclass_dump.ace

Find 2_point_data *
show -a -T -f $dump_dir/2_pt_dump.ace

Find Multi_pt_data *
show -a -T -f $dump_dir/multi_pt_dump.ace

Find Pos_neg_data *
show -a -T -f $dump_dir/posneg_dump.ace

Find Laboratory *
show -a -T -f $dump_dir/lab_dump.ace

quit
END

  open (DUMP, "| $tace $database") || warn "Failed to connect to Geneace";
  print DUMP $command;
  close DUMP;

  `cat $dump_dir/gene_dump.ace $dump_dir/variation_dump.ace $dump_dir/strain_dump.ace $dump_dir/2_pt_dump.ace $dump_dir/multi_pt_dump.ace $dump_dir/posneg_dump.ace $dump_dir/geneclass_dump.ace $dump_dir/lab_dump.ace > $dump_dir/class_dump.ace`;

  my @dumps = qw (gene_dump.ace allele_dump.ace strain_dump.ace geneclass_dump.ace 2_pt_dump.ace multi_pt_dump.ace posneg_dump.ace class_dump.ace lab_dump.ace);

  foreach (@dumps){system ("chmod 777 $dump_dir/$_")}

  open(IN, "$dump_dir/class_dump.ace") || warn $!;

# look for person/author names that needs to be converted in 8 classes (regardless of which tag as the scripts greps from string pattern

  my $evid_errors = 0;
  my $updates = 0;
  my $info_num = 0;
  my (@counters, $class_obj, $class, $obj, $tag, $ori, $b4_evi, $name, $paper, $author, $last_name, $initials);

  while (<IN>){
    chomp;
    if ($_ =~ /^(Gene) : \"(.+)\" -O .+/){$class = $1; $obj = $2}
    if ($_ =~ /^(Variation) : \"(.+)\" -O .+/){$class = $1; $obj = $2}
    if ($_ =~ /^(Strain) : \"(.+)\" -O .+/){$class = $1; $obj = $2}
    if ($_ =~ /^(Gene_class) : \"(.+)\" -O .+/){$class = $1; $obj = $2}
    if ($_ =~ /^(2_point_data) : \"(.+)\" -O .+/){$class = $1; $obj = $2}
    if ($_ =~ /^(Multi_pt_data) : \"(.+)\" -O .+/){$class = $1; $obj = $2}
    if ($_ =~ /^(Pos_neg_data) : \"(.+)\" -O .+/){$class = $1; $obj = $2}
    if ($_ =~ /^(Laboratory) : \"(.+)\" -O .+/){$class = $1; $obj = $2}

    if ($_ =~ /((\w+)\s+.+)Person_evidence -O .+\"(.+)\" -O.+/){
      $ori = $_;
      $b4_evi = $1;	
      $tag = $2;
      $name = $3;
      if ($name !~ /WBPerson\d+/){
	#$evid_errors++;
	print LOG "\nERROR: $class $obj has non-Person $name under main tag $tag\n";

	@counters = get_WBPerson_ID($name, $class, $obj, $tag, $ori, $b4_evi, "PtoA");
	$evid_errors += $counters[0];
	$updates += $counters[1];
	$info_num += $counters[2];
      }
    }
    if ($_ =~ /(\w+)\s+.+Paper_evidence -O .+\"(.+)\" -O.+/){
      $ori = $_;
      $tag = $1;
      $paper = $2;
      if ($paper !~ /WBPaper\d+/){
        $evid_errors++;
	print LOG "\nERROR: $class $obj has Paper $paper under main tag $tag\n";
      }
    }

    if ($_ =~ /((\w+)\s+.+)Author_evidence -O .+\"(.+)\" -O.+/){
      $ori = $_;
      $b4_evi = $1;
      $tag = $2;
      $author = $3;

      @counters = get_WBPerson_ID($author, $class, $obj, $tag, $ori, $b4_evi);
      $evid_errors += $counters[0];
      $updates += $counters[1];
      $info_num += $counters[2];
    }

    if ($_ =~ /(Unregistered_lab_members)\s+-O\s\"\d+.+_\w+\" \"(.+)\" -O .+/){
      $ori = $_;
      $tag = $1;
      $author = $2;
      #print $author, "#\n";
      # print "$author, $class, $obj, $tag, $ori@\n";
      @counters = get_WBPerson_ID($author, $class, $obj, $tag, $ori);
      $evid_errors += $counters[0];
      $updates += $counters[1];
      $info_num += $counters[2];
    }

  }
  print LOG "\n\nThere are $evid_errors Evidence errors in 8 classes checked\n";
  print LOG "\n$updates Authors can be converted to Persons\n";
  print LOG "\n$info_num Authors are not Persons\n" if $verbose;
  system ("rm -f $dump_dir/*_dump.ace");

 #########################################################################
  # functions called by check_evidence() to do Author to Person convertion
  #########################################################################

  # make hashes for last name and its corresponding first/middle name and WBPersonID
  # these hashes are used by get_WBPerson_ID subroutine

  sub process_WBPerson_names {
    my ($def, $db)=@_;
    my ($WBPerson, $F_name, $M_name, $L_name, $F_char, $M_char);
    open (FH, "echo '$def' | $tace $db | ") || warn "Couldn't access current_DB\n";
    while (<FH>){
      chomp($_);
      if ($_ =~ /^\"(WBPerson\d+)\"\s+\"(\w+)\"\s+\"(\w+|\w+.)\"\s+\"(\w+|\w+-\w+)\"$/){
	$WBPerson = $1;
	$F_name = $2;
	$M_name = $3;
	$L_name = $4;
	$F_char = substr ($F_name, 0, 1);
	$M_char = substr ($M_name, 0, 1);
	push (@{$L_name_F_M_WBP{$L_name}}, $F_char.$M_char, $F_char, $WBPerson);
      }
      if ($_ =~ /^\"(WBPerson\d+)\"\s+\"(\w+)\"\s+\"(\w+|\w+-\w+)\"$/){
	$WBPerson = $1;
	$F_name = $2;
	$L_name = $3;
	$F_char = substr ($F_name, 0, 1);
	push (@{$L_name_F_WBP{$L_name}}, $F_char, $WBPerson);
      }
    }
    close FH;
  }

  # last name (key) is checked against hashes created from process_WBPerson_names subroutine
  # the values of the key (first name, middle name) are checked to assign identify

  # convert authors that have WBPersonID to person 
  # (or unregistered_lab_members to registered_lab_members for laboratory class)

  # move names in person_evidence that are not person to author_evidence

  sub get_WBPerson_ID {
    my ($name, $class, $obj, $tag, $ori, $b4_evi, $conversion) = @_;
    my ($last_name, $initials, $num);
    my $class_obj = $class." : "."\"$obj\"";
    my $convert = 0;
    my $info_count = 0;
    my $evid_errors = 0;

    if ($name =~ /\w+\,\w+/ || $name =~ /\w+\,\s+\w+/){
      $evid_errors++;
      if (!defined $conversion){
	($last_name, $initials) = split(/ /, $name);
	$last_name =~ s/,//;
	print LOG "\nERROR(a): $class $obj has $name (Author) under main tag $tag\n";
	print LOG"=====>Correct $name as $last_name $initials\n";
	if ($ace){
	  print ACE "\n$class_obj\n";
	  print ACE "-D $ori\n\n";
	  print ACE "\n$class_obj\n";
	  if (defined $b4_evi){
	    print ACE "$b4_evi Author_evidence \"$last_name $initials\"\n";
	  }
	  else {
	    print ACE "$tag \"$last_name $initials\"\n";
	  }
	}
      }
    }
    else {
      ($last_name, $initials) = split(/ /, $name);
    }

    if (!exists $L_name_F_WBP{$last_name} && !exists $L_name_F_M_WBP{$last_name}){
      if (defined $conversion){
	print LOG "=====>Move $name under Author_evidence as NO corresponding WBPersonID exists\n";
      }
      if ($verbose && !defined $conversion){
	$info_count++; 
	print LOG "\nINFO: $class $obj has $name (Author under $tag tag): NOT yet a Person\n";
      }
      if ($ace && defined $conversion){
	print ACE "\n$class_obj\n";
	print ACE "-D $ori\n\n";
	print ACE "\n$class_obj\n";
	if (defined $b4_evi){
	  print ACE "$b4_evi Author_evidence \"$last_name $initials\"\n";
	}
      }
    }

    if (exists $L_name_F_WBP{$last_name}){
      $num = scalar @{$L_name_F_WBP{$last_name}};

      for (my $i=0; $i< $num; $i=$i+2){
	if ($initials eq ${@{$L_name_F_WBP{$last_name}}}[$i]){
	  $convert++;
	  if (!defined $conversion){
	    print LOG "\nUPDT(a): $class $obj has $name (Author) under $tag tag\n";
          }
	  if ($num == 2){
	    print LOG "=====>$name can now be Person ${@{$L_name_F_WBP{$last_name}}}[$i+1]\n";
	    if ($ace){
              print ACE "\n$class_obj\n";
              print ACE "-D $ori\n\n";
	      print ACE "\n$class_obj\n";
	      if (defined $b4_evi){
		print ACE "$b4_evi"." Person_evidence \"${@{$L_name_F_WBP{$last_name}}}[$i+1]\"\n";
	      }
	      else {
		print ACE "Registered_lab_members \"${@{$L_name_F_WBP{$last_name}}}[$i+1]\"\n";
	      }	
            }
	  }
	  else {
	    print LOG "=====>$name might be Person ${@{$L_name_F_WBP{$last_name}}}[$i+1]\n";
	    if ($ace){
              print ACE "\n$class_obj\n";
              print ACE "-D $ori\n\n";
	      print ACE "\n$class_obj\n";
	      if (defined $b4_evi){
		print ACE "$b4_evi"." Person_evidence \"${@{$L_name_F_WBP{$last_name}}}[$i+1]\"\n";
	      }
	      else {
		print ACE "Registered_lab_members \"${@{$L_name_F_WBP{$last_name}}}[$i+1]\"\n";
	      }	
            }
	  }
	}
      }
    }

    if (exists $L_name_F_M_WBP{$last_name}){
      $num = scalar @{$L_name_F_M_WBP{$last_name}};

      for (my $i=0; $i< $num; $i=$i+3){
	if ($initials eq ${@{$L_name_F_M_WBP{$last_name}}}[$i] ||
            $initials eq ${@{$L_name_F_M_WBP{$last_name}}}[$i+1] ){
          $convert++;
          if (!defined $conversion){
	    print LOG "\nUPDT(a): $class $obj has $name (Author) under $tag tag\n";
          }
       	  if ($num == 3){
	    print LOG "=====>$name can now be Person ${@{$L_name_F_M_WBP{$last_name}}}[$i+2]\n";
            if ($ace){
              print ACE "\n$class_obj\n";
              print ACE "-D $ori\n\n";
	      print ACE "\n$class_obj\n";
	      if (defined $b4_evi){
		print ACE "$b4_evi"." Person_evidence \"${@{$L_name_F_M_WBP{$last_name}}}[$i+2]\"\n";
	      }
	      else {
		print ACE "Registered_lab_members \"${@{$L_name_F_M_WBP{$last_name}}}[$i+2]\"\n";
	      }
            }
	  }
	  else {
	    print LOG "=====>$name might be Person ${@{$L_name_F_M_WBP{$last_name}}}[$i+2]\n";
            if ($ace){
              print ACE "\n$class_obj\n";
              print ACE "-D $ori\n\n";
	      print ACE "\n$class_obj\n";
	      if (defined $b4_evi){
		print ACE "$b4_evi"." Person_evidence \"${@{$L_name_F_M_WBP{$last_name}}}[$i+2]\"\n";
	      }
	      else {
		print ACE "Registered_lab_members \"${@{$L_name_F_M_WBP{$last_name}}}[$i+2]\"\n";
	      }
            }
	  }
	}
      }
    }
    $conversion =();
    return $evid_errors, $convert, $info_count;
  }
}	

                          ###############################################################
                          #         SUBROUTINES FOR -class laboratory option            #
                          ###############################################################


sub process_laboratory_class{

  print "\n\nChecking Laboratory class for errors:\n";
  print LOG "\n\nChecking Laboratory class for errors:\n";
  print LOG "-------------------------------------\n";

  my @labs = $db->fetch(-class => 'Laboratory',
		        -name  => '*');
  my $e = 0;
  # test for Allele_designation and Representative tags
  foreach my $lab (@labs){
    $e = 0;
    if ((defined $lab->Remark) && ($lab->Remark =~ /closed/)) {
      print LOG "INFO: $lab ".$lab->Remark."\n" if ($verbose);
    }
    else {
      if ( $verbose && (!defined $lab -> Allele_designation) && $lab =~ /[A-Z]{3,3}/ ){
	print LOG "INFO: $lab has no Allele_designation (exception)\n";
	$e++;
      }
      if ( (!defined $lab -> Allele_designation) && $lab !~ /[A-Z]{3,3}/ ){
	print LOG "ERROR: $lab has no Allele_designation\n";
	$e++;
      }
      if( (!defined $lab->Representative) && $lab !~ /DR/ && $lab !~ /XA|CGC/ ){
	print LOG "ERROR: $lab has no Representative\n";
	$e++;
      }
      if( defined $lab->at('Laboratory.CGC.Representative') && !defined $lab->Representative ){
	print LOG "ERROR: $lab has Representative tag without value\n";
	$e++;
      }
      if ( defined $lab->Representative && $lab -> Representative !~ /^WBPerson\d+/ ){
	print LOG "$lab representative (", $lab->Representative, ") is not using WBPerson# format\n";
	$e++;
      }
      print LOG $lab->Remark.": $lab Remark\n\n" if (($e ne 0) && (defined $lab->Remark));
    }
  }
}


                          ###########################################################
                          #         SUBROUTINES FOR -class allele option            #
                          ###########################################################


sub process_allele_class{

  print"\n\nChecking Allele class for errors:\n";
  print LOG "\n\nChecking Allele class for errors:\n";
  print LOG "---------------------------------\n";

  # make hash of allele to lab connections
  my %allele2lab;
  my $def="Table-maker -p \"$def_dir/allele_designation_to_LAB.def\"\nquit\n";

  open (FH, "echo '$def' | $tace $database | ") || warn "Couldn't access $db\n";
  while (<FH>) {
      print if ($debug);
      chomp;
      if ($_ =~ /^\"(.+)\"\s+\"(.+)\"/) {
	  %allele2lab->{$2}=$1;	# $2 is allele_designation $1 is lab designation	
      }
  }

  my $query = "find Variation Allele";
  #my $query = "find Variation WBVar00088961";

  if ($excludeprojects) {@skip_methods = qw(NBP_knockout_allele KO_consortium_allele NemaGENETAG_consortium_allele Million_mutation SNP Mos_insertion Transposon_insertion CGH_allele WGS_Hobert WGS_Rose WGS_Jarriault WGS_McGrath WGS_Flibotte)}
  foreach my $meth (@skip_methods) {
    $query .= " AND Method != \"$meth\"";
  }
  if (defined $one_variation) {
    $query = "find Variation $one_variation";
  }

  my $alleles_it = $db->fetch_many(-query => "$query");

  # now loop through all alleles looking for problems
  while( my $allele = $alleles_it->next){
    print "$allele\n" if ($verbose);

    # check allele has no laboratory tag
    if (!defined $allele->Laboratory ) {
      if (!$ace) {
	print LOG "ERROR: $allele has no Laboratory tag\n";
      } else {
	print LOG "ERROR(a): $allele has no Laboratory tag\n";
	if ($allele->Public_name->name =~ /^WBVar/) {
	    print LOG "\tWARNING: WBVar named allele cannot fix Laboratory\n";
	    next;
	}
	# try to find lab designation for alleles with CGC-names (1 or 2 letters and then numbers)
	if ($allele->Public_name->name =~ /^([a-z]+)\d+$/) {
	  my $allele_name = $1;	  
	  my $allele_ID = $allele->name;
	  if (defined $allele2lab{$allele_name}){
	      print ACE "\nVariation : \"$allele_ID\"\n";
	      print ACE "-D Laboratory\n";
	      print ACE "\nVariation : \"$allele_ID\"\n";
	      print ACE "Laboratory \"$allele2lab{$allele_name}\"\n";
	      print LOG "\tSUCCESS: $allele2lab{$allele_name} found for $allele_ID\n";
	      next;
	  }
	  else {
	      print LOG "\tCHECK: $allele_name is not a recognised allele designation ($allele_ID), cannot fix Laboratory\n";
	      next;
	  }
	}	
      }
    }
    

    # check that Lab designation (when present) is correct (for CGC named alleles)
    if (defined $allele->Laboratory) {
      my $lab = $allele->Laboratory;
      if ($allele =~ m/^([a-z]{1,2})\d+$/) {
	my $allele_prefix = $1;
	if ( $allele2lab{$allele_prefix} ne $lab ) {
	  print LOG "ERROR: $allele is connected to $lab lab which is incorrect according to $allele2lab{$allele_prefix} lab info\n";
	}
      }
    }

    # Check for Public_name tag missing
    if (!defined($allele->Public_name)) {
	print LOG "ERROR: $allele has no Public_name tag\n";
    }

    # Check for Status tag missing
    if (!defined($allele->Status)) {
      print LOG "ERROR: $allele has no Status tag\n";
    }

    # Check for Variation_type
    if (!defined($allele->Variation_type )) {
	print LOG "ERROR: $allele has no Variation_type tag\n";
    }

    # Check for SeqStatus tag missing
    if (!defined($allele->SeqStatus)) {
	print LOG "ERROR: $allele has no SeqStatus tag\n" if ($allele->Status eq "Live");
    }

    # Check for method tag missing
    if (!defined($allele->Method)) {
      print LOG "ERROR: $allele has no Method tag\n";
    }

    # check for Method tag which has no value
    if (!defined $allele->Method) {
      my $method = $allele->at('Method');
      if (defined($method)) {
	print LOG "ERROR: $allele has a Method tag but not associated value\n";
      }
    }

    # check that Allele_type is set for all alleles with flanking sequences
    if (defined($allele->Flanking_sequences) && !defined($allele->Type_of_mutation)) {
      print LOG "ERROR: $allele has flanking sequences but no 'Type_of_mutation' tag\n";  
    }
    

    # test allele has method that matches Allele_type, eg, Deletion tag -> Deletion_allele
    if (defined $allele->Type_of_mutation) {
      my $expected_method; 
      my $observed_method;


      if (defined $allele->Method) {
	$observed_method = $allele->Method;
      } else {
	$observed_method = '';
      }


      my @mut_type = $allele->Type_of_mutation;
      if ( scalar @mut_type == 1 ) {

	# It is acceptable to have Method 'NemaGENETAG_consortium_allele' when we have a Transposon and so might expect the Method to be 'Transposon_insertion'
	if ($mut_type[0] eq "Insertion" && defined $allele->Transposon_insertion && $observed_method eq 'NemaGENETAG_consortium_allele') {
	  $expected_method = 'NemaGENETAG_consortium_allele';

	} elsif (defined $allele->Production_method && $allele->Production_method eq 'CRISPR_Cas9') {
	  $expected_method = 'Engineered_allele';

	} elsif ($mut_type[0] eq "Deletion" ) {
	  $expected_method = "Deletion_allele";

	} elsif ($mut_type[0] eq "Insertion" && !defined $allele->Transposon_insertion) {
	  $expected_method = "Insertion_allele";

	} elsif (defined $allele->Transposon_insertion && $allele->Transposon_insertion && $allele->at('Sequence_details.Type_of_mutation.Insertion')) {
	  $expected_method = "Transposon_insertion";

	} elsif (defined $allele->Type_of_mutation && $allele->Type_of_mutation eq "Substitution" ) {
	  $expected_method = "Substitution_allele";
	}

      } elsif (scalar @mut_type > 1) {
	
	if ( grep(/Deletion/, @mut_type) and grep(/Insertion/, @mut_type) ) {
	  $expected_method = "Deletion_and_insertion_allele";
	} else {
	  print LOG "ERROR: $allele has multiple Type_of_mutation @mut_type with method $observed_method. Unsure how to check this!\n";
	}
	
      }

      # does $observed method tag agree with expected method tag (based on Type_of_mutation tag)?
      if ($expected_method ne $observed_method) {
	if ($ace) {
	  print LOG "ERROR(a): $allele has wrong method ($observed_method): change to $expected_method\n";
	  print ACE "\nAllele : \"$allele\"\n";
	  print ACE "-D Method\n";
	  print ACE "\nAllele : \"$allele\"\n";
	  print ACE "Method \"$expected_method\"\n";
	} else {
	  print LOG "ERROR: $allele has method $observed_method which might need to be $expected_method\n";
	}
      }

      # find alleles that have flanking_seqs but no SMAPPED sequence
      if ( $allele ->Flanking_sequences && ! defined $allele ->Mapping_target ) {
	print LOG "ERROR: Allele $allele has Flanking_sequences tag but has no Mapping_target tag\n";
      }
    }  
  }
}
                          ###########################################################
                          #         SUBROUTINES FOR -class feature option           #
                          ###########################################################

sub process_feature_class{

  print "\n\nChecking Feature class for errors:\n";
  print LOG "\n\nChecking Feature class for errors:\n";

  # find features that have Mapping_target but no flanking_seqs
  #   if ( $feature ->Sequence && ! defined $feature ->Flanking_sequences ) {
  #       print LOG "WARNING: Feature $features has a Mapping_target tag but not Flanking_sequences\n";
  #   }
  # }   

  # find features that have flanking_seqs but no Mapping_target
  foreach my $feature ($db->fetch(-query=>'Find Feature WHERE Flanking_sequences AND NOT Mapping_target')){
    print LOG "ERROR: $feature has Flanking_sequences tag but has no Mapping_target tag\n"; 
  }

  foreach my $feature ($db->fetch(-query=>'Find Feature WHERE Sequence AND NOT Flanking_sequences')){
    print LOG "WARNING: $feature has Mapping_target but no Flanking_sequences tag\n"; 
  }

}




                          ###########################################################
                          #         SUBROUTINES FOR -class strain option            #
                          ###########################################################


sub process_strain_class {

  # Check if sequence name of a strain genotype is now connected to a CGC_name
  print"\n\nChecking Strain class for errors:\n";
  print LOG "\n\nChecking Strain class for errors:\n";
  print LOG "---------------------------------\n";
  print LOG "Loads of alleles are still connected to lin-15\n";
  my $strain;
  my @strains = $db->fetch('Strain','*');
  my $straincount = scalar(@strains);
  my $laststrain = $strains[-1]->name;
  unless ($laststrain =~ "WBStrain000$straincount") {print LOG "The last WBStrain ID $laststrain is not equal the strain class count $straincount (assuming 3x0 padding).\n";}

  # Accepted mutagen values
  #HCHO | Formaldehyde
  #ICR | Umbrell for these compounds?  ICR 191, 170, 292, 372, 191-OH, and 170-OH.
  #DEO | 1,2:7,8-Diepoxyoctane
  #DMH | 1,6-dimethoxyhexane (DMH)
  #60Co Irradiation | 60Co Irradiation
  #Acetaldehyde | Acetaldehyde
  #CRE | Cre-mediated Mutagenesis
  #CRISPR_Cas9 | CRISPR_Cas9
  #Cs137 Irradiation | Cs137 Irradiation
  #DEB | DiEpoxyButane
  #Diethoxybutane | Diethoxybutane
  #DES | DiEthyl Sulphate
  #EMS | Ethyl MethaneSulfonate (EMS)
  #Gamma Irradiation | gamma irradiation
  #Heat shock | Heat shock
  #IEF | Isoelectric focusing (IEF)
  #MMS | MethylMethane Sulphonate (MMS)
  #Microinjection | Microinjection
  #Microparticle bombardment | Microparticle bombardment
  #mini-Mos | mini-Mos
  #Mos1 transposon | Mos1 transposon
  #MosDel | MosDel
  #MosSCI | MosSCI
  #ENU | N-Ethyl-N-nitrosoUrea (ENU)
  #Nitrosoguanidine | Nitrosoguanidine
  #Phage transduction | Phage transduction
  #32P Irradiation | Phosphorus-32
  #Spontaneous | Spontaneous
  #TALEN | Transcription activator-like effector nuclease
  #Tc1 | Tc1
  #Tc4 | Tc4
  #Tc5 | Tc5
  #TMP | trimethylpsoralen (TMP)
  #UV | Ultra Violet
  #X-ray | X-ray
  #UNKNOWN | UNKNOWN
  #EMS+ENU | Mix
  #EMS+Formaldehyde | Mix
  #EMS+Tc3 | Mix
  #EMS+Tc4 | Mix
  #IEF+EMS | Mix
  #Microinjection+Gamma Irradiation | Mix
  #Microparticle bombardment+MosSCI | Mix
  #Microparticle bombardment+X-ray | Mix
  #TMP+UV+Gamma Irradiation | Mix
  #UV+Formaldehyde | Mix
  #UV+TMP | Mix
  my @mutagens = ('HCHO','ICR','DEO','DMH','60Co Irradiation','Acetaldehyde','CRE','CRISPR_Cas9','Cs137 Irradiation','DEB','Diethoxybutane','DES','EMS','Gamma Irradiation','Heat shock','IEF','MMS','Microinjection','Microparticle bombardment','mini-Mos','Mos1 transposon','MosDel','MosSCI','ENU','Nitrosoguanidine','Phage transduction','32P Irradiation','Spontaneous','TALEN','Tc1','Tc4','Tc5','TMP','UV','X-ray','UNKNOWN','EMS+ENU','EMS+Formaldehyde','EMS+Tc3','EMS+Tc4','IEF+EMS','Microinjection+Gamma Irradiation','Microparticle bombardment+MosSCI','Microparticle bombardment+X-ray','TMP+UV+Gamma Irradiation','UV+Formaldehyde','UV+TMP');
  
  
  foreach $strain (@strains){
    #N2 should not have Variants associated with it.
    if ($strain eq "WBStrain00000001") {
	  if (defined $strain->Contains) {
	      print LOG "WARNING: The N2 Strain ($strain) has been connected to variants in error!!!!!!\n";
	  }
	  else {
	      print LOG "RESULT: The N2 Strain ($strain) is clear of Variation data....which is good!\n" if ($debug);
	  }
    }
    if (defined $strain->Mutagen){
	my $mutagen = $strain->Mutagen->name;
	unless ($mutagen ~~ @mutagens){
	    print LOG "WARNING: Mutagen $mutagen is not on the known mutagen list, please investigate $strain\n";
	}
    }
    unless ($strain =~ /WBStrain\d{8}/){
	print LOG "WARNING: $strain has not been accessioned or merged into the WBStrain record.\n";
    }
    if (!$strain->Location){
      print LOG "WARNING(a): Strain $strain has no location tag\n";
      if ($ace){
	  if ($strain->Public_name->name =~ /([A-Z]+)\d+/){
	      print ACE "\n\nStrain : \"$strain\"\n";
	      print ACE "Location \"$1\"\n";
	  }
      }
    }
    else {
      my $cgc=$strain->Location;
      if ($strain->Genotype){
	my $genotype = $strain->Genotype;
	my $extract = $genotype;
	   $extract =~ s/\(|\)|\/|\+|;|\?|\{|\}|,|\=|\.$/ /g;
	   $extract =~ s/ I | II | III | IV | V | X / /g;
	   $extract =~ s/ I | II | III | IV | V | X | f / /g; # further tidying up of chromosomes
	   $extract =~ s/^\s|\w{3}-\s| f | A //g;
	   $extract =~ s/\s{1,}/ /g;
	my @items=split(/ /,$extract);

	foreach (@items){
	  if( my $wbgeneID = $Gene_info{$_} ){
	    if( $Gene_info{"$wbgeneID"}{'CGC_name'}) {
	      if ($cgc eq "CGC"){	
		print LOG "WARNING: CGC Strain $strain has sequence_name $_ in Genotype, which can now become $Gene_info{$wbgeneID}{'CGC_name'}\n";
	      }
	      else {
		print LOG "WARNING: Non_CGC Strain $strain has sequence_name $_ in Genotype, which can now become $Gene_info{$wbgeneID}{'CGC_name'}\n";
	      }
	    }
	  }
	}
      }
    }
  }
 


  my ($locus, %allele_locus, %strain_genotype, $cds, %locus_cds, $main, $allele);

  my $geneaceWqueryDir=$ENV{CVS_DIR}.'/../wquery/geneace';

  my $get_genotype_in_strain="Table-maker -p \"$def_dir/strain_genotype.def\"\nquit\n";
  my $allele_to_locus="Table-maker -p \"$def_dir/allele_to_locus.def\"\nquit\n";

  open (FH1, "echo '$get_genotype_in_strain' | $tace $database | ") || warn $!;
  open (FH2, "echo '$allele_to_locus' | $tace $database | ") || warn $!;

  # looks like:
  # "AA1"	"daf-12(rh257) X."
  while(<FH1>){
    chomp;
     if ($_ =~ /\"(.+)\"\s+\"(.+)\"/){
      $strain_genotype{$1} = $2;
    }
  }
  close FH1;

  # looks like:
  # "ad418"	"WBGene00000086"
  while (<FH2>){ # only genes with Alleles
    chomp $_;
    if ($_ =~ /\"(.+)\"\s+\"(.+)\"/){
      push(@{$allele_locus{$1}}, $Gene_info{$2}{'CGC_name'}) if $Gene_info{$2}{'CGC_name'};
    }
  }
  close FH2;

  # test with: "daf-12(rh257) X."
  foreach my $strain (keys %strain_genotype){
    my @matches = ($strain_genotype{$strain} =~  /((C[brjn]-\w{3,4}-\d+\w{0,1}|\w{3,4}-\d+\w{0,1})\(\w+\d+\))/g); # daf-12(rh257)
    foreach (@matches){
      my @la = split(/\s+/, $_);
      foreach (@la){
	if ($_ =~ /(C[brjn]-\w{3,4}-\d+\w{0,1}|\w{3,4}-\d+\w{0,1})\((\w+\d+)\)/){
	  $locus = $1; $allele = $2;

	  # diff allele->locus link in geneace and allele->locus in strain genotype
	  # if diff, print error LOG
	
	  if ( exists $allele_locus{$allele} ) {
	    my @LOci = @{$allele_locus{$allele}};
	    my %LOci;
	    map {$LOci{$_}++} @LOci;


	    unless ($LOci{$locus} ){
		my @obj;
		@obj = $db->fetch(-query => "FIND Strain $strain");
		my $public_name = $obj[0]->Public_name->name;
		print LOG "ERROR: Strain $strain\($public_name\) has $locus($allele) in genotype: ";
		print LOG "change each $locus to @{$allele_locus{$allele}}\n";
	    }
	  }
	}
      }
    }
  }
}

                          ##################################################################
                          #         SUBROUTINES FOR -class rearrangement option            #
                          ##################################################################


sub process_rearrangement {

  print"\n\nChecking Rearrangement class for errors:\n";
  print LOG "\n\nChecking Rearrangement class for errors:\n";
  print LOG "----------------------------------------\n";
  # checks presence of non-rearrangement object 
  # as objects of Rearrangement class

  my @rearr;
  my $count = 0;
  @rearr = $db -> fetch('Rearrangement','*');
  foreach (@rearr){
    if ($_ !~/\w+(Df|Dp|Ex|T|In|C|D)\d*/){
      $count++;
      print LOG "WARNING: $_ is NOT an object of Rearrangement\n";
    }
  }
  print LOG "No errors found\n" if $count == 0;
}


                          ############################################################
                          #         SUBROUTINES FOR -class paper option              #
                          ############################################################

sub process_paper_class {


    print "\nChecking Paper class for errors:\n";
    print LOG "\nChecking Paper class for errors:\n--------------------------------\n";

    #fetches the data from Geneace
    my @Papers =$db->fetch (-query => 'FIND Paper');
    #iterates through @Papers
    foreach my $paper_id ( @Papers ) {
	unless ($paper_id->name =~ (/WBPaper[0-9]{8}$/)) {
	    print LOG "Error: $paper_id is not formatted correctly\n"; # prints an error if the Paper ID does not match WBPaper12345678
	    print "Error: $paper_id is not formatted correctly\n" if ($debug); 
	} 
    }
}
  
                          ############################################################
                          #         SUBROUTINES FOR -class mapping option            #
                          ############################################################


sub check_genetics_coords_mapping {

  print "\nChecking discrepancies in genetics/coords mapping:\n";
  print LOG "\nChecking discrepancies in genetics/coords mapping:\n";
  print LOG "--------------------------------------------------\n";
  print JAHLOG "\nChecking discrepancies in genetics/coords mapping:\n";
  print JAHLOG "--------------------------------------------------\n";
  system ("$ENV{'CVS_DIR'}/GENEACE/get_interpolated_gmap.pl -database $database -diff");

  my $map_diff = "/nfs/wormpub/logs/mapping_diff.".$rundate;
  open(IN, $map_diff) || die $!;
  while(<IN>){
    print LOG $_;
    print JAHLOG $_;
  }
}




                          ############################################################
                          #      SUBROUTINES FOR -class multipoint option            #
                          ############################################################


sub check_dubious_multipt_gene_connections {

  my $header = "\n\nChecking dubious multi-pt <-> Gene association\n";
  my $sep = "-----------------------------------------------------------------------------------\n";
  print $header;
  print LOG $header, $sep;

  system("perl5.8.0 $ENV{'CVS_DIR'}/GENEACE/check_dubious_multi_pt_2_locus.pl");
  print LOG `cat $log_dir/dubious_multiPt_2_locus.$rundate`;
  print JAHLOG `cat $log_dir/dubious_multiPt_2_locus.$rundate`;
  print ACE `cat $database/CHECKS/multiPt_2_locus.ace` if $ace;
  #`rm -f /nfs/wormpub/DATABASES/geneace/CHECKS/multiPt_2_locus.ace`;
}

sub usage {
  my $error = shift;
  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

sub create_log_files{

  $log = "$log_dir/geneace_check.$rundate.$$";

  open (LOG, ">$log") or warn "cant open $log";
  print LOG "geneace_check\n";
  print LOG "Started at ",`date`,"\n";
  print LOG "=============================================\n";
  print LOG "\n";


  if($ace){
    print LOG "The (a) following ERROR, or UPDT, eg, indicates ace output \n$acefile for direct upload to correct problems.\n\n";
  }
  
}



sub check_operons {
  my @dead_genes = $db->fetch(-query=>'Find Gene; dead; Contained_in_operon');
  foreach my $gene (@dead_genes){
	print LOG $gene->name." dead but in operon\n";
  }
}

__END__

=head2 NAME - geneace_check.pl

=head3 <USAGE>


=head2 Options: [help] [debug] [class] [ace]


B<-help:>
            Displays usage of the script: surprised?

B<-debug:>
            Debug mode, follow by user name, eg. -d ck1 or -debug ck1

            This allows checking results mailed only to the user,
            otherwise results would email to the gangs of Sanger
            Wormbase group and Caltech (Erich)

B<-class:>
            Allows checking only specified classes in Geneace.
            All classes will be checked without this option.
            Choosing multiple classes is supported.
            Class names are case insensitive.
            Current valid options to check these classes are:

               gene
               allele
               laboratory
               strain
               rearrangement
               evidence
               mapping
               gmap
               multipoint

            For example: -class allele -class gene OR -class sequence -class rearrangment

B<-database:>
            Allows specifying path to a specific database.
            Default database path is /nfs/wormpub/DATABASES/geneace without this option.

            For example: -database /wormsrv2/autoace or -database /wormsrv1/someone/Test_DB


B<-ace:>
            Allows generating ace file for fixing erros spotted by
            this checking script.
            Default location and filename of ace file:
            /nfs/wormpub/DATABASES/geneace/CHECKS/geneace_check.rundate.processid.ace
            For example: -ace


B<-verbose:>
            Toggles on extra output.   Useful when running on command line and not on crob job
            For the locus class it will display each locus name as it is processes it and show
            (next to the locus name) a dot for each error it had found in that locus


=head3 <RUN geneace_check.pl>

            Geneace check is now set to run on Sundays.
            ##Run Geneace check over weekend
            20 7 * * 0 $ENV{'CVS_DIR'}/GENEACE/geneace_check.pl -a &
