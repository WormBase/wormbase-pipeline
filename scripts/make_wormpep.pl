#!/software/bin/perl -w
#
# make_wormpep
# 
# Usage : make_wormpep.pl 
#
# Builds a wormpep data set from the current autoace database
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2015-05-06 12:03:48 $

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use IO::Handle;
use Ace;
use File::Copy;
use Bio::SeqIO;

##############################
# command-line options       #
##############################

my ($help, $debug, $test, $verbose, $store, $wormbase);

my $initial;             # run script in initial mode at start of build
my $final;               # for full run at end of build
my $species;# = 'elegans';
my ($table,$fasta,$accession,$history,$pepfile,$dna,$pepace,$additional,$pid, $all);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "test"       => \$test,
            "verbose"    => \$verbose,
            "store:s"    => \$store,
            "species:s"  => \$species,
            "initial"    => \$initial,
            "final"      => \$final,
            "table"      => \$table,
            "fasta"	 => \$fasta,
            "accession"  => \$accession,
            "history"    => \$history,
            "pepfile"    => \$pepfile,
            "dna"	 => \$dna,
            "pepace"	 => \$pepace,
            "additional" => \$additional,
            "pid"	 => \$pid,
            "all"	 => \$all,	    	
           );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism=> $species
			     );
}

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);
$log->write_to("Building protein set for ".$wormbase->full_name('-short' => 1)."\n");


##########################################
# Set up database paths                  #
##########################################

# Set up top level base directory which is different if in test mode
# Make all other directories relative to this
my $basedir   = $wormbase->basedir;
my $dbdir     = $wormbase->orgdb;
my $new_wpdir = $wormbase->wormpep;
my $release   = $wormbase->get_wormbase_version;
my $tace      = $wormbase->tace;

my $PEP_FILE_STEM = $wormbase->pepdir_prefix . "pep";
my $PEP_PREFIX    = $wormbase->pep_prefix;

# need to get previous build WORMPEP
my ($pepdir, $old_release ) = $new_wpdir =~ /^(\S+)\/${PEP_FILE_STEM}(\d+)/;
my ($wpdir, $old_wpdir_test_file);
do {
  $old_release--;
  $wpdir = "$pepdir/${PEP_FILE_STEM}${old_release}";
  $old_wpdir_test_file = "$wpdir/${PEP_FILE_STEM}.fasta${old_release}";
} while (not -e $old_wpdir_test_file);


# Make new directory for current release         
if (-e $new_wpdir){
  $log->write_to("$new_wpdir already exists\n");
}
else {
  mkdir ("$new_wpdir" , 0755) || $log->log_and_die("=> Failed to create a new directory for wormpep release wormpep$release\n\n");
  $log->write_to($wormbase->runtime.": making ${PEP_FILE_STEM}$release\n\n");
}

#GLOBAL VARS
my $TOP_ID =0;
my %cds2aa;
my %cds2id;
my %cds2gene_id;
my %id2cds;
my %aa_seq2id;

if ($initial) {
  $fasta=$accession=$history=$pepfile=$dna=$pepace=$additional=1 if $all;
  &write_dna 			if $dna;	#writes dna file and pops %cds2aa [DataDumps]
  &read_old_build;					#reads old fasta & pops %aa_seq2id. 
  &assign_new_ids;					#iterates over %cds2aa assigns id. pops %id2cds and fills in %aa_seq2id for new
  &write_wormpep_history_and_diff if $history;	#reads old history &updates based on %id2cds
  &write_fasta 		if $fasta;
  &run_pepace			if $pepace; #req history and fasta files
  &write_accession 	if $accession;
  &write_blast_pep 	if $pepfile;
}
elsif( $final ) {
  $pid=$pepfile=$table=1 if $all;
  &write_final_pep 	if $pepfile;
  &write_table 		if $table; #normally run within final pepfile
}

$log->mail;

sub write_dna {  
  my $pep_dna_file = $new_wpdir."/".$wormbase->pepdir_prefix."pep.dna".${release};
  my $pep_aa_file  = $new_wpdir."/".$wormbase->pepdir_prefix."pep${release}.pep";
  my $query = "query find CDS where method = \"curated\" AND species = \"".$wormbase->full_name."\"";
  
  $log->write_to("dumping dna file\n");
  my $command =<<EOF;
$query 
dna -f "$pep_dna_file"
quit
EOF

  open (TACE, "echo '$command' | tace $dbdir | ") or die "Failed to open $!\n";
  while (<TACE>){ print }
  close TACE;
  
  $log->write_to("dumping peptide file\n");
  $command =<<EOF;
$query 
Peptide -f "$pep_aa_file"
quit
EOF

  open (TACE, "echo '$command' | tace $dbdir | ") or die "Failed to open $!\n";
  while (<TACE>){ 
    print;
  }
  close TACE;
  
  my $cdspeps = Bio::SeqIO->new(-file => $pep_aa_file, '-format' => 'Fasta');
  while(my $cdspep = $cdspeps->next_seq) {
    $cds2aa{$cdspep->id} = $cdspep->seq;
  }
  
  open (G2P,">".$wormbase->common_data."/cds2aa.dat");
  print G2P Data::Dumper->Dump([\%cds2aa]);	
  close G2P;
}

# read in the old wormpep.fasta file and use this to create new wormpep entries for any
# cds's with novel amino acid sequences.
sub read_old_build {
  # read in the wp.fasta file, contains all protein sequences ever assigned
  $log->write_to($wormbase->runtime." : initiate parsing of fasta file\n");	
  my $fasta = $wpdir."/$PEP_FILE_STEM.fasta".$old_release;
  my $old_fasta = Bio::SeqIO->new(-file => $fasta, '-format' => 'Fasta');
  while(my $seq = $old_fasta->next_seq) {
    my $id = $seq->id;
    $id =~ s/$PEP_PREFIX//;
    $aa_seq2id{$seq->seq} = $id unless $aa_seq2id{$seq->seq};
    my ($id_num) = $seq->id =~/(\d+)/;
    $TOP_ID = $id_num if($TOP_ID < $id_num);
  }
}

sub assign_new_ids {

  unless (%cds2aa) {
    %cds2aa = $wormbase->FetchData('cds2aa'); #cr01.sctg0.wum.3.1 => MKIFLAILTLVA....
  }
  
  foreach my $cds (keys %cds2aa) {	
    if($aa_seq2id{$cds2aa{$cds}}) {       #if peptide already exists
      $cds2id{$cds}= $aa_seq2id{$cds2aa{$cds}}; #fill cds2id hash
      push(@{$id2cds{$aa_seq2id{$cds2aa{$cds}}}},$cds);
    }
    else {
      # this peptide sequence doesn't exist
      my $new_id = ++$TOP_ID;
      $aa_seq2id{$cds2aa{$cds}} = $new_id;
      $cds2id{$cds}= $new_id; # store only the number, not prefix (eg CE)#fill cds2id hash
      push(@{$id2cds{$new_id}},$cds);
    }
  }
  #dump hash
  open (C2P,">".$wormbase->common_data."/cds2pepid.dat");
  print C2P Data::Dumper->Dump([\%cds2id]);
  close C2P;	
  
  #dump hash
  open (AA2ID,">".$wormbase->common_data."/aa2pepid.dat");
  print AA2ID Data::Dumper->Dump([\%aa_seq2id]);
  close AA2ID;
}

sub write_wormpep_history_and_diff {
  open (OLDHISTORY, "<$wpdir/${PEP_FILE_STEM}.history$old_release")  
      or $log->log_and_die("cannot open $wpdir/${PEP_FILE_STEM}.history$old_release: $!\n");
  open (HISTORY,    ">$new_wpdir/${PEP_FILE_STEM}.history$release") 
      or $log->log_and_die("cannot create $wpdir/${PEP_FILE_STEM}.history$release: $!\n");
  open (DIFF,       ">$new_wpdir/${PEP_FILE_STEM}.diff$release")    
      or $log->log_and_die("cannot create $wpdir/${PEP_FILE_STEM}.diff$release: $!\n");
  open (NEW,        ">$new_wpdir/new_entries.WS$release")
      or $log->log_and_die("cannot create new entries file :$!\n");
  unless (%cds2id){
    %cds2id = $wormbase->FetchData('cds2pepid');
  }
  unless (%cds2aa){
    %cds2aa = $wormbase->FetchData('cds2aa');
  }
  my %line;
  my %live_last_build;
  while (my $line = <OLDHISTORY>) {
    chomp $line;
    my ($cds , $wpid , $start , $end) = split (/\t+/ , $line);
    $wpid =~ /${PEP_PREFIX}(\d+)/ ; my $num = $1;
    $line{$cds} = $line; # track which CDSs have an entry
    
    if ($end){ #stuff dead in last release
      if($cds2id{$cds}){ #now coding
        if($cds2id{$cds} eq $num)  { #reappeared
          print DIFF "reappeared:\t$cds\t${PEP_PREFIX}".&pad($cds2id{$cds})."\n";
        }
        else { #no change - coding for other peptide (or the same one reinstated)
          #B0303.7 CE00538 8       17
          #B0303.7 CE00538 18      77
          #B0303.7 CE30412 77      123
          print HISTORY "$cds\t$wpid\t$start\t$end\n";
        }
      }
      else { #no change
        print HISTORY "$cds\t$wpid\t$start\t$end\n";
      }
    }
    else { #live in last release
      $live_last_build{$cds} =1;
      if($cds2id{$cds}){ #still coding
        if("${PEP_PREFIX}$cds2id{$cds}" eq $wpid) {#same peptide
          print HISTORY "$cds\t$wpid\t$start\n"; #no change
        }
        else{ #cds coding for something else
          print HISTORY "$cds\t$wpid\t$start\t$release\n";
          print HISTORY "$cds\t${PEP_PREFIX}".&pad($cds2id{$cds})."\t$release\n";
          print DIFF "changed:\t$cds\t$wpid --> ${PEP_PREFIX}".&pad($cds2id{$cds})."\n";
          my $pad = &pad($cds2id{$cds});
          print NEW ">$PEP_PREFIX".$pad."\n".$wormbase->format_sequence($cds2aa{$cds})."\n";
        }
      }
      else { #now dead
        print HISTORY "$cds\t$wpid\t$start\t$release\n";
        print DIFF "lost:\t$cds\t$wpid\n";
      }
    }
  }
  close OLDHISTORY;
  #catch resurrections
  foreach my $cds (keys %cds2id){
    unless ($live_last_build{$cds}) {
      my $pad = &pad($cds2id{$cds});
      print HISTORY "$cds\t${PEP_PREFIX}$pad\t$release\n";
      print DIFF "new:\t$cds\t${PEP_PREFIX}$pad\n";
      print NEW ">$PEP_PREFIX".$pad."\n".$wormbase->format_sequence($cds2aa{$cds})."\n";
    }
  }
  
  close HISTORY;
  close DIFF;
  close NEW;
}

sub pad {
  my $num = shift;
  return sprintf "%05d" , $num;
}

sub write_fasta {
  #>CE00001
  #MLRYIFFAAILFFVFPNTESAGFLKFELTADRDCLLHLEHSSTYSETVRLLAYESRPLEI
  #YTQGSINEIPVHFQLLHHFSGKALSEAKFQIFQLKNNGLWDSKVIDTDKVILSVRSTFYC
  #ENGYFGPICDRRSRTFAPKSDIQTSTPGYQTQVLKFDFKISDDIIIYSSLAFFVLLLIIF

  unless (%aa_seq2id) {
    %aa_seq2id = $wormbase->FetchData('aa2pepid');
  }
  open (FASTA,">$new_wpdir/${PEP_FILE_STEM}.fasta$release")
      or $log->log_and_die("cannot create $wpdir/${PEP_FILE_STEM}.fa$release\n");
  foreach my $pep (sort { $aa_seq2id{$a} cmp $aa_seq2id{$b} } keys %aa_seq2id) {
    print FASTA ">$PEP_PREFIX".&pad($aa_seq2id{$pep})."\n".$wormbase->format_sequence($pep)."\n";
  }
  close FASTA;
}

sub write_blast_pep {
  #>2L52.1 CE32090 
  #MSMVRNVSNQSEKLEILSCKWVGCLKSTEVFKTVEKLLDHVTADHIPEVIVNDDGSEEVV
  #CQWDCCEMGASRGNLQKKKEWMENHFKTRHVRKAKIFKCLIEDCPVVKSSSQEIETHLRI

  #load data from previous subs if not populated.
  unless (%cds2id){
    %cds2id = $wormbase->FetchData('cds2pepid');
  }
  unless (%cds2aa){
    %cds2aa = $wormbase->FetchData('cds2aa');
  }
	my $pep_file = $wormbase->wormpep."/$PEP_FILE_STEM".$release;
  open (PEP,">$pep_file") or $log->log_and_die("cant write $pep_file : $!\n");
  my ($cds,$id);
  while(($cds,$id) = each (%cds2id)){
    print PEP ">$cds $PEP_PREFIX".&pad($id)."\n";
    print PEP $wormbase->format_sequence($cds2aa{$cds})."\n";
  }
  close PEP;
}

sub write_final_pep {
  #>2L52.1 CE32090 WBGene00007063 Zinc finger, C2H2 type status:Partially_confirmed TR:A4F336 protein_id:ABO33278.1
  #MSMVRNVSNQSEKLEILSCKWVGCLKSTEVFKTVEKLLDHVTADHIPEVIVNDDGSEEVV
  #CQWDCCEMGASRGNLQKKKEWMENHFKTRHVRKAKIFKCLIEDCPVVKSSSQEIETHLRI
  $log->write_to("writing final peptide file\n");
  
  unless (%cds2id){
    %cds2id = $wormbase->FetchData('cds2pepid');
  }
  #load data from previous subs if not populated.
  unless (%cds2aa){
    %cds2aa = $wormbase->FetchData('cds2aa');
  }
  
  my %cds_info;
  my $def = $wormbase->autoace."/wquery/SCRIPT:make_wormpep.def";
  print STDERR "parsing $def\n";
  my $tm_query = $wormbase->table_maker_query($dbdir,$def);
  while(<$tm_query>) {
    next if (/>/ or /\/\//);
    s/\"//g;#"
    s/\n//g;
    my($cds,$pid,$pidver,$unip,$unip_ac,$gene,$cgc,$status,$brief) = split(/\t/,$_);
    next unless ($cds and $gene);
    $cds_info{$cds}->{'gene'} = $gene;
    $cds_info{$cds}->{'pid'} = "$pid.$pidver" if ($pid and $pidver);
    $cds_info{$cds}->{'unip'} = $unip_ac if $unip_ac;
    $cds_info{$cds}->{'cgc'} = $cgc if $cgc;
    $cds_info{$cds}->{'status'} = $status if $status;
    $cds_info{$cds}->{'brief'} = $brief if $brief;
  }
  
  my $pep_file = $wormbase->wormpep."/$PEP_FILE_STEM".$release;
  open (PEP,">$pep_file") or $log->log_and_die("cant write $pep_file : $!\n");
  foreach my $cds (sort keys %cds_info){	
    unless (defined $cds2aa{$cds}) {
      $log->write_to("no sequence for $cds\n");
      next;
    }
    my @desc = (sprintf("wormpep=%s%s", $PEP_PREFIX, &pad($cds2id{$cds})));
    push @desc, sprintf("gene=%s",  $cds_info{$cds}->{'gene'})        if $cds_info{$cds}->{'gene'};
    push @desc, sprintf("locus=%s", $cds_info{$cds}->{'cgc'})         if $cds_info{$cds}->{'cgc'};
    push @desc, sprintf("status=%s", $cds_info{$cds}->{'status'})     if $cds_info{$cds}->{'status'};
    push @desc, sprintf("uniprot=%s", $cds_info{$cds}->{'unip'})      if $cds_info{$cds}->{'unip'};
    push @desc, sprintf("insdc=%s", $cds_info{$cds}->{'pid'})         if $cds_info{$cds}->{'pid'} ;
    push @desc, sprintf("product=\"%s\"", $cds_info{$cds}->{'brief'}) if $cds_info{$cds}->{'brief'};

    print "$cds\n";
    printf PEP ">%s %s\n", $cds, join(" ", @desc);
    print PEP $wormbase->format_sequence( $cds2aa{$cds} )."\n";
  }
  close(PEP) or $log->log_and_die("Could not successfully write wormpep file; probably disk full\n");
  
  &write_table;
  &count_isoforms if ($wormbase->species eq 'elegans');
}

sub write_table {
  #>2L52.1 CE32090         Zinc finger, C2H2 type  Partially_confirmed     TR:A4F336       ABO33278.1
  my $pep_file = $wormbase->wormpep."/$PEP_FILE_STEM".$release;
  my $table_file = $wormbase->wormpep."/$PEP_FILE_STEM.table".$release;
  if( -e $pep_file) {
    $wormbase->run_command("grep '>' $pep_file > $table_file", $log);
  }
  else {
    $log->log_and_die("$pep_file does not exist\n");
  }
}

sub write_accession {
  #CE32090 2L52.1 
  $log->write_to("writing accession file\n");
  my $acc_file = $wormbase->wormpep."/$PEP_FILE_STEM.accession".$release;
  open (PEP,">$acc_file") or $log->log_and_die("cant write $acc_file : $!\n");
  foreach my $id (sort keys %id2cds) {
    my $cds = $id2cds{$id};
    print PEP "$PEP_PREFIX$id\t";
    print PEP join("\t",@$cds) if $cds;
    print PEP "\n";
  }
  close PEP;		
}


sub run_pepace {
  # build database history info 
  my $history_file = "$new_wpdir/${PEP_FILE_STEM}.history$release";
  if( -e $history_file ){
    $wormbase->run_script("build_pepace.pl", $log);
  }
  else {
    $log->error("cant build pepace as $history file no good : $!\n");
  }
}

sub get_embl_data {
  my $pid_mail = $wormbase->autoace."/protein_ID.mail";
  $log->log_and_die("no mail $pid_mail: !$\n") unless (-e $pid_mail);
  $wormbase->run_script("get_EMBL_data.pl", $log);
}


sub count_isoforms{
  my @wormpep_CDSs;
  open (FILE , "$new_wpdir/wormpep$release") or $log->log_and_die("cannot open the wormpep$release file\n");
  while (<FILE>) {
    if (/^\>(\S+)\s+\S+/) {
      my $cds = $1;
      push (@wormpep_CDSs , $cds);
    }
  }
  close (FILE);
  
  my $total_cds_count  = 0;

  ##########################################################################################################
  # WARNING: there are 3 pairs of genes in C. elegans that share a Sequence ID (B0564.1 Y105E8A.7 ZC416.8) #
  # leading to an off by 3 calculation in the $no_isoform_count so could do the following:                 #
  ##########################################################################################################
  # my $no_isoform_count;
  # if ($species eq 'elegans') {
  #   $no_isoform_count = 3;
  # }
  # else {
  #   $no_isoform_count = 0;
  # }

  my $no_isoform_count = 0;
  my $isoform_count    = 0;
  my %new_name2x;
  
  foreach (@wormpep_CDSs) {
    $total_cds_count++;
    /^([A-Z0-9_]+)\.(.*)$/i;
    my $cds_prefix = $1;
    my $cds_suffix = $2;
    if ($cds_suffix =~ /^[1-9][0-9]{0,4}$/) {
      $no_isoform_count++;
    }
    elsif ($cds_suffix =~ /(^[1-9][0-9]{0,4})([a-z])/) {
      my $number = $1;
      my $letter = $2;
      my $new_name = $cds_prefix."_".$number;
      unless (exists ($new_name2x{$new_name})) {
        $no_isoform_count++;
      } 
      else {
        $isoform_count++;
      }
      $new_name2x{$new_name} = "x";
    } 
    else {
      $log->write_to("$_ has a non\-acceptable name in wormpep$release \(has not been counted\)\n");
      next;
    }
  }
  $log->write_to("\n\nThis release of Wormpep derives from $total_cds_count CDSs (from $no_isoform_count protein-coding genes)\n");
  
  # write the release letter (also does some checks)
  $wormbase->release_wormpep($no_isoform_count,$total_cds_count,$isoform_count, $log);
  chmod (0444 , "$new_wpdir/*") or $log->write_to("cannot chmod $new_wpdir/ files\n");
}
