#!/usr/local/bin/perl5.8.0 -w

# incorporate_EMBL_names.pl


# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-05-27 13:16:35 $ 

use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use strict;
use Getopt::Long;
use GENEACE::Geneace;

my $version;
GetOptions ("version|v=s"  => \$version);

# warning
if (!$version){
  print "You need to specify version like -v R75 to proceed!\n";
  exit(0);
}

#################

# dowload EMBL AC

#################

my $start = &runtime;

# directory, files used for EMBL update
#my $script_dir = "/wormsrv2/scripts/GENEACE/";
my $script_dir = "/nfs/team71/worm/ck1/WORMBASE_CVS/scripts/GENEACE/";

my $output_dir = "/wormsrv1/geneace/EMBL_GENE_NAMES/";
my $temp_dir   = "/wormsrv1/geneace/EMBL_GENE_NAMES/TEMP/";

my $download = "$temp_dir/embl";
my $flatfile= "$temp_dir/$version"."_EMBL_entry_no_seqs";

print "\nDownloading EMBL Release $version . . .\n";

# query parameters
`getz -e "((([emblrelease-Organism:caenorhabditis*]&[emblrelease-Division:inv]))>[emblrelease-FtQualifier:gene]) > parent" > $download`;

print "\nDeleting protein and DNA sequences from EMBL flatfile . . . ";

open(IN, $download) || die "Can't read in file!";
open(EMBL, ">$flatfile") || die "Can't write to file!";

## fetch only ID, AC, SV, KW, De, OS, FT RT and no DNA/protein seq
while(<IN>){
  chomp;
  if ($_ =~ /(^ID.+)|(^AC.+)|(^SV.+)|(^KW.+)|(^DE.+)|(^OS.+)|(^RT.+)|(^\/\/)/){
    print EMBL $_,"\n";
  } 
  if (($_ =~ /^FT[\w\W.]+/)&&($_ !~ /^FT\s+\/translation[\w\W.]+/) &&($_ !~ /^FT\s+[A-Z]+/)){
    print EMBL $_,"\n";
  }
}
print "Done\n";

close EMBL;
close IN;

system("rm -f $download");

############################

# Processing EMBL entry data

############################

##################
# global variables
##################

my (%AC_SV_ID_DE, %AC_ID_Protid_DE, %Non_cgc_style_name_Protid_AC_OS, %Non_cgc_style_name_Protid_AC_OS_SV,
    %ProtID_AC_SV_Pver_OS_Gene, %ProtID_AC_SV_Pver_OS_Gene_Chksum);

my ($AC, $ID, $DE, $DEs, $seq_version, $OS, $cgc_style_name, $non_cgc_style_name, $product, $proteinid);
my $AC_line = 0; my $omit = 0;

open(IN, $flatfile) || die "Can't open the file!";

while(<IN>){

  chomp;
  my $each_line = $_;

  #############################################################
  # get AC:  get only 1st accession of each AC line if multiple
  #############################################################
  if ($each_line =~ /^AC\s+(\w+);$/ || $each_line =~ /^AC\s+(\w+);.+$/){
    $AC_line++;
    $AC = $1 if $AC_line == 1; # some AC may have > 1 AC lines, read in only the first AC line
    next;
  }

  ########
  # get ID
  ########
  if ($each_line =~ /^ID\s+(\w+)\s+.+/ ){
    $ID = $1;
    next;
  }

  ########
  # get DE
  ######## 
  if ($each_line =~ /^DE\s{3}(.+)/){
    $DE=$1;
    $DEs .= $DE;  # concat DE that has two lines
    next;
  }

  #################
  # get seq version
  #################
  if ($each_line =~ /^SV.+\.(\d)+/){
    $seq_version = $1;
    next;	
  }

  ###################
  # get organism name
  ###################
  if ($each_line =~ /^OS\s{3}(.+)/){
    $OS = $1 if $1 eq "Caenorhabditis elegans"; # do only C. elegans for now
    next;
  }

  #######################
  # omit AC from WormBase
  #######################

  if ($each_line =~ /RT\s+\"Genome sequence of the nematode C/i ){
    $AC_line = 0; $omit = 1;
    next;
  }


  ###########################################################
  # get gene name / sequence name / standard name (seq. name)
  ###########################################################
  if ($omit == 0){
    if ($each_line =~ /^FT\s+\/gene=\"(.+)\"/ ){
      $cgc_style_name =();
      my $name = $1;
  #    if ($name =~  /.+\..+/ )
      if ($name =~ /^[a-z]{3,3}-\d+$/ || $name =~ /^\w{3,3}-\d+\.[a-zA-Z0-9]+$/ ){
	$cgc_style_name = $name;
      }
      else {
	$non_cgc_style_name = $name;
      }
      next;
    }
  }

  #############
  # get product
  #############
  if ($each_line =~ /^FT\s+\/product=\"(^\w{3,3}-\d+$)\"/ && $omit == 0){
    $product =();
    # grep only 3-letter style product name
    $product = $1;
    next;
  }

  ################
  # get protein_id
  ################
  if ($each_line =~ /^FT\s+\/protein_id=\"(.+)\.(\d+)\"/ && $omit == 0){
    my $prot_id = $1;
    my $pid_ver = $2;
    $proteinid = ();
    $proteinid =  $prot_id.".".$pid_ver;

    if ( $cgc_style_name ){
      # associate AC, SV, OS, Gene name with a prot_id and avoid duplications
      if (!exists $ProtID_AC_SV_Pver_OS_Gene{$prot_id}){
         push(@{$ProtID_AC_SV_Pver_OS_Gene{$prot_id}}, $AC, $seq_version, $pid_ver, $OS, $cgc_style_name);
      }
    }
    # at absence of gene name, 3 letter-number style product becomes gene name, if available, otherwise ignore
    if ( $product && !$cgc_style_name ){
      # associate AC, SV, OS, product with a prot_id
      push(@{$ProtID_AC_SV_Pver_OS_Gene{$prot_id}}, $AC, $seq_version, $pid_ver, $OS, $product);  # $product might be upper case

    }
    # omit grappig $non_cgc_style_name if $cgc_style_name already is avalable (no need for duplication)
    if ( $non_cgc_style_name && !$cgc_style_name ){
      push(@{$Non_cgc_style_name_Protid_AC_OS_SV{$non_cgc_style_name}}, $prot_id.".".$pid_ver, $AC, $OS, $seq_version);
    }
    $cgc_style_name =();  $product =();  $non_cgc_style_name =();
    next;	
  }

  ############################################
  # creating further hashes at end of an entry
  ############################################
  if ($each_line =~ /^\/\//){

    $proteinid = "NA" if !$proteinid;
    push(@{$AC_ID_Protid_DE{$AC}}, $ID, $proteinid, $DEs);

    # reinitialization at end of each accession or screwed up
    $cgc_style_name = ();  $product = ();  $non_cgc_style_name =();  $AC_line = 0; $DEs =(); $DE=(); $omit = 0;
  }
}
close IN;


#####################################################################################
#   assign protein_id chksum from SWALL to names from EMBL /gene tag with protein_id
#####################################################################################

# output swall entries to a file, parse it and prepare a hash: key(prot_id), value (checksum)
# for now fetch only C. elegans entries

print "\nDownloading SWALL flatfile via getz . . .";
`getz -e "[swall-org:Caenorhabditis elegans]" > $temp_dir/SWALL.$version`;
print "Done\n";

open(SWALL, "$temp_dir/SWALL.$version") || die $!;
my (%pid_ver, $protid, $pid_ver, %ProtID_Chksum_Ver, $checksum);
while (<SWALL>){
  chomp;
  if ($_ =~ /^DR\s{3}EMBL;\s.+;\s(.+)\.(\d+);.+/){ # the line: DR   EMBL; L26546; AAA20077\.1; -\.

    $pid_ver{$1} = $2;
  }
  if ($_ =~ /^SQ.+MW;\s+(.+)\s+.+/){$checksum = $1}
  if ($_ =~ /\/\//){
    foreach (keys %pid_ver){
      push(@{$ProtID_Chksum_Ver{$_}}, $checksum, $pid_ver{$_});
    }
    %pid_ver =();
  }
}

# remove temp swall file
#system("rm -f $temp_dir/SWALL.$version");

# key: prot_id, values: gene name, AC, SV, OS and Gene, checksum (appended)
# checksum is for C.elegans only, for the time being

foreach (keys %ProtID_AC_SV_Pver_OS_Gene){  # EMBL
  if (exists $ProtID_Chksum_Ver{$_}){	            # SWALL
    push(@{$ProtID_AC_SV_Pver_OS_Gene{$_}}, $ProtID_Chksum_Ver{$_}->[0]);
  }
  else {
    push(@{$ProtID_AC_SV_Pver_OS_Gene{$_}}, "NA");
  }
}

my %Non_cgc_style_name_Protid_AC_OS_SV_Chksum;

foreach (keys %Non_cgc_style_name_Protid_AC_OS_SV){  # EMBL
  for (my $i = 0; $i < scalar @{$Non_cgc_style_name_Protid_AC_OS_SV{$_}}; $i=$i+4 ){
    my $pid = $Non_cgc_style_name_Protid_AC_OS_SV{$_}->[$i];
    $pid =~ /(.+)\.(\d+)/;
    $pid = $1;
    my $version = $2;
    my $ac = $Non_cgc_style_name_Protid_AC_OS_SV{$_}->[$i+1];	
    my $os = $Non_cgc_style_name_Protid_AC_OS_SV{$_}->[$i+2];
    my $sv = $Non_cgc_style_name_Protid_AC_OS_SV{$_}->[$i+3];

    if ( exists $ProtID_Chksum_Ver{$pid} ){	 # SWALL
      push(@{$Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}}, $pid, $ac, $os, $sv, $ProtID_Chksum_Ver{$pid}->[0]);
    }
    else {
      push(@{$Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}}, $pid, $ac, $os, $sv, "NA");
    }
  }
}

%ProtID_AC_SV_Pver_OS_Gene_Chksum = %ProtID_AC_SV_Pver_OS_Gene; # rename
%ProtID_AC_SV_Pver_OS_Gene = ();                                # undef
%Non_cgc_style_name_Protid_AC_OS_SV = ();                       # undef

open(GN, ">$temp_dir/EMBL_cgc_style_name_dataset.$version") || die $!;
foreach (sort keys %ProtID_AC_SV_Pver_OS_Gene_Chksum){
  print GN "$_ -> @{$ProtID_AC_SV_Pver_OS_Gene_Chksum{$_}}\n";
}
open(OUT, ">$temp_dir/EMBL_non_cgc_style_name_dataset.$version") || die $!;
foreach (keys %Non_cgc_style_name_Protid_AC_OS_SV_Chksum){
  print OUT "$_\t@{$Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}}\n";
}

##############################
#  prepare geneace dataset
##############################

# ----- associate Geneace cgc names with protein id, chksum, other_names, sequence names 
# ----- via wormpep.tableXXX and SWALL file

my $ga = init Geneace();
my $ga_dir   = $ga->geneace();

my $db = Ace->connect(-path =>$ga_dir) || die "Connection failure: ",Ace->error;

my (@wormpep_tbl, %GA_info);

my $wormpep_ver = get_wormbase_version() - 1;
@wormpep_tbl = `cut -f 1,3,7 /wormsrv2/WORMPEP/wormpep$wormpep_ver/wormpep.table$wormpep_ver`;
if (!@wormpep_tbl){
  print "wormpep.table$wormpep_ver not available . . . cannot proceed, script ends!!\n";
  exit(0)
}

my %wormpep_seqs; # keys are all the valid seq. in wormpep (without isoforms suffixes)

foreach (@wormpep_tbl){
  chomp;
  my($seq, $locus, $prot_id) = split(/\t/, $_);

  $seq =~ s/>//;
  if ( $seq =~ /(.+\.\d+)\w+/ ){$seq = $1};  # make all isoforms the format of Sequence_name, ie, B0034.1a becomes B0034.1
  $wormpep_seqs{$seq}++;

  $locus = "NA" if !$locus;
  $prot_id =~ s/\.\d+$//; # drop suffix

  if (exists $ProtID_Chksum_Ver{$prot_id}) {
    push(@{$GA_info{$ProtID_Chksum_Ver{$prot_id}->[0]}}, $locus, $seq); # key is checksum, $locus is the CGC name and sequence name (no isoforms) in wormpep.table
                                                                        # the values maybe several CGC name and seq. name 
  }
}

open(DEBUG, ">$temp_dir/Debug_info.$version") || die $!;
foreach (keys %GA_info){
  print DEBUG "\n(1) $_ -> @{$GA_info{$_}}\n";
}

# -- hash of CGC_name->Other_name
my %main_other = $ga->other_name($db,"main_other");
my @other_names = $ga->other_name($db,"other");

my %other;
foreach (@other_names){$other{$_}++};  # put all other_names as keys of %other


# -- hash for Gene id <-> CGC name conversion
# eg, $Gene_info{'WBG0000001'}{'CGC_name'} returns CGC name of the gene id , similar syntax for Sequence name / Other_name / Public_name
#     $Gene_info{'CGC_name'}{'Gene'} to get gene id, eg, WBG0000001 of a CGC_name

my %Gene_info = $ga -> gene_info();


# get SV of existing other_sequence in Geneace from ?Accession_number
my %AC_SV;
my $SV_query = "find Accession_number *.*";
my @AC_SV_in_GA = $db->find($SV_query);

foreach (@AC_SV_in_GA){
  $_  =~ /(.+)\.(\d+)/;   #eg. AF03445.1
  $AC_SV{$1} = $2;        # AF03445 is key, 1 is version
}

$db->close;

#######################################################################################
# compare checksum of EMBL /gene and Geneace info : look only for identical checksums
#######################################################################################

open(PARSE, ">$temp_dir/EMBL_to_Geneace_pre_longtext.$version") || die $!;
# cols of the file EMBL_to_Geneace_pre_longtext.$version:
#     AC       gene id     ID       SV prot_id    DE                                                      OS
#--------------------------------------------------------------------------------------------------------------------------------
# eg: AB112928 WBG0000001  AB112928 1  BAD07033.1 Caenorhabditis elegans rab3 mRNA for Rab3, complete cds.Caenorhabditis elegans

open(NOMATCH, ">$temp_dir/EMBL_names_with_diff_chksum_to_wormpep.$version") || die $!;
open(UPDT, ">$output_dir/EMBL_update.$version.ace") || die $!; # load this before longtext
open(INFO, ">$output_dir/EMBL_name_incorporation_stats.$version") || die $!;

my ($pepace_chksum, $ga_name, $ga_cds);

# update counter info
my $new_other_name =0;
my $updated_AC = 0;
my $new_AC = 0;
my @new_ACs = (); my @updated_ACs = ();
my @new_Other_names = ();
my @new_ACs1 = ();
my @new_Other_names1 = ();
my %new_AC = ();
my %updated_AC = ();
my %new_other_name = ();
my $no_match = 0;


#--- deal with standard names from EMBL /gene tag

foreach (keys %ProtID_AC_SV_Pver_OS_Gene_Chksum){ # EMBL dataset: key is prefix of a protein_id only

  my $prot_prefix = $_;
  my $embl_ac     = $ProtID_AC_SV_Pver_OS_Gene_Chksum{$_}->[0];
  my $embl_sv     = $ProtID_AC_SV_Pver_OS_Gene_Chksum{$_}->[1];
  my $pver        = $ProtID_AC_SV_Pver_OS_Gene_Chksum{$_}->[2];
  my $OS          = $ProtID_AC_SV_Pver_OS_Gene_Chksum{$_}->[3];
  my $embl_name   = $ProtID_AC_SV_Pver_OS_Gene_Chksum{$_}->[4];
  my $embl_chksum = $ProtID_AC_SV_Pver_OS_Gene_Chksum{$_}->[5];
  my $prot_id = $_.".".$pver;

  #--- EMBL prot chksum and wormpep chksum is identical
  if ( exists $GA_info{$embl_chksum} ){
    for ( my $i = 0; $i < scalar @{$GA_info{$embl_chksum}}; $i = $i+2 ){	
      check_other_name_and_seq_ver($GA_info{$embl_chksum}->[$i], $embl_name, $embl_ac, $embl_sv, $OS, $GA_info{$embl_chksum}->[$i+1]);
    }
  }
  #--- EMBL prot chksum and wormpep chksum is different: requires hand check
  else {
    $no_match++;
    print NOMATCH "DIFF chksum: EMBL: $embl_ac ($embl_name) : Prot_id ($prot_prefix) : chksum: $embl_chksum\n";
  }
}

#--- deal with non-standard names from EMBL /gene tag

foreach (keys %Non_cgc_style_name_Protid_AC_OS_SV_Chksum){  # key is non_CGC style name from EMBL

  for (my $i = 0; $i < scalar @{$Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}}; $i=$i+5 ){
    my $embl_name = $_;
    my $pid = $Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}->[$i]; # prefix only
    my $embl_ac = $Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}->[$i+1];
    my $OS = $Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}->[$i+2];
    my $embl_sv = $Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}->[$i+3];
    my $embl_chksum = $Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}->[$i+4];

    # --- EMBL prot chksum and wormpep chksum is identical
  #  if ( exists $ProtID_Chksum_Ver{$pid} && $ProtID_Chksum_Ver{$pid}->[0] eq $chksum){
    if ( exists $GA_info{$embl_chksum} ){
      for ( my $i = 0; $i < scalar @{$GA_info{$embl_chksum}}; $i = $i+2 ){	
	check_other_name_and_seq_ver($GA_info{$embl_chksum}->[$i], $embl_name, $embl_ac, $embl_sv, $OS, $GA_info{$embl_chksum}->[$i+1]);
      }
    }
    #--- EMBL prot chksum and wormpep chksum is different: requires hand check
    else {
      $no_match++;
      print NOMATCH "DIFF chksum: EMBL: $embl_ac ($embl_name) : Prot_id ($pid) : chksum: $embl_chksum\n";
    }
  }
}

# --- print out update stats
print INFO "Incorporating /gene names from EMBL Release $version\n\n";
print INFO "Update overview:\n";
print INFO "----------------\n";
print INFO "Updated AC due to version change: ", scalar @updated_ACs, "\n\n";
print INFO "New Other_name found: ", scalar @new_Other_names1 + scalar @new_Other_names, "\n";
print INFO "New Other_name to add: ", scalar @new_Other_names, "\n";
print INFO "New AC found as Other_sequence: ", scalar @new_ACs1 + scalar @new_ACs, "\n";
print INFO "New AC to add as Other_sequence: ", scalar @new_ACs, "\n";
print INFO "EMBL /gene name cannot be linked to a CGC name due to different protein checksum: $no_match\n";

print INFO "\nUpdate details:\n";
print INFO "---------------\n";
print INFO "\nUpdated AC based on SV: ", scalar @updated_ACs, "\n@updated_ACs";
print INFO "\nNew Other_names to add: ", scalar @new_Other_names,"\n@new_Other_names";
print INFO "\nNew Other_names that cannot be added: ", scalar @new_Other_names1,"\n@new_Other_names1";
print INFO "\nNew AC to add: ", scalar @new_ACs,"\n@new_ACs";
print INFO "\nNew AC that cannot be added: ", scalar @new_ACs1,"\n@new_ACs1";


sub check_other_name_and_seq_ver {
  my ($locus, $embl_name, $embl_ac, $embl_sv, $OS, $seq_name) = @_;

  #--- EMBL prot chksum and wormpep chksum is identical and EMBL /gene tag = CGC_name
  if ( $locus eq $embl_name || $locus eq lc($embl_name) ) {

    #--- check if AC version needs to be updated
    if ( exists $AC_SV{$embl_ac} && $embl_sv != $AC_SV{$embl_ac} ){

      # for use as stats in the EMBL_name_incorporation_stats.$version file
      # $updated_AC{$embl_ac}++;
      #remove duplicates, eg, AF320903 have multiple /gene
      if ( $updated_AC{$embl_ac}++ == 1 ){
        push (@updated_ACs, $embl_ac." to $locus ($Gene_info{$locus}{'Gene'})\n");

        print DEBUG "\n(A) Update this EMBL AC $embl_ac ($embl_name, SV: $embl_sv) for CGC $locus ($Gene_info{$locus}{'Gene'}, SV: $AC_SV{$embl_ac})\n"; # prints gene id as well
        print PARSE "\n\/\/Update this EMBL AC $embl_ac ($embl_name) for CGC $locus ($Gene_info{$locus}{'Gene'})\n";
        print UPDT "-D Sequence : \"$embl_ac\"\n";
        # cols:      EMBL_AC   gene id                     EMBL_ID                          SV        prot_id                          DE                               OS
        print PARSE "$embl_ac\t$Gene_info{$locus}{'Gene'}\t$AC_ID_Protid_DE{$embl_ac}->[0]\t$embl_sv\t$AC_ID_Protid_DE{$embl_ac}->[1]\t$AC_ID_Protid_DE{$embl_ac}->[2]\t$OS\n";
      }
    }
    #--- add EMBL AC as Other_sequence if not yet linked to a CGC_name in geneace
    if ( !exists $AC_SV{$embl_ac} ){

      $new_AC{$embl_ac}++;
      # for use as stats in the EMBL_name_incorporation_stats.$version file
      # remove duplicates, eg, AF320903 have multiple /gene
       if ( $new_AC{$embl_ac} == 1 ){
	 push (@new_ACs, $embl_ac." ($embl_name) to $locus ($Gene_info{$locus}{'Gene'})\n");
	 print DEBUG "\n(B) Add this EMBL AC $embl_ac ($embl_name) as an other_sequence for $locus ($Gene_info{$locus}{'Gene'})\n";
	 print PARSE "\n\/\/Add this EMBL AC $embl_ac ($embl_name) as an other_sequence for $locus ($Gene_info{$locus}{'Gene'})\n";
	 # cols:      EMBL_AC   gene id                     EMBL_ID                          SV        prot_id                          DE                               OS
	 print PARSE "$embl_ac\t$Gene_info{$locus}{'Gene'}\t$AC_ID_Protid_DE{$embl_ac}->[0]\t$embl_sv\t$AC_ID_Protid_DE{$embl_ac}->[1]\t$AC_ID_Protid_DE{$embl_ac}->[2]\t$OS\n";
       }
    }

    #--- EMBL AC already in Geneace and SV is up-to-date: do nothing
    if ( exists $AC_SV{$embl_ac} && $embl_sv == $AC_SV{$embl_ac} ){
      print DEBUG "\n(C) DO nothing: this EMBL $embl_name is same as a CGC_name and its AC $embl_ac already exists as an other_sequence (same version) in geneace\n";
    }
  }

  #--- EMBL prot chksum and wormpep chksum is identical but EMBL /gene tag != GA CGC name
  #    check if EMBL /gene exists as an Other_name of a CGC_name, if yes, check if AC already is an Other_sequence and if SV needs update
  #                                                               if no, make this /gene an Other_name and AC an Other_sequence
  if ( $locus ne lc($embl_name) && $locus ne "NA" ) {
    &assign_other_name_other_sequence($locus, $embl_name, $embl_ac, $embl_sv, $OS);
  }

  if ( ($locus ne lc($embl_name) || $locus ne $embl_name) && $locus eq "NA") {
    &assign_other_name_other_sequence($locus, $embl_name, $embl_ac, $embl_sv, $OS, $seq_name);
  }

  sub assign_other_name_other_sequence {

    my ($locus, $embl_name, $embl_ac, $embl_sv, $OS, $seq_name) = @_;

    $locus = $seq_name if $locus eq "NA"; # use sequence name if no CGC_name is available for a gene id

    if ( !exists $other{lc($embl_name)} && !exists $other{$embl_name} ) {
      # make this EMBL /gene as an Other_name, and this AC an Other_sequence

      # for use as stats in the EMBL_name_incorporation_stats.$version file
      # remove duplicates, eg, AF320903 have multiple /gene and EMBL /gene which is same as Wormpep sequence name
      $new_AC{$embl_ac}++;

      if ( $new_AC{$embl_ac} == 1 && ( !exists $wormpep_seqs{$embl_name} || exists $main_other{$locus} ) ){
	if ( !exists $Gene_info{$locus}{'Gene'} ){
	  push (@new_ACs1, $embl_ac." ($embl_name) to $locus (Gene id not available: NO Other_sequence assignment)\n");
	  push (@new_Other_names1, $embl_name." ($embl_ac) to $locus (Gene id not available: NO Other_name assignment)\n");
	}
	else {
	  print DEBUG "\n(D) Add EMBL $embl_name as an Other_name & the AC $embl_ac as an Other_sequence for $locus ($Gene_info{$locus}{'Gene'})\n";
	
	  print UPDT "\nGene : \"$Gene_info{$locus}{'Gene'}\"\n";
	  print UPDT "Other_name \"$embl_name\" Accession_evidence \"EMBL\" \"$embl_ac\"\n";
	  push (@new_ACs, $embl_ac." ($embl_name) to $locus ($Gene_info{$locus}{'Gene'})\n");
	  push (@new_Other_names, $embl_name." ($embl_ac) to $locus ($Gene_info{$locus}{'Gene'})\n");
	
	  print PARSE "\/\/Add EMBL $embl_name as an Other_name & $embl_ac as an Other_sequence for $locus ($Gene_info{$locus}{'Gene'})\n";
	  # cols:      EMBL_AC   gene id                     EMBL_ID                          SV        prot_id                          DE                               OS
	  print PARSE "$embl_ac\t$Gene_info{$locus}{'Gene'}\t$AC_ID_Protid_DE{$embl_ac}->[0]\t$embl_sv\t$AC_ID_Protid_DE{$embl_ac}->[1]\t$AC_ID_Protid_DE{$embl_ac}->[2]\t$OS\n";
	}
      }	
    }

    if ( exists $other{$embl_name} ||  exists $other{lc($embl_name)} ) {
	
      #--- check if AC version needs to be updated
      if ( exists $AC_SV{$embl_ac} && $embl_sv != $AC_SV{$embl_ac} ) {
	
	# for use as stats in the EMBL_name_incorporation_stats.$version file
	# remove duplicates, eg, AF320903 have multiple /gene and EMBL /gene which is same as Wormpep sequence name
	$updated_AC{$embl_ac}++;

	if ( $updated_AC{$embl_ac} == 1 ){
	  push (@updated_ACs, $embl_ac." ($embl_name) to $locus ($Gene_info{$locus}{'Gene'})\n");
	  print DEBUG "\n(E) Update this EMBL AC $embl_ac ($embl_name, SV: $embl_sv) for $locus ($Gene_info{$locus}{'Gene'}, SV: $AC_SV{$embl_ac})\n"; # prints gene id as well
	  print PARSE "\n\/\/Update this EMBL AC $embl_ac ($embl_name) for $locus ($Gene_info{$locus}{'Gene'})\n";
	  print UPDT "\nGene : \"$Gene_info{$locus}{'Gene'}\"\n";
	  print UPDT "-D Sequence : \"$embl_ac\"\n";
	  # cols:      EMBL_AC   gene id                     EMBL_ID                          SV        prot_id                          DE                               OS
	  print PARSE "$embl_ac\t$Gene_info{$locus}{'Gene'}\t$AC_ID_Protid_DE{$embl_ac}->[0]\t$embl_sv\t$AC_ID_Protid_DE{$embl_ac}->[1]\t$AC_ID_Protid_DE{$embl_ac}->[2]\t$OS\n";
	}
      }
    }
  }
}

# write longtext
system("perl $script_dir/get_EMBL_longtext.pl -i $temp_dir/EMBL_to_Geneace_pre_longtext.$version -o $output_dir/EMBL_to_Geneace_longtext.$version.ace");

my $end = &runtime;
print "\n$0 started at $start, finished at $end\n";



__END__

NOTE; gene names are lower-cases for string comparison

Deal with standard names from EMBL /gene tag

  (1) EMBL prot checksum and wormpep checksum is identical
      (A) EMBL /gene tag is the SAME as geneace CGC_name
          => check if AC version needs to be updated
          => add EMBL AC as Other_sequence if not yet linked to a CGC_name in geneace

      (B) EMBL /gene tag is NOT the SAME as geneace CGC name
          => check if EMBL /gene exists as an Other_name of a CGC_name,
          => if yes, check if AC already is an Other_sequence and if SV needs update
          => if no, make this /gene an Other_name and AC an Other_sequencet of the CGC locus

   (2) EMBL Checksum and wormpep checksum is different: requires hand check

Deal with non-standard names from EMBL /gene tag

   (1) EMBL prot checksum and wormpep checksum is identical
       (A)check if EMBL /gene exists as an Other_name of a CGC_name,
          => if yes, (a) check if AC version needs to be updated
                     (b) add EMBL AC as Other_sequence if not yet linked to a CGC_name in geneace

          => if no, make this /gene an Other_name and AC an Other_sequencet of the CGC locus
