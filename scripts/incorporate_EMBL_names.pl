#!/usr/local/bin/perl5.8.0 -w

# incorporate_EMBL_names.pl

# What it does: Parse EMBL entry with /gene tags to 
#   1. Make EMBL gene name an other-name of an existing CGC locus, if different from CGC locus
#   2. OR, if EMBL name has corresponding protein (100% match), but not linked to a CGC locus, submit to CGC
#      for approving names 
#   3. Make the AC becomes an other-sequence of corresponding locus


# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-03-16 17:36:02 $ 

use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use strict;
use Getopt::Long;
#use lib "/home/ck1/SCRIPTS";
#use lib "/nfs/team71/worm/ck1/WORMBASE_CVS/scripts/";
use Geneace;

my ($version, $verbose, $debug);
#my $start = &runtime;

GetOptions ("version|v=s"  => \$version,
            "verbose|vb"   => \$verbose, 
	    "debug|d"      => \$debug,
           );

# warning
print "You need to specify version like -v R75 to proceed!\n" if !$version;

#################

# dowload EMBL AC

#################

my $start = &runtime;
my $script_dir = "/nfs/team71/worm/ck1/WORMBASE_CVS/scripts";

# directory, files used for EMBL update
my $output_dir = "/nfs/team71/worm/ck1/EMBL";
my $download = "$output_dir/embl";
my $flatfile= "$output_dir/$version"."_EMBL_entry_no_seqs";

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
    $OS = $1;
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

`getz -e "[swall-org:Caenorhabditis elegans]" > $output_dir/SWALL.$version`;

open(SWALL, "$output_dir/SWALL.$version") || die $!;
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
#system("rm -f $output_dir/SWALL.$version");

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

    if ( exists $ProtID_Chksum_Ver{$pid} ){	            # SWALL
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

##############################################################
#  in verbose mode, output list of EMBL /gene infos
##############################################################

if ($verbose){
  open(GN, ">$output_dir/EMBL_cgc_style_name_dataset.$version") || die $!;
  open(GN_CHK, ">$output_dir/EMBL_cgc_style_name_Non_elegans.$version") || die $!;
  foreach (sort keys %ProtID_AC_SV_Pver_OS_Gene_Chksum){
    print GN "$_ -> @{$ProtID_AC_SV_Pver_OS_Gene_Chksum{$_}}\n";
    print  GN_CHK "$_ -> @{$ProtID_AC_SV_Pver_OS_Gene_Chksum{$_}}\n" if $ProtID_AC_SV_Pver_OS_Gene_Chksum{$_}->[5] eq "NA";
  }
}

if ($verbose){
  open(OUT, ">$output_dir/EMBL_non_cgc_style_name_dataset.$version") || die $!;
  foreach (keys %Non_cgc_style_name_Protid_AC_OS_SV_Chksum){
    print OUT "$_\t@{$Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}}\n";
  }
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

if (!$debug){
  my $wormpep_ver = get_wormbase_version() - 1;
  @wormpep_tbl = `cut -f 1,3,7 /wormsrv2/WORMPEP/wormpep$wormpep_ver/wormpep.table$wormpep_ver`;
}

if ($debug){
  @wormpep_tbl = `cut -f 1,3,7 $output_dir/wormpep.table120`;
}

foreach (@wormpep_tbl){
  chomp;
  my($seq, $locus, $prot_id);
  my @items = split(/\t/, $_);
  if (scalar @items == 3){
    ($seq, $locus, $prot_id) = split(/\t/, $_);
    $prot_id =~ s/\..+//; # drop suffix
  }
  else {
    ($seq, $prot_id) = split(/\t/, $_);
    $locus = "NA";
    $prot_id =~ s/\..+//; # drop suffix
  }
  if (exists $ProtID_Chksum_Ver{$prot_id}) {
    push(@{$GA_info{$ProtID_Chksum_Ver{$prot_id}->[0]}}, $locus, $seq); # key is checksum 
  }
}

open(DEBUG, ">$output_dir/Debug_info.$version") || die $!;
foreach (keys %GA_info){
  print DEBUG "$_ -> @{$GA_info{$_}}\n";
}

# -- here code for fetch other_names of a cgc name
my %main_other = $ga->other_name($db,"main_other");

# get SV of existing other_sequence in Geneace
my %AC_SV;
my $SV_query = "find Accession_number *.*";
my @AC_SV_in_GA = $db->find($SV_query);

foreach (@AC_SV_in_GA){
  $_  =~ /(.+)\.(\d+)/;   #eg. AF03445.1
  $AC_SV{$1} = $2;
}

$db->close;

#######################################################################################
# compare checksum of EMBL /gene and Geneace info : look only for identical checksums
#######################################################################################

open(CGC, ">$output_dir/EMBL_name_to_submit_to_CGC.$version") || die $!;
open(PARSE, ">$output_dir/EMBL_to_Geneace_pre_longtext.$version") || die $!;
open(NOMATCH, ">$output_dir/EMBL_cgc_style_name_with_diff_chksum.$version") || die $!;
open(UPDT, ">$output_dir/EMBL_update.$version.ace") || die $!; # load this before longtext
open(INFO, ">$output_dir/EMBL_name_incorporation_stats.$version") || die $!; 

my ($pepace_chksum, $ga_name, $ga_cds);

# update counter info
my $new_other_name =0;
my $updated_AC = 0;
my $new_AC = 0;

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

  #--- identical Chksum: EMBL /gene tag = GA CGC name 
  if ( exists $GA_info{$embl_chksum} && $embl_name eq $GA_info{$embl_chksum}->[0] ){
    print DEBUG "EMBL $embl_ac ($embl_name) : $prot_prefix -> other_sequence\n";
    check_other_name_and_seq_ver("CGC", $embl_name, $embl_ac, $embl_sv, $OS);
  }

  #--- identical Chksum: EMBL /gene tag != GA CGC name : 
  #    check if EMBL /gene is alreayd a CGC other_name or if to make it an other-name
  if ( exists $GA_info{$embl_chksum} && $embl_name ne $GA_info{$embl_chksum}->[0] && $GA_info{$embl_chksum}->[0] ne "NA" ){

    check_other_name_and_seq_ver($GA_info{$embl_chksum}->[0], $embl_name, $embl_ac, $embl_sv, $OS);
  }
  
  #--- identical Chksum: GA has no CGC name 
  #    ask CGC if EMBL name should be made as CGC name
  if ( exists $GA_info{$embl_chksum} && $GA_info{$embl_chksum}->[0] eq "NA" ){
    print CGC "$embl_name has identical chksum to $GA_info{$embl_chksum}->[1]: Make it CGC name?\n";
  }

  #--- requires hand check
  if ( !exists $GA_info{$embl_chksum} ) {
    print NOMATCH "DIFF chksum:  ($embl_name) : $prot_prefix : $embl_ac : $embl_chksum\n";
  }
}

#--- deal with non-standard names from EMBL /gene tag

foreach (keys %Non_cgc_style_name_Protid_AC_OS_SV_Chksum){  # EMBL
#  print "$_ => @{$Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}}##########\n";

  for (my $i = 0; $i < scalar @{$Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}}; $i=$i+5 ){
    my $pid = $Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}->[$i]; # prefix only
    
    my $ac = $Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}->[$i+1];
    my $os = $Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}->[$i+2];
    my $embl_sv = $Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}->[$i+3];
    my $chksum = $Non_cgc_style_name_Protid_AC_OS_SV_Chksum{$_}->[$i+4];

    # check chksum of non-standard name equals a SWALL checksum
    if ( exists $ProtID_Chksum_Ver{$pid} && $ProtID_Chksum_Ver{$pid}->[0] eq $chksum){
      # is the chksum linked to a GA CGC name?
      if ( exists $GA_info{$chksum}->[0] && exists $main_other{$GA_info{$chksum}->[0]} ){ 
	check_other_name_and_seq_ver($GA_info{$chksum}->[0], $_, $ac, $embl_sv, $OS);
      }
      else {
	print NOMATCH "DIFF chksum: cannot associate EMBL $_ ($ac) with a CGC_name\n"; 
      }
    }
  }
} 

print INFO "\nNew Other_name: $new_other_name\nUpdated AC based on SV: $updated_AC\nNew AC: $new_AC\n";

sub check_other_name_and_seq_ver {
  my ($locus, $embl_name, $embl_ac, $embl_sv, $OS) = @_;
  my %other;

  #  print "$locus(1), $embl_name(2), $embl_ac(3) ======\n";

  if ( exists $main_other{$locus} ){
    my @other_names = @{$main_other{$locus}};
    
    foreach (@other_names){$other{$_}++};
    
    if (!exists $other{$embl_name} ){
      $new_other_name++; $new_AC++;
      print DEBUG "EMBL $embl_name as an Other_name & $embl_ac as an Other_sequence for CGC $locus\n";
      print PARSE "\/\/EMBL $embl_name as an Other_name & $embl_ac as an Other_sequence for CGC $locus\n";
      print UPDT "\nLocus : \"$locus\"\n";
      print UPDT "Other_name \"$embl_name\" Accession_evidence \"EMBL\" $embl_ac\"\n";
      print PARSE "$embl_ac\t$locus\t$AC_ID_Protid_DE{$embl_ac}->[0]\t$embl_sv\t$AC_ID_Protid_DE{$embl_ac}->[1]\t$AC_ID_Protid_DE{$embl_ac}->[2]\t$OS\n";
    }
    if (exists $other{$embl_name} or $locus eq "CGC") {
      #--- check if AC version needs to be updated
      if ( exists $AC_SV{$embl_ac} && $embl_sv != $AC_SV{$embl_ac} ){
        $updated_AC++;
	print DEBUG "Update this EMBL AC $embl_ac for CGC $locus\n";
        print PARSE "\/\/Update this EMBL AC $embl_ac for CGC $locus\n";
        print UPDT "-D Sequence : \"$embl_ac\"\n";
        print PARSE "$embl_ac\t$locus\t$AC_ID_Protid_DE{$embl_ac}->[0]\t$embl_sv\t$AC_ID_Protid_DE{$embl_ac}->[1]\t$AC_ID_Protid_DE{$embl_ac}->[2]\t$OS\n";
      }
      #--- add EMBL AC as Other_sequence
      if ( !exists $AC_SV{$embl_ac} ){
        $new_AC++;
     	print DEBUG "Make EMBL $embl_name an other_name and $embl_ac an other_sequence for CGC $locus\n";
        print PARSE "\/\/Make EMBL $embl_name an other_name and $embl_ac an other_sequence of CGC $locus\n";
        print PARSE "$embl_ac\t$locus\t$AC_ID_Protid_DE{$embl_ac}->[0]\t$embl_sv\t$AC_ID_Protid_DE{$embl_ac}->[1]\t$AC_ID_Protid_DE{$embl_ac}->[2]\t$OS\n";
      }
      #--- EMBL AC already in Geneace and SV is up-to-date: do nothing
      if ( exists $AC_SV{$embl_ac} && $embl_sv == $AC_SV{$embl_ac} ){
	print DEBUG "EMBL name and AC already exists as an other_sequence with same version\n";
      }
    }
  }
}


# write longtext 
system("perl $script_dir/get_EMBL_longtext.pl -i $output_dir/EMBL_to_Geneace_pre_longtext.$version -o $output_dir/EMBL_to_Geneace_longtext.$version.ace");

my $end = &runtime;

print "\n$0 started at $start, finished at $end\n";



__END__

In verbose mode, the script writes two EMBL datasets file:
  EMBL_cgc_style_name: fetch EMBL locus name info
  EMBL_non_cgc_style_name_dataset:  if an EMBL sequence name matches 
                          the sequence of a geneace locus, 
                          that locus gets assigned EMBL AC 
                          as an Other_name
  
and one Wormpep dataset file:
  Wormpep_dataset: assign SWALL checksum to Wormpep seq, protid, locus
