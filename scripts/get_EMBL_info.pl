#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2003-11-28 11:32:20 $ 

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'}; 
use Wormbase;
use Ace;
use Carp;
use Getopt::Long;
use EMBL_feature_parser;
use Coords_converter;
use Feature_mapper;

##############################
#     lexical variables
##############################

my ($version, $feature, $name);

GetOptions ("version|v=s"  => \$version,
            "feature|f"    => \$feature,
            "name|n"       => \$name,
           );

###############
#   warning
###############

if (!$version){
  print "\n-----------------------------------------------------------------------------\n";
  print " You need to specify EMBL version and task(s) to proceed, like:\n";
  print " -v R75 -n <to do gene name incorporation> AND/OR -f <to fetch feature data>\n";
  print "-----------------------------------------------------------------------------\n";
  exit(0);
}

my $start = &runtime;

#####################
#  initializing obj
#####################

 my $embl = init EMBL_feature_parser($version);
 my $out_dir = "/nfs/team71/worm/ck1/EMBL";
 my $embl_query = "$out_dir/embl_60bp_L_flank";   # input FASTA format file for blat
 my $database = "/nfs/disk100/wormpub/DATABASES/current_DB/";

#############################################################
# tasks of get_EMBL_info.pl: 
#      -download  download EMBL release flatfile first 
#      -name:     incorporate EMBL gene name into geneace  
#      -feature:  fetch EMBL feature into camace/Stlace
############################################################# 

# download before doing -f or -n
$embl->download_EMBL_release();


# do incorporating EMBL gene name stuff
if ($name) {
print "To be continued...";
}


# do incorporating EMBL feature stuff
if ($feature){

  ################################################################################
  #  ANY intersted EMBL features AND/OR their qualifiers can be passed in to the 
  #  get_feature_info method 
  #  Note: be sure to match EMBL terms exactly
  #################################################################################
  my %FEATURE = $embl->get_feature_info(
					polyA_site   => ["evidence", "note"],
					polyA_signal => ["evidence", "note"],
                                       );
 
  open(FH, ">$embl_query") || carp "Can't open $!\n";

  ############################################## 
  #   BLAT EMBL feature coord to chrom coord   
  ##############################################
  
  my ($AC, %seen, @coords, $coord, $chrom, $qualifier, $remark, $Feature, $L, $R, $psl, @blat);
  my @Unique_coords = ();
  my ($ac, $number, $num_a, $num_b, @ac_feature, %ac_feature, $psl_file);

  foreach (keys %FEATURE){

    $num_a = $num_b = 0;
    my %Coord_F_Q_R =();   # Feature, Qualifier, Remark
    my %Coord_F_Q_R_1 =(); # temp hash for dealing with duplication

    $AC = $_;
    foreach (my $i =0; $i < scalar @{$FEATURE{$_}}; $i=$i+6){
      $chrom    = $FEATURE{$_}[$i+1];
      $Feature  = $FEATURE{$_}[$i+2];
      $coord    = $FEATURE{$_}[$i+3];
      $qualifier = $FEATURE{$_}[$i+4];
      $remark     = $FEATURE{$_}[$i+5];
       
      # remove duplicates resulted in null qualifiers (evidence/note) and its remark 
      if ($qualifier ne "NA"){
	push(@{$Coord_F_Q_R{$coord}}, $Feature, $qualifier, $remark);
      }
      # $qualifier is "NA"
      else {
	push(@{$Coord_F_Q_R_1{$coord}}, $Feature, "NA", "NA");
      }
    }  
    # merge coords that have NA as qualifier only into %Coord_F_Q_R
    foreach my $c (keys %Coord_F_Q_R_1){
      if (!exists $Coord_F_Q_R{$c}){
	push(@{$Coord_F_Q_R{$c}}, $Coord_F_Q_R_1{$c}->[0], "NA", "NA");
      }
    }
    
    # pass EMBL feature coords to subroutine to generate an input file for BLATting
    foreach my $c (keys %Coord_F_Q_R){
      if ( $c =~ /(\d+)\.\.(\d+)$/ && $Coord_F_Q_R{$c}->[0] eq "polyA_signal" ){
	$L = $1; $R = $2;
	$R = ($2-$L);  # number of bp downstream of the first coord of polyA_signal 
	$number = get_60bp_left_to_EMBL_coord($AC, $L, $R, *FH, $chrom) if $R < 10; # filter out signal seqs > 10 bp apart
	$num_a += $number;
      }
      if ( $c =~ /(\d+)\.\.(\d+)$/ && $Coord_F_Q_R{$c}->[0] eq "polyA_site" ){
        $L = $1; $R = $2;

        # when polyA_site coords look like x..y and y-x =1
        if ($R - $L == 1){
          $R = "NA";  # second coord of polya_site is assigned "NA"
          $number = get_60bp_left_to_EMBL_coord($AC, $L, $R, *FH, $chrom);
          $num_b += $number;
        }
      }
      if ( $c =~ /^(\d+)$/ && $Coord_F_Q_R{$c}->[0] eq "polyA_site" ){
	$L = $1; $R = "NA";  # seconde coord of polyA_site is assigned "NA"
	$number = get_60bp_left_to_EMBL_coord($AC, $L, $R, *FH, $chrom);
	$num_b += $number;
      }
    }
    # hash to record number of featrue obj. per AC
    $ac_feature{$AC} = $num_a + $num_b;
  }

  #################################
  #   BLATting and parse psl file
  #################################
 
  $psl_file = blat();                # BLAT output in $psl_file

  @blat = parse_psl($psl_file);
  open(BLAT, ">$out_dir/Blat_columns_output") || die $!;
  foreach (@blat){print BLAT $_, "\n";}
 
  ####################################################################
  #   get clone name and 30 bp flank seqs for blatted chrom coords   
  ####################################################################

  my %blat_hits = get_30bp_flank_seq(@blat);

  # compare number of features per EMBL acc with blat hits of each acc
  #if the number differs, then there are multiple blat hits   
  open(EXCP, ">$out_dir/exception_ACs") || die $!;  
  foreach (sort keys %ac_feature){
    if (exists $blat_hits{$_}){
      print EXCP "$_ -> $ac_feature{$_} : ", scalar @{$blat_hits{$_}}, "\n" if $ac_feature{$_} != scalar @{$blat_hits{$_}};
    }
  }
}

print "Script started at $start, finished at ", &runtime, "\n";
exit(0);

####################################################
#            s u b r o u t i n e s
####################################################


sub get_60bp_left_to_EMBL_coord {
  ######################################################################################
  # this generates a file (specified by $fh) with left flanking seq of the left coord 
  # for mapping EMBL seq coords, which spcify features in submitted seq, to chrom coords
  ######################################################################################

  my ($ac, $L, $R, $fh, $chrom) = @_;
  my ($database, @lines, $line, @seqs, $DNA, $embl_query, $L_flank);

  my $ac_count = 0; 

  #print "$ac, $L, $R, $fh, $chrom -----\n";
  $line = ();
  @seqs = `pfetch -q $ac`;
  foreach (@seqs){chomp; $line .= $_}

  if ($R =~ /\d+/){
    $DNA = substr($line, $L-61, 60);        # for blatting for polyA_signal, the EMBK flanking seq. are 60 bp left to the 1st coord
    $L_flank = substr($line, $L-31, 30);    # for used as verification later
    print $fh ">$ac.$L.+$R.$chrom.$L_flank\n";
    $ac_count++;
    print $fh lc($DNA);
    print $fh "\n";  
  }

  else {
    $DNA = substr($line, $L-61, 60);         # for polyA_site, the EMBK flanking seq. are 60 bp left to the 1st coord
    $L_flank = substr($line, $L-31, 30);     # for used as verification later
    print $fh ">$ac.$L.$chrom.$L_flank\n";
    $ac_count++;
    print $fh lc($DNA);  
    print $fh "\n";  
  }
  # return number of feature objs per AC
  return $ac_count;
}

sub blat {
  ############
  # BLATting
  ############
 
  my $psl_file = "$out_dir/embl.psl";                                   # specify BLAT psl output file 
  my $blat = "/nfs/disk100/wormpub/blat/blat";                          # blat binary
  my $DNAs = "/nfs/team71/worm/ck1/SEQUENCES/CHROMOSOME_all.dna";       # sequence input for blat

  # run blat and output a psl file as $output
  # BLAT query file is $embl_query generated by get_60bp_left_to_EMBL_coord routine
  `$blat -noHead $DNAs $embl_query $psl_file`;

  return $psl_file;
}

sub parse_psl {
  ############################################################################
  # get these columns from psl file
  # 2: mismatch, 9: strand, 10: query seq name, 14: chrom, 17: end chrom coord
  ############################################################################

  my $psl = shift;

  # get these columns:
  # miss_match(2), strand(9), AC(10), chrom(14), coord of start bp matched (16), coord of end bp matched (17)
  my @blat = `cut -f 2,9,10,14,16,17 $psl`;

  return @blat;	  
}

sub get_30bp_flank_seq {
  #######################################################################
  # retrieve 30 bp flank seqs and clone name from parsed BLAT psl columns  
  #######################################################################
  
  my $mapper = Feature_mapper->new( $database);

  my @blat = @_;
  my $chrom_dir = "$database/CHROMOSOMES/";
  my (%blat_hits, $L_flank);
  my $blat_count = 0;
  my $wb = 13999;  # feature obj starting at number 14000

  # acefile output for uploading to db
  open (ACE, ">$out_dir/feature_ace.ace") || die $!;
  foreach (@blat){
    chomp;
    $blat_count = 1;
    # $coord_S (start), $coord_E (end) of BLAT coord
    my ($miss_match, $strand, $ac, $chrom, $coord_S, $coord_E)= split(/\s+/, $_);

    $chrom =~ s/CHROMOSOME_//;   # chrom returned from BLAT
    # $ac is of this format : AF326938.1885.+5.NA (AC.feature coord.+5(plus L coord = R coord).Chrom(NA in this case)
    # OR                    : AY052771.1521.III
    my ($acc, $embl_chrom, $R_end, $feature); 
      
    if ($ac =~ /^(.+)\.(.+)\.(.+)\.(.+)/){
      
      $acc = $1;
      $embl_chrom = $3;         # chrom returned from parsing EMBL AC
      $L_flank = $4;
      $feature = "polyA_site";
    }
      
    if ($ac =~ /^(.+)\.(.+)\.\+(\d+)\.(.+)\.(.+)/){
      $acc = $1;
      $R_end = $3;	# true coord needs to add $R_end, see later
 #     print "$R_end (C)\n";
      $embl_chrom = $4;
      $L_flank = $5;
      $feature = "polyA_signal_sequence";   
    }
  
    # criteria: if BLAT mismatch under 3 bp and chrom from EMBL and BLAT is identical, then generate a feature with flank seqs
    #-----------------------------------------------------------------------------------------------------------------------
    if ($miss_match == 0 ){    
#      print "$acc ###\n";
        
      # count blat hits of each EMBL acc (to be compare with the expected one from above)
      push(@{$blat_hits{$acc}}, $blat_count);

      my $dna_file = "$chrom_dir/CHROMOSOME_".$chrom.".dna";
      
      my @line = `egrep "[atcg]" $dna_file`;
      my $line;
      foreach (@line){chomp; $line .= $_}
      
      my ($L_30, $R_30, $DNA_L, $DNA_R);

      # take $L_flank as EMBL left flank 
      $DNA_L = lc($L_flank);  

      # 30 bp flanking seq. on both side      
      if ($strand eq "-"){
        #$DNA_L = substr($line, $coord_S, 30);
	#$DNA_L = reverse $DNA_L; $DNA_L =~ tr/atcg/tagc/;

	$DNA_R = substr($line, $coord_S-30-$R_end-1, 30) if $feature eq "polyA_signal_sequence";
        $DNA_R = substr($line, $coord_S-30-2, 30)        if $feature eq "polyA_site";   # 2 bp followd by right flank
	$DNA_R = reverse $DNA_R; $DNA_R =~ tr/atcg/tagc/; 
      }
      if ($strand eq "+"){
	#$DNA_L = substr($line, $coord_E-30, 30); 
	
	$DNA_R = substr($line, $R_end+$coord_E+1, 30) if $feature eq "polyA_signal_sequence";
        $DNA_R = substr($line, $coord_E+2, 30)        if $feature eq "polyA_site"; # 2 bp following by right flank
      }

      # get clone info
      my $coords = Coords_converter->invoke("$database");
      my @clone = $coords->LocateSpan("CHROMOSOME_$chrom",$coord_S, $coord_E);
      my @V = $mapper->_check_for_match($line, $DNA_L, $DNA_R);
      my $verify = $V[1]-$V[0] if @V;
      #print "$acc   : $verify -> $V[0] -- $V[1]\n";

      if ( ($feature eq "polyA_site" && $verify == 1) || ($feature eq "polyA_signal_sequence" && $verify == $R_end) ){
	$wb++;
	
	print ACE "\nFeature : \"WBsf0$wb\"\n";
	print ACE "Sequence \"$clone[0]\"\n";
	print ACE "Flanking_sequences\t\"$clone[0]\"\t\"$DNA_L\"\t\"$DNA_R\"\n";
	print ACE "Defined_by_sequence\t\"$acc\"\n";
	print ACE "Method\t\"$feature\"\n";
      }
    }
  }
  return %blat_hits;
}



__END__
