#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2003-10-23 09:23:09 $ 

use strict;
#use lib -e "/wormsrv2/scripts/"  ? "/wormsrv2/scripts/" : "/nfs/team71/worm/ck1/WORMBASE_CVS/scripts/";
#use lib "/wormsrv2/scripts/"; 

use Wormbase;
use Ace;
use Carp;
use Getopt::Long;
use EMBL_feature_parser;

`date`;
##############################
#     lexical variables
##############################

my ($version, $feature, $name);
my $start = `date`;

GetOptions ("version|v=s"  => \$version,
            "feature|f"    => \$feature,
            "name|n"       => \$name,
           );

###############
#   warning
###############

if (!$version){
  print "\n------------------------------------------------------------------------\n";
  print " You need to specify EMBL version and task(s) to proceed, like:\n";
  print " -v R75 -n <do gene name incorporation> AND/OR -f <fetch feature data>\n";
  print "------------------------------------------------------------------------\n";
  exit(0);
}

#####################
#  initializing obj
#####################

 my $embl = init EMBL_feature_parser($version);

#############################################################
# tasks of get_EMBL_info.pl: 
#      -name:     incorporate EMBL gene name into geneace  
#      -feature:  fetch EMBL feature into camace/Stlace
############################################################# 

# do incorporating EMBL gene name stuff
if ($name) {
print "OKKKK";
}


# do incorporating EMBL feature stuff
if ($feature){
 
  $embl->download_EMBL_release();

  ################################################################################
  #  ANY intersted EMBL features AND/OR their qualifiers can be passed in to the 
  #  get_feature_info method 
  #  Note: be sure to match EMBL terms exactly
  #################################################################################

  my %FEATURE = $embl->get_feature_info(
					polyA_site   => ["evidence", "note"],
					polyA_signal => ["evidence", "note"],
                                       );

  open(FH, ">embl_query") || carp "Can't open $!\n";

  ############################################## 
  #   BLAT EMBL feature coord to chrom coord   
  ##############################################

  my ($AC, %seen, @coords, $coord, $chrom, $L, $R, $psl, @blat);
  my @Unique_coords = ();
  my ($ac, $number, $num_a, $num_b, @ac_feature, %ac_feature);

  foreach (keys %FEATURE){
    #print "$_ -> @{$FEATURE{$_}} #\n";
    $AC = $_;
    foreach (my $i =0; $i < scalar @{$FEATURE{$_}}; $i=$i+6){
      $chrom   = $FEATURE{$_}[$i+1];
      $feature = $FEATURE{$_}[$i+2];
      $coord   = $FEATURE{$_}[$i+3];
      push(@coords, $coord);
    }  
    %seen = ();
    foreach my $e (@coords){push(@Unique_coords, $e) unless $seen{$e}++}
   # print "@Unique_coords\n";
   
    foreach (@Unique_coords){
      if ($_ =~ /(\d+)\.\.(\d+)$/ && $feature ne "polyA_site"){
	$L = $1; $R = $2;
	#print "$AC, $L, $R, $chrom @@\n";
	$number = get_100bp_left_to_EMBL_coord($AC, $L, $R, *FH, $chrom);
        $num_a += $number;
      }
      if ($_ =~ /^(\d+)$/){
	$L = $1; $R = "NA";
	#print "$AC, $L, $R, $chrom @@\n";
	$number = get_100bp_left_to_EMBL_coord($AC, $L, $R, *FH, $chrom);
        $num_b += $number;
      }
    } 
    $ac_feature{$AC} = $num_a + $num_b;
    $num_a = $num_b = 0;
    #print "$AC, $L, $R, $chrom @@\n";
    @Unique_coords = (); 
    @coords = ();
  }

  # check number of features per EMBL acc
  #foreach (sort keys %ac_feature){
  #  print "$_ -> $ac_feature{$_}\n";	
  #}

  #################################
  #   BLATTING and parse psl file
  #################################
  blat("/nfs/team71/worm/ck1/out.psl");
  @blat = parse_psl("/nfs/team71/worm/ck1/out.psl");
 # foreach (@blat){print $_, "\n";}
 
  ####################################################################
  #   get clone name and 30 bp flank seqs for blatted chrom coords   
  ####################################################################
  my %blat_hits = get_30bp_flank_seq(@blat);
   
  # compare number of features per EMBL acc with blat hits of each acc
  # if the number differs, then there are multiple blat hits   
  foreach (sort keys %ac_feature){
    if (exists $blat_hits{$_}){
      print "$_ -> $ac_feature{$_} : ", scalar @{$blat_hits{$_}}, "\n" if $ac_feature{$_} != scalar @{$blat_hits{$_}};
    }
  }
}

`date`;

sub get_100bp_left_to_EMBL_coord {
  ######################################################################################
  # for mapping EMBL seq coords, which spcify features in submitted seq, to chrom coords
  ######################################################################################

  my ($ac, $L, $R, $fh, $chrom) = @_;
  my ($database, @lines, $line, @seqs, $DNA, $embl_query);

  my $ac_count = 0; 

 # print "$ac, $L, $R, $fh, $chrom --\n";
 # print "$L", "\n";
  $line = ();
  @seqs = `pfetch -q $ac`;
  foreach (@seqs){chomp; $line .= $_}
  
  $DNA = substr($line, $L-60, 60);
  print $fh ">$ac.$L.$chrom \n";
  $ac_count++;
  print $fh lc($DNA);
  print $fh "\n";  
  

  if ($R ne "NA"){
    $DNA = substr($line, $R-60, 60);
    print $fh ">$ac.$R.$chrom\n";
    $ac_count++;
    print $fh lc($DNA);  
    print $fh "\n";  
  }
  return $ac_count;
}

sub blat {

  my $output = shift;
  my $blat = "/nfs/disk100/wormpub/blat/blat";
  my $database = "/nfs/team71/worm/ck1/SEQUENCES/CHROMOSOME_all.dna";

  # run blat
  `rm -f $output`;
  `$blat -noHead $database embl_query $output`;
}

sub parse_psl {
  ############################################################################
  # get these columns from psl file
  # 2: mismatch, 9: strand, 10: query seq name, 14: chrom, 17: end chrom coord
  ############################################################################

  my $psl = shift;
  my @blat = `cut -f 2,9,10,14,17 $psl`;
  return @blat;	  
}

sub get_30bp_flank_seq {
  #######################################################################
  # retrieve 30 bp flank seqs and clone name from parsed BLAT psl columns  
  #######################################################################

  my @blat = @_;
  my $chrom_dir = "/nfs/disk100/wormpub/DATABASES/current_DB/CHROMOSOMES/";
  my %blat_hits;
  my $blat_count = 0;
  
  foreach (@blat){
    chomp;
    $blat_count = 1;
    my ($miss_match, $strand, $seq, $chrom, $coord)= split(/\s+/, $_);

   # print "$miss_match, $strand, $seq, $chrom, $coord  1#\n";
 
    my $embl_chrom = $seq;
    $seq =~ s/\..+\..+//;
    print $seq, "\n";

    # count blat hits of each EMBL acc (to be compare with the expected one from above)
    push(@{$blat_hits{$seq}}, $blat_count);
 
    $embl_chrom  =~ s/.+\.//;
   # print $embl_chrom, "  2#\n";

    $chrom =~ s/CHROMOSOME_//;

    my $dna_file = "$chrom_dir/CHROMOSOME_".$chrom.".dna";
    
    my @line = `egrep "[atcg]" $dna_file`;
    my $line;
    foreach (@line){chomp; $line .= $_}
    
    my $L_30 = $coord -31;
    my $R_30 = $coord;
    my ($DNA_L, $DNA_R);
    
    if ($strand eq "-"){
      $DNA_L = substr($line, $L_30, 30);
      $DNA_L = reverse $DNA_L; $DNA_L =~ tr/atcg/tagc/;
      $DNA_R = substr($line, $R_30, 30);
      $DNA_R = reverse $DNA_R; $DNA_R =~ tr/atcg/tagc/;
    }
    if ($strand eq "+"){
      $DNA_L = substr($line, $L_30, 30);
      $DNA_R = substr($line, $R_30, 30);
    }
    my $coords = Coords_converter->invoke("/nfs/disk100/wormpub/DATABASES/current_DB");
    my @clone = $coords->LocateSpan("CHROMOSOME_$chrom",$coord, $coord);

    if ($chrom eq $embl_chrom){
      print "$DNA_L  $DNA_R -> $clone[0] -> same CHROM verified -> $miss_match mismatch\n";
    }
    else {
      print "$DNA_L  $DNA_R -> $clone[0] -> CHROM NOT verified -> $miss_match mismatch\n";
    }
  }
  return %blat_hits;
}




__END__
