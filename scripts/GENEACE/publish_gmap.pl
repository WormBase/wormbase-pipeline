#!/usr/local/bin/perl5.8.0 -w
# 
# gene ace_check.el
#
# Chao-Kung Chen
#
# Script to prepare for 2003 C. elegans genetic map publication
#
# Last updated by: $Author: ck1 $
# Last updated on: $Date: 2004-03-19 11:59:04 $


use strict;
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'}; 
use Wormbase;
use Ace;
use Spreadsheet::WriteExcel;


#########################
#   global variables
#########################

my $tace = &tace;     
my $def_dir = "/wormsrv1/geneace/wquery";
my $geneace = "/wormsrv1/geneace";
my $file_dir= "/nfs/team71/worm/ck1/CGC";

#################################
# table maker definitions files
#################################

my $table_one=<<EOF;
  Table-maker -p "$def_dir/CGC_table-1.def" quit 
EOF

my $locus_to_seq=<<EOF;
  Table-maker -p "$def_dir/locus_to_CDS.def" quit 
EOF

my $table_three_map=<<EOF;
  Table-maker -p "$def_dir/CGC_table-3_map.def" quit 
EOF
my $table_three_intmap=<<EOF;
  Table-maker -p "$def_dir/CGC_table-3_intmap.def" quit 
EOF

#######################################
# database connection and fetch data
#######################################

my $db = Ace->connect(-path  => $geneace,
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};
__END__
my @loci = $db->fetch(-class => 'Locus',
		      -name  => '*');

my (%locus_seqs, $main_name, @seqs);

#foreach (@loci){
#  @seqs =();	
#  if ( defined $_ -> Other_name_for(1) ){
#    $main_name = $_ -> Other_name_for(1);
#  }
#  if ( defined $_ -> Genomic_sequence(1) ){
#    push(@seqs, $_ -> Genomic_sequence(1) );
#  }
#  if ( defined $_ -> Transcript(1) ){
#    push(@seqs, $_ -> Transcript(1) );
#  }
#  if ( defined $_ -> Pseudogene(1) ){
#    push(@seqs, $_ -> Pseudogene(1) );
#  }
#  @seqs = "NA" if ! @seqs;	  
#  push(@{$locus_seqs{$_}}, @seqs);
#}

#foreach (keys %locus_seqs){
#	print "$_ ->@{$locus_seqs{$_}}\n";
#}


#################################################
# define excel table paremeters for output format
#################################################

my $rows_per_page = 65;
my($left_header, $right_header, $left_cells, $right_cells);
my ($wbk, $wksheet, @data, @data_1);

#######################
#    start working
#######################

#&table_one;

&table_two;

#&table_three;

#&table_four;

#&other_name_table;

#&marker_table;

#&int_map_to_map_loci;  # intepolated loci to display on genetics map

#&make_multipt_obj;
#&get_LG_from_strain_in_Locus_obj;

###########################
# s u b r o u t i n e s
###########################

sub get_LG_from_strain_in_Locus_obj {
  my (@loci, @strains);
  my $query = "Find Locus * where CGC_approved & !Map & !Interpolated_map_position & Strain & Species = \"*elegans\"";
  my $strain_query = "Find Locus * where !Map & !Interpolated_map_position & Strain & Species = \"*elegans\"; >strain";
  push(@loci, $db->find( $query ) );
  push(@strains, $db->find( $strain_query ) );

  foreach (sort @loci){
    my @ST_symbol = $_->Strain(1);
    foreach my $e (@ST_symbol){
      my $info = $e->Genotype(1);
      print $info, "##\n";
      if ($info =~ /$_(\(.+\)).+(I|II|III|IV|V|X)/){
	print "$_-> $1 -> $2\n";
      }
    }
  }
}


sub make_multipt_obj {
  
  $file_dir= "/nfs/team71/worm/ck1/CGC";

  # concat lists of gene summary data from gmap and parse fields 1,4,5 from resulting file
  my @gene_boundaries = `cat $file_dir/gene_summary_* | cut -f 1,4,5`;
 
  my %gene_boundaries;
  foreach (@gene_boundaries){
    chomp;
    my($locus, $LLM, $RLM) = split(/\s+/, $_);
    push(@{$gene_boundaries{$locus}}, $LLM, $RLM);
  }

  # make a list of loci that have map+allele+sequence but no mapping_data
  my $multi_query  = "find Locus * where map & allele & (CDS|transcript|pseudogene) & !mapping_data";
  push(my @loci_to_make_multi_pt_obj, $db->find($multi_query) );


  my $multi = 4134;  # current last number of multi_pt obj in db 

  open(ACE, ">$file_dir/inferred_multi_pt_obj_with_allele") || die $!;
  foreach my $locus (@loci_to_make_multi_pt_obj){
    foreach my $e ( $locus-> Allele(1) ){ 
      $multi++;
      print ACE "\n\nLocus : \"$locus\"\n";
      print ACE "Multi_point $multi\n";
      print ACE "\n\nMulti_pt_data : $multi\n";
      print ACE "Locus_A \"$locus\" \"$e\"\n";
      print ACE "Locus \"$locus\" \"$e\"\n";
      print ACE "Combined Locus \"$gene_boundaries{$locus}->[0]\" 1 Locus \"$locus\" 1 Locus \"$gene_boundaries{$locus}->[1]\"\n";
      print ACE "Remark \"Data inferred from $e, sequence of $locus and interpolated map position (which became genetics map)\" Inferred_automatically\n";
      
    }
  }
} 


sub table_four {

  @data = `cat $file_dir/gene_summary_ALL`;

  init("Table_four.xls");

  my $num = (scalar @data + 2) / (2 * $rows_per_page); 
  $num =~ /(\d+)\.(\d+)/;        
  my $page = $1 + 1  if $2;
  $page = $num if !$2;

  my @titles = qw(GENE LG POS LFG RFG LLIM RLIM CLONE);   
  my $line = my $increase = my $header = 0;

  for(my $i = 1; $i < $page + 1; $i++){
    print "page $i => $left_header : $left_cells : $right_header : $right_cells\n";
    # header
    my $row = 0;
    foreach my $e (@titles){
      $wksheet->write($header, $row , $e); 
      $row++;
    }

    my $cell = $header + 1;
    
    for (my $line = $increase ; $line < $left_cells - 1; $line++){
      
      my ($cell_1, $cell_2, $cell_3, $cell_4, $cell_5, $cell_6, $cell_7, $cell_8) = split(/\t/, $data[$line]);
    #  print "$cell_1, $cell_2, $cell_3, $cell_4, $cell_5, $cell_6, $cell_7, $cell_8\n";
      chomp $cell_8;
      
      $wksheet->write($cell, 0, $cell_1);
      $wksheet->write($cell, 1, $cell_2);
      $wksheet->write($cell, 2, $cell_3);
      $wksheet->write($cell, 3, $cell_4);
      $wksheet->write($cell, 4, $cell_5);
      $wksheet->write($cell, 5, $cell_6);
      $wksheet->write($cell, 6, $cell_7);
      $wksheet->write($cell, 7, $cell_8);
      
      $cell++;
    }  
    $increase = $left_cells - 2;     
    $header = $header + 65;
    $left_cells  = 65 + $left_cells; 
  }

}

sub other_name_table {

  init("Other_name_table.xls");

 # my $other_name = "find Gene_name * where *-* AND !(Cb-* OR Cr-*)"; 
#  my ($main_name, @seqs, @other); 
#  open(ON, ">/nfs/team71/worm/ck1/CGC/other_name_table") || die $!;
 @data=(); 

#  @other = $db->find($other_name);
#  foreach (sort @other){	
#    if ( defined $_ -> Other_name_for(1) ){	
#      $main_name = $_ -> Other_name_for(1);
#      if (scalar @{$locus_seqs{$main_name}} == 1){
#	print ON "$_\t$main_name\t@{$locus_seqs{$main_name}}\n";
#	push(@data, "$_\t$main_name\t@{$locus_seqs{$main_name}}\n");
#      }
#      else {
#	my @SEQ = sort @{$locus_seqs{$main_name}};
	
#	print ON "$_\t$main_name\t$SEQ[0]\n";
#        push(@data, "$_\t$main_name\t$SEQ[0]\n");

#	shift @SEQ;
#	foreach my $e (@SEQ){
#	  print ON "\t\t$e\n";
#          push(@data, "\t\t$e\n");
#	}
#      }
#    }
#  }
  
  @data = `less /nfs/team71/worm/ck1/CGC/other_name_table`;
  print scalar @data, "\n";
  my $num = (scalar @data + 2) / (2 * $rows_per_page); # plus 1 for header cell row, so each 59 lines fit into one spreadsheet page
  $num =~ /(\d+)\.(\d+)/;
  my $page = $1 + 1  if $2;
  $page = $num if !$2;
  print "$page pages##########\n";
  my $titles = "Other-name\tMain name\tSequence";  

  my $line = my $increase = my $header = 0;

  for(my $i = 1; $i < $page + 1; $i++){
    print "page $i => $left_header : $left_cells : $right_header : $right_cells\n";

    excel_panel_3($titles, $increase, $header, @data);

    $increase = $right_cells - 2;     
    $header = $header + 65;
    $left_cells  = 65 + $right_cells;  $right_cells  = 130 + $right_cells;
  }
}

sub marker_table {

  init("Marker_table.xls");

  my $marker_query  = "find Locus * where map & mapping_data & polymorphism";
  my (%marker, @I, @II, @III, @IV, @V, @X); 
  #  open(MARKER,   ">/nfs/team71/worm/ck1/CGC/marker_table") || die $!;
  
  push( my @markers, $db->find($marker_query) );
  foreach (sort @markers){
    my $LG    = $_ -> Map(1);
    my $pos   = $_ -> Map(3);
   # $pos = sprintf("%4.2f", "$pos");  
    my $clone = $_ -> Positive(2);
    push(@I,   "$_\t$LG\t$pos\t$clone\n") if $LG eq "I";
    push(@II,  "$_\t$LG\t$pos\t$clone\n") if $LG eq "II";
    push(@III, "$_\t$LG\t$pos\t$clone\n") if $LG eq "III";
    push(@IV,  "$_\t$LG\t$pos\t$clone\n") if $LG eq "IV";
    push(@V,   "$_\t$LG\t$pos\t$clone\n") if $LG eq "V";
    push(@X,   "$_\t$LG\t$pos\t$clone\n") if $LG eq "X";
  }
  
  @data = (@I,@II,@III,@IV,@V,@X);
  my $num = (scalar @data + 2) / (2 * $rows_per_page); 
  $num =~ /(\d+)\.(\d+)/;        
  my $page = $1 + 1  if $2;
  $page = $num if !$2;

  my $titles = "MARKER\tLG\tPOS\tCLONE";   
  my $line = my $increase = my $header = 0;

  for(my $i = 1; $i < $page + 1; $i++){
    print "page $i => $left_header : $left_cells : $right_header : $right_cells\n";
    
    excel_panel_4($titles, $increase, $header, @data);
    
    $increase = $right_cells - 2;     
    $header = $header + 65;
    $left_cells  = 65 + $right_cells;  $right_cells  = 130 + $right_cells;
  }
}

sub table_one {
 
  open (FH, "echo '$table_one' | $tace $geneace | ") || die $!;
  open(TBL1, ">/nfs/team71/worm/ck1/CGC/table_1") || die $!;
  open(TBL1_1, ">/nfs/team71/worm/ck1/CGC/table_1-1") || die $!;
  my ($lab, $allele, $pi, $location);
  @data= (); 

  while(<FH>){
    chomp;
    #print $_, "\n";
    if ($_ =~ /^\"([A-Z]{2,2})\"\s+\"(.+)\"\s+\"(.+)\"\s+\"(.+)\"\s+\"(.+)\"/){
      $lab      = $1;
      $allele   = $2;
      $pi       = $4; #$pi = sprintf("%-30s", $pi);
      $location = $5; 
      $location =~ s/\\//g;
      print TBL1 "$lab\t$allele\t$pi\t$location\n";
    }
    elsif ($_ =~ /^\"([A-Z]{3,3})\"\s+\"(.+)\"\s+\"(.+)\"\s+\"(.+)\"/){
      my $person;
      $lab      = $1;
      $pi       = $3;
      $location = $4;
      print TBL1_1"$lab\t$pi\t$location\n";
    }
  } 
  close FH; 
  close TBL1; close TBL1_1;
}

sub table_two {

  open(TBL2,   ">/nfs/team71/worm/ck1/CGC/table_2") || die $!;
  open(TBL2_1, ">/nfs/team71/worm/ck1/CGC/table_2-1") || die $!;

  my ($des, $pheno, $lab);
  my @GC = $db->fetch(-class => 'Gene_class',
		      -name  => '*');

  foreach (@GC){
    $des =();
    $pheno =();
    $lab =();

    if ( defined $_ -> Designating_laboratory(1) ){
      $lab = $_ -> Designating_laboratory(1);
    }
    if ( defined $_ -> Description(1) ){
      $des = $_ -> Description(1);
    }
    if ( defined $_ -> Phenotype(1) ){
      $pheno = $_ -> Phenotype(1);
    }
    if ($lab  && $des){
      print TBL2  "$_\t$lab\t$des\n";
    }
    if (!$lab && $des){ 
      print TBL2  "$_\t-\t$des\n";
    }
    if (!$lab && !$des && $pheno){
      print TBL2_1 "$_\t-\t$pheno\n";
    }
  }
  close TBL2; close TBL2_1;
}

sub table_three {

  my ($pos, $pos_i, $LG, $LG_i, $error, %vague_map_loci); 
  open(TBL3,   ">$file_dir/table_3") || die $!;

 # # sorting all loci that have map position by map position 
#  $map_pos_query =  "Find locus * where MAP AND NEXT AND NEXT =\"*\"";
#  my @maps = $db->find( $map_pos_query );
#  my %Locus_map;
#  foreach (@maps){
#    $Locus_map{$_} = $_->Map(3);
    
#  }

  # get list of genes with mapping ambiguity > 3 map units
  my @gene_boundaries = `cut -f 1,6,7 $file_dir/gene_summary_ALL`;
  foreach (@gene_boundaries){
    chomp;
    my $line = $_;
    $line =~ s/\*//g;
    my ($gene, $LLIM, $RLIM) = split(/\t+/, $line);

    if ($LLIM <  0 && $RLIM <  0) { $error = -$LLIM + $RLIM }
    if ($LLIM >  0 && $RLIM >  0) { $error =  $RLIM - $LLIM }
    if ($LLIM <= 0 && $RLIM >= 0) { $error =  $RLIM - $LLIM }
    print "$gene -> $error\n";
    $vague_map_loci{$gene} = $error if $error >= 3;
  }
    
  my $C = "Find Locus * where CGC_approved & Map & (CDS|Transcript|Pseudogene)";
  my $P = "Find Locus * where CGC_approved & Interpolated_map_position & (CDS|Transcript|Pseudogene) & !Allele";
  my $N = "Find locus * where CGC_approved & !(cds|transcript|pseudogene) & map & !mapping_data";
  my $G = "Find Locus * where CGC_approved & Map &  Mapping_data & !(CDS|Transcript|Pseudogene)";


  my (%cat_C, %cat_P, %cat_G, %cat_N, $e);	
  push( my @cat_C, $db->find($C) ); foreach $e (@cat_C){$cat_C{$e}++}
  push( my @cat_P, $db->find($P) ); foreach $e (@cat_P){$cat_P{$e}++} 
  push( my @cat_N, $db->find($N) ); foreach $e (@cat_N){$cat_N{$e}++}
  push( my @cat_G, $db->find($G) ); foreach $e (@cat_G){$cat_G{$e}++}
  
  foreach (@loci){
    $LG  = $_ -> Map(1) if defined $_->Map(1);
    $pos = $_ -> Map(3) if defined $_->Map(3);
    $LG_i  = $_ -> Interpolated_map_position(1) if defined $_->Interpolated_map_position(1);
    $pos_i = $_ -> Interpolated_map_position(2) if defined $_->Interpolated_map_position(2); 
    $pos = "" if !defined $_ -> Map(3) && defined $LG;
    print TBL3 "$_\t$LG\t$pos\tC\n"     if $cat_C{$_};
    print TBL3 "$_\t$LG_i\t$pos_i\tP\n" if $cat_P{$_};
 #   print TBL3 "$_\t$LG\t$pos\tN\n" if $LG && $pos ne "" && $cat_N{$_};
    print TBL3 "$_\t$LG\t$pos\tN\n" if $cat_N{$_};
    print TBL3 "$_\t$LG\t$pos\tG\n" if !exists $vague_map_loci{$_} && $cat_G{$_}; 
    print TBL3 "$_\t$LG\t$pos\tN\n" if exists $vague_map_loci{$_} && $cat_G{$_};
 }   
}
__END__
=start
  foreach (@cat_C){
    $LG  = $_ -> Map(1);
    $pos = $_ -> Map(3); 
    print TBL3 "$_\t$LG\t$pos\tC\n";
    push(@data, "$_\t$LG\t$pos\tC\n");
  }
  foreach (@cat_P){
    $LG  = $_ -> Interpolated_map_position(1);
    $pos = $_ -> Interpolated_map_position(2); 
    print TBL3 "$_\t$LG\t$pos\tP\n";
    push(@data, "$_\t$LG\t$pos\tP\n");
  }

  # define vague loci as those having map units > 3 
  foreach (@cat_N){
    $LG  = $_ -> Map(1);
    $pos = $_ -> Map(3) if defined $_ -> Map(3);
    $pos = "" if !defined $_ -> Map(3);
    print TBL3 "$_\t$LG\t$pos\tN\n" if exists $vague_map_loci{$_};
    print TBL3 "$_\t$LG\t$pos\tN\n" if $pos eq "";
    push(@data, "$_\t$LG\t$pos\tN\n") if exists $vague_map_loci{$_};
    push(@data, "$_\t$LG\t$pos\tN\n") if $pos eq "";

  }
  # loci of this category have map unit ambiguity < 3
  foreach (@cat_G){
    $LG  = $_ -> Map(1);
    if ( defined $_ -> Map(3) ){
      $pos = $_ -> Map(3); 
      print TBL3 "$_\t$LG\t$pos\tG\n" if !exists $vague_map_loci{$_};
      print TBL3 "$_\t$LG\t$pos\tN\n" if exists $vague_map_loci{$_};
      push(@data, "$_\t$LG\t$pos\tG\n") if !exists $vague_map_loci{$_};
      push(@data, "$_\t$LG\t$pos\tN\n") if exists $vague_map_loci{$_};
    }
  }
  open(TBL_3, ">$file_dir/table_3_sorted") || die $!;
  @data = sort {$a<=>$b} @data;
  foreach (@data){
    print TBL_3 $_;
  }
=end
=cut
}

sub int_map_to_map_loci {

  # display some loci without map and mapping_data but has allele and seq. connection and interpolated_map_position
  # by converting their Interpolated_map to Map

  my $int_loci  = "find Locus * where !map & !mapping_data & allele & (CDS|transcript|pseudogene) & Interpolated_map_position & species =\"*elegans\"";
  my %INT_loci;

  push( my @int_loci, $db->find($int_loci) );
  foreach (@int_loci){
    my $int_map = $_ -> Interpolated_map_position(1);
    my $int_pos = $_ -> Interpolated_map_position(2);
    print "\nLocus : \"$_\"\n";
    print "-D Interpolated_map_position \n";
    print "\nLocus : \"$_\"\n";
    print "Map \"$int_map\" Position $int_pos\n";
  }
}
    
sub init {

  my $table = shift;
  # allocated lines from data to current spreadsheet page
  $left_header = 1;  $right_header = 66;      
  $left_cells = 65;  $right_cells  = 130;  

  # create an Excel workbook wbk and initilize parameters for Excel row/column

  $wbk = Spreadsheet::WriteExcel->new($table);
  $wksheet = $wbk ->addworksheet();
}  

