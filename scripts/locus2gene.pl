#!/usr/local/bin/perl5.8.0 -w
#
# locus2gene.pl
# 
# by Keith Bradnam, aged 12 and half
#
# A script to convert ?Locus objects into the new ?Gene object
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2004-04-07 10:34:08 $   

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Ace;
use Getopt::Long;


################################################################
# set up input/output streams
################################################################

my $file; # where input file is located
my $fix_locus; # option to fix locus names created by Hides_under tags
my $fix_2point; # option for fixing one way links to Locus from 2_point_data
my $fix_pos_neg; # option for fixing one way links to Locus class from Pos_neg_data
my $fix_multi_pt; # option for fixing one way links to Locus class from multi_pt_data and from Results hash of multi_pt_data
my $fix_map; # option for fixing links to Locus class from Main_marker and Locus tags in Map class
my $fix_rearrangement; # option for fixing links to Locus class from Locus_inside and Locus_outside tags of Rearrangement class

GetOptions ("file=s"       => \$file,
	    "fix_locus"    => \$fix_locus,
	    "fix_2point"   => \$fix_2point,
	    "fix_pos_neg"  => \$fix_pos_neg,
	    "fix_multi_pt" => \$fix_multi_pt,
	    "fix_map"      => \$fix_map,
	    "fix_rearrangement" => \$fix_rearrangement);

my $db_path = "/wormsrv1/geneace2";

&process_locus_dump if ($file);

&fix_locus_class if($fix_locus);

&fix_two_point_data_class if ($fix_2point);

&fix_pos_neg_data_class if ($fix_pos_neg);

&fix_multi_pt_data_class if ($fix_multi_pt);

&fix_map_class if ($fix_map);

&fix_rearrangement_class if ($fix_rearrangement);



################################################################
# main conversion of input file
################################################################

sub process_locus_dump{

  open(IN,"<$file") || die "Can't open input file\n";
  open(OUT,">$file.gene") || die "Couldn't open output file\n";
  
  # pattern to match timestamps
  my $ts = "-O \"(\\d{4}\-\\d{2}\-\\d{2}_\\d{2}:\\d{2}:\\d{2}\.?\\d?_\\S+|original)\"";
  my $id = 1; # the numerical part of the gene ID
  my $date = "2004-02-04";
  $/ = ""; # reset input line separator



# process input file, changing lines as appropriate

  while(<IN>){
    
    my $id_padded = sprintf "%08d" , $id;
    my $name = "WBGene$id_padded";  
    my ($public_name, $cgc_name, $non_cgc_name, $sequence_name, $other_name);
    
    
    
    # Change Non_CGC_name field to Other_name
    
#  s/Non_CGC_name\s+$ts\s+(\"[\w\d]+\.[\d\w:]+\")\s+$ts\s+(\d+)\s+$ts\s+(\d+)\s+$ts/CDS $3 $5 $7/g;
    if(m/Non_CGC_name/){
      s/Non_CGC_name\s+$ts\s+(\S+)\s+$ts/Other_name -O \"$1\" $2 -O \"$3\"/;
      $non_cgc_name = $2;
    }
    # Change Transcript_name or Pseudogene_name to Sequence_name
    s/Transcript_name/Sequence_name/g if (m/Transcript_name/);
    s/Pseudogene_name/Sequence_name/g if (m/Pseudogene_name/);
    
    # Remove Sequence_name splice variant suffix
    s/Sequence_name\s+$ts\s+(\"[A-Z0-9]*?\.\d+)[a-z]\"\s+/Sequence_name -O \"$1\" $2\" /g;
    
    # grab other names when present: e.g. cgc_name, sequence_name or other_name
    if(m/\s+CGC_name\s+$ts\s+(\S+)\s+$ts/){
      $cgc_name = $2;
    }
    if(m/\s+Sequence_name\s+$ts\s+(\S+)\s+$ts/){
      $sequence_name = $2;
    }
    if(m/\s+Other_name\s+$ts\s+(\S+)\s+$ts/){
      $other_name = $2;
    }
    
    
    # Assign a 'public_name' based on a priority CGC_name -> Sequence_name -> Other_name
    # For non-elegans genes, maybe not use sequence name (where it exists) for public_name
    if($cgc_name){
      $public_name = $cgc_name;
    }
    elsif($sequence_name){
      $public_name = $sequence_name;
    }
    elsif($non_cgc_name){
      $public_name = $non_cgc_name;
    }
    else{
      $public_name = $other_name;
    }
    
    
    # Main conversion to Gene object with basic history/identity information
    s/^Locus :.*/Gene : \"$name\"\nVersion 1\nVersion_change 1 now \"WBPerson1971\" Imported \"Initial conversion from geneace\"\nPublic_name $public_name\nLive/;
    
    # Remove CGC_approved, now assumed if CGC_name is present
    s/Type\s+$ts\s+Gene\s+$ts\s+CGC_approved\s+$ts/Type -O \"$1\" Gene -O \"$2\"/;
    
  # Remove Type.Gene lines which have no more data (i.e. they are the last line of the object)
    s/Type\s+$ts\s+Gene\s+$ts\n$/\n/;  
    
    # prune Type.Gene part of tree in other cases where tags follow Gene tag
    s/Type\s+$ts\s+Gene\s+$ts\s+//g;  
    
    # change tag names to new shorter formats
    s/Gene_information\s+/Gene_info /g;
    s/Molecular_information\s+/Molecular_info /g;
    
    # fix other tag name changes
    s/Contained_in\s+/In_cluster /g;
    
    # Negative_clone now uses Evidence hash and needs Author_evidence rather than text
    s/Negative\s+$ts\s+Negative_clone\s+$ts\s+(\S+)\s+$ts\s+(\S+? \S+?)\s+$ts/Negative_clone -O \"$2\" $3 -O \"$4\" Author_evidence $5 -O \"$6\"/g;
    
    # Same for Inside_rearr / Outside_rearr
    s/Positive\s+$ts\s+Inside_rearr\s+$ts\s+(\S+)\s+$ts\s+(\S+? \S+?)\s+$ts/Inside_rearr -O \"$2\" $3 -O \"$4\" Author_evidence $5 -O \"$6\"/g;
    s/Negative\s+$ts\s+Outside_rearr\s+$ts\s+(\S+)\s+$ts\s+(\S+? \S+?)\s+$ts/Outside_rearr -O \"$2\" $3 -O \"$4\" Author_evidence $5 -O \"$6\"/g;
    
    print OUT;
    
    $id++;
#  last if $id>205;
  }
  
  
# tidy up
  close(IN);
  close(OUT);
}

exit(0);


 

###############################################################################

sub fix_locus_class{
  # This is just to deal with Gene objects created by Hide_under/Representative_for tags.  These are
  # objects with CGC style names, so just need to look up what the correct gene ID is in the Gene_name
  # class and make a simple fix file

  # After loading the fix file, you still then need to go into geneace and remove any Gene objects
  # that have CGC-style names.  Believe me, if you have loaded the fix file you can then safely
  # remove any Gene object which isn't called WBGene*

  open(OUT, ">fix_locus.ace") || die "Can't create output file\n";

  my $db = Ace->connect(-path  => $db_path) || do { print "Connection failure: ",Ace->error; die();};

  # This query gets all genes which link to other genes using old CGC style names,
  # which independently creates them via XREF
  my $query = "Find Gene WBGene* WHERE Representative_for";
  push(my @genes, $db->find($query));
  foreach my $gene (@genes){
    my @loci = $gene->Representative_for;
    foreach my $locus (@loci){

      my ($obj) = $db->fetch(-class=>'Gene_name',-name=>"$locus");
      my $gene_id = $obj->Public_name_for;      
      print OUT "Gene : \"$gene\"\n-D Representative_for $locus\n\nGene : $gene\nRepresentative_for $gene_id\n\n";
    }


  }

  $db->close;  
  close(OUT);
}


############################################################################



sub fix_two_point_data_class{
  # This is just to deal with Gene objects created by one-way links to Locus from the 2_point_data clas.  
  # These need to be changed to use the Gene_1 and Gene_2 tags rather than Locus_1 and Locus_2 and also to use 
  # correct Gene ID

  # After loading the fix file, you still then need to go into geneace and remove any Gene objects
  # that have CGC-style names.  Believe me, if you have loaded the fix file you can then safely
  # remove any Gene object which isn't called WBGene*

  open(OUT, ">fix_2_pt_data.ace") || die "Can't create output file\n";

  my $db = Ace->connect(-path  => $db_path) || do { print "Connection failure: ",Ace->error; die();};

  # This query gets all 2_point_data objects which link to Locus objects using old CGC style names,
  my $query = "Find 2_point_data * WHERE (Locus_1 OR Locus_2)";
  push(my @twopoints, $db->find($query));
  foreach my $twopoint (@twopoints){
    print "$twopoint\n";
    my $locus_1 = $twopoint->Locus_1;
    my $locus_2 = $twopoint->Locus_2;

    my $extra_1 = $twopoint->Locus_1(2);
    my $extra_2 = $twopoint->Locus_2(2);

    # write ace file according to whether one or two loci were present in the 2_point object
    if ($locus_1){
      next if (($locus_1 eq "act-123") || ($locus_1 =~ m/\w+P\d+/)); # exeception for Polymorphisms and Gene_clusters

      my ($obj) = $db->fetch(-class=>'Gene_name',-name=>"$locus_1");
      my $gene_id = $obj->Public_name_for;      
      print OUT "2_point_data : \"$twopoint\"\n-D Locus_1 $locus_1\n\n";
      if($extra_1){
	print OUT "2_point_data : \"$twopoint\"\nGene_1 $gene_id $extra_1\n\n";
      }
      else{
	print OUT "2_point_data : \"$twopoint\"\nGene_1 $gene_id\n\n";
      }
    }
    if ($locus_2){
      next if (($locus_2 eq "act-123") || ($locus_2 =~ m/\w+P\d+/)); # exeception for Polymorphisms and Gene_clusters

      my ($obj) = $db->fetch(-class=>'Gene_name',-name=>"$locus_2");
      my $gene_id = $obj->Public_name_for;      
      if($extra_2){
	print OUT "2_point_data : \"$twopoint\"\nGene_2 $gene_id $extra_2\n\n";
      }
      else{
	print OUT "2_point_data : \"$twopoint\"\nGene_2 $gene_id\n\n";
      }
      print OUT "2_point_data : \"$twopoint\"\n-D Locus_2 $locus_2\n\n";
      print OUT "2_point_data : \"$twopoint\"\nGene_2 $gene_id\n\n";
    }

  }

  $db->close;  
  close(OUT);
}

############################################################################



sub fix_pos_neg_data_class{
  # This is just to deal with Gene objects created by one-way links to Locus from the Pos_neg_data class.  
  # These need to be changed to use the Gene_1 and Gene_2 tags rather than Locus_1 and Locus_2 and also to use 
  # correct Gene ID

  # After loading the fix file, you still then need to go into geneace and remove any Gene objects
  # that have CGC-style names.  Believe me, if you have loaded the fix file you can then safely
  # remove any Gene object which isn't called WBGene*

  open(OUT, ">fix_pos_neg_data.ace") || die "Can't create output file\n";

  my $db = Ace->connect(-path  => $db_path) || do { print "Connection failure: ",Ace->error; die();};

  # This query gets all 2_point_data objects which link to Locus objects using old CGC style names,
  my $query = "Find Pos_neg_data * WHERE (Locus_1 OR Locus_2)";
  push(my @pos_negs, $db->find($query));
  foreach my $pos_neg (@pos_negs){
    print "$pos_neg\n";
    my $locus_1 = $pos_neg->Locus_1;
    my $locus_2 = $pos_neg->Locus_2;

    my $extra_1 = $pos_neg->Locus_1(2);
    my $extra_2 = $pos_neg->Locus_2(2);

    # write ace file according to whether one or two loci were present in the 2_point object
    if ($locus_1){
      next if (($locus_1 eq "act-123") || ($locus_1 =~ m/\w+P\d+/)); # exeception for Polymorphisms and Gene_clusters

      my ($obj) = $db->fetch(-class=>'Gene_name',-name=>"$locus_1");
      my $gene_id = $obj->Public_name_for;      
      print OUT "Pos_neg_data : \"$pos_neg\"\n-D Locus_1 $locus_1\n\n";
      if($extra_1){
	print OUT "Pos_neg_data : \"$pos_neg\"\nGene_1 $gene_id $extra_1\n\n";
      }
      else{
	print OUT "Pos_neg_data : \"$pos_neg\"\nGene_1 $gene_id\n\n";
      }
    }
    if ($locus_2){
      next if (($locus_2 eq "act-123") || ($locus_2 =~ m/\w+P\d+/)); # exeception for Polymorphisms and Gene_clusters

      my ($obj) = $db->fetch(-class=>'Gene_name',-name=>"$locus_2");
      my $gene_id = $obj->Public_name_for;      
      if($extra_2){
	print OUT "Pos_neg_data : \"$pos_neg\"\nGene_2 $gene_id $extra_2\n\n";
      }
      else{
	print OUT "Pos_neg_data : \"$pos_neg\"\nGene_2 $gene_id\n\n";
      }
      print OUT "Pos_neg_data : \"$pos_neg\"\n-D Locus_2 $locus_2\n\n";
      print OUT "Pos_neg_data : \"$pos_neg\"\nGene_2 $gene_id\n\n";
    }

  }

  $db->close;  
  close(OUT);
}


############################################################################

sub fix_multi_pt_data_class{
  # This is just to deal with Gene objects created by one-way links to Locus from the Multi_pt_data class.  
  # These need to be changed to use the Gene_A and Gene_B tags rather than Locus_A and Locus_B and also to use 
  # correct Gene ID.  Also there is a separate link to Locus objects using the Locus tag

  # After loading the fix file, you still then need to go into geneace and remove any Gene objects
  # that have CGC-style names.  Believe me, if you have loaded the fix file you can then safely
  # remove any Gene object which isn't called WBGene*

  open(OUT, ">fix_multi_pt_data.ace") || die "Can't create output file\n";

  my $db = Ace->connect(-path  => $db_path) || do { print "Connection failure: ",Ace->error; die();};

  # This query gets all Multi_pt_data objects which link to Locus objects using old CGC style names,
  my $query = "Find Multi_pt_data * WHERE (Locus OR Locus_1 OR Locus_2)";
  push(my @multipoints, $db->find($query));
  foreach my $multipoint (@multipoints){
    print "$multipoint\n";
    my $locus_A = $multipoint->Locus_A;
    my $locus_B = $multipoint->Locus_B;
    my $locus   = $multipoint->Locus;

    my $extra_1 = $multipoint->Locus_A(2);
    my $extra_2 = $multipoint->Locus_B(2);
    my $extra_3 = $multipoint->Locus(2);

    # write ace file according to whether Locus, Locus_A, and Locus_B tags  were present in the 
    # Multi_pt_data object

    if ($locus_A){
      next if (($locus_A eq "act-123") || ($locus_A =~ m/\w+P\d+/)); # exeception for Polymorphisms and Gene clusters

      my ($obj) = $db->fetch(-class=>'Gene_name',-name=>"$locus_A");
      my $gene_id = $obj->Public_name_for;      

      # remove old data
      print OUT "Multi_pt_data : \"$multipoint\"\n-D Locus_A $locus_A\n\n";

      # add new data
      if($extra_1){
	print OUT "Multi_pt_data : \"$multipoint\"\nGene_A $gene_id $extra_1\n\n";
      }
      else{
	print OUT "Multi_pt_data : \"$multipoint\"\nGene_A $gene_id\n\n";
      }
    }

    if ($locus_B){
      next if (($locus_B eq "act-123") || ($locus_B =~ m/\w+P\d+/)); # exeception for Polymorphisms and Gene clusters

      my ($obj) = $db->fetch(-class=>'Gene_name',-name=>"$locus_B");
      my $gene_id = $obj->Public_name_for;      

      # delete old data
      print OUT "Multi_pt_data : \"$multipoint\"\n-D Locus_B $locus_B\n\n";

      # add new data
      if($extra_2){
	print OUT "Multi_pt_data : \"$multipoint\"\nGene_B $gene_id $extra_2\n\n";
      }
      else{
	print OUT "Multi_pt_data : \"$multipoint\"\nGene_B $gene_id\n\n";
      }
    }

    if ($locus){
      next if (($locus eq "act-123") || ($locus =~ m/\w+P\d+/) || ($locus eq "RW#L17")); # exeception as not a Gene
      next if ($locus eq "TCPAR1"); # exeception as not a Gene

      my ($obj) = $db->fetch(-class=>'Gene_name',-name=>"$locus");
      my $gene_id = $obj->Public_name_for;      
      # delete old link
      print OUT "Multi_pt_data : \"$multipoint\"\n-D Locus $locus\n\n";

      # add new data
      if($extra_3){
	print OUT "Multi_pt_data : \"$multipoint\"\nGene $gene_id $extra_3\n\n";
      }
      else{
	print OUT "Multi_pt_data : \"$multipoint\"\nGene $gene_id\n\n";
      }
    }
  }

  print OUT "\n\/\/ New fix\n\n";
  
  # This is just to deal with Gene objects created by one-way links to Locus from the Multi_pt_data class.  
  # This is to look up Locus names that might be in the Multi_counts hash under the Results tag

  # After loading the fix file, you still then need to go into geneace and remove any Gene objects
  # that have CGC-style names.  Believe me, if you have loaded the fix file you can then safely
  
  # remove any Gene object which isn't called WBGene*

  
  # This query gets all Multi_pt_data objects which link to Locus objects using old CGC style names,
  my $query2 = "Find Multi_pt_data * WHERE Results";
  push(my @multipoints2, $db->find($query2));

  foreach my $multipoint (@multipoints2){
    print "$multipoint\n";

    # Store results in @row array
    my @rows;

    # grab row data if it exists for each of the three different types
    (push(@rows, $multipoint->at("Results.Combined")->row)) if ($multipoint->at("Results.Combined"));
    (push(@rows, $multipoint->at("Results.A_non_B")->row)) if ($multipoint->at("Results.A_non_B"));
    (push(@rows, $multipoint->at("Results.B_non_A")->row)) if ($multipoint->at("Results.B_non_A"));

    # Loop through loop elements changing Locus for Gene and substituting Gene IDs for locus names
    foreach (my $i=0; $i<@rows;$i++){

      # need to add new line when there are multiple rows of results
      ($rows[$i] = "\n\nMulti_pt_data : \"$multipoint\"\n"."$rows[$i]") if ($rows[$i] eq "A_non_B");
      ($rows[$i] = "\n\nMulti_pt_data : \"$multipoint\"\n"."$rows[$i]") if ($rows[$i] eq "B_non_A");
      ($rows[$i] = "\n\nMulti_pt_data : \"$multipoint\"\n"."$rows[$i]") if ($rows[$i] eq "Combined");

      if(($rows[$i] eq "Locus") && ($rows[$i+1] eq "act-123")){
	$i += 2;
	next;
      }
      $rows[$i] = "Gene" if ($rows[$i] eq "Locus");
      if($rows[$i] =~ m/\w+\-.*/){
	my ($obj) = $db->fetch(-class=>'Gene_name',-name=>"$rows[$i]");
	my $gene_id = $obj->Public_name_for;
	$rows[$i] = $gene_id;
      }
    }
    if(@rows){
      print OUT "\nMulti_pt_data : \"$multipoint\"\n";
      print OUT "-D Results";
      
      print OUT "@rows\n\n";
    }
  }

  $db->close;  
  close(OUT);
}



###################################################################################################

sub fix_map_class{

  open(OUT, ">fix_map.ace") || die "Can't create output file\n";

  my $db = Ace->connect(-path  => $db_path) || do { print "Connection failure: ",Ace->error; die();};

  my $query = "Find Genetic_map *";

  push(my @maps, $db->find($query));
  foreach my $map (@maps){
    print "\n\n$map\n";

    # Two types of Locus name to correct
    my @main_loci;
    push(@main_loci, $map->at("Main_Marker.Main_Locus")->col);

    foreach my $main_locus (@main_loci){
      my ($obj) = $db->fetch(-class=>'Gene_name',-name=>"$main_locus");
      my $gene_id = $obj->Public_name_for;
      print OUT "\nMap : \"$map\"\n";
      print OUT "-D Main_Locus $main_locus\n\n";
      print OUT "\nMap : \"$map\"\n";
      print OUT "Main_Locus $gene_id\n\n";
    }

    
    my @loci;
    push(@loci, $map->at("Contains.Locus")->col);

    foreach my $locus (@loci){
      print "$locus\n";
      next if (($locus eq "act-123") || ($locus eq "rDNA_cluster") || ($locus =~ m/\w+P\d+/)); # exeception for Polymorphisms and Gene clusters
      my ($obj) = $db->fetch(-class=>'Gene_name',-name=>"$locus");
      my $gene_id = $obj->Public_name_for;
      print OUT "\nMap : \"$map\"\n";
      print OUT "-D Locus $locus\n\n";
      print OUT "\nMap : \"$map\"\n";
      print OUT "Gene $gene_id\n\n";
    }
  }  



  $db->close;  
  close(OUT);

}


###################################################################################################

sub fix_rearrangement_class{

  open(OUT, ">fix_rearrangement.ace") || die "Can't create output file\n";

  my $db = Ace->connect(-path  => $db_path) || do { print "Connection failure: ",Ace->error; die();};

  my $query = "Find Rearrangement * WHERE Positive OR Negative";

  push(my @rearrangements, $db->find($query));
  foreach my $rearrangement (@rearrangements){
    print "$rearrangement\n";

    # Two types of Locus name to correct...
    my @inside;
    push(@inside, $rearrangement->at("Positive.Locus_inside")->col) if ($rearrangement->at("Positive.Locus_inside"));

    foreach my $locus (@inside){
      next if (($locus eq "act-123") || ($locus eq "rDNA_cluster") || ($locus =~ m/\w+P\d+/)); # exeception for Polymorphisms and Gene clusters

      my ($obj) = $db->fetch(-class=>'Gene_name',-name=>"$locus");
      my $gene_id = $obj->Public_name_for;
      print OUT "\nRearrangement : \"$rearrangement\"\n";
      print OUT "-D Locus_inside $locus\n\n";
      print OUT "\nRearrangement : \"$rearrangement\"\n";
      print OUT "Gene_inside $gene_id\n\n";
    }

    # Two types of Locus name to correct
    my @outside;
    push(@outside, $rearrangement->at("Negative.Locus_outside")->col) if ($rearrangement->at("Negative.Locus_outside"));

    foreach my $locus (@outside){
      next if (($locus eq "act-123") || ($locus eq "rDNA_cluster") || ($locus =~ m/\w+P\d+/)); # exeception for Polymorphisms and Gene clusters

      my ($obj) = $db->fetch(-class=>'Gene_name',-name=>"$locus");
      my $gene_id = $obj->Public_name_for;
      print OUT "\nRearrangement : \"$rearrangement\"\n";
      print OUT "-D Locus_outside $locus\n\n";
      print OUT "\nRearrangement : \"$rearrangement\"\n";
      print OUT "Gene_outside $gene_id\n\n";
    }


  }  



  $db->close;  
  close(OUT);

}
