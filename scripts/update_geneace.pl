#! /usr/bin/perl -w

use Tk;
use strict;

#######################
# Create a main windows
#######################

my $mw=MainWindow->new();
$mw->configure (title => "Geneace Update Accelerator\t\tby Chao-Kung Chen\t\tVersion 1.0 (2003/01)",
                background => "grey",
	       );
#$mw->minsize(qw(800 750));#slow
#$mw->geometry('+620+770');

######################################
# Instantiate widgets and arrange them
######################################

#######################################################
# create one frame for menu, another one for work bench
#######################################################

my $btn_frame = $mw ->Frame(relief => 'groove', borderwidth => 2)
                    ->pack(side => 'top', anchor => 'n', expand => 1, fill => 'x');

my $txt_frame = $mw->Frame(relief => 'groove', borderwidth => 2,
                           width => 735, height => 750,
			   #label => 'Genetics update'
			   )
                    ->pack(in => $btn_frame,
			   after => $btn_frame,
			   side => 'top',
			   anchor => 'n',
			   expand => 1, fill => 'x'
			  );


#my @menus;
#foreach (qw/File Help/){
 # push(@menus, $frame->Menubutton(text => $_));
#}

my $menu_File=$btn_frame->Menubutton(text => 'File')->pack (side => 'left', anchor => 'n', fill => 'x');
$menu_File->AddItems(
		     ["command" => "Open Update file", command => \&load_file_input],
		     ["command" => "Save Update file", command => \&save_file_input],
                     ["command" => "Open ace file", command => \&load_file_ace],
		     ["command" => "Save ace file", command => \&save_file_ace],
		     "-",
		     ["command" => "Load ace file to Geneace", command => \&upload_ace],
                     "-",
                     ["command" => "Exit", command => sub{exit}]
                    );
my $menu_Option=$btn_frame->Menubutton(text => 'Options')->pack (side => 'left', anchor => 'n', fill => 'x');
$menu_Option->AddItems(
		       ["command" => "Gene_class & Other_name -> .ace", command => \&geneclass_loci_other_name],
		       ["command" => "Gene mapping -> .ace", command => \&gene_mapping],
                       ["command" => "Laboratory -> .ace", command => \&update_lab],
		       "-",
	               ["command" => "CGC strain -> .ace", command => \&update_strain],
		       "-",
                       ["command" => "EMBL Release to Geneace -> .ace", command => \&update_EMBL_geneace],
		       "-",
                       ["command" => "Caltech function annotation -> .ace", command => \&update_caltech],
		      
                      );
my $btn_ClearBench=$btn_frame->Button(text    => 'Clear_Bench',
				      command => \&clear_bench
				     )
                             ->pack (side => 'left', anchor => 'n', fill => 'x');

my($filename, $info, $btn_save);

$btn_frame->Label(text => 'Filename')->pack(side => 'left', anchor => 'w');
$btn_frame->Entry(textvariable => \$filename)->pack(side => 'left', anchor => 'w', fill => 'x', expand => 1);
$btn_frame->Button(text => 'Load Updt', command => \&load_file_input)->pack(side => 'left', anchor => 'w');
$btn_frame->Button(text => 'Save Updt', command => \&save_file_input)->pack(side => 'left', anchor => 'w');
$btn_frame->Button(text => 'Load ace', command => \&load_file_ace)->pack(side => 'left', anchor => 'w');
$btn_frame->Button(text => 'Save ace', command => \&save_file_ace)->pack(side => 'left', anchor => 'w');
$mw->Label(textvariable => '\$info', -relief => 'ridge')-> pack(side => 'bottom', fill => 'x');

#my $btn_UpdateFrame=$btn_frame->Button(text => 'Update_Input',
	#				   command => \&update_input
		#			  )->pack (side => 'left', anchor => 'n', fill => 'x');

my ($input_window_lbl, $input_window, $ace_window_lbl, $ace_window, $scroll, $input, @content, $update);
&create_windows_and_lbl;


#open (IN, $filename) || die "Can't read in file!";
#chomp (@content = <IN>);
#&geneclass_loci_other_name;


MainLoop();

#############
# subroutines
#############

sub clear_bench {
  $input_window_lbl -> destroy();
  $input_window -> destroy();
  $ace_window_lbl -> destroy();
  $ace_window -> destroy();
  &create_windows_and_lbl;
}

sub load_file_input {
  $info = "Loading file $filename";
  $input_window->delete ('1.0', 'end');
  if (!open(IN, "$filename")){
    $input_window->insert('end', "ERROR: $!\n");
    return;
  }
  while (<IN>){
    $input_window->insert('end', $_);
  }
  close IN;
  $info = "File $filename loaded";
  return;
}

sub save_file_input {
  $info="Saving $filename";
  open(IN, ">$filename");
  print IN $input_window->get('1.0', 'end');
  $info = "File saved";
}

sub load_file_ace {
  
  $info = "Loading file $filename";
  $ace_window->delete ('1.0', 'end');
  if (!open(IN, "$filename")){
    $ace_window->insert('end', "ERROR: $!\n");
    return;
  }
  while (<IN>){
    $ace_window->insert('end', $_);
  }
  close IN;
  $info = "File $filename loaded";
  return;
}

sub save_file_ace {
  $info="Saving $filename";
  open(IN, ">$filename");
  print IN $ace_window->get('1.0', 'end');
  $info = "File saved";
}

sub upload_ace {
 #$filename="/home/ck02/Dump/dump_2003-01-04_A.1.ace";
 my $command=<<END;
pparse $filename
save
quit
END

  my $test_db_dir="/wormsrv1/chaokung/wormdata/CK1_GENEACE/";
  open (Load_testGA,"| tace $test_db_dir ") || die "Failed to upload to test_Geneace";
  print Load_testGA $command;
  close Load_testGA;
}


sub update_input {

  open(OUT, ">>/home/ck02/srg") || die "Can't update file!";
  #$input= "a line of text";
  #$input_window -> insert('end',"$input\n");  # next line is at end of previous one 
  $update=$input_window -> get('1.0', 'end');
  $ace_window -> insert('end', "$update\n");
  print OUT $update;
}

sub geneclass_loci_other_name {
  open (IN, $filename) || die "Can't read in file!";
  #open (IN, "/wormsrv1/chaokung/wormdata/JAH/$filename") || die "Can't read in file!"; 

  my ($gene_class, @Update, $locus, $seq, $rest, @parts, $num_parts, 
      $head, $tail, @variants, $i, $person, $pmid, $evidence, @persons);

  while(<IN>){
    if ($_ =~ ""){}
    if ($_ =~ /^(\w{3,3})\s{2,}(\w+)\s{2,}(.+)/) {
      $gene_class = $1;
      push(@Update,"\n\nGene_Class : \"$gene_class\"\n");
      #@parts = split(/\s{2,}/, $2);
      #print "Description\t\"$parts[1]\"\n";
      push(@Update,"Description\t\"$3\"\n");
      push(@Update,"Designating_laboratory\t\"$2\"\n");
      push(@Update,"CGC_approved\n");
      #print "Designating_laboratory\t\"$parts[0]\"\n";
    }

    if ($_ =~ /^((\w{3,3})-\d+)\s+=\s+(\w+\..+)/) {
      push(@Update, "\n\nLocus : \"$1\"\n");
      push(@Update, "Gene\n");
      print $2, "\n";
      push(@Update, "Gene_Class\t\"$2\"\n");
      push(@Update, "CGC_approved\n");
      push(@Update, "Species\t\"Caenorhabditis elegans\"\n");

      $locus = $1;
      $rest = $3;
      @parts = split(/\s{2,}/, $rest);
      $seq = $parts[0];
      #print "$seq\n";
      $num_parts = scalar @parts;
      if ($num_parts == 1){
	$seq = uc($seq);
	push(@Update, "Genomic_sequence\t\"$seq\"\n"); 
	print "Genomic_sequence\t\"$seq\"\n";  
      }
      if ($num_parts > 1){
	$head = $seq;
	$tail = $seq;
	$head =~ /[\w\d]+\.\d+/;
	$head = $&;
	$tail =~ s/$head//;
	$head = uc($head);
	$seq=$head.$tail;
	if($tail ne ""){
	  @variants=split(/,/, $tail);
	  #print "@variants\n";
	  foreach (@variants){
	    $seq = $head.$_;
	    for ($i = 1; $i < $num_parts; $i++){
	      #print $parts[$i],"\n";
	      if ($parts[$i] =~ /(Paper_evidence)\s(.+)/){
		#print $1, "\n";
		if ($2 =~ /PMID.+/){
		  $pmid = $&;
		  $pmid =~ s/PMID:\s|\[|\]|//;
		  $pmid =~ s/\]//;		
		  print $pmid, "\n";
		  push(@Update, "Genomic_sequence\t\"$seq\"\tPMID_evidence\t\"$pmid\"\n");
		  push(@Update, "Evidence\tPMID_evidence\t\"$pmid\"\n");
		  #print "Genomic_sequence\t\"$seq\"\t$1\t\"$pmid\"\n";	
		}
		else {
		  push(@Update, "Genomic_sequence\t\"$seq\"\tPaper_evidence\t\"$2\"\n");
		  push(@Update, "Evidence\tPaper_evidence\t\"$2\"\n");
		}
	      }
	      if ($parts[$i] =~ /(Person_evidence)\s(.+)/){
		#print $1, "\n";
		$evidence = $1;
		$person = $2;
		$person =~ s/\[|\]//g;
		@persons = split(/,/, $person);
		foreach (@persons){
		  $_ =~ s/^\s//;
		  push(@Update, "Genomic_sequence\t\"$seq\"\tPerson_evidence\t\"$_\"\n");
		  push(@Update, "Evidence\tPerson_evidence\t\"$_\"\n");
		  #print "Genomic_sequence\t\"$seq\"\tPerson_evidence\t\"$_\"\n";
                }
	      }
	    }
	  }
	}
	if ($tail eq "") {
	  for ($i = 1; $i < $num_parts; $i++){
	    print $parts[$i],"\n";
	    if ($parts[$i] =~ /(Paper_evidence)\s(.+)/){
	      #print $1, "\n";
	      if ($2 =~ /PMID.+/){
		  $pmid = $&;
		  $pmid =~ s/PMID:\s|\[|\]|//;
		  $pmid =~ s/\]//;		
		  print $pmid, "\n";
		  push(@Update, "Genomic_sequence\t\"$seq\"\tPMID_evidence\t\"$pmid\"\n");
		  push(@Update, "Evidence\tPMID_evidence\t\"$pmid\"\n");
		  #print "Genomic_sequence\t\"$seq\"\t$1\t\"$pmid\"\n";	
	      }
	      else {
		  push(@Update, "Genomic_sequence\t\"$seq\"\tPaper_evidence\t\"$2\"\n");
		  push(@Update, "Evidence\tPaper_evidence\t\"$2\"\n");
	      }
	    }
	    if ($parts[$i] =~ /(Person_evidence)\s(.+)/){
	      push(@Update, $1);
	      $person = $2;
	      $person =~ s/\[|\]//g;
	      @persons = split(/,/, $person);
	      foreach (@persons){
		$_ =~ s/^\s//;
		push(@Update, "Genomic_sequence\t\"$seq\"\tPerson_evidence\t\"$_\"\n");
		push(@Update, "Evidence\tPerson_evidence\t\"$_\"\n");
		#print "Genomic_sequence\t\"$seq\"\tPerson_evidence\t\"$_\"\n";
              }
	    }
	  }
	}
      }
    }
  }
  foreach (@Update){
    print $_,;
    $ace_window -> insert('end',"$_");
  }
}

sub gene_mapping {

  open (IN, $filename) || die "Can't read in file!";
  #open (IN, "/wormsrv1/chaokung/wormdata/JAH/$filename") || die "Can't read in file!";  
  my (%obj_info, $info, $obj, $remark, $combined, $mapper, @mappers, $date,
      @update, $seq, @yr_mon, $locus, $provider, @providers, $person, $temp,
      @clones, $clone, $allele, @alleles, @distances, $A_non_B, $B_non_A, $ga,
      @gene_allele, @genes, @tested, $gene_class, $rearr);

  ############################################
  # fetch the numbers of the last obj of 
  # ?multi_pt_data ?2_point_data ?pos_neg_data
  ############################################
  
  my $multi_pt_count=4065;
  my $two_pt_count=7142;
  my $pos_neg_count=10683;
  
  
  
  ############################
  # date of update: for remark
  ############################
  
  $date = `date +%y%m%d`; chomp $date;
  
  #####################################
  # parses updates of gene mapping data
  #####################################
  
  while (<IN>){
    chomp;
    
    ######################

    # parse multip_pt_data
    ######################
    
    if ($_ =~ /^\/\/Multi[_-]pt.+/){
      $multi_pt_count++;
      write_ace($multi_pt_count, "\n\nMulti_pt_data : ", $multi_pt_count);
    }
    if ($_ =~ /^Multi[_-]pt_[cC]omment[:|\s:]\s(.+)/){
      $remark = $1; $remark =~ s/^\s//; $remark .= " [$date ck1]";
      write_ace($multi_pt_count, "Remark", $remark);
    }
    if ($_ =~ /^Multi[_-]pt_[cC]ombined_results[:|\s:]\s(.+)/){
      @distances = split(/ /, $1);
      $combined = "Locus \"$distances[0]\" $distances[1] Locus \"$distances[2]\" $distances[3] Locus \"$distances[4]\"";
      write_ace($multi_pt_count, "Combined", $combined)}
    if ($_ =~ /^Multi[_-]pt_[gG]enotype[:|\s:]\s(.+)/){write_ace($multi_pt_count, "Genotype", $1)}
    if ( ($_ =~ /^Multi[_-]pt_[mM]apper[:|\s:]\s(.+)/) || ($_ =~ /Multi_pt_Data_provider[:|\s:]\s(.+)/) ){
      @mappers = split(/,/, $1);
      foreach $mapper (@mappers){
	$mapper =~ s/^\s//;
	write_ace($multi_pt_count, "Mapper", $mapper);
      }
    }
    if ($_ =~ /^Multi[_-]pt_[dD]ate[:|\s:]\s(.+)/){
      @yr_mon = split(/\//, $1);
      write_ace($multi_pt_count, "Date", "$yr_mon[1]-$yr_mon[0]");
    }
    if ($_ =~ /^Multi[_-]pt_GeneA[:|\s:]\s(.+)/){write_ace($multi_pt_count, "Locus_A", $1)}
    if ($_ =~ /^Multi[_-]pt_GeneB[:|\s:]\s(.+)/){write_ace($multi_pt_count, "Locus_B", $1)}
    if ($_ =~ /^Multi[_-]pt_A_non_B_results[:|\s:]\s(.+)/){
      @distances = split(/ /, $1);
      $A_non_B = "Locus \"$distances[0]\" $distances[1] Locus \"$distances[2]\" $distances[3] Locus \"$distances[4]\"";
      write_ace($multi_pt_count, "A_non_B", $A_non_B);
    }
    if ($_ =~ /^Multi[_-]pt_B_non_A_results[:|\s:]\s(.+)/){
      @distances = split(/ /, $1);
      $B_non_A = "Locus \"$distances[0]\" $distances[1] Locus \"$distances[2]\" $distances[3] Locus \"$distances[4]\"";
      write_ace($multi_pt_count, "B_non_A", $B_non_A);
    }
    
    #################
    # parse 2_pt_data
    #################
    
    if ($_ =~ /^\/\/2_[pP]oint.+/){
      $two_pt_count++;
      write_ace($two_pt_count, "\n\n2_point_data : ", $two_pt_count);
    }
    if ($_ =~ /^2_point_[tT]emp[:|\s:]\s(.+)/){$temp = $1; $temp =~ s/°//; write_ace($two_pt_count, "Temperature", $temp)}
    if ($_ =~ /^2_point_[mM]apped_genes[:|\s:]\s(.+)/){
      @genes = split(/,/, $1);
      write_ace($two_pt_count, "Locus_1", $genes[0]);
      write_ace($two_pt_count, "Locus_2", $genes[1]);
    }
    if ($_ =~ /^2_point_[gG]enotype[:|\s:]\s(.+)/){write_ace($two_pt_count, "Genotype", $1)}
    if ($_ =~ /^2_point_[rR]esults[:|\s:]\s(.+)/){write_ace($two_pt_count, "Results", $1)}
    if ($_ =~ /^2_point_[mM]apper[:|\s:](.+)/){
      @mappers = split(/,/, $1);
      foreach (@mappers){write_ace($two_pt_count, "Mapper", $_)}
    }
    if ($_ =~ /^2_point_[dD]ate[:|\s:]\s(.+)/){
      @yr_mon = split(/\//, $1); write_ace($two_pt_count, "Date", "$yr_mon[1]-$yr_mon[0]");
    }
    if ($_ =~ /^2_point_[cC]alculation[:|\s:]\s(.+)/){
      @tested = split(/ /, $1);
      write_ace($two_pt_count, "Tested", "$tested[1] $tested[2]");
    }
    if ($_ =~ /^2_point_[cC]omment[:|\s:]\s(.+)/){
      $remark = $1; $remark =~ s/^\s//; $remark .= " [$date ck1]";
      write_ace($two_pt_count, "Remark", $remark);
    }
    
    ####################
    # parse pos_neg_data
    ####################
    
    if ($_ =~ /^\/\/Df[\/|_]Dp.+/){
      $pos_neg_count++; write_ace($pos_neg_count, "\n\nPos_neg_data : ", $pos_neg_count);
    }
    if ($_ =~ /^Df[\/|_]Dp_[rR]esults[:|\s:]\s(.+)/){write_ace($pos_neg_count, "Results", $1)}
    if ($_ =~ /^Df[\/|_]Dp_[cC]omment[:|\s:]\s(.+)/){
      $remark = $1; $remark =~ s/^\s//; $remark .= " [$date ck1]";
      write_ace($pos_neg_count, "Remark", $remark);
    }
    if ($_ =~ /^Df[\/|_]Dp_[rR]earrangement[:|\s:]\s(.+)/){write_ace($pos_neg_count, "Rearrangement_2", $1)}
    if ($_ =~ /^Df[\/|_]Dp_[gG]enotype[:|\s:]\s(.+)/){write_ace($pos_neg_count, "Genotype", $1)}
    if ( ($_ =~ /^Df[\/|_]Dp_[dD]ata_[pP]rovider[:|\s:]\s(.+)/) || ($_ =~ /^Df[\/|_]Dp_[mM]apper[:|\s:]\s(.+)/) ) {
      @providers = split(/,/, $1);
      foreach (@providers){write_ace($pos_neg_count, "Mapper", $_);}
    }
    if ($_ =~ /^Df[\/|_]Dp_[dD]ate[:|\s:]\s(.+)/){
      @yr_mon = split(/\//, $1); write_ace($pos_neg_count, "Date", "$yr_mon[1]-$yr_mon[0]");
    }
    if ($_ =~ /^Df[\/|_]Dp_[gG]ene[:|\s:]\s(.+)/){
      $ga = $1; $ga =~ s/\)$//;
      @gene_allele = split(/\(/, $ga);
      if ($gene_allele[1] ne ""){
	write_ace($pos_neg_count, "Locus_1", "$gene_allele[0]\" \"$gene_allele[1]"); # w/ allel info
      }
      else {write_ace($pos_neg_count, "Locus_1", "$gene_allele[0]")}; # w/o allele info
    }
    
    
    ##################
    # parse locus data
    ##################
    
    if ( ($_ =~ /^Locus_[lL]ocus_[nN]ame[:|\s:]\s(.+)/) ||
	 ($_ =~ /^Locus_[gG]ene_[nN]ame[:|\s:]\s(.+)/) ){
      $locus = $1;
      $gene_class = $1;
      $gene_class =~ s/-.+//;
      write_ace($locus, "\n\nLocus : ", $locus);
      write_ace($locus, "CGC_approved");
      write_ace($locus, "Gene_Class", $gene_class);
    }
    if ($_ =~ /^Locus_[sS]equence[:|\s:]\s(.+)/){$seq = uc($1); write_ace($locus, "Genomic_sequence", $seq)}
    if ($_ =~ /^Locus_[nN]egative_[mM]ethod[:|\s:]\s(.+)/){
      $remark = $1;
      $remark = "Negative method is ".$remark;
      write_ace($locus, "Remark", $remark);
    }
    if ($_ =~ /^Locus_[pP]ositive_[mM]ethod[:|\s:]\s(.+)/){
      $remark = $1;
      $remark = "Postitive method is ".$remark;
      write_ace($locus, "Remark", $remark);
    }
    if ($_ =~ /^Locus_[cC]omment[:|\s:]\s(.+)/){
      $remark = $1; $remark =~ s/^\s//; $remark .= " [$date ck1]";
      write_ace($locus, "Remark", $remark);
    }
    if ($_ =~ /^Locus_[nN]egative_[cC]lone[:|\s:]\s(.+)/){
      @clones = split(/,|and/, $1);
      foreach $clone (@clones){
	$clone =~ s/^\s//;
	write_ace($locus, "Negative_clone", $clone);
      }
    }
    if ($_ =~ /^Locus_[pP]ositive_[cC]lone[:|\s:]\s(.+)/){
      @clones = split(/,|and/, $1);
      foreach $clone (@clones){
	$clone =~ s/^\s//;
	write_ace($locus, "Positive_clone", $clone);
      }
    }
    if ($_ =~ /^Locus_[dD]escription[:|\s:]\s(.+)/){write_ace($locus, "Description", $1)}
    if ($_ =~ /^Locus_[cC]hromosome[:|\s:]\s(.+)/){write_ace($locus, "Map", $1)}
=start
      if ($_ =~ /^Locus_[dD]ata_provider[:|\s:]\s(.+)/){
	@providers = split(/,|and/, $1);
	foreach $person (@providers){
	  $person =~ s/^\s|\s$//;
	}
      }
=end
=cut
    if ($_ =~ /^Locus_[aA]llele[:|\s:]\s(.+)/){
      @alleles = split(/,/, $1);
      foreach (@alleles){
        $allele = $_;
	$allele =~ s/^\s|\s$//;
	write_ace($locus, "Allele", $allele);
      }
    }
    if ($_ =~ /^Locus_[gG]ene_[pP]roduct[:|\s:]\s(.+)/){write_ace($locus, "Product", $1)}
    
    ###############################
    # parse breakpt (rearrangement)
    ###############################
    
    if ($_ =~ /^Breakpt_Rearrangement[:|\s:]\s(.+)/){$rearr = $1; write_ace($1, "\n\nRearrangement : ", $1)}
    if ($_ =~ /^Breakpt_Comment[:|\s:]\s(.+)/){
      $remark = $1; $remark =~ s/^\s//; $remark .= " [$date ck1]";
      write_ace($rearr, "Remark", $remark);
    }
    if ($_ =~ /^Breakpt_Negative_clone[:|\s:]\s(.+)/){
      @clones = split(/,/, $1);
      foreach $clone (@clones){
	$clone =~ s/^\s//;
	write_ace($rearr, "Clone_outside", $clone);
      }
    }
    if ($_ =~ /^Breakpt_Positive_clone[:|\s:]\s(.+)/){
      @clones = split(/,/, $1);
      foreach $clone (@clones){
	$clone =~ s/^\s//;
	write_ace($rearr, "Clone_inside", $clone);
      }
    }
    if ($_ =~ /^Breakpt_Data_provider[:|\s:]\s(.+)/){
      @providers = split(/,/, $1);
      foreach $person (@providers){
	$person =~ s/^\s//;
	write_ace($rearr, "Author", $person);
      }
    }	

    ##################
    # parse gene_class
    ##################
    
    if ($_ =~ /GeneClass_[gG]ene_name[:|\s:]\s(.+)/){
      $gene_class = $1; 
      write_ace($1, "\n\nGene_Class : ", $1);
      write_ace($1, "CGC_approved");
    }
    if ($_ =~ /GeneClass_[pP]henotype[:|\s:]\s(.+)/){write_ace($gene_class, "Phenotype", $1)}
    if ($_ =~ /GeneClass_[lL]aboratory[:|\s:]\s(.+)/){write_ace($gene_class, "Designating_laboratory", $1)}
    if ($_ =~ /GeneClass_[cC]omment[:|\s:]\s(.+)/){
      $remark = $1; $remark =~ s/^\s//; $remark .= " [$date ck1]";
      write_ace($gene_class, "Remark", $remark);
    }
  }
  
  sub write_ace {
    my($obj, $tag, $value)=@_;
    my(%obj_info, $info, $gene);
    if ( ($tag eq "A_non_B") || ($tag eq "B_non_A") || ($tag eq "Combined") || ($tag eq "Tested") ){
      push(@{$obj_info{$obj}},  "$tag\t$value");
    }
    elsif ($tag eq "CGC_approved"){
      push(@{$obj_info{$obj}},  "$tag");
    }
    else {
      push(@{$obj_info{$obj}},  "$tag\t\"$value\""); # quotation different
    }
    
    foreach (keys %obj_info){
      #print @{$obj_info{$obj}};
      $ace_window -> insert('end',"@{$obj_info{$_}");  # next line is at end of previous one
    }
  }
}

sub create_windows_and_lbl {

  $input_window_lbl = $txt_frame -> Label(
		    text       => 'Genetics Update',   # property of label widget
		    anchor     => 'n', # anchor text to north within the defined area of label by width/height
                    #side       => 'top',
		    relief     => 'groove',         # border style
		    foreground => 'white',
		    background => 'blue',
		    #image     => '/home/ck02/warning.gif',
		    width      => 30,               # unit is char
		    height     => 1
		    ) -> pack(side => 'top', anchor => 'n'); # give it a default place within the main windows

  ##################################################
  # create a scrollable text window for input window
  ##################################################

  $input_window =$txt_frame -> Scrolled(
					"Text",
					width => 100,
					height => 30,
					
				       )-> pack(side => 'top', anchor => 'n', fill => 'x'); # takes "Text" inside () 

  $ace_window_lbl = $txt_frame -> Label(
		    text       => '.ace file',   # property of label widget
		    anchor     => 'n',              # anchor text to north within the defined area of label by width/height
		    relief     => 'groove',         # border style
		    foreground => 'white',
		    background => 'blue',
		    #image      => '/home/ck02/warning.gif',
		    width      => 30,               # unit is char
		    height     => 1
		    ) -> pack(side =>'top', anchor => 'n'); # give it a default place within the main windows

  ################################################
  # create a scrollable text window for ace window
  ################################################
  $ace_window =$txt_frame -> Scrolled(
                                      "Text",
				      width => 100,
				      height => 30
				     )-> pack(side => 'top', anchor => 'n', fill => 'x'); # takes "Text" inside () 

}}


__END__
=start
$menus[0]->pack(side => 'left');
$menus[0]->AddItems(["command" => "Open", command => \$open_file],
                    ["command" => "Make ace file", command => \$make_ace_file],
                    ["command" => "Upload ace file", command => \$upload_ace_file]
                   )

$menus[1]->pack(side => 'left');
=end
=cut
$scroll = $txt_frame -> Scrollbar (orient => 'horizontal');
  $ace_window = $txt_frame -> Entry(relief => 'sunken',
				    width => 100, 
				    #height => 30
				    xscrollcommand => ['set' => $scroll]				   )						
			   ->pack (-side =>'top', anchor => 'n', fill => 'x');      # takes "Text" inside ()
