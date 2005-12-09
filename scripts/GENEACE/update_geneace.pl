#!/usr/local/bin/perl5.8.0 -w

use Tk;
use strict;
use Cwd;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use GENEACE::Geneace;

################
# database paths
################

my $tace = &tace;
my $ga_dir = "/nfs/disk100/wormpub/DATABASES/geneace";

#######################
# Create a main windows
#######################

my $mw=MainWindow->new();
$mw->configure (title => "Geneace Update Accelerator   by Chao-Kung Chen   Version 1.0 (2003/01)",
                background => "white",
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

my $input_frame = $mw ->Frame(relief => 'groove', borderwidth => 2)
                    ->pack(after => $btn_frame, side => 'top', anchor => 'n', expand => 1, fill => 'x');


my $txt_frame = $mw->Frame(relief => 'groove', borderwidth => 2,
                           width => 735, height => 750,
			   #label => 'Genetics update'
			   )
                    ->pack(#in => $btn_frame,
			   after => $input_frame,
			   side => 'top',
			   anchor => 'n',
			   expand => 1, fill => 'x'
			  );


#my @menus;
#foreach (qw/File Help/){
 # push(@menus, $frame->Menubutton(text => $_));
#}

my $menu_File=$btn_frame->Menubutton(text => 'File', relief => 'groove', borderwidth => 2)->pack (side => 'left', anchor => 'n', fill => 'x');
$menu_File->AddItems(
		     ["command" => "Open Update file", "accelerator" => "Ctrl-o", command => \&open_update_file],
		     ["command" => "Save Update file", command => \&save_update_file],
                     ["command" => "Open ace file", "accelerator" => "Ctrl-a", command => \&open_ace_file],
		     ["command" => "Save ace file", command => \&save_ace_file],
		     "-",
		     ["command" => "Load ace file to Test_DB", "accelerator" => "Ctrl-t", command => \&upload_ace_test],
		     ["command" => "Load ace file to Geneace", "accelerator" => "Ctrl-g", command => \&upload_ace_GA],
                     "-",
                     ["command" => "View data format this script parses", "accelerator" => "Ctrl-r", command => \&show_format],
                     "-",
                     ["command" => "Quit", "accelerator" => "Ctrl-q", command => sub{exit}]
                    );

$mw->bind('<Control-Key-o>' => \&open_update_file);
$mw->bind('<Control-Key-a>' => \&open_ace_file);
$mw->bind('<Control-Key-s>' => \&save_ace_file);
$mw->bind('<Control-Key-t>' => \&upload_ace_test);
$mw->bind('<Control-Key-g>' => \&upload_ace_GA);
$mw->bind('<Control-Key-r>' => \&show_format);
$mw->bind('<Control-Key-q>' => sub{exit});

my $clear_txt = $btn_frame->Menubutton(text => 'Window', relief => 'groove', borderwidth => 2)->pack (side => 'left', anchor => 'n', fill => 'x');
$clear_txt->AddItems(["command" => "Clear upper window", "accelerator" => "Ctrl-u", command => \&clear_up],
                       ["command" => "Clear lower window", "accelerator" => "Ctrl-l", command => \&clear_down],
                      );

$mw->bind('<Control-Key-u>' => \&clear_up);
$mw->bind('<Control-Key-l>' => \&clear_down);


#my $menu_GA=$btn_frame->Menubutton(text => 'Correct Geneace_check', relief => 'groove')->pack (side => 'left', anchor => 'n', fill => 'x');
#$menu_GA->AddItems(
#		     ["command" => "Add loci to Gene_class obj", command => \&add_loci_to_geneclass],
#		     ["command" => "Add location to Allele obj", command => \&add_location_to_allele]
#		  );


my $menu_Option=$btn_frame->Menubutton(text => 'Options', relief => 'groove')->pack (side => 'left', anchor => 'n', fill => 'x');
$menu_Option->AddItems(
		       ["command" => "Gene_class/Gene/Other_name assignment -> .ace", command => \&geneclass_loci_other_name],
		       ["command" => "Genetics mapping -> .ace", command => \&gene_mapping],
                       ["command" => "Lab / Allele / PI -> .ace", command => \&update_lab],
		      );

my $clear = $btn_frame->Photo(-file =>"/wormsrv1/chaokung/DOCS/empty.gif");
my $btn_ClearBench=$btn_frame->Button(text    => 'Clear_Bench',
				      image   => $clear,
				      command => \&clear_up_down_window
				     )
                             ->pack (side => 'left', anchor => 'n', fill => 'x');

my($filename, $info, $btn_save, $mesg);

$input_frame->Label(text => 'Filename')->pack(side => 'left', anchor => 'w');
$input_frame->Entry(textvariable => \$filename)->pack(side => 'left', anchor => 'w', fill => 'x', expand => 1);

my $open_update = $btn_frame->Photo(-file =>"/wormsrv1/chaokung/DOCS/open_update.gif");
$btn_frame->Button(text => 'Open Updt', 
		   image => $open_update, command => \&type_update_file)->pack(side => 'left', anchor => 'w');

$btn_frame->Button(text => 'Save Updt', command => \&save_update_file)->pack(side => 'left', anchor => 'w');

my $open_ace = $btn_frame->Photo(-file =>"/wormsrv1/chaokung/DOCS/open_ace.gif");
$btn_frame->Button(text => 'Open ace',
		   image => $open_ace, command => \&type_ace_file)->pack(side => 'left', anchor => 'w');

$btn_frame->Button(text => 'Save ace', command => \&save_ace_file)->pack(side => 'left', anchor => 'w');


my ($input_window_lbl, $input_window, $ace_window_lbl, $ace_window, $scroll, $input, @content, $update);

&create_windows_and_lbl;

my $gaobj = init Geneace();
my %Gene_info = $gaobj -> gene_info();


#open (IN, $filename) || die "Can't read in file!";
#chomp (@content = <IN>);
#&geneclass_loci_other_name;

MainLoop();

#############
# subroutines
#############

sub clear_up {
  $input_window->delete('1.0', 'end')
}

sub clear_down {
  $ace_window->delete('1.0', 'end')
}

sub clear_up_down_window {
   $input_window->delete('1.0', 'end');
   $ace_window->delete('1.0', 'end')
}

sub show_format {
  $filename = "/wormsrv1/chaokung/JAH/JAH_data_format";

  $input_window->delete ('1.0', 'end');
  if (!open(IN, "$filename")){
    $input_window->insert('end', "ERROR: $!\n");
    return;
  }
  while (<IN>){
    $input_window->insert('end', $_);
  }
  close IN;
  return;
}

sub open_update_file {
  $filename=$btn_frame->getOpenFile(
			  -filetypes        =>
			  [
			     ['All Files',      '*'              ],
			     ['ace Files',      '.ace'           ],
			     ['Text Files',     ['.txt', '.text']],
			     ['Perl Scripts',   '.pl'            ],
			  ],
			 # -initialdir       => Cwd::cwd(),
			  -initialdir       => chdir "/wormsrv1/chaokung/wormdata/ACE/",
			  #-initialfile      => "getopenfile",
			  -title            => "Open update file",
			 ),

  #$info = "Loading file $filename";

  $input_window->delete ('1.0', 'end');
  if (!open(IN, "$filename")){
    $input_window->insert('end', "ERROR: $!\n");
    return;
  }
  while (<IN>){
    $input_window->insert('end', $_);
  }
  close IN;
  #$info = "File $filename loaded";
  return;
}


sub type_update_file {
  #$info = "Loading file $filename";

  $input_window->delete ('1.0', 'end');
  if (!open(IN, "$filename")){
    $input_window->insert('end', "ERROR: $!\n");
    return;
  }
  while (<IN>){
    $input_window->insert('end', $_);
  }
  close IN;
  #$info = "File $filename loaded";

  return;
}

sub save_update_file {
  $info="Saving $filename";
  open(IN, ">$filename");
  print IN $input_window->get('1.0', 'end');
  $info = "File saved";
}

sub open_ace_file {
  my $init_dir = chdir "/wormsrv1/chaokung/DOCS"; 
  $filename=$btn_frame->getOpenFile(-defaultextension =>  ".ace",
			  -filetypes        =>
			  [
			   ['ace Files',        '.ace'],
			   ['All Files',        '*'],
                          ],
			
			  -initialdir       => $init_dir,
			  #-initialfile      => "getopenfile",
			  -title            => "Open ace file",
			 ),

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

sub type_ace_file {

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

sub save_ace_file {
  $info="Saving $filename";
  open(IN, ">$filename");
  print IN $ace_window->get('1.0', 'end');
  $info = "File saved";
}

sub upload_ace_test {

  if (!$filename){
    my $command="pparse $filename\nsave\nquit\n";
    my $test_db_dir="/nfs/disk100/wormpub/DATABASES/TEST_DBs/CK1TEST/";

    open (Load_testGA,"| tace $test_db_dir ") || die "Failed to upload to test_Geneace";
    print Load_testGA $command;
    close Load_testGA;
    print "Done\n";
  }
  else {
    print "Save your ace file first before upload\n";
  }
}

sub upload_ace_GA {
  my $command="pparse $filename\nsave\nquit\n";
  my $db_dir="/nfs/disk100/wormpub/DATABASES/geneace/";
  open (Load_GA,"| tace $db_dir ") || die "Failed to upload to test_Geneace";
  print Load_GA $command;
  close Load_GA;
  print "Done\n";
}

sub geneclass_loci_other_name {
  open (IN, $filename) || die "Can't read in file!";

  my ($gene_class, @Update, $locus, $seq, $rest, @parts, $num_parts,$cgc_paper, $paper,
      $head, $tail, @variants, $i, $person, $pmid, $other_name, @evidence, $evidence, @persons);

  while(<IN>){
    if ($_ =~ ""){}
    #if ($_ =~ /^(\w{3,3})\s{2,}(\w+)\s{2,}(.+)/) {
    if ($_ =~ /^(\w{3,3})\s+\[phenotype: (.+)\]/ ) {
      push(@Update,"\n\nGene_class : \"$1\"\n");
      push(@Update,"CGC_approved\n");
      push(@Update,"Phenotype \"$2\"\n");
    }

    if ($_ =~ /^(\w{3,3})\s+([A-Z]+)\s+(.+)/ && $_ !~ /^New|NEW/ ) {
      $gene_class = $1;
      push(@Update,"\n\nGene_class : \"$gene_class\"\n");
      #@parts = split(/\s{2,}/, $2);
      #print "Description\t\"$parts[1]\"\n";
      push(@Update,"Description\t\"$3\"\n");
      push(@Update,"Designating_laboratory\t\"$2\"\n");
      #print "Designating_laboratory\t\"$parts[0]\"\n";
    }

    if ($_ =~ /^(\w{3,3})\s+see\s+(.+)/ && $_ !~ /^New|NEW/ ) {
      $gene_class = $1;
      push(@Update,"\n\nGene_Class : \"$1\"\n"); 
      push(@Update,"Description\t\"See $2\"\n");
    }

    if ($_ =~ /^((\w{3,4})-\d+)\s+=\s+(\w+\..+)/) {
      push(@Update, "\n\nGene : \"$Gene_info{$1}{'Gene'}\"\n");
      push(@Update, "Public_name\t\"$1\"\n");
      push(@Update, "CGC_name\t\"$1\"\n");
      push(@Update, "Gene_class\t\"$2\"\n");
      push(@Update, "Species\t\"Caenorhabditis elegans\"\n");
      push(@Update, "Live\n");
      push(@Update, "Version  1\n");
      push(@Update, "Method Gene\n");
      push(@Update, "Version_change   1 now \"WBPerson1971\" Created\n");

      $locus = $1;
      $rest = $3;
      @parts = split(/\s{2,}/, $rest);
      $seq = $parts[0];
      #print "$seq\n";
      $num_parts = scalar @parts;
      if ($num_parts == 1){
	$seq = uc($seq);
	push(@Update, "Sequence_name\t\"$seq\"\n");
	print "Sequence_name\t\"$seq\"\n";  
      }
      if ($num_parts > 1){
	$head = $seq;
	$tail = $seq;
	$head =~ /[\w\d]+\.\d+/; # get seq. names like C18H9.7 or M01E5.5 from M01E5.5a,b
	$head = $&;	$tail =~ s/$head//;  # would retrieve a,b from M01E5.5a,b
	$head = uc($head);
	$seq=$head.$tail;
	#print $head, "##1##\n";
	#print $seq, "\##2##\n";
	if($tail ne ""){
	  @variants=split(/,/, $tail);
	  #print "@variants\n";
	  foreach (@variants){
	    $seq = $head.$_;
	    for ($i = 1; $i < $num_parts; $i++){
	      #print $parts[$i],"\n";
	      if ($parts[$i] =~ /(Paper_evidence)\s(.+)/){
		#print $1, "\n";
		$cgc_paper=$2;
		$cgc_paper =~ s/\[|\]//g;
		push(@Update, "Sequence_name\t\"$seq\"\n");
		push(@Update, "Evidence\tPaper_evidence\t\"WBPaper\"\n");	      
	      }
	      if ($parts[$i] =~ /(Person_evidence)\s(.+)/){
		
		$evidence = $1;
		$person = $2;
		$person =~ s/\[|\]//g;
		@persons = split(/,/, $person);
		foreach (@persons){
		  $_ =~ s/^\s//;
		  push(@Update, "Sequence_name\t\"$seq\"\n");   
		  push(@Update, "Evidence\tAuthor_evidence\t\"$_\"\n");
                }
	      }
	    }
	  }
	}
	if ($tail eq "") {
	  for ($i = 1; $i < $num_parts; $i++){
	    if ($parts[$i] =~ /(Paper_evidence)\s(.+)/){
	      #print $1, "\n";
	      $cgc_paper=$2;
	      $cgc_paper =~ s/\[|\]//g;
	      push(@Update, "Sequence_name\t\"$seq\"\n");
	      push(@Update, "Evidence\tPaper_evidence\t\"WBPaper\"\n");
	      
	    }
	    if ($parts[$i] =~ /Person_evidence\s(.+)/){
	      #push(@Update, $2);
	      $person = $1;
	      $person =~ s/\[|\]//g;
	      @persons = split(/,/, $person);
	      foreach (@persons){
		$_ =~ s/^\s//;
		push(@Update, "Sequence_name\t\"$seq\"\n");
		push(@Update, "Evidence\tAuthor_evidence\t\"$_\"\n");
              }
	    }
	  }
	}
      }
    }
    if($_ =~ /(\w{3,4}-\d+)\s+\(common name\)\s+=\s+(\w{3,4}-\d+)\s+\(other name\)\s+(.+)/){
      my $main = $1;
      $other_name = $2;
      print "$main : $other_name : $3\n";
      my $evidence = $3;
      push(@Update, "\nLocus : \"$2\"\n");
      push(@Update, "\n-R Locus : \"$2\" \"$1\"\n");
      push(@Update, "\nLocus : \"$1\"\n");
      #@evidence = split(/\s+/, $evidence);

      if ($evidence =~ /Paper_evidence\s+(.+)/){
	#print $evidence[1], "\n";
	$paper = $1;
	$paper =~ s/\[|://g;
	$paper =~ s/\[|\]//g;
	push(@Update, "Other_name\t\"$other_name\"\tPaper_evidence\t\"WBPaper\"\n");
	
      }
      if ($evidence =~ /Person_evidence\s+(.+)/){
	
	$person = $1;
	$person =~ s/\[|\]//g;
	@persons = split(/,/, $person);
	foreach (@persons){
	  $person = $_;
	  $person =~ s/^\s//;
	  push(@Update, "Other_name\t\"$other_name\"\tAuthor_evidence\t\"$person\"\n");
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
  print $filename, "\n";
  #open (IN, "/wormsrv1/chaokung/wormdata/JAH/$filename") || die "Can't read in file!";  
  my (%obj_info, $info, $obj, $remark, $combined, $mapper, @mappers, $date,
      @update, $seq, @yr_mon, $locus, $provider, @providers, $person, $temp,
      @clones, $clone, $allele, @alleles, @distances, $A_non_B, $B_non_A, $ga,
      @gene_allele, @genes, @tested, $gene_class, $rearr);

  ############################################
  # fetch the numbers of the last obj of 
  # ?multi_pt_data ?2_point_data ?pos_neg_data
  ############################################

  my ($multi_pt_count, $two_pt_count, $pos_neg_count);
  
  my $dialog =  $mw -> DialogBox(-title   => "Get next number...",
				 -buttons => ["Close"]);
  
  $dialog->geometry("220x260");
  $dialog->add('Label',
	       -anchor => 'n',
	       -justify => 'left',
	       -text => "Enter next number of\npos-neg-data\nmulti-pt-data\n2-pt-data objects!\n")
         ->pack();

  $dialog->Label(text => 'Last multi-point #')->pack(side => 'top', anchor => 'w');
  $dialog->Entry(textvariable => \$multi_pt_count)->pack(side => 'top', anchor => 'w', fill => 'x', expand => 1);

  $dialog->Label(text => 'Last 2-point #')->pack(side => 'top', anchor => 'w');
  $dialog->Entry(textvariable => \$two_pt_count)->pack(side => 'top', anchor => 'w', fill => 'x', expand => 1);

  $dialog->Label(text => 'Last pos-neg-data #')->pack(side => 'top', anchor => 'w');
  $dialog->Entry(textvariable => \$pos_neg_count)->pack(side => 'top', anchor => 'w', fill => 'x', expand => 1);

  $dialog->Show();

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
    if ($_ =~ /^Multi[_-]pt_[cC]omment(:|\s:)\s(.+)/){
      $remark = $2; $remark =~ s/^\s//; #$remark .= " [$date ck1]";
      write_ace($multi_pt_count, "Remark", "$remark\"\t\"CGC_data_submission");
    }
    if ($_ =~ /^Multi[_-]pt_[cC]ombined_results(:|\s:)\s+(.+)/){
      @distances = split(/\s+/, $2);
      $combined = "Gene \"$Gene_info{$distances[0]}{'Gene'}\" $distances[1] Gene \"$Gene_info{$distances[2]}{'Gene'}\" $distances[3] Gene \"$Gene_info{$distances[4]}{'Gene'}\"";
      write_ace($multi_pt_count, "Combined", $combined)}
    if ($_ =~ /^Multi[_-]pt_[gG]enotype(:|\s:)\s+(.+)/){write_ace($multi_pt_count, "Genotype", $2)}
    if ( ($_ =~ /^Multi[_-]pt_[mM]apper(:|\s:)\s+(.+)/) || ($_ =~ /Multi_pt_Data_provider(:|\s:)\s+(.+)/) ){
      @mappers = split(/,/, $2);
      foreach $mapper (@mappers){
	$mapper =~ s/^\s//;
	write_ace($multi_pt_count, "Mapper", $mapper);
      }
    }
    if ($_ =~ /^Multi[_-]pt_[dD]ate(:|\s:)\s+(.+)/){
      @yr_mon = split(/\//, $2);
      write_ace($multi_pt_count, "Date", "$yr_mon[1]-$yr_mon[0]");
    }
    if ($_ =~ /^Multi[_-]pt_GeneA(:|\s:)\s+(.+)/){write_ace($multi_pt_count, "Gene_A", $2)}
    if ($_ =~ /^Multi[_-]pt_GeneB(:|\s:)\s+(.+)/){write_ace($multi_pt_count, "Gene_B", $2)}

    if ($_ =~ /^Multi[_-]pt_A_non_B_results(:|\s:)\s+(.+)/){
      @distances = split(/\s+/, $2);
      $A_non_B = "Gene \"$Gene_info{$distances[0]}{'Gene'}\" $distances[1] Gene \"$Gene_info{$distances[2]}{'Gene'}\" $distances[3] Gene \"$Gene_info{$distances[4]}{'Gene'}\"";
      write_ace($multi_pt_count, "A_non_B", $A_non_B);
    }
    if ($_ =~ /^Multi[_-]pt_B_non_A_results(:|\s:)\s+(.+)/){
      @distances = split(/\s+/, $2);
      $B_non_A = "Genes \"$Gene_info{$distances[0]}{'Gene'}\" $distances[1] Gene \"$Gene_info{$distances[2]}{'Gene'}\" $distances[3] Gene \"$Gene_info{$distances[4]}{'Gene'}\"";
      write_ace($multi_pt_count, "B_non_A", $B_non_A);
    }

    #################
    # parse 2_pt_data
    #################
 
    if ($_ =~ /^\/\/2_[pP]oint.+/){
      $two_pt_count++;
      write_ace($two_pt_count, "\n\n2_point_data : ", $two_pt_count);
    }
    if ($_ =~ /^2_point_[tT]emp(:|\s:)\s+(.+)/){$temp = $2; $temp =~ s/°//; write_ace($two_pt_count, "Temperature", $temp)}
    if ($_ =~ /^2_point_[mM]apped_genes(:|\s:)\s+(.+)/){
      @genes = split(/,/, $2);
      write_ace($two_pt_count, "Gene_1", $Gene_info{$genes[0]}{'Gene'});
      write_ace($two_pt_count, "Gene_2", $Gene_info{$genes[1]}{'Gene'});
    }
    if ($_ =~ /^2_point_[gG]enotype(:|\s:)\s+(.+)/){write_ace($two_pt_count, "Genotype", $2)}
    if ($_ =~ /^2_point_[rR]esults(:|\s:)\s+(.+)/){write_ace($two_pt_count, "Results", $2)}
    if ($_ =~ /^2_point_[mM]apper(:|\s:)(.+)/){
      @mappers = split(/,/, $2);
      foreach (@mappers){write_ace($two_pt_count, "Mapper", $_)}
    }
    if ($_ =~ /^2_point_[dD]ate(:|\s:)\s+(.+)/){
      @yr_mon = split(/\//, $2); write_ace($two_pt_count, "Date", "$yr_mon[1]-$yr_mon[0]");
    }
    if ($_ =~ /^2_point_[cC]alculation(:|\s:)\s+(.+)/){
      @tested = split(/ /, $2);
      write_ace($two_pt_count, "Tested", "$tested[1] $tested[2]");
    }
    if ($_ =~ /^2_point_[cC]omment(:|\s:)\s+(.+)/){
      $remark = $2; $remark =~ s/^\s//; #$remark .= " [$date ck1]";
      write_ace($two_pt_count, "Remark", "$remark\"\t\"CGC_data_submission");
    }

    ####################
    # parse pos_neg_data
    ####################

    if ($_ =~ /^\/\/Df[\/|_]Dp.+/){
      $pos_neg_count++; write_ace($pos_neg_count, "\n\nPos_neg_data : ", $pos_neg_count);
    }
    if ($_ =~ /^Df[\/|_]Dp_[rR]esults(:|\s:)\s+(.+)/){write_ace($pos_neg_count, "Results", $2)}
    if ($_ =~ /^Df[\/|_]Dp_[cC]omment(:|\s:)\s+(.+)/){
      $remark = $2; $remark =~ s/^\s//; #$remark .= " [$date ck1]";
      write_ace($pos_neg_count, "Remark", "$remark\"\t\"CGC_data_submission");
    }
    if ($_ =~ /^Df[\/|_]Dp_[rR]earrangement(:|\s:)\s+(.+)/){write_ace($pos_neg_count, "Rearrangement_2", $2)}
    if ($_ =~ /^Df[\/|_]Dp_[gG]enotype(:|\s:)\s+(.+)/){write_ace($pos_neg_count, "Genotype", $2)}
    if ( ($_ =~ /^Df[\/|_]Dp_[dD]ata_[pP]rovider(:|\s:)\s+(.+)/) || ($_ =~ /^Df[\/|_]Dp_[mM]apper(:|\s:)\s+(.+)/) ) {
      @providers = split(/,/, $2);
      foreach (@providers){write_ace($pos_neg_count, "Mapper", $_);}
    }
    if ($_ =~ /^Df[\/|_]Dp_[dD]ate(:|\s:)\s+(.+)/){
      @yr_mon = split(/\//, $2); write_ace($pos_neg_count, "Date", "$yr_mon[1]-$yr_mon[0]");
    }
    if ($_ =~ /^Df[\/|_]Dp_[gG]ene(:|\s:)\s+(.+)/){
      $ga = $2; $ga =~ s/\)$//;
      @gene_allele = split(/\(/, $ga);
      if ($gene_allele[1] ne ""){
	write_ace($pos_neg_count, "Gene_1", "$Gene_info{$gene_allele[0]}{'Gene'}", "$gene_allele[1]"); # w/ allel info
      }
      else {write_ace($pos_neg_count, "Gene_1", "$Gene_info{$gene_allele[0]}{'Gene'}")}; # w/o allele info
    }


    ##################
    # parse Gene data
    ##################

    if ( ($_ =~ /^Locus_[lL]ocus_[nN]ame(:|\s:)\s+(.+)/) ||
	 ($_ =~ /^Locus_[gG]ene_[nN]ame(:|\s:)\s+(.+)/) ){
      $locus = $2;
      $gene_class = $2;
      $gene_class =~ s/-.+//;
      write_ace($locus, "\n\nGene : ", $Gene_info{$locus}{'Gene'});
      write_ace($locus, "Gene_class", $gene_class);
    }
    if ($_ =~ /^Locus_[sS]equence(:|\s:)\s+(.+)/){$seq = uc($2); write_ace($locus, "CDS", $seq)}
    if ($_ =~ /^Locus_[nN]egative_[mM]ethod(:|\s:)\s+(.+)/){
      $remark = $2;
      $remark = "Negative method is ".$remark;
      write_ace($locus, "Remark", "$remark\"\t\"CGC_data_submission");
    }
    if ($_ =~ /^Locus_[pP]ositive_[mM]ethod(:|\s:)\s+(.+)/){
      $remark = $2;
      $remark = "Postitive method is ".$remark;
      write_ace($locus, "Remark", "$remark\"\t\"CGC_data_submission");
    }
    if ($_ =~ /^Locus_[cC]omment(:|\s:)\s+(.+)/){
      $remark = $2; $remark =~ s/^\s//; #$remark .= " [$date ck1]";
      write_ace($locus, "Remark", "$remark\"\t\"CGC_data_submission");
    }
    if ($_ =~ /^Locus_[nN]egative_[cC]lone(:|\s:)\s+(.+)/){
      @clones = split(/,|and/, $2);
      foreach $clone (@clones){
	$clone =~ s/^\s//;
	write_ace($locus, "Negative_clone", $clone);
      }
    }
    if ($_ =~ /^Locus_[pP]ositive_[cC]lone(:|\s:)\s+(.+)/){
      @clones = split(/,|and/, $2);
      foreach $clone (@clones){
	$clone =~ s/^\s//;
	write_ace($locus, "Positive_clone", $clone);
      }
    }

    if ($_ =~ /^Locus_(d|D)ata_provider(:|\s:)\s+(.+)/){
      @providers = split(/,|and/, $3);
    }


    if ($_ =~ /^Locus_[cC]hromosome(:|\s:)\s+(.+)/){write_ace($locus, "Map", $2)}

    if ($_ =~ /^Locus_[aA]llele(:|\s:)\s+(.+)/){
      @alleles = split(/,/, $2);
      foreach (@alleles){
        $allele = $_;
	$allele =~ s/^\s|\s$//;
	write_ace($locus, "Allele", $allele);
      }
    }

    ###############################
    # parse breakpt (rearrangement)
    ###############################

    if ($_ =~ /^Breakpt_Rearrangement(:|\s:)\s+(.+)/){$rearr = $2; write_ace($2, "\n\nRearrangement : ", $2)}
    if ($_ =~ /^Breakpt_Comment(:|\s:)\s+(.+)/){
      $remark = $2; $remark =~ s/^\s//; #$remark .= " [$date ck1]";
      write_ace($rearr, "Remark", "$remark\"\tCGC_data_submission");
    }
    if ($_ =~ /^Breakpt_Negative_clone(:|\s:)\s+(.+)/){
      @clones = split(/,/, $2);
      foreach $clone (@clones){
	$clone =~ s/^\s//;
	write_ace($rearr, "Clone_outside", $clone);
      }
    }
    if ($_ =~ /^Breakpt_Positive_clone(:|\s:)\s+(.+)/){
      @clones = split(/,/, $2);
      foreach $clone (@clones){
	$clone =~ s/^\s//;
	write_ace($rearr, "Clone_inside", $clone);
      }
    }
    if ($_ =~ /^Breakpt_Data_provider(:|\s:)\s+(.+)/){
      @providers = split(/,/, $2);
      foreach $person (@providers){
	$person =~ s/^\s//;
	write_ace($rearr, "Author", $person);
      }
    }	

    ##################
    # parse gene_class
    ##################

    if ($_ =~ /GeneClass_[gG]ene_name(:|\s:)\s+(.+)/){
      $gene_class = $1; 
      write_ace($1, "\n\nGene_class : ", $1);
      write_ace($1, "CGC_approved");
    }
    if ($_ =~ /GeneClass_[pP]henotype(:|\s:)\s+(.+)/){write_ace($gene_class, "Phenotype", $2)}
    if ($_ =~ /GeneClass_[lL]aboratory(:|\s:)\s+(.+)/){write_ace($gene_class, "Designating_laboratory", $2)}
    if ($_ =~ /GeneClass_[cC]omment(:|\s:)\s+(.+)/){
      $remark = $2; $remark =~ s/^\s//; #$remark .= " [$date ck1]";
      write_ace($gene_class, "Remark", "$remark\"\tCGC_data_submission");
    }
  }

  sub update_lab {
    open (IN, $filename) || die "Can't read in file!";

    my %person_id = $gaobj->get_WBPersonID();

    my ($gene_class, @Update, $locus, $seq, $rest, @parts, $num_parts,$cgc_paper, $paper,
	$head, $tail, @variants, $i, $person, $pmid, $other_name, @evidence, $evidence, @persons);

    while(<IN>){
      if ($_ =~ ""){}

      # 2-letter lab
      if ( $_ =~ /^([A-Z]{2,2})\s+(\w{2,2}|--|\w{1,1})\s+(\w{1,} \w{1,},\s\w{1,}|\w{1,}-\w{1,},\s\w{1,}|\w{1,}-\w{1,},\s\w{1,}-\w{1,}|\w{1,},\s\w{1,}|\w{1,},\s\w{1,}-\w{1,})\s+(.+)$/ ) { 

     #  if ( $_ =~ /^([A-Z]{2,2})\s+.+/ ) { 	
	push(@Update,"\n\nLaboratory : \"$1\"\n"); 
	push(@Update,"Allele_designation \"$2\"\n");
	push(@Update,"Representative \"$person_id{$3}\"\n") if $person_id{$3};
	push(@Update,"Mail\t\"$4\"\n");
      }
      if ( $_ =~ /^([A-Z]{2,2})\s+(\w{2,2}|--)\s+\[(.+)\]\s+\[(.+)\]/){
	push(@Update,"\n\nLaboratory : \"$1\"\n"); 
	push(@Update,"Allele_designation \"$2\"\n");
	push(@Update,"Representative \"$person_id{$3}\"\n") if $person_id{$3};
	push(@Update,"Mail\t\"$4\"\n");
      }
      # 3-letter lab 
      if ( $_ =~ /^([A-Z]{3,3})\s+(\w{1,} \w{1,},\s\w{1,}|\w{1,}-\w{1,},\s\w{1,}|\w{1,}-\w{1,},\s\w{1,}-\w{1,}|\w{1,},\s\w{1,}|\w{1,},\s\w{1,}-\w{1,})\s+(.+)$/ ) { 
	push(@Update,"\n\nLaboratory : \"$1\"\n");
	push(@Update,"Representative \"$person_id{$2}\"\n") if $person_id{$2};
	push(@Update,"Mail\t\"$3\"\n");
      }
    }
    foreach (@Update){
      print $_,;
      $ace_window -> insert('end',"$_");
    }
  }

  sub write_ace {
    my($obj, $tag, $value)=@_;
    my(%obj_info, $info, $gene);
    if ( ($tag eq "A_non_B") || ($tag eq "B_non_A") || ($tag eq "Combined") || ($tag eq "Tested") ){
      push(@{$obj_info{$obj}},  "$tag\t$value");
    }
    elsif ($tag eq "CGC_approved"){
      push(@{$obj_info{$obj}}, "$tag");
    }

    else {
      push(@{$obj_info{$obj}}, "$tag\t\"$value\""); # quotation different
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

}
}




__END__



