#! /usr/bin/perl5.6.1 -w

# get_allele_flanking_seqs_TK.pl

# by Chao-Kung Chen [030625]

# Last updated on: $Date: 2003-09-15 15:44:50 $
# Last updated by: $Author: ck1 $

use Tk;
use strict;
use Cwd;
use Term::ANSIColor;
use lib "/wormsrv2/scripts";
use Wormbase;
use TextANSIColor;
use Tk::DialogBox;


#######################
# Create a main windows
#######################

#----------- top level frame ----------

my $mw=MainWindow->new();
$mw->configure (title => "Allele Flanking Sequences Retriever   by Chao-Kung Chen   Version 1.0 (2003/06/20)",
                background => "white",
               );

$mw->geometry("750x900+0+0");
#$mw->resizable(0,0);

######################################
# Instantiate widgets and arrange them
######################################

#----------- menu frame ----------


my $menu_frame = $mw ->Frame(relief => 'groove', borderwidth => 2)
                    ->pack(side => 'top', anchor => 'n', expand => 1, fill => 'x');

my $menu_File=$menu_frame->Menubutton(text => 'File')->pack (side => 'left', anchor => 'w');
$menu_File->AddItems(
		     ["command" => "Open ace file", "accelerator" => "Ctrl-o", command => \&open_ace_file],
		     ["command" => "Save ace file", "accelerator" => "Ctrl-s", command => \&save_ace_file],
		     "-",
		     ["command" => "Load ace file to TEST_DB", "accelerator" => "Ctrl-t", command => \&upload_ace_TEST],
		     ["command" => "Load ace file to Geneace", "accelerator" => "Ctrl-g", command => \&upload_ace_GA],
                     "-",
                     ["command" => "Quit", "accelerator" => "Ctrl-q", command => sub{exit}]
                    );

$mw->bind('<Control-Key-o>' => \&open_ace_file);
$mw->bind('<Control-Key-s>' => \&save_ace_file);
$mw->bind('<Control-Key-t>' => \&upload_ace_TEST);
$mw->bind('<Control-Key-g>' => \&upload_ace_GA);
$mw->bind('<Control-Key-q>' => sub{exit});

my $menu_Window=$menu_frame->Menubutton(text => 'Window')->pack (side => 'left', anchor => 'w');
$menu_Window->AddItems(["command" => "Clear upper window", "accelerator" => "Ctrl-u", command => \&clear_up],
                       ["command" => "Clear lower window", "accelerator" => "Ctrl-l", command => \&clear_down],
                      );

$mw->bind('<Control-Key-u>' => \&clear_up);
$mw->bind('<Control-Key-l>' => \&clear_down);


my $menu_Help=$menu_frame->Menubutton(text => 'Help')->pack (side => 'left', anchor => 'w');
$menu_Help->AddItems(["command" => "About allele flanking sequence retriever", "accelerator" => "Ctrl-h", command => \&read_doc]);

$mw->bind('<Control-Key-h>' => \&read_doc);

#----------- parameter example label ----------

my $param_frame = $mw ->Frame(relief => 'groove', borderwidth => 2)
                    ->pack(side => 'top', anchor => 'n', after => $menu_frame, expand => 1, fill => 'x');

$param_frame -> Label(text=>"Query parameters: ", fg=>"blue") # create a horizontal space 
           -> pack(side => "left");

$param_frame -> Label(text => "Eg: 4R79.1 -aa 324y ok12 abc-1 pmid/[cgcxxx]   OR   4R79.1 -dna 1324t ok12 abc-1 pmid/[cgcxxx]", fg => "black")
                -> pack(side => "left");

#----------- entry box ----------

my $params;
my @info;

my $entry_frame = $mw ->Frame(relief => 'groove', borderwidth => 2)
                    ->pack(side => 'top', anchor => 'n', after => $param_frame, expand => 1, fill => 'x');

my $entry = $entry_frame->Entry(textvariable => \$params, bg => "white", fg => "black", width => 60)->pack(side =>"left");

$entry_frame->Button(text => "RUN", command => \&run)
          -> pack(side => "left", fill => "x", expand => 1);

$entry_frame->Button(text => "Reset", command => \&Reset)
          -> pack(side => "left", fill => "x", expand => 1);


#----------- button frame ----------

my $btn_frame = $mw ->Frame(relief => 'groove', borderwidth => 2)
                    ->pack(side => 'top', anchor => 'n', after => $entry_frame, expand => 1, fill => 'x');

my $label_1 = $btn_frame->Label(text => ""); # not displayed
my $label_2 = $btn_frame->Label(text => ""); # not displayed
my $label_3 = $btn_frame->Label(text => ""); # not displayed
my $label_4 = $btn_frame->Label(text => ""); # not displayed


my $btn_1 = $btn_frame->Button(text => "1st site mutation", command => \&get_1_site_flanks)
                      -> pack(side => "left",  fill => "x", expand => 1);

my $btn_2 = $btn_frame->Button(text => "2nd site mutation", command => \&get_2_site_flanks)
                      -> pack(side => "left", fill => "x", expand => 1);

my $btn_3 = $btn_frame->Button(text => "3rd site mutation", command => \&get_3_site_flanks)
                      -> pack(side => "left", fill => "x", expand => 1);

my $btn_4 = $btn_frame->Button(text => "Flanking codon", command => \&get_codon_flanks)
                      -> pack(side => "left", fill => "x", expand => 1);


#----------- color coding message frame ----------

my $msg_frame_1 = $mw ->Frame(height => 30)
                    ->pack(side => 'top', anchor => 'n', after => $btn_frame, expand => 1, fill => 'x');

$msg_frame_1 -> Label(text => "Codon triplet color coding: mutated site in")
             -> pack(side => "left", anchor => "n");
$msg_frame_1 -> Label(text => "RED", fg => "red")
             -> pack(side => "left", anchor => "n");
$msg_frame_1 -> Label(text => ", the other two sites behind or before it in the triplet in")
             -> pack(side => "left", anchor => "n");
$msg_frame_1 -> Label(text => "BLUE", fg => "blue")
             -> pack(side => "left", anchor => "n");

my $msg_frame_2 = $mw ->Frame(height => 30)
                    ->pack(side => 'top', anchor => 'n', after => $msg_frame_1, expand => 1, fill => 'x');

$msg_frame_2 -> Label(text => "Flanking bp coding at frame shift:  upsteam of mutated site in")
             -> pack(side => "left");
$msg_frame_2 -> Label(text => "MAGENTA", fg => "magenta")
             -> pack(side => "left");
$msg_frame_2 -> Label(text => ", downstream in")
             -> pack(side => "left");
$msg_frame_2 -> Label(text => "GREEN", fg => "green")
             -> pack(side => "left");

my $msg_frame_3 = $mw ->Frame(height => 30)
                    ->pack(side => 'top', anchor => 'n', after => $msg_frame_2, expand => 1, fill => 'x');

$msg_frame_3 -> Label(text => "For DNA query, the mutated nt is in")
             -> pack(side => "left");
$msg_frame_3 -> Label(text => "RED", fg => "red")
             -> pack(side => "left");
$msg_frame_3 -> Label(text => ". This color coding does not apply to multiple mutation sites.")
             -> pack(side => "left");

#----------- upper window frame ----------

my $run_window = $mw ->TextANSIColor(height => 23)
                  ->pack(side => 'top', anchor => 'n', after => $msg_frame_3, expand => 1, fill => 'x');

#----------- lower window frame ----------

my $ace_frame = $mw->Frame(relief => 'groove', borderwidth => 2)
		   ->pack(after => $run_window, side => 'top', anchor => 'n', expand => 1, fill => 'x');

my $ace_window =$ace_frame -> Scrolled("Text",  -scrollbars => "ow", height => 400)
                           -> pack(side => 'top', anchor => 'n', fill => 'x');

#---------- using ANSIColor in PERL/TK ----------

my $red = color('red');  # Retrieve color codes
my $green = color('green');
my $magenta = color('magenta');
my $blue = color('blue');
my $bold = color('bold');
my $black = color('black');
my $reset = color('reset');


MainLoop();

#----------- END OF WIDGET LAYOUT ----------

##################
# global variables
##################

my ($cds_or_locus, $aa_or_dna, $mutation, $allele, $locus, $position, $cds, $current, @archive, $archive,
    @output, @ace, @out, @DNA, @prot, %code, @CDS_coords, $molecule, $filename);


######################################################
#              s u b r o u t i n e s
######################################################

sub clear_up {
  $run_window->delete('1.0', 'end')
}

sub clear_down {
  $ace_window->delete('1.0', 'end')
}

sub read_doc{
  my $dialog =  $mw -> DialogBox(-title   => "About Allele Flanking Sequences Retriever",
                                 -buttons => ["Close" ]);

  $dialog->geometry("940x400");
  $dialog->resizable(0,0);

  my @doc = `perldoc /wormsrv1/chaokung/my-scripts/$0`;
  my $doc;
  foreach(@doc){
    if ($_ eq "" ){}
    elsif ($_ =~ /=(.+)/){$doc .= $1}
    else {$doc .= $_}
  }
  my $txt=$dialog->Scrolled("Text",  -scrollbars=>"ow", height=>60, width=> 130)->pack(side => "left", anchor => "w");
  $txt -> insert('end', "USAGE:\n\nQuery parameters are in the order of\n\n(1) CDS/Transcript/Pseudogene\n(2) -aa (amino acid) or (-dna) nucleotide\n(3) nucleotide or amino acid position followed by mutation in single letter (or three-letter code for amino acid)\n(4) allele name \n(5) locus name \n\nseparated by space (all parameters are case insensitive)\n\n\n---------------------------------------------------------------------------------------------------------------------------------\n                                            Detailed description of this program\n---------------------------------------------------------------------------------------------------------------------------------");  
  $txt -> insert('end', "$doc");
  $dialog->Show();

}

sub open_ace_file{

  #$ace_window->delete ('1.0', 'end'); # maybe useful, leave here for the moment
  my $filetype = [["ace file", [".ace"]]];
  my $filename = $menu_frame->getOpenFile(-title =>"Select file to open",
			       -filetypes => $filetype,
			       -defaultextension => ".ace");

  open(IN, $filename);
  while(<IN>){
    $ace_window->insert('end', "$_");
  }
   close IN;
}

sub save_ace_file{

  my $filetype = [["ace file", [".ace"]]];
  my $filename = $menu_frame->getSaveFile(-title =>"Select file to save to",
			       -filetypes => $filetype,
			       -defaultextension => ".ace");

# maybe useful, leave here for the moment

#  my $fail_dialog= $mw->DialogBox(-title          => "Allele Flanking Sequences Retriever",
#                                  -default_button => "Exit",
#                                  -text           => "No filename selected!");

#  return $fail_dialog->Show() if !$filename;

  $run_window->insert('end', "$filename");
  my @out = $ace_window->get('1.0', 'end');
  if (open(SAVE, ">$filename")){
   print SAVE @out;
   close SAVE;
  }
}

sub upload_ace_TEST{

  my $log = "/wormsrv2/tmp/load_allele.log";
  system("chmod 777 $log");

  my @out = $ace_window->get('1.0', 'end');
  my $filename = "/wormsrv2/tmp/tmp.ace";
  system("chmod 777 $filename");

  open(SAVE, ">$filename");
  print SAVE @out;
  close SAVE;

  my $command=<<EOF;
pparse $filename
save
quit
EOF

  my $test_db_dir="/nfs/disk100/wormpub/DATABASES/TEST_DBs/CK1TEST/";
  open (Load_testGA,"| tace $test_db_dir > $log") || die "Failed to upload to test_Geneace";
  print Load_testGA $command;
  close Load_testGA;

  my $dialog2 =  $mw -> DialogBox(-title   => "Uploading ace file . . .",
                                 -buttons => ["Close" ]);

  $dialog2->geometry("800x500");
  #$dialog2->resizable(0,0);
  my $txt=$dialog2->Scrolled("Text",  -scrollbars=>"ow", height=>60, width=> 170)->pack(side => "left", anchor => "w");
  open(LOG, "$log");
  while(<LOG>){$txt -> insert('end', "$_")}
  $dialog2->Show();
}

sub upload_ace_GA{

  my $user = `whoami`;
  chomp $user;
  if ($user ne "wormpub"){
    my $dialog3 =  $mw -> DialogBox(-title   => "ERROR loading ace file . . .",
                                    -buttons => ["Close" ]);

    $dialog3->geometry("250x100");
    $dialog3->add('Label',
		-anchor => 'n',
		-justify => 'left',
		-text => "You have no write access to Geneace!")
            ->pack();
    $dialog3->Show();
  }
  if ($user eq "wormpub") {
    my $log = "/wormsrv2/tmp/load_allele.log";
    system("chmod 777 $log");
    
    my @out = $ace_window->get('1.0', 'end');
    my $filename = "/wormsrv2/tmp/tmp.ace";
    system("chmod 777 $filename");
    
    open(SAVE, ">$filename");
    print SAVE @out;
    close SAVE;

    my $command=<<EOF;
pparse $filename
save
quit
EOF

    my $db_dir="/wormsrv1/geneace/";
    open (Load_GA,"| tace $db_dir > $log") || die "Failed to upload to test_Geneace";
    print Load_GA $command;
    close Load_GA;
    
    my $dialog4 =  $mw -> DialogBox(-title   => "Uploading ace file to Geneace . . .",
				    -buttons => ["Close" ]);
    
    $dialog4->geometry("800x500");
    #$dialog4->resizable(0,0);
    my $txt=$dialog4->Scrolled("Text",  -scrollbars=>"ow", height=>60, width=> 170)->pack(side => "left", anchor => "w");
    open(LOG, "$log");
    while(<LOG>){$txt -> insert('end', "$_")}
    $dialog4->Show();
  }
}

sub Reset{
   $entry -> delete('0.1', 'end');
}

sub run {

  my $text = $entry -> cget(-textvariable);

  # look up genetic code of aa
  %code = (
	    'A' => ["gct", "gcc", "gca", "gcg"],
	    'B' => ["aat", "aac", "gat", "gac"],
	    'C' => ["tgt", "tgc"],
	    'D' => ["gat", "gac"], 
	    'E' => ["gaa", "gag"], 
	    'F' => ["ttt", "ttc"],
	    'G' => ["ggt", "ggc", "gga", "ggg"],
	    'H' => ["cat", "cac"],
	    'I' => ["att", "atc", "ata"],
	    'K' => ["aaa", "aag"],
	    'L' => ["tta", "ttg", "ctt", "ctc", "cta", "ctg"],
	    'M' => ["atg"],
	    'N' => ["aat", "aac"], 
	    'P' => ["cct", "ccc", "cca", "ccg"],
	    'Q' => ["caa", "cag"],
	    'R' => ["cgt", "cgc", "cga", "cgg", "aga", "agg"],
	    'S' => ["tct", "tcc", "tca", "tcg", "agt", "agc"],
	    'T' => ["act", "acc", "aca", "acg"],
	    'V' => ["gtt", "gtc", "gta", "gtg"],
	    'W' => ["tgg"], 
	    'Y' => ["tat", "tac"], 
	    'Z' => ["gaa", "gag", "caa", "cag"],
	    'X' => ["taa", "tag", "tga"], #stop codons
	     );

  # inverting %code to %three_ltr_code
  my %three_ltr_code;

  foreach my $ea (keys %code){
    foreach (@{$code{$ea}}){
      $three_ltr_code{$_} = $ea;
    }
  }

  #####################################
  # retrieve info from query parameters
  #####################################

  if (!$text){
    $run_window -> insert('end',"NO parameter!");
  }

  else {
    my $ref;
    @info=();
    ($cds_or_locus, $aa_or_dna, $mutation, $allele, $locus, $ref) = split(/\s+/, $params);
     $run_window -> insert('end', "$cds_or_locus, $aa_or_dna, $mutation, $allele, $locus, $ref");
    $position = $mutation;
    $position =~ s/\D//g;
  
    $mutation =~ s/\d+//; 
    if ($mutation =~ /\w{3,3}/){$mutation = $three_ltr_code{lc($mutation)}}
    $mutation = uc($mutation);
  
    push(@info, uc($mutation), $ref);

    if ($aa_or_dna eq "-aa" || $aa_or_dna eq "-AA"){$molecule = "aa"} else {$molecule = "DNA"}

    if ($cds_or_locus =~ /(.+\.\d+)(\w)/){
      my $variant = $2; my $seq = uc($1); 
      $cds = $seq.$variant; 
    }
    else {
      $cds = uc($cds_or_locus);
    }

  }

  #####################################
  # check for latest exon table version
  #####################################

  my $tace = &tace;
  my $curr_db = "/nfs/disk100/wormpub/DATABASES/current_DB";
  my $exon_tbl = "/wormsrv1/geneace/ALLELE_DATA/EXON_TABLES";

  $current = `grep "NAME WS" $curr_db/wspec/database.wrm`; $current =~ s/NAME //; chomp $current;
  @archive = `ls $exon_tbl`; foreach (@archive){chomp; if ($_ =~ /CDS.+(WS\d+)/){$archive = $1}} 

  if ("$current" ne "$archive") {	
    push(@output, "New WS release available . .\nFetching latest source exons of all CDS/Transcripts and 6 chromosomal DNA sequencs\n\n");
    system ("rm -f $exon_tbl/* ");
    `echo "table-maker -o $exon_tbl/CDS_exons_$current -p /wormsrv1/geneace/wquery/get_CDS_source_exons.def" | $tace $curr_db`;  
    system ("chmod 775 $exon_tbl/* ");
  }

  ####################################################
  # get DNA sequence (exon/intron) of a CDS/Transcript
  ####################################################

  my ($DNA, @coords, $chrom, $left, $right, $strand, $CDS);

  chdir "/wormsrv2/autoace/GFF_SPLITS/$current/";

  @CDS_coords = `grep $cds *.genes.gff *.rna.gff | cut -f 1,4,5,7,9`; 
 
  foreach (@CDS_coords){
    chomp;
    ($chrom, $left, $right, $strand, $CDS)= split(/\s+/, $_);
    $chrom =~ s/.+CHROMOSOME_//;
    push(@coords, $left, $right);
    @coords = sort {$a <=> $b} @coords;
  }

 # # 30 bp extension beyond 1st/last nucleotide
  $left = $coords[0] - 30;
  $right = $coords[-1] + 30;

  my $dna_file = "$curr_db/CHROMOSOMES/CHROMOSOME_".$chrom.".dna";
 
  my @line = `egrep "[atcg]" $dna_file`;
  my $line;

  foreach (@line){chomp; $line .= $_}

  if ($strand eq "-"){
    $DNA = substr($line, $left-1, $right-$left+1);
    $DNA = reverse $DNA; $DNA =~ tr/atcg/tagc/;
  }
  if ($strand eq "+"){
    $DNA = substr($line, $left-1, $right-$left+1);
  }

  @DNA = split(//, $DNA);

  if ($molecule eq "DNA"){

    $cds =~ /(.+)\..+/;
    my $seq = $1;

    my @dna_L = @DNA[$position-1..$position+28];
    my @dna_R = @DNA[$position+30..$position+59];
    my $dna_L = join('', @DNA[$position-1..$position+28]);
    my $dna_R = join('', @DNA[$position+30..$position+59]);

    $run_window -> delete('1.0', 'end');
    $run_window -> insert('end', "$current\n$cds\n\n");
    $run_window -> insert('end', "$DNA[$position+29] ($position) [Full-length DNA sequence: ");
    my $length = scalar @DNA-60;
    $run_window -> insert('end', "$length]\n\n");
    $run_window -> insert('end', "$dna_L");
    $run_window -> insert('end', "$red $DNA[$position-1]");
    $run_window -> insert('end', " $dna_R");

    $ace_window->insert('end', "\nLocus : \"$locus\"\n");
    $ace_window->insert('end', "\/\/Allele \"$allele\" Paper_evidence \"\"\n");
    $ace_window->insert('end', "\/\/Allele \"$allele\" PMID_evidence \"\"\n");
    $ace_window->insert('end', "\nAllele : \"$allele\"\n");
    $ace_window->insert('end', "\/\/Evidence Paper_evidence \"\"\n");
    $ace_window->insert('end', "\/\/Evidence PMID_evidence \"\"\n");
    $ace_window->insert('end', "Sequence \"$seq\"\n");
    $ace_window->insert('end', "\/\/Substitution \"[\/]\"\n");
    $ace_window->insert('end', "\/\/Deletion \"[\/]\"\n");
    $ace_window->insert('end', "Flanking_sequences \"$dna_L\" \"$dna_R\"\n");
    $ace_window->insert('end', "Species \"Caenorhabditis elegans\"\n");
    $ace_window->insert('end', "\/\/Remark \"\" Paper_evidence \"\"\n");
    $ace_window->insert('end', "\/\/Remark \"\" PMID_evidence \"\"\n");
    $ace_window->insert('end', "\/\/Remark \"\" Curator_confirmed \"WBPerson1845\"\n");
    $ace_window->insert('end', "\/\/Method \"Allele\"\n");
    $ace_window->insert('end', "\/\/Method \"Deletion_allele\"\n");
  }

  else {

    #############################################################################
    # get protein sequence (exon/intron) of a CDS/Transcript wormpep.current file
    #############################################################################

    open(IN1, "/nfs/disk100/wormpub/WORMPEP/wormpep_current") || die $!;
 
    my ($prot_seq, $DNA_seq);  

    $prot_seq = get_seq($cds, *IN1);
    #print "$cds\n$prot_seq\n\n";   
    @prot = split(//, $prot_seq);
    #print "aa $position = $prot[$position-1] [length = ", scalar @prot, "]\n";

    ########################################
    # fetch source exons of a CDS/Transcript
    ########################################

    my @exons = `grep $cds $exon_tbl/CDS*`;

    ########################################################################################
    # retrieving flank seq of a specified codon or mutation site via exons_to_codons routine
    ########################################################################################

    exons_to_codons($cds, \@exons, \@DNA, \@prot, $mutation, $position, $allele, $locus);
  }
}

sub get_seq {

  # get protein seq via wormpep.current file
  my ($cds, $FH) = @_;
  my $seq =();
  my $line_count = 0;

  while(<$FH>){
    chomp;
    if ($_ =~ /^>$cds/){$line_count++}
    if ($_ !~ /^>/ && $line_count > 0){$seq .= $_; $line_count++}
    if ($_ =~ /^>/ && $line_count > 1){last}	
  }
  return $seq;
}


#############################################################################################
# This chunk does several things:
# 1. process soruce exon coods to figure out frame shift 
# 2. the result of 1 is passed into codon_to_seq routine to retrieves 30 bp flanks of a codon
#############################################################################################

sub exons_to_codons {

  my ($cds, $exons, $DNA, $prot, $mutation, $position) = @_;
  my ($i, $j, $start, $end, @exon_start_end, $num_aa, $remainder, %codon_seq, $codon_seq, $ref, %remainder_hash, $start_bp, @return, @all);
  my $position_length = 0;
  my $position_codon = 0;


  foreach (@$exons){
    chomp;
    my ($cds, $start, $end) = split (/\s+/, $_);
    $start += 30; $end += 30; # 30 bp extension to get flank seq. of 1st/last amino acid
    push (@exon_start_end, $start, $end);
    $start=(); $end=();
  }

  for ($i = 0; $i < scalar @exon_start_end; $i=$i+2){
    if ($i == 0 || ($i > 0 && $remainder_hash{$i-2} == 0)){
      $remainder = ($exon_start_end[$i+1]-$exon_start_end[$i]+1)%3;
      $remainder_hash{$i} = $remainder;
    }
    if ($i > 0 && $remainder_hash{$i-2} == 1){
      $remainder = ($exon_start_end[$i+1]-$exon_start_end[$i]-1)%3;
      $remainder_hash{$i} = $remainder;
    }
    if ($i > 0 && $remainder_hash{$i-2} == 2){
      $remainder = ($exon_start_end[$i+1]-$exon_start_end[$i])%3;
      $remainder_hash{$i} = $remainder;
    }

    $start_bp = $exon_start_end[$i];

    ########################################################################################################################
    if ($remainder == 0){
      if ($i == 0 || ($i > 0 && $remainder_hash{$i-2} == 0) ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]+1;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "A", $position_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $position_codon = $return[1];
        push(@all, %codon_seq); next;
      }
      if ($i > 0 && $remainder_hash{$i-2} == 1 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-1;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "B", $position_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $position_codon = $return[1];
        push(@all, %codon_seq); next;
      }
      if ($i > 0 && $remainder_hash{$i-2} == 2 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i];	
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "C", $position_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $position_codon = $return[1];
        push(@all, %codon_seq); next;
      }   
    }
    #########################################################################################################################
    if ($remainder == 1){
      if ($i == 0 || ($i > 0 && $remainder_hash{$i-2} == 0)){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i];
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "D", $position_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $position_codon = $return[1];
        push(@all, %codon_seq); next;
      }
      if ($i > 0 && $remainder_hash{$i-2} == 1 ){
	my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-2;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "H", $position_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $position_codon = $return[1];
        push(@all, %codon_seq); next;
      }
      if ($i > 0 && $remainder_hash{$i-2} == 2 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-1;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "I", $position_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $position_codon = $return[1];
        push(@all, %codon_seq); next;
      }
    }
    ######################################################################################################################
    if ($remainder == 2){
      if ($i == 0 || ($i > 0 && $remainder_hash{$i-2} == 0) ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-1;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "E", $position_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $position_codon = $return[1];
        push(@all, %codon_seq); next;
      }
      if ($i > 0 && $remainder_hash{$i-2} == 1 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-3;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "F", $position_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $position_codon = $return[1];
        push(@all, %codon_seq); next;
      }
      if ($i > 0 && $remainder_hash{$i-2} == 2 ){
     	my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-2;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "G", $position_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $position_codon = $return[1];
        push(@all, %codon_seq); next;
      }
    }
    ######################################################################################################################
  }

  my $codons = scalar @all;
  for($i=0; $i < scalar @all; $i=$i+2){
    push(@{$codon_seq{$all[$i]}}, ${@{$all[$i+1]}}[0], ${@{$all[$i+1]}}[1], ${@{$all[$i+1]}}[2], ${@{$all[$i+1]}}[3], ${@{$all[$i+1]}}[4], ${@{$all[$i+1]}}[5]);
#     print "${@{$all[$i+1]}}[0], ${@{$all[$i+1]}}[1], ${@{$all[$i+1]}}[2], ${@{$all[$i+1]}}[3], ${@{$all[$i+1]}}[4], ${@{$all[$i+1]}}[5]\n";
  }  


  push(@output, "\n$prot[$position-1]($position) = "." @{$codon_seq{$position}}->[0] (". (@{$codon_seq{$position}}->[1]-30) .") @{$codon_seq{$position}}->[2] (". (@{$codon_seq{$position}}->[3]-30) . ") @{$codon_seq{$position}}->[4] (". (@{$codon_seq{$position}}->[5]-30) . ") [full-length aa of this gene: ". scalar @prot. "]\n");

  push(@info, "$prot[$position-1]($position)");

  my $codon = "@{$codon_seq{$position}}->[0]"."@{$codon_seq{$position}}->[2]"."@{$codon_seq{$position}}->[4]";


  ################################
  # output 30 bp flanks of a codon
  ################################

  push(@output, "-------------------------------------\n");
  push(@output, "   	Codon      ($prot[$position-1]):\t\t$codon\n");
  for ($i=0; $i < scalar @{$code{$mutation}}; $i++){
    push(@output, "	Mutated to \($mutation\):\t\t@{$code{$mutation}}->[$i-1]\n") if $mutation ne "X";
    push(@output, "	Mutated to \(STOP\):\t@{$code{$mutation}}->[$i-1]\n") if $mutation eq "X";
  }
  push(@output, "-------------------------------------");    

  my ($first_bp, $second_bp, $third_bp, $first_site, $second_site, $third_site);
  $first_bp = @{$codon_seq{$position}}->[1];
  $second_bp = @{$codon_seq{$position}}->[3];
  $third_bp = @{$codon_seq{$position}}->[5];
  $first_site = @{$codon_seq{$position}}->[0];
  $second_site = @{$codon_seq{$position}}->[2];
  $third_site = @{$codon_seq{$position}}->[4];

  #######################################################################
  # output 30 bp flank seq under frame shift or no frame shift situations
  #######################################################################

  my ($dna_L, $dna_R, $bp_num, @ace, $dna_Lf, $dna_Rf);

  @ace =();

  ################
  # no frame shift
  ################

  if ($first_bp == $second_bp-1 && $second_bp == $third_bp-1){
    for($i=1; $i<4; $i++){
      push(@output, "\n\n# 1st site mutation:\n") if $i == 1;
      push(@output,"\n\n# 2nd site mutation:\n") if $i == 2;
      push(@output, "\n\n# 3rd site mutation:\n") if $i == 3;
      $dna_L = join('', @DNA[$first_bp-32+$i..$first_bp-2]);
      push(@output, "black $dna_L");
      if ($i == 1){push(@output, "red $first_site"); push(@output, "blue $second_site"); push(@output, "blue $third_site"); $dna_Lf = $dna_L}
      if ($i == 2){push(@output, "blue $first_site"); push(@output, "red $second_site"); push(@output, "blue $third_site")}
      if ($i == 3){push(@output, "blue $first_site"); push(@output, "blue $second_site"); push(@output, "red $third_site")}
      $dna_R = join('', @DNA[$third_bp..$third_bp+26+$i]);
       push(@output, "black $dna_R");
      if ($i == 1){push(@ace, "1st $dna_L $second_site$third_site$dna_R\n")}
      if ($i == 2){push(@ace, "2nd $dna_L$first_site $third_site$dna_R\n")}
      if ($i == 3){push(@ace, "3rd $dna_L$first_site$second_site $dna_R\n"); $dna_Rf = $dna_R}
    }
    push(@ace, "codon $dna_Lf $dna_Rf\n");
  }

  ######################
  # second frame shifted
  ######################

  if ($first_bp != $second_bp-1 && $second_bp == $third_bp-1){

    push(@output,  "\n\n# 1st site mutation:\n");
    $bp_num = $first_bp+1-30; 
    push(@output,  "green_G $DNA[$first_bp]:  $bp_num");
    $dna_L = join('', @DNA[$first_bp-31..$first_bp-2]);
    push(@output, "black $dna_L");
    push(@output,  "red $first_site");
    push(@output,  "green $DNA[$first_bp]");
    $dna_R = join('', @DNA[$first_bp+1..$first_bp+29]);
    push(@output, "black $dna_R");
    push(@ace, "1st $dna_L $DNA[$first_bp]$dna_R\n");
    $dna_Lf = $dna_L;

    push(@output,  "\n\n# 2nd site mutation:\n"); 
    $bp_num = $second_bp-1-30;
    push(@output,  "magenta_M $DNA[$second_bp-2]:  $bp_num");
    $dna_L = join('', @DNA[$second_bp-31..$second_bp-3]);
    push(@output, "black $dna_L");
    push(@output,  "magenta $DNA[$second_bp-2]");
    push(@output,  "red $second_site");
    push(@output,  "blue $third_site");
    $dna_R = join('', @DNA[$third_bp..$third_bp+28]);
    push(@output, "black $dna_R");
    push(@ace, "2nd $dna_L$DNA[$second_bp-2] $third_site$dna_R\n");

    push(@output,  "\n\n# 3rd site mutation:\n");
    $bp_num = $second_bp-1-30;
    push(@output,  "magenta_M $DNA[$second_bp-2]: $bp_num");
    $dna_L = join('',  @DNA[$second_bp-30..$second_bp-3]);
    push(@output,  "black $dna_L");
    push(@output,  "magenta $DNA[$second_bp-2]");
    push(@output,  "blue $second_site");
    push(@output,  "red $third_site");
    $dna_R = join('', @DNA[$third_bp..$third_bp+29]);
    push(@output, "black $dna_R");
    push(@ace, "3rd $dna_L$DNA[$second_bp-2]$second_site $dna_R\n");
    $dna_Rf = $dna_R;

     push(@ace, "codon $dna_Lf $dna_Rf\n");
  }

  #####################
  # third frame shifted
  #####################

  if ($first_bp == $second_bp-1 && $second_bp != $third_bp-1){

    push(@output,  "\n\n# 1st site mutation:\n");
    $bp_num = $second_bp+1-30;
    push(@output,  "green_G $DNA[$second_bp]: $bp_num");
    $dna_L = join('',  @DNA[$first_bp-31..$first_bp-2]);
    push(@output,  "black $dna_L");
    push(@output,  "red $first_site");
    push(@output,  "blue $second_site");
    push(@output,  "green $DNA[$second_bp]");
    $dna_R = join('', @DNA[$second_bp+1..$second_bp+28]);
    push(@output, "black $dna_R");
    push(@ace, "1st $dna_L $second_site$DNA[$second_bp]$dna_R\n");
    $dna_Lf = $dna_L;

    push(@output,  "\n\n# 2nd site mutation:\n");
    $bp_num = $second_bp+1-30;
    push(@output,  "green_G $DNA[$second_bp]: $bp_num");
    $dna_L = join('',  @DNA[$first_bp-30..$first_bp-2]);
    push(@output,  "black $dna_L");
    push(@output,  "blue $first_site");
    push(@output,  "red $second_site"); 
    push(@output,  "green $DNA[$second_bp]");
    $dna_R = join('',  @DNA[$second_bp+1..$second_bp+29]);
    push(@output,  "black $dna_R");
    push(@ace, "2nd $dna_L$first_site $DNA[$second_bp]$dna_R\n");

    push(@output,  "\n\n# 3rd site mutation:\n");
    $bp_num = $third_bp-1-30;
    push(@output,  "magenta_M $DNA[$third_bp-2]: $bp_num");
    $dna_L = join('',  @DNA[$third_bp-31..$third_bp-3]);
    push(@output,  "black $dna_L");
    push(@output,  "magenta $DNA[$third_bp-2]");
    push(@output,  "red $third_site");
    $dna_R = join('', @DNA[$third_bp..$third_bp+29]);
    push(@output, "black $dna_R");
    push(@ace, "3rd $dna_L$DNA[$third_bp-2] $dna_R\n");
    $dna_Rf = $dna_R;

    push(@ace, "codon $dna_Lf $dna_Rf\n");
  }
  push(@output,  "\n");

  $run_window -> delete('1.0', 'end');
  $run_window -> insert('end', "$archive\n$cds\n");

  foreach (@output){
    if ($_ =~ /(red) (.+)/ || $_ =~ /(blue) (.+)/ || $_ =~ /(black) (.+)/ || $_ =~ /(magenta) (.+)/ || 
	$_ =~ /(green) (.+)/ || $_ =~ /(magenta_M) (.+)/ || $_ =~ /(green_G) (.+)/){
      my $color = $1; my $nt = $2;

      $run_window -> insert('end', "$red $nt") if $color eq "red";
      $run_window -> insert('end', "$black $nt") if $color eq "black";
      $run_window -> insert('end', "$blue $nt") if $color eq "blue";
      $run_window -> insert('end', "$green $nt\n") if $color eq "green_G";
      $run_window -> insert('end', "$magenta $nt\n") if $color eq "magenta_M";
      $run_window -> insert('end', "$green $nt") if $color eq "green";
      $run_window -> insert('end', "$magenta $nt") if $color eq "magenta";
    }
    else {
      $run_window -> insert('end',"$_");
    }
  }

   $label_1 -> configure(text => "@ace");
   $label_2 -> configure(text => "$allele");
   $label_3 -> configure(text => "$locus");

   $cds =~ /(.+)\..+/;
   my $parent = $1;
   $label_4 -> configure(text => "$parent");

   @output=();
}

sub get_1_site_flanks{
  my @ace = $label_1->cget("text");
  my $allele = $label_2->cget("text");
  my $locus = $label_3->cget("text");
  my $seq = $label_4->cget("text");

  write_ace(@ace);
}
  
sub get_2_site_flanks{
  my @ace = $label_1->cget("text");
  my $allele = $label_2->cget("text");
  my $locus = $label_3->cget("text");
  my $seq = $label_4->cget("text");

  write_ace(@ace);
}
  

sub get_3_site_flanks{
  my @ace = $label_1->cget("text");
  my $allele = $label_2->cget("text");
  my $locus = $label_3->cget("text");
  my $seq = $label_4->cget("text");
  write_ace(@ace);
}
  

sub get_codon_flanks{
  my @ace = $label_1->cget("text");
  my $allele = $label_2->cget("text");
  my $locus = $label_3->cget("text");
  my $seq = $label_4->cget("text");

  write_ace(@ace);
}

# get DNA triplet of a specified amino acid based on source exons processed in the above routine


sub codon_to_seq {
  my ($start_bp, $num_bp, $DNA, $option, $position_codon, $i, @exon_start_end) = @_;
  my %codon_seq;

 # print "\$i = $i\n";
 # print "start: $start_bp end: $exon_start_end[$i+1] = ", $exon_start_end[$i+1] - $exon_start_end[$i]+1, " Bp = $num_bp","\n";

  for (my $j=0; $j < $num_bp; $j=$j+3){
    my $pos = $j + $start_bp;
    $position_codon++;
    if ($option eq "A" || $option eq "D" || $option eq "E"){
      #print "Codon $position_codon(ADE): ${@$DNA}[$pos-1] ", $pos-30, " ${@$DNA}[$pos] ", $pos-30+1, " ${@$DNA}[$pos+1] ", $pos-30+2,"\n";
      push(@{$codon_seq{$position_codon}}, ${@$DNA}[$pos-1], $pos, ${@$DNA}[$pos], $pos+1, ${@$DNA}[$pos+1], $pos+2);
    }
    if ($option eq "B" || $option eq "H" || $option eq "F"){
      #print "Codon $position_codon(FH): ${@$DNA}[$pos+1] ", $pos-30+2, " ${@$DNA}[$pos+2] ", $pos-30+3, " ${@$DNA}[$pos+3] ", $pos-30+4,"\n";
      push(@{$codon_seq{$position_codon}}, ${@$DNA}[$pos+1], $pos+2, ${@$DNA}[$pos+2], $pos+3, ${@$DNA}[$pos+3], $pos+4);
    }
    if ($option eq "C" || $option eq "G" ||  $option eq "I"){
      #print "Codon $position_codon(CG): ${@$DNA}[$pos] ", $pos-30+1, " ${@$DNA}[$pos+1] ", $pos-30+2, " ${@$DNA}[$pos+2] ", $pos-30+3,"\n";
      push(@{$codon_seq{$position_codon}}, ${@$DNA}[$pos], $pos+1, ${@$DNA}[$pos+1], $pos+2, ${@$DNA}[$pos+2], $pos+3);
    }
  }

  if ($option eq "D" || $option eq "H"){
    $position_codon++;
    my $pos_up = $exon_start_end[$i+1];
    my $pos_down = $exon_start_end[$i+2];
    # print "\$i = $i\n";
    #print "Codon $position_codon(CDH): ${@$DNA}[$pos_up-1], ", $pos_up-30, " ${@$DNA}[$pos_down-1], ", $pos_down-30, " ${@$DNA}[$pos_down], ", $pos_down-30+1, "\n";
    push(@{$codon_seq{$position_codon}}, ${@$DNA}[$pos_up-1], $pos_up, ${@$DNA}[$pos_down-1], $pos_down, ${@$DNA}[$pos_down], $pos_down+1);
  }
  if ($option eq "E" || $option eq "F" || $option eq "G"){
    $position_codon++;
    my $pos_up = $exon_start_end[$i+1];
    my $pos_down = $exon_start_end[$i+2];
    #print "Codon $position_codon(EFG): ${@$DNA}[$pos_up-2], ", $pos_up-30-1, " ${@$DNA}[$pos_up-1], ", $pos_up-30, " ${@$DNA}[$pos_down-1], ", $pos_down-30, "\n";
    push(@{$codon_seq{$position_codon}}, ${@$DNA}[$pos_up-2], $pos_up-1, ${@$DNA}[$pos_up-1], $pos_up, ${@$DNA}[$pos_down-1], $pos_down);
  }
  if ($option eq "I"){
    $position_codon++;
    my $pos_up = $exon_start_end[$i+1];
    my $pos_down = $exon_start_end[$i+2];
    #print "Codon $position_codon(I): ${@$DNA}[$pos_up-1] ", $pos_up-30, " ${@$DNA}[$pos_down-1], ", $pos_down-30, " ${@$DNA}[$pos_down], ", #$pos_down-30+1,"\n";
    push(@{$codon_seq{$position_codon}}, ${@$DNA}[$pos_up-1], $pos_up, ${@$DNA}[$pos_down-1], $pos_down, ${@$DNA}[$pos_down], $pos_down+1);
  }

  return \%codon_seq, $position_codon;
}

sub write_ace {

  my @ace = @_;
  my $seq;

  foreach (@ace){
    if ($_ =~ /1st (\w+) (\w+)/){
      $ace_window->insert('end', "\nLocus : \"$locus\"\n");
      if ($info[1] =~ /\[.+\]/){ 
	$ace_window->insert('end', "Allele \"$allele\" Paper_evidence \"$info[1]\"\n");
      }
      else{
	$ace_window->insert('end', "Allele \"$allele\" PMID_evidence \"$info[1]\"\n");
      }
      $ace_window->insert('end', "\nAllele : \"$allele\"\n");
      if ($info[1] =~ /\[.+\]/){ 
	$ace_window->insert('end', "Evidence Paper_evidence \"$info[1]\"\n");
      }
      else {
	$ace_window->insert('end', "Evidence PMID_evidence \"$info[1]\"\n");
      }
      $ace_window->insert('end', "Sequence \"$seq\"\n");
      $ace_window->insert('end', "\/\/Substitution \"[\/]\"\n");
      $ace_window->insert('end', "\/\/Deletion \n");
      $ace_window->insert('end', "\/\/Insertion\n");
      $ace_window->insert('end', "\/\/Deletion_with_insertion\n");
      $ace_window->insert('end', "Flanking_sequences \"$1\" \"$2\"\n");
      $ace_window->insert('end', "Gene  \"$locus\"\n");
      $ace_window->insert('end', "Predicted_gene  \"$cds\"\n");
      $ace_window->insert('end', "Species \"Caenorhabditis elegans\"\n");
      $ace_window->insert('end', "\/\/Reference \"\"\n");
      if ($info[0] eq "X"){
	$ace_window->insert('end', "\/\/Nonsense \"Amber_UAG\"\n");
	$ace_window->insert('end', "\/\/Nonsense \"Ochre_UAA\"\n");
	$ace_window->insert('end', "\/\/Nonsense \"Opal_UGA\"\n");
      }
      else {
	$ace_window->insert('end', "Missense\n");
      }
      if ($info[1] =~ /\[.+\]/){ 
	$ace_window->insert('end', "Remark \"$allele is a $info[2] to $info[0] mutation\" Paper_evidence \"$info[1]\"\n");
      }
      else {
	$ace_window->insert('end', "Remark \"$allele is a $info[2] to $info[0] mutation\" PMID_evidence \"$info[1]\"\n");
      }
      $ace_window->insert('end', "\/\/Remark \"\" Curator_confirmed \"WBPerson1845\"\n");
      $ace_window->insert('end', "\/\/Method \"Allele\"\n");
      $ace_window->insert('end', "\/\/Method \"Deletion_allele\"\n");
      $ace_window->insert('end', "\/\/Method \"Transposon_insertion\"\n");
    }
  }
}

__END__

#Output: G(45 )= g (590) g (591) g (640) [full-length aa of this gene = 255]

=head2 NAME - get_allele_flank_seq.pl  

DESCRIPTION

  First thing first: perl path is perl5.6.1, 5.8.0 won't run

  This script is suitable for curating allele flanking sequences described in paper such as "
  .... amino acid Q at position 235 is mutated to G ..." or "...nucleotide t in position 1234 of 4R79.1 is mutated  to g ...".

  This script can output 30 bp flanking DNA seqs on two sides of any mutated site of a triplet of a defined gene (including the 
  first and last codon, although such alleles rarely occur). 

  You need to supply a CDS/Transcript name (case-insensitive) and amino acid coordinate (or nucleotide coordinate) followed by 
  a mutation in single- or three-letter code (case-insensitive) as arguments (see USAGE) to retrieve the flanking sequences.  

SCENARIO A: DNA coordinate

  Nucleotide 1234 is mutated from t to g. You want to retrieve the 30 bp flanking sequences on two sides of mutation side.

  OUTPUT:

  The output is fairly simple and is similar to:

  (1) t (123) [Full-length DNA sequence: 2017]

  (2) acagactacttaaacattgtaaaaggatat g ggtaagaatatatatcttatacaaccctta

  Comments:

    (1) Verification. If the nt specified is identical, then it is a good sign. 

    (2) The nucleotide g in the middle separated by spaces is the mutation.
        ace file template with minimal information is generated without clicking on the mutation site buttons in the lower
        window where you can modify and finish according to some extra information described in paper.


SCENARIO B: amino acide coordinate

  An amino acid G has DNA triplet ggg at positions 590, 591 and 640 (3rd frame shifted),
  and has a missense mutation to E. 
  You want to retrieve the flanking seq. of the mutated sites in a triplet.

  OUTPUT:

  (1) G(45): g (590) g (591) g (640) [full-length aa of this gene: 255]

  (2) Triplet color coding: mutated site in RED, the other two sites behind or before it in the triplet in BLUE
      Flanking bp coding when frame shift:  upsteam of mutated site in MAGENTA, downstream in GREEN
      This color coding does not apply to dinucleotide mutation sites. 

  (3)
     ------------------------------
        Codon       (G):    ggg
        Mutation to (E):    gag
        Mutation to (E):    gaa
     ------------------------------

  (4)

     # 1st site mutation:
     g: 592
     acagactacttaaacattgtaaaaggatat g ggtaagaatatatatcttatacaaccctta

     # 2nd site mutation:
     g: 592
     cagactacttaaacattgtaaaaggatatg g gtaagaatatatatcttatacaacccttac

     # 3rd site mutation:
     g: 639 t: 641
     tacaacccttactgaattttaatttttcag g tgctactctcaagttggacgaactggagga


     Comments on output (generated in the upper window, which is overwritten by a new query result)
	
     (1) Verification. This simply tells you that G(45) described in paper (scenario) is the same as current WS dataset. 
         So should be OK to run script for this allele, . . . usually.  

     (2) Color coding (cannot be seen here, but when you run the script) helps you quickly identify the flanking 
         sequences (4) of a mutated site (the one in between spaces), especially in cases where frame shift occur so that 
         the immediate flanking nucleotide maybe in the intron between two sites of a codon.  

     (3) As three potential single-site mutation can occur in a codon, (3) gives you genetic codes of the amino acid 
         resulted in mutation and allows you a quick look up of bp substitution. This table is also fine 
         for a dinucleotide mutation.
         In the scenario, eg, if codon (G) has a missense mutation and changed to E, the codon table conveniently tells you 
         that ggg has been mutated to gag. So, this would be a 2nd site mutation or a [g/a] substitution.
         You should then choose the matching flanking sequences in (4).       

     (4) By clicking on the 1 site mutation, 2 site mutation ...etc, an ace file template with minimal information is 
         generated in the lower window where you can modify and finish according to some extra information described in paper.

UPLOADING ACE FILE TO GENEACE:

     This requires write access to the database.
#####


foreach (@ace){
    if ($_ =~ /2nd (\w+) (\w+)/){
      $ace_window->insert('end', "\nLocus : \"$locus\"\n");
      $ace_window->insert('end', "\/\/Allele \"$allele\" Paper_evidence \"\"\n");
      $ace_window->insert('end', "\/\/Allele \"$allele\" PMID_evidence \"\"\n"); 
      $ace_window->insert('end', "\nAllele : \"$allele\"\n");
      $ace_window->insert('end', "\/\/Evidence Paper_evidence \"\"\n");
      $ace_window->insert('end', "\/\/Evidence PMID_evidence \"\"\n");
      $ace_window->insert('end', "Sequence \"$seq\"\n");
      $ace_window->insert('end', "\/\/Substitution \"[\/]\"\n");
      $ace_window->insert('end', "\/\/Deletion\n");
      $ace_window->insert('end', "\/\/Insertion\n");
      $ace_window->insert('end', "\/\/Deletion_with_insertion\n");
      $ace_window->insert('end', "Flanking_sequences \"$1\" \"$2\"\n");
      $ace_window->insert('end', "Gene  \"$locus\"\n");
      $ace_window->insert('end', "Predicted_gene  \"$cds\"\n");
      $ace_window->insert('end', "Species \"Caenorhabditis elegans\"\n");
      $ace_window->insert('end', "\/\/Reference \"\"\n");
      if ($info[0] eq "X"){
	$ace_window->insert('end', "\/\/Nonsense \"Amber_UAG\"\n");
	$ace_window->insert('end', "\/\/Nonsense \"Ochre_UAA\"\n");
	$ace_window->insert('end', "\/\/Nonsense \"Opal_UGA\"\n");
      }
      $ace_window->insert('end', "Missense\n");
      $ace_window->insert('end', "\/\/Remark \"$allele is a $info[1] to $info[0] mutation\" Paper_evidence \"\"\n");
      $ace_window->insert('end', "\/\/Remark \"$allele is a $info[1] to $info[0] mutation\" PMID_evidence \"\"\n");
      $ace_window->insert('end', "\/\/Remark \"\" Curator_confirmed \"WBPerson1845\"\n");
      $ace_window->insert('end', "\/\/Method \"Allele\"\n");
      $ace_window->insert('end', "\/\/Method \"Deletion_allele\"\n");
      $ace_window->insert('end', "\/\/Method \"Transposon_insertion\"\n");
    }
  }
}
###
foreach (@ace){
    if ($_ =~ /3rd (\w+) (\w+)/){
      $ace_window->insert('end', "\nLocus : \"$locus\"\n");
      $ace_window->insert('end', "\/\/Allele \"$allele\" Paper_evidence \"\"\n");
      $ace_window->insert('end', "\/\/Allele \"$allele\" PMID_evidence \"\"\n"); 
      $ace_window->insert('end', "\nAllele : \"$allele\"\n");
      $ace_window->insert('end', "\/\/Evidence Paper_evidence \"\"\n");
      $ace_window->insert('end', "\/\/Evidence PMID_evidence \"\"\n");
      $ace_window->insert('end', "Sequence \"$seq\"\n");
      $ace_window->insert('end', "\/\/Substitution \"[\/]\"\n");
      $ace_window->insert('end', "\/\/Deletion \n");
      $ace_window->insert('end', "\/\/Insertion\n");
      $ace_window->insert('end', "\/\/Deletion_with_insertion\n");
      $ace_window->insert('end', "Flanking_sequences \"$1\" \"$2\"\n");
      $ace_window->insert('end', "Gene  \"$locus\"\n");
      $ace_window->insert('end', "Predicted_gene  \"$cds\"\n");
      $ace_window->insert('end', "Species \"Caenorhabditis elegans\"\n");
      $ace_window->insert('end', "\/\/Reference \"\"\n");
      if ($info[0] eq "X"){
	$ace_window->insert('end', "\/\/Nonsense \"Amber_UAG\"\n");
	$ace_window->insert('end', "\/\/Nonsense \"Ochre_UAA\"\n");
	$ace_window->insert('end', "\/\/Nonsense \"Opal_UGA\"\n");
      }
      $ace_window->insert('end', "Missense\n");
      $ace_window->insert('end', "\/\/Remark \"$allele is a $info[1] to $info[0] mutation\" Paper_evidence \"\"\n");
      $ace_window->insert('end', "\/\/Remark \"$allele is a $info[1] to $info[0] mutation\" PMID_evidence \"\"\n");
      $ace_window->insert('end', "\/\/Remark \"\" Curator_confirmed \"WBPerson1845\"\n");
      $ace_window->insert('end', "\/\/Method \"Allele\"\n");
      $ace_window->insert('end', "\/\/Method \"Deletion_allele\"\n");
      $ace_window->insert('end', "\/\/Method \"Transposon_insertion\"\n");
    }
  }
}
##

  foreach (@ace){
    if ($_ =~ /codon (\w+) (\w+)/){
      $ace_window->insert('end', "\nLocus : \"$locus\"\n");
      $ace_window->insert('end', "\/\/Allele \"$allele\" Paper_evidence \"\"\n");
      $ace_window->insert('end', "\/\/Allele \"$allele\" PMID_evidence \"\"\n");
      $ace_window->insert('end', "\nAllele : \"$allele\"\n");
      $ace_window->insert('end', "\/\/Evidence Paper_evidence \"\"\n");
      $ace_window->insert('end', "\/\/Evidence PMID_evidence \"\"\n");
      $ace_window->insert('end', "Sequence \"$seq\"\n");
      $ace_window->insert('end', "\/\/Substitution \"[\/]\"\n");
      $ace_window->insert('end', "\/\/Deletion \n");
      $ace_window->insert('end', "\/\/Insertion\n");
      $ace_window->insert('end', "\/\/Deletion_with_insertion\n");
      $ace_window->insert('end', "Flanking_sequences \"$1\" \"$2\"\n");
      $ace_window->insert('end', "Gene  \"$locus\"\n");
      $ace_window->insert('end', "Predicted_gene  \"$cds\"\n");
      $ace_window->insert('end', "Species \"Caenorhabditis elegans\"\n");
      $ace_window->insert('end', "\/\/Reference \"\"\n");
      if ($info[0] eq "X"){
	$ace_window->insert('end', "\/\/Nonsense \"Amber_UAG\"\n");
	$ace_window->insert('end', "\/\/Nonsense \"Ochre_UAA\"\n");
	$ace_window->insert('end', "\/\/Nonsense \"Opal_UGA\"\n");
      }
      $ace_window->insert('end', "Missense\n");
      $ace_window->insert('end', "\/\/Remark \"$allele is a $info[1] to $info[0] mutation\" Paper_evidence \"\"\n");
      $ace_window->insert('end', "\/\/Remark \"$allele is a $info[1] to $info[0] mutation\" PMID_evidence \"\"\n");
      $ace_window->insert('end', "\/\/Remark \"\" Curator_confirmed \"WBPerson1845\"\n");
      $ace_window->insert('end', "\/\/Method \"Allele\"\n");
      $ace_window->insert('end', "\/\/Method \"Deletion_allele\"\n");
      $ace_window->insert('end', "\/\/Method \"Transposon_insertion\"\n");
    }
  }
}
