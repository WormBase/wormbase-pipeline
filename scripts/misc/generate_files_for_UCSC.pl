#!/software/bin/perl -w

#This is the details that Jim Kent originally sent on how he'd like the files to load WormBase genes in to 
# UCSC genome browser.


#1) The genome sequence in a fasta file with one sequence per
#      chromosome.  The fasta headers should be:
#           >chrI
#            sequenc here
#            >chrII
#  and so forth.  That is, just with "chr" in front of the traditional upper case
#  chromosome name.  Please include the mitochondria as chrM.
#2) An AGP file that describes how the genome is built from clones.
#   It looks like Sanger might have worm AGP files already as described
#   at http://www.sanger.ac.uk/Projects/C_elegans/DOCS/agp_files.shtml.
#   NCBI has an AGP file format spec at
#    http://www.ncbi.nlm.nih.gov/genome/guide/Assembly/AGP_Specification.html
#   It works out easiest for us if the first column of the AGP has "chr" in front of
#   the chromosomes too.
#3) The WormBase gene models in GFF format. Please include *both* exon
#   and CDS lines for the coding portion of coding exons, and just exon for the
#   UTR. These should be in chromosome coordinates, and again the chromosomes
#   with chr in front of them.  For the column containing the gene name, please
#   use transcript ID (something like ZC101.2e).  It's ok to keep the GFF file very simple in
#   that final field.  If you want to be more complex, please use the GTF subset of
#   GFF.  See http://genome.ucsc.edu/FAQ/FAQformat#format3 and
#   http://genome.ucsc.edu/FAQ/FAQformat#format4 for our interpretation of GFF
#   and GTF.
#4) A tab-separate two column table with first column transcript ID and second column
#     gene name (things like unc-47).
#5) A tab-separated file with first column transcript ID and second column a descriptive
#     sentence or so.
#6) A tab-separated file with first column transcript ID and second column a descriptive
#     paragraph or two.
#7) A tab-separated file with first column transcript ID and second column a URL to
#     link into your web works.
#8) A tab-separated file with first column transcript ID and second column a synonym for
#     the gene.
#9) A tab-separated file with first column transcript ID and second column UniProt ID.

#see wiki for details on editing the storable file if you need to do this outside of the build
#http://scratchy.internal.sanger.ac.uk/wiki/index.php/UsefulTitBits

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use File::Path;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($version, $dna, $agp, $data, $gff, $chroms, $package, $ftp, $all);

GetOptions (	"help"          => \$help,
            	"debug=s"       => \$debug,
		"test"          => \$test,
		"verbose"       => \$verbose,
		"store:s"       => \$store,
		"version:i"     => \$version,
		"dna"	        => \$dna,
		"agp"           => \$agp,
		"data"          => \$data,
		"gff"           => \$gff,
		"package"       => \$package,
		"all"           => \$all,
		"chromosomes:s" => \$chroms,
		"ftp"           => \$ftp,
	   );

($dna=$agp=$data=$gff=$package= $all) if $all;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

my $log = Log_files->make_build_log($wormbase);
$log->write_to("options used : ",join("\n",@ARGV),"\n");
$log->log_and_die("Please enter a WS version number : -version WS170") unless $version;

$version = 'WS'.$version unless ($version =~ /\WS/);
my $ftp_dir = "/nfs/disk69/ftp/pub2/wormbase/";
my $base_dir = "/nfs/wormpub/DATABASES";
my $database;
if (-e "${base_dir}/$version") {
  $database = $wormbase->database($version);
}
elsif (-e "${base_dir}/frozen_release/$version") {
  $database = "${base_dir}/frozen_release/$version";
  $log->write_to("This is a frozen release we might need to unarchive some data\n------------------------------------------\n");
  unless (-e $database."/CHROMOSOMES/CHROMOSOME_I.gff") {
    $log->write_to("Copying GFF data.....\n");
    mkpath $database."/CHROMOSOMES";
    $wormbase->run_command("scp ${ftp_dir}WS205/genomes/c_elegans/genome_feature_tables/GFF2/*.gff $database/CHROMOSOMES/",$log);
  }
  unless (-e $database."/CHROMOSOMES/CHROMOSOME_I.dna") {
    $log->write_to("Copying AGP data.....\n");
    $wormbase->run_command("${ftp_dir}WS205/genomes/c_elegans/sequences/dna/*.agp $database/CHROMOSOMES/",$log);
    $log->write_to("Copying DNA data.....\n");
    $wormbase->run_command("${ftp_dir}WS205/genomes/c_elegans/sequences/dna/*.dna $database/CHROMOSOMES/",$log);
    $log->write_to("Unzipping data.....\n");
    $wormbase->run_command("gunzip $database/CHROMOSOMES/*.gz",$log);
  }
}
else {
  $log->log_and_die("Database specified: $version could not be located.\n");
}

my @chromosomes;
defined $chroms ? (@chromosomes = split(/\,/,$chroms)) : (@chromosomes = qw(I II III IV V X));

$log->write_to("\nCreating UCSC genome browser files for chromosomes ".join(" ",@chromosomes)." in $database\n");
my $output_dir = $database."/UCSC";
mkpath $output_dir unless -e $output_dir;

&dna_files if ($dna or $all);
&agp_files if ($agp or $all);
&gff_files if ($gff or $all);
&transcript_data_files if ($data or $all);

#tar and zip the data ready for the ftp site.
my $tmp_dir = "/nfs/wormpub/tmp/";
my $outfile = "wormbase_${version}_UCS.tar.gz";

$wormbase->run_command("tar -C $output_dir -zcf ${tmp_dir}$outfile .",$log) if ($package or $all);
#mv tar back to UCSC dir.
$wormbase->run_command("cp ${tmp_dir}$outfile ${output_dir}",$log) if ($package or $all);
$wormbase->run_command("rm ${tmp_dir}$outfile",$log) if ($package or $all);
# create md5sum file.
$wormbase->run_command("md5sum ${output_dir}/$outfile > $output_dir/${version}_md5sum",$log) if ($package or $all);
#copy the data to the ftp dir.
$wormbase->run_command("scp ${output_dir}/$outfile ${ftp_dir}data/UCSC/",$log) if ($ftp or $all);
$wormbase->run_command("scp ${output_dir}/${version}_md5sum ${ftp_dir}data/UCSC/",$log) if ($ftp or $all);

$log->mail;
exit;


#################
#  Subroutines  #
#################

sub gff_files {
  $log->write_to("Creating GFF files . . . ");
  foreach my $chrom (@chromosomes) {
    my $gff = "$database/UCSC/chr${chrom}.gff";
    if (-e $database."/GFF_SPLITS/CHROMOSOME_${chrom}_Coding_transcript.gff") {
      $wormbase->run_command("cat $database/GFF_SPLITS/CHROMOSOME_${chrom}_Coding_transcript.gff | grep exon | sed s/exon/CDS/  | sed s/CHROMOSOME_/chr/ > $gff", $log);
      $wormbase->run_command("cat $database/GFF_SPLITS/CHROMOSOME_${chrom}_UTR.gff | grep coding_exon | sed s/coding_exon/exon/ | sed s/CHROMOSOME_/chr/ >> $gff", $log);
      $wormbase->run_command("cat $database/GFF_SPLITS/CHROMOSOME_${chrom}_*RNA*.gff | grep exon | sed s/CHROMOSOME_/chr/ >> $gff", $log);
    }
    elsif (-e $database."/CHROMOSOMES/CHROMOSOME_${chrom}.gff") {
      $wormbase->run_command("cat $database/CHROMOSOMES/CHROMOSOME_${chrom}.gff | grep Coding_transcript | grep exon | sed s/exon/CDS/  | sed s/CHROMOSOME_/chr/ > $gff", $log);
      $wormbase->run_command("cat $database/CHROMOSOMES/CHROMOSOME_${chrom}.gff | grep Coding_transcript | grep coding_exon | grep CDS | sed s/coding_exon/exon/ | sed s/CHROMOSOME_/chr/ >> $gff", $log);
      $wormbase->run_command("cat $database/CHROMOSOMES/CHROMOSOME_${chrom}.gff | grep exon | grep RNA | grep -v RNASEQ | grep -v history | grep -v jigsaw | sed s/CHROMOSOME_/chr/ >> $gff", $log);
    }
    else {
      $log->log_and_die("GFF data for Coding_transcript/UTR/RNA exons is missing.\n");
    }
    `cat $gff | sort -u | sed s/\\"//g | sed s/'Transcript '// | cut -f1 -d';' > tmp_file`;
    $wormbase->run_command("mv -f tmp_file $gff",$log);
  }
  $log->write_to("Done \n");
}

sub agp_files {
  $log->write_to("Creating AGP files . . . ");
  #	III     1       3906    1       F       AL031226.2      1       3906    +
  # just need to change III -> chrIII and accn to clone name
  my %acc2clone = $wormbase->FetchData('accession2clone', undef, "$database/COMMON_DATA");
  foreach my $chrom (@chromosomes) {
    
    #Find the agp files
    my $agpf;
    if (-e $database."/yellow_brick_road/CHROMOSOME_${chrom}.agp") {
      $agpf = "$database/yellow_brick_road/CHROMOSOME_${chrom}.agp";
    }
    elsif (-e $database."/CHROMOSOMES/CHROMOSOME_${chrom}.agp") {
      $agpf = "$database/CHROMOSOMES/CHROMOSOME_${chrom}.agp";
    }
    else {
      $log->log_and_die("Can't locate AGP files for $database, do you need to unarchive them?\n");
    }
    open(AGP,"<$agpf") or die "cant open agp for $chrom\n";
    open(NEW,">$output_dir/chr${chrom}.agp")or die "cant open agp for $chrom\n";
    while(<AGP>) {
      next if (/\/\//);
      next if (/acedb>/);
      my @data = split;
      $data[0] = 'chr'.$data[0];
      my ($acc) = $data[5] =~ m/(\w+)\.\d/;
      die "bad clone\n" unless ($acc and $acc2clone{$acc});
      $data[5] = $acc;#$acc2clone{$acc};
      print NEW join("\t",@data), "\n";
    }
  }
  $log->write_to("Done \n");
}


sub dna_files {
  $log->write_to("Creating DNA files . . . ");
  foreach my $chrom (@chromosomes) {
    unless (-e "${database}/CHROMOSOMES/CHROMOSOME_${chrom}.dna") {
      $log->log_and_die("The $database/CHROMOSOMES/*.dna files are missing, do you need to unarchive them?\n");
    }
    $wormbase->run_command("cat $database/CHROMOSOMES/CHROMOSOME_${chrom}.dna | sed s/CHROMOSOME_/chr/ > $database/UCSC/chr${chrom}.dna",$log);
  }
  $log->write_to("Done \n");
}


sub transcript_data_files {
  $log->write_to("Creating transcript data files . . . ");
  my $def;
  if (-e $database."/wquery/SCRIPT:UCSC_browser_files.def") {
    $def = $database."/wquery/SCRIPT:UCSC_browser_files.def";
  }
  else {
    $def = $wormbase->autoace."/wquery/SCRIPT:UCSC_browser_files.def";
  }
  $log->error("ERROR: $def does not exist\n") unless -e $def;
  my $data = $wormbase->table_maker_query($database,$def);
  
  #open files to write
  open (PUB_NAME,">$database/UCSC/gene_name.txt") or die ("cant open $database/UCSC/gene_name.txt");
  open (SENTENCE,">$database/UCSC/sentence.txt")  or die ("cant open $database/UCSC/sentence.txt");
  open (PARA,    ">$database/UCSC/paragraph.txt") or die ("cant open $database/UCSC/paragraph.txt");
  open (URL,     ">$database/UCSC/url.txt")       or die ("cant open $database/UCSC/url.txt");
  open (UNIPROT,">$database/UCSC/uniprot.txt")    or die ("cant open $database/UCSC/uniprot.txt");
  
  my $error6;
  while(<$data>) {
    s/\"//g;#"
    chomp;
    next if (/\/\//);
    next if (/acedb>/);
    my @info = split("\t",$_);
    if (!defined $info[6]) {
      $log->write_to("$info[0] Lacks UniProt data\n") if ($verbose);
      $error6++;
      next;
    }
    if ($info[6] ne 'UniProtAcc'){
      next;
    }
    #		next unless (/SwissProt_AC/ or /TrEMBL_AC/)
    print PUB_NAME "$info[0]\t$info[2]\n";
    print SENTENCE "$info[0]\t$info[3]\n";
    #paragraph is DB_remark too unless concise description.		
    my $paragraph = ($info[4] eq "") ? $info[3] : $info[4];
    print PARA     "$info[0]\t$paragraph\n";
    print URL      "$info[0]\thttp://wormbase.org/db/gene/gene?name=$info[1];class=Gene\n";
    print UNIPROT  "$info[0]\t$info[7]\n";
  }
  $log->write_to("\n$error6 objects missing UniProtAcc numbers.\n");
  close $data;
  
  # and non-coding genes
  my $nc_def;
  if (-e $database."/wquery/SCRIPT:UCSC_browser_files_nc.def") {
    $nc_def = $database."/wquery/SCRIPT:UCSC_browser_files_nc.def";
  }
  else {
    $nc_def = $wormbase->autoace."/wquery/SCRIPT:UCSC_browser_files_nc.def";
  }
  my $nc_data = $wormbase->table_maker_query($database,$nc_def);
  my $linecount = "0";
  my $error2 = "0";
  my $error3 = "0";
  my $error4 = "0";

  while(<$nc_data>) { 
    next if (/\/\//);
    next if (/acedb>/);
    s/\"//g;#"
    chomp;
    my @info = split("\t",$_);
    $linecount++;
    next unless defined $info[2];
   #paragraph is DB_remark too unless concise description.
    if (!defined $info[2]) {
      $info[2] = "";
      $log->write_to("Line: $0 Lacks data in element 2\n")if ($verbose);
      $error2++;
    }
    if (!defined $info[3]) {
      $info[3] = "";
      $log->write_to("Line: $0 Lacks data in element 3\n")if ($verbose);
      $error3++;
    }
    if (!defined $info[4]) {
      $info[4] = "";
      $log->write_to("Line: $0 Lacks data in element 4\n") if ($verbose);
      $error4++;
    }
    print PUB_NAME "$info[0]\t$info[2]\n" if (defined $info[2]);
    print SENTENCE "$info[0]\t$info[3]\n" if (defined $info[3]);
    my $paragraph = ($info[4] eq "") ? $info[3] : $info[4];
    print PARA     "$info[0]\t$paragraph\n";
    print URL      "$info[0]\thttp://wormbase.org/db/gene/gene?name=$info[1];class=Gene\n";
  }
   $log->write_to("\nError Summary:\n--------------\n\@info Lacking element 2 = $error2\n\@info Lacking element 3 = $error3\n\@info Lacking element 4 = $error4\n------------------------------\nLines Processed = $linecount \n\n");
  close $nc_data;
  close PUB_NAME;
  close SENTENCE ;
  close PARA ;
  close URL;
  close UNIPROT;
  $log->write_to("Done \n");
}
