#!/usr/local/bin/perl5.8.0 -w
#
# find_utrs.pl
#
# usage find_utrs.pl -database <database> -output_dir <output dir>
#
# Ashwin Hajarnavis ah3@sanger.ac.uk  August 2002
#
# Last updated by: $Author: krb $                 
# Last updated on: $Date: 2003-12-08 16:40:33 $   


use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Data::Dumper;

#################################
# Command-line options          #
#################################

my ($debug, $database, $output_dir, $test);

GetOptions ("debug=s"      => \$debug,
            "database=s"   => \$database,
            "output_dir=s" => \$output_dir,
            "test"         => \$test);


my $usage = "find_utrs.pl -database [WS\#\#] -output_dir [path to output directory] (for Ace file, gff file, error file)\n";
die "$usage" unless ($database);


# database/file paths and locations
my $basedir  = "/wormsrv2";
$basedir     = glob("~wormpub")."/TEST_BUILD" if ($test); 
my $gff_dir  = "$basedir/autoace/GFF_SPLITS/$database";
my $dbdir    = "$basedir/$database/";

# touch logfile for run details
$0 =~ m/\/*([^\/]+)$/; system("touch $basedir/logs/history/$1.`date +%y%m%d`");

my $tace      = &tace;


##############################
# Paths etc                  #
##############################


if ($database eq "current_DB") {
    my $WS = &get_wormbase_version($dbdir);
    $gff_dir = "$basedir/autoace/GFF_SPLITS/WS$WS";
}
elsif ($database eq "autoace") {
    my $WS = &get_wormbase_version($dbdir);
    $WS = "WS666" if ($test);
    $gff_dir = "$basedir/autoace/GFF_SPLITS/GFF_SPLITS";
}



my $query_def = "$basedir/autoace/wquery/cDNA_CDS_EST5_EST3_METHOD.def";
my $genes     = "genes.gff";
my $est_file  = "BLAT_EST_BEST.gff";
my $mrna_file = "BLAT_mRNA_BEST.gff";

&log_file($output_dir, $dbdir) if $output_dir;

my %link_coordinate_start;
my %link_coordinate_end;
my %link_chrom;
my $chrom2parse;

# Link data
open (LINK, "<$basedir/autoace/BLAT/chromosome.ace") || die "whoops $!";
while (<LINK>) {
  if (/^Sequence \: \"CHROMOSOME_(\S+)\"/) {
    $chrom2parse = $1;
    next;
  }
  
  if (/Subsequence\s+\"(\S+)\"\s+(\d+)\s+(\d+)/) {
    $link_coordinate_start{$1} = $2;
    $link_coordinate_end{$1}   = $3;
    $link_chrom{$1}            = $chrom2parse;
  }
}
close LINK;


# check to see if EST hash data exists
# make it via tablemaker queries if absent

my %ESTacc2ESTname;
my %EST_name;
my %EST_dir;

unless (-e "$basedir/autoace/BLAT/EST.dat") {
    (%EST_name,%EST_dir) = &make_EST_hash;
}
# else read it into memory
else {
    print "Read hash from file EST.dat\n";
    open (FH, "<$basedir/autoace/BLAT/EST.dat") or die "EST.dat : $!\n";
    undef $/;
    my $data = <FH>;
    eval $data;
    die if $@;
    $/ = "\n";
    close FH;
}

print "Reverse key/value for EST_name hash .....";
%ESTacc2ESTname = reverse %EST_name;
print "OK.\n";

#
# Create hashes for gene positions and transcript mappings

my ($gene_est_map, $est_data) = parse_table($dbdir, $query_def);

# write gene gff positions
my $gene_data = read_gff($gff_dir, $genes);

# write EST gff positions
$est_data = read_gff($gff_dir, $est_file, $est_data);

# write mRNA gff positions
$est_data = read_gff($gff_dir, $mrna_file, $est_data);

 ##############################
 # Write output               #
 ##############################

if ($output_dir)  {
  open(OUTFILE, ">$output_dir/ALL.UTRs.tmp")|| die "Cannot open outfile $output_dir\n";
  open(ERRORFILE, ">$output_dir/data_errors.txt")|| die "Cannot dump errors\n";
  open(ACEFILE, ">$output_dir/UTRs.ace")||die "Cannot open ace file\n";
}

&find_transcript($gene_est_map, $gene_data, $est_data);

if ($output_dir) {
  close(OUTFILE);
  close(ERRORFILE);
  close(ACEFILE);
  &sort_gff();
}
 

exit(0);

#####################################################################################################
#                                                                                                   #
# Subroutines:                                                                                      #
#                                                                                                   #
#####################################################################################################



sub find_transcript {
    my ($gene_est_map,$gene_data,$est_data) = @_;
    my $a = "a";
    my ($gene_id, $gene_chrom, $gene_start, $gene_end, $gene_strand, $three_start, $three_end, $five_start, $five_end);
    
    my $ests;
    while (($gene_id, $ests) = each (%{$gene_est_map})) {
	unless (defined($gene_data->{$gene_id})) {
	    print ERRORFILE "$gene_id is not in GFF\n" if $output_dir;  ### gene_id appears in tablemaker query but not in GFF
	    next;
	}
	
	if ($gene_id =~ /(.*)([a-z])$/) {
	    if ($2 ne "a") {
		unless (defined($gene_data->{$1.$a})) {
		    print ERRORFILE "$gene_id exists but $1"."a does not\n" if $output_dir;  ### for cases where there appears to be a .b variant, but no .a variant in GFF.
		    next;
		}
	    }
	}
	
	$gene_chrom  = $gene_data->{$gene_id}->[0];
	$gene_start  = $gene_data->{$gene_id}->[1];
	$gene_end    = $gene_data->{$gene_id}->[2];
	$gene_strand = $gene_data->{$gene_id}->[3];
	
	
	my ($transcript_start, $transcript_end, $min_cdna, $max_cdna) = find_extent($ests, $est_data, $gene_chrom, $gene_start, $gene_end, $gene_id, $gene_strand);
	my $comments;
	if (($transcript_start == 1.e12) || ($transcript_end == -1.e12)) {
	    print ERRORFILE "cDNAs associated with $gene_id not found in BEST.GFF\n" if $output_dir; ### if none of the cDNAs associated with a CDS are BEST matches.
	    next;
	}
		
	if(($transcript_start > $gene_start)&&($transcript_end < $gene_end)) {
	    $comments = "cds $gene_id spans transcript";
	}
	else {
	    $comments = "$gene_id";
	}

	if ($gene_strand eq "+") {
	    $five_start  = $transcript_start;
	    $five_end    = $gene_start - 1;
	    $three_start = $gene_end +1;
	    $three_end   = $transcript_end;
	    if ($five_start < $gene_start) {
		print OUTFILE "CHROMOSOME_$gene_chrom\tfind_utrs\t5_UTR\t$five_start\t$five_end\t.\t$gene_strand\t.\t$comments\n" if $output_dir;
		print_ace("5_UTR", $gene_id, $min_cdna, $gene_chrom, $five_start, $five_end) if $output_dir;
	    }
	    if ($three_end > $gene_end) {
		print OUTFILE "CHROMOSOME_$gene_chrom\tfind_utrs\t3_UTR\t$three_start\t$three_end\t.\t$gene_strand\t.\t$comments\n" if $output_dir;
		print_ace("3_UTR", $gene_id, $max_cdna, $gene_chrom, $three_start, $three_end) if $output_dir;
	    }
	}   
	
	if ($gene_strand eq "-") {
	    $five_start  = $gene_end+1;
	    $five_end    = $transcript_end;
	    $three_start = $transcript_start;
	    $three_end   = $gene_start -1;
	    if ($five_end > $gene_end) {
		print OUTFILE "CHROMOSOME_$gene_chrom\tfind_utrs\t5_UTR\t$five_start\t$five_end\t.\t$gene_strand\t.\t$comments\n" if $output_dir;
		print_ace("5_UTR", $gene_id, $max_cdna, $gene_chrom, $five_end, $five_start) if $output_dir;
	    }
	    if ($three_start < $gene_start) {
		print OUTFILE "CHROMOSOME_$gene_chrom\tfind_utrs\t3_UTR\t$three_start\t$three_end\t.\t$gene_strand\t.\t$comments\n" if $output_dir;
		print_ace("3_UTR", $gene_id, $min_cdna, $gene_chrom, $three_end, $three_start) if $output_dir;
	    }
	}   
    }
}

sub print_ace {
   my ($type,$gene_id,$evidence,$gene_chrom,$five_end,$three_end) = @_;

   print ACEFILE "\nUTR :\t\"$type:$gene_id\"\n";
   print ACEFILE "Matching_CDS\t$gene_id\n";
   if ($ESTacc2ESTname{$evidence} ne "") {
       print ACEFILE "$type\tAccession_evidence\t\"EMBL\"\t\"$ESTacc2ESTname{$evidence}\"\n";
   }
   else {
       print ACEFILE "$type\n";
   }
   print ACEFILE "Species\t\"Caenorhabditis elegans\"\n";
   print ACEFILE "Method\t\"UTR\"\n\n";

   print "CDS $gene_id on chromosome $gene_chrom [$five_end - $three_end]\n\n";

   my $link; my $rel_five_end; my $rel_three_end;
   foreach $link (sort {$link_coordinate_start{$a} <=> $link_coordinate_start{$b}} keys %link_coordinate_start) {
       print "$link is on chromosome $link_chrom{$link} from $link_coordinate_start{$link} to $link_coordinate_end{$link}\n";
       next unless ($link_chrom{$link} eq $gene_chrom); 

       print "Checking $link [$link_chrom{$link} $link_coordinate_start{$link} - $link_coordinate_end{$link}]\n\n";


       next unless ( ($five_end  > $link_coordinate_start{$link}) && ($three_end < $link_coordinate_end{$link}) ); 

       $rel_five_end  = $five_end  - $link_coordinate_start{$link} + 1;
       $rel_three_end = $three_end - $link_coordinate_start{$link} + 1;

       print ACEFILE "Sequence :\t\"$link\"\n";
       print ACEFILE "UTR\t\"$type:$gene_id\"\t$rel_five_end\t$rel_three_end\n\n";
       print ACEFILE "//Found match - UTR maps to $link at $rel_five_end - $rel_three_end [abs: $five_end $three_end]\n\n";

       print "Found match - UTR maps to $link at $rel_five_end - $rel_three_end\n\n";
       last;
   }
} 
      

sub find_extent {
    my ($ests,$est_data,$gene_chrom,$gene_start,$gene_end,$gene_id,$gene_strand) = @_;
    my $current_start = 1.e12;
    my $current_end = -1.e12;
    my $min_cdna;
    my $max_cdna;
    my $i;
    foreach (@{$ests}) {
	unless (defined($est_data->{$_}->[0])) {
	    print ERRORFILE"$_ is not found in BLAT_EST_BEST or BLAT_mRNA_BEST\n" if $output_dir;
	    next;
	}
	unless(defined($est_data->{$_}->[3])) {
	    print ERRORFILE "$_ has not been allocated a method\n" if ($output_dir);
	    next;
	}
	unless ($gene_chrom eq $est_data->{$_}->[0]) {
	    print ERRORFILE "$gene_id on chromosome $gene_chrom, best $_ match is on  $est_data->{$_}->[0]\n" if $output_dir;
	    next;
	}
	 
	if ($est_data->{$_}->[3] eq "EST_3") {
	    if($est_data->{$_}->[2] eq  $gene_strand) {
		print ERRORFILE "$_ is an $est_data->{$_}->[3], ($est_data->{$_}->[2]) but appears on the same strand as $gene_id ($gene_strand)\n" if ($output_dir);
		next;
	    }
	}
	elsif ($est_data->{$_}->[3] eq "EST_5") {
	    if($est_data->{$_}->[2] ne  $gene_strand) {
		print  ERRORFILE "$_ is an $est_data->{$_}->[3], ($est_data->{$_}->[2]) but appears on the other strand to $gene_id ($gene_strand)\n" if ($output_dir);
		next;
	    }
	}
	elsif ($est_data->{$_}->[3] eq "NDB") {
	    if($est_data->{$_}->[2] ne $gene_strand) {
		print ERRORFILE "$_ is an $est_data->{$_}->[3], ($est_data->{$_}->[2]) but appears on the other strand to $gene_id ($gene_strand)\n" if ($output_dir);
		next;
	    }
	} 
	
	for($i=0; $i<scalar(@{$est_data->{$_}->[1]}); $i++) {
	#    print "$_\t $est_data->{$_}->[0] $est_data->{$_}->[1]->[$i]->[0]\t $est_data->{$_}->[1]->[$i]->[1]\t $i\n";
	    if ((($gene_start - $est_data->{$_}->[1]->[$i]->[0]) > 10000) || ($est_data->{$_}->[1]->[$i]->[1] - $gene_end) >10000) {
		print ERRORFILE "cDNA $_\t$est_data->{$_}->[0]\t$est_data->{$_}->[1]->[$i]->[0]\t$est_data->{$_}->[1]->[$i]->[1]\tis matched to $gene_id\t$gene_start\t$gene_end\n" if $output_dir;
		next;
	    }
	    
	    if ($current_start > $est_data->{$_}->[1]->[$i]->[0]) {
		$current_start = $est_data->{$_}->[1]->[$i]->[0];
		$min_cdna = $_;
	    }
	    if ($current_end < $est_data->{$_}->[1]->[$i]->[1]) {
		$current_end = $est_data->{$_}->[1]->[$i]->[1];
		$max_cdna = $_;
	    }
	}
    }
    return ($current_start, $current_end, $min_cdna, $max_cdna);
}

sub read_gff {
    my $gff_dir = shift;
    my $gff_type = shift;
    my @chromosomes = qw(I II III IV V X);
    my $rh;
    if (@_) {
	print "second pass\n";
	$rh = shift;
    }
    else {
	print "first pass\n";
    }
	
    print "Reading $gff_type.\n";
    foreach(@chromosomes) {
	open (GFFFILE, "<$gff_dir/CHROMOSOME_$_\.$gff_type") || die "Cannot open $gff_type file $gff_dir/CHROMOSOME_$_\.$gff_type\n";
	while (<GFFFILE>) {
	    chomp;
	    next if /^\#/;
	    my ($chrom, $curated, $sequence, $start, $end, $score, $strand, $frame, $comments) = split /\t/;
	    
	    if ($curated eq "curated") {
		$comments =~ s/.*"(.*)".*/$1/;
	
		$chrom =~ s/CHROMOSOME_//;
		
		$rh->{$comments}->[0] = $chrom;
		$rh->{$comments}->[1] = $start;
		$rh->{$comments}->[2] = $end;
		$rh->{$comments}->[3] = $strand;
	    }
	    elsif (($curated eq "BLAT_EST_BEST")||($curated eq "BLAT_mRNA_BEST")) {
		$comments =~ s/.*"(.*)".*/$1/;
		$comments =~ s/.*://;
		$chrom =~ s/CHROMOSOME_//;
		$rh->{$comments}->[0] = $chrom;
		$rh->{$comments}->[2] = $strand;
		
		if (defined(@{$rh->{$comments}->[1]})) {
		    my @coords = ($start, $end);
		    push (@{$rh->{$comments}->[1]}, \@coords);
		}
		else {
		    $rh->{$comments}->[1]->[0]->[0] =$start;
		    $rh->{$comments}->[1]->[0]->[1] =$end;
		}
	    }
	    else {
		#	print "is this a trna $curated $comments\n";
		next;
	    }
	}
	close (GFFFILE);
    }
    print "Finished reading $gff_type.\n";
    return $rh;
}



    

sub parse_table {
    $ENV{'ACEDB'} = shift;
    my $query_def = shift;
    my %cds_est = ();
    my $cds_est = \%cds_est;
    my $est_data;
    my $q = 0;
    my $command = "Table-maker -p $query_def\nquit\n";
    open (TACE, "echo '$command' | $tace | ") || die "Cannot query acedb. $command  $tace\n";
    while (<TACE>) {
	chomp;
	s/acedb\>//g;
	next if ($_ eq "");
	next if (/\/\//);
	next if (/^\*/); ###removes lines starting with * 
	next if (/@/);   ###removes lines with email addresses (@)
	next if !(/\d/); ###ACHTUNG - this removes all lines with no numbers in them.
	s/\"//g;
	(/^(\S+)\s/);
	my ($est, $cds, $five_tag, $three_tag, $method) = split /\t/;
	
	if (defined(@{$cds_est->{$cds}})) {
	    push (@{$cds_est->{$cds}}, $est);
	}
	else {
	    $cds_est->{$cds}->[0] = $est;
	}
	if ($method eq "NDB") {
	    $est_data->{$est}->[3] = $method;
	}
	elsif ($method eq "EST_elegans") {
	    $est_data->{$est}->[3] = $five_tag if ($five_tag); 
	    $est_data->{$est}->[3] = $three_tag if ($three_tag); 
	    
	    unless(($five_tag)xor($three_tag)) {
		print ERRORFILE  "Check 3 and 5 tags on $est\n";
		next;
	    }
	}
    }
    print "Tablemaker query complete.\n";
    close(TACE);
    return (\%cds_est, $est_data);
}

sub printtest {
    my $cds_est = shift;
    my ($cds, $est, $i);
    while (($cds, $est) = each(%{$cds_est})) {
	print "$cds";
	if (scalar(@{$cds_est->{$cds}})>1) {
	    for ( $i=0; $i<scalar(@{$cds_est->{$cds}}); $i++) {
		print  "\t$cds_est->{$cds}->[$i]\n";
	    }
	}
	else {
	    print "\t$cds_est->{$cds}->[0]\n";
	}
    }
}
    

sub log_file {
    my ($file,$dbdir) = @_;
    my $WS = get_wormbase_version($dbdir);
    $WS = "WS666" if ($test);
    open (LOG, ">>$file/resultslog.txt")|| die "Cannot open LOG.txt\n";
    print LOG "\n\n\n##################################################\n";
    print LOG "$file/CHROMOSOME_n_.transcripts.gff  CREATED BY $0 \n\n";
    print LOG "Connecting to WS$WS\nperl $0 -database $database -output_dir $output_dir\n\nON ";
    print LOG  `date`;
    print LOG "BY ";
    print LOG `whoami`;
    close(LOG);
}



sub sort_gff {
    
    my (@a, @n, @s, @e, @st, @f, $f, $i, $fh);
    
    open(OUTFILE, "<$output_dir/ALL.UTRs.tmp")|| die "Cannot open outfile $output_dir\n";

    open(ONEFILE,   ">$output_dir/CHROMOSOME_I.UTRs.gff")   || die "Cannot open final file in  $output_dir\n";
    open(TWOFILE,   ">$output_dir/CHROMOSOME_II.UTRs.gff")  || die "Cannot open final file in  $output_dir\n";
    open(THREEFILE, ">$output_dir/CHROMOSOME_III.UTRs.gff") || die "Cannot open final file in  $output_dir\n";
    open(FOURFILE,  ">$output_dir/CHROMOSOME_IV.UTRs.gff")  || die "Cannot open final file in  $output_dir\n";
    open(FIVEFILE,  ">$output_dir/CHROMOSOME_V.UTRs.gff")   || die "Cannot open final file in  $output_dir\n";
    open(XFILE,     ">$output_dir/CHROMOSOME_X.UTRs.gff")   || die "Cannot open final file in  $output_dir\n";
    
    while (<OUTFILE>) {
	s/\#.*//;
	next unless /\S/;
	@f = split /\t/;
	push @a,  $_;
	push @n,  $f[0];
	push @s,  $f[3];
	push @e,  $f[4];
	push @st, $f[6];
    }
    
    foreach $i (sort { $n[$a] cmp $n[$b] or  $st[$a] cmp $st[$b] or $s[$a] <=> $s[$b] or $e[$a] <=> $e[$b] } 0..$#a) {
	print ONEFILE   "$a[$i]" if ($a[$i] =~ /CHROMOSOME_I\t/);
	print TWOFILE   "$a[$i]" if ($a[$i] =~ /CHROMOSOME_II\t/);
	print THREEFILE "$a[$i]" if ($a[$i] =~ /CHROMOSOME_III\t/);
	print FOURFILE  "$a[$i]" if ($a[$i] =~ /CHROMOSOME_IV\t/);
	print FIVEFILE  "$a[$i]" if ($a[$i] =~ /CHROMOSOME_V\t/);
	print XFILE     "$a[$i]" if ($a[$i] =~ /CHROMOSOME_X\t/);
    }     
    close(OUTFILE);
    close(ONEFILE);
    close(TWOFILE);
    close(THREEFILE);
    close(FOURFILE);
    close(FIVEFILE);
    close(XFILE);
}

sub get_wormbase_version {
  my $dbdir = shift;
  my $WS_version = `grep "NAME WS" $dbdir/wspec/database.wrm`;
  chomp($WS_version);
  $WS_version =~ s/.*WS//;    
  return($WS_version);
}






__END__;


// Spread sheet definition for the ACeDB software 
// User: ah3
// Date: 2002-08-13_14:36:19

// %n (%%n in the graphic) are parameter to be given on the command line in tace
// or by default by the Parameters given in this file
// \%n (%n in the graphic) are substituted by the value of column n at run time
// Line starting with // are ignored, starting with # are comments

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Sequence 
From 1 
Condition (Method = "EST_elegans") OR (Method = "NDB")
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class Sequence 
From 1 
Tag Matching_CDS  
 
Colonne 3 
Width 12 
Optional 
Visible 
Show_Tag 
From 1 
Tag EST_5  
 
Colonne 4 
Width 12 
Optional 
Visible 
Show_Tag 
From 1 
Tag EST_3  
 
Colonne 5 
Width 12 
Optional 
Visible 
Class 
Class Method 
From 1 
Tag Method 
 
 

// End of these definitions

