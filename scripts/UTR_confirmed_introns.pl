#!/usr/local/bin/perl5.6.1 -w

use strict;
use lib glob("~ar2/wormbase/scripts/");
use Ace;
use Getopt::Long;
use Data::Dumper;

my ($debug,$update,$track,$datadump, $database, $gff_dir, $specified_genes);

GetOptions ( "debug" => \$debug,
	     "track" => \$track,
	     "datadump" => \$datadump,
	     "update" => \$update,
	     "database:s" => \$database,
	     "gff_dir:s"  => \$gff_dir,
	     "specify"    => \$specified_genes
	   );

my @chromosomes = qw(CHROMOSOME_I CHROMOSOME_II CHROMOSOME_III CHROMOSOME_IV CHROMOSOME_V CHROMOSOME_X);

$database = glob("~wormpub/DATABASES/WS98_allele_fix") unless $database;
$gff_dir = glob("~wormpub/DATABASES/WS98_allele_fix") unless $gff_dir;
my $acefile = "$database/utr.ace";

open (ACE,">$acefile") or die "cant write $acefile\n";

# this bit is to separate wrongly attached cDNAs to centre specific acefile output
my $STL_bad_mcDNAs = "$database/STL_bad_matching_cDNA.ace";
my $CAM_bad_mcDNAs = "$database/CAM_bad_matching_cDNA.ace";
open (STLBAD,">$STL_bad_mcDNAs") or die "cant open $STL_bad_mcDNAs : $!\n";
open (CAMBAD,">$CAM_bad_mcDNAs") or die "cant open $CAM_bad_mcDNAs : $!\n";

my %BAD_FH = ( 'HX' => *CAMBAD,
	       'RW' => *STLBAD
	     );

my %gene2centre;
&reFetchData("gene2centre",\%gene2centre);

print "\nAce->connect $database\t" if $track;
my $db = Ace->connect(-path => "$database") or die Ace->error,"\n";
print "DONE\n" if $track;

# sample GFF
# C43G2   curated         Sequence     10841   13282   .       +       .       Sequence "C43G2.4"
# C43G2   BLAT_EST_BEST   similarity   8877    9029    100     +       .       Target "Sequence:yk1125e01.5" 37 189
# C43G2   *UNKNOWN*       intron       3582    3934    .       -       .       Confirmed_by_EST ; Confirmed_in_UTR

# get genes to output from command line or entire list
my @genes_to_dump;
if ($specified_genes ) {
  my $g;
  while ($g = shift) { push(@genes_to_dump,$g) };
}

my %genes_span;
my %cDNA;
my @conf_intr;
my $superlink;

foreach ( @chromosomes ) {

  my $gff = "$gff_dir/$_.gff";
  $superlink = $_;
  %genes_span = ();
  %cDNA = ();
  @conf_intr = ();
  @genes_to_dump = () unless $specified_genes;

  # parse GFF file to get Sequence, exon and cDNA info
  if ($update ) {
    open (GFF,"<$gff") or die "gff\n";
    print "updating data . . \n" if $track;
    while (<GFF>) {
      next if (/\#\#/);
      my @data = split(/\t/,$_);
      #print " . . \t@data\n" ;
      if ( (defined $data[1]) and ($data[1] ne ".")) {
	if ($data[1] eq "curated") {
	  if( defined $data[8] ) {
	    $data[8] =~ /Sequence \"(.*)\"/;
	    $data[8] = $1;
	  }
	  if( $data[2] eq "Sequence" ) {
	    $genes_span{$data[8]} = [($data[3], $data[4], $data[6])];
	  }
	}
	elsif( ($data[1] eq "BLAT_EST_BEST" ) or ($data[1] eq "BLAT_mRNA_BEST" ) ){
	  $data[8] =~ /Target \"Sequence:(.*)\"/ ;
	  $data[8] = $1;
	  $cDNA{$data[8]}[2] = $data[6]; # strand
	  if( !(defined $cDNA{$data[8]}[0]) or ($cDNA{$data[8]}[0] > $data[3]) ) {
	    $cDNA{$data[8]}[0] = $data[3];
	  }
	  if( !(defined $cDNA{$data[8]}[1]) or ($cDNA{$data[8]}[0] < $data[4]) ) {
	    $cDNA{$data[8]}[1] =  $data[4];
	  }
	}
      }
      elsif( ($data[2] eq "intron") and ( $data[8] =~ /Confirmed_in_UTR/ ) ) {
	push( @conf_intr, [($data[3], $data[4]) ] );
      }
    }

    close GFF;

    if( $datadump ) {
      print "writing data  . . \n" unless $track;
      open (FH,">$database/cdna.dat") or die;
      print FH Data::Dumper->Dump([\%cDNA]);
      close FH;

      open (FH,">$database/genes.dat") or die;
      print FH Data::Dumper->Dump([\%genes_span]);
      close FH;

      open (FH,">$database/conf_int.dat") or die;
      print FH Data::Dumper->Dump([\@conf_intr]);
      close FH;
    }
  }
  else {
    print "loading data  . . \n" unless $track;
    &reFetchData("cdna",\%cDNA);
    &reFetchData("genes",\%genes_span);
    &reFetchArrayData("conf_int",\@conf_intr);
  }


  # sort the array of cIntrons
  print "sorting confirmed introns\t" if $track;
  my @sorted_cis = sort { $a->[0] <=> $b->[0] } @conf_intr;
  print "DONE\n" if $track;
  # get matching_cDNA s for genes

  # process each cDNA
  my $gene_count;

  print "\n\nStarting Genes loop\n" if $track;

  # if no specific gene to dump do all of them
  unless ( defined $specified_genes ) {
    foreach (keys %genes_span ) {
      push (@genes_to_dump,$_);
    }
  }

  foreach my $gene (@genes_to_dump ) {
    print "doing gene $gene\n" if $debug;
    my $utr5start;
    my $utr3end;
    my (%utr5, %utr3);

    my $seq = $db->fetch(-name    => "$gene",
			 -class   => 'Sequence'
			);
    if( $seq ) {
      my @matching_cDNA = $seq->Matching_cDNA;
      my @extremes = ($genes_span{$gene}[0],$genes_span{$gene}[1] );

      # determine extremes of transcript based on mcDNAs
      foreach my $mcDNA (@matching_cDNA) {
	if ( $cDNA{$mcDNA} ) {
	  &FindExtremes(\@extremes, \$mcDNA, \$gene);
	}
	elsif($debug) { print $gene," -> ",$mcDNA->name," not in GFF file\n"; }
      }

      # find overlapping introns
      foreach my $ci (@sorted_cis) {
	last if (&Check_Introns(\$genes_span{$gene}[2],
		       \$gene ,
		       \@extremes ,
		       \$ci ,
		       \%utr5 ,
		       \%utr3
		      ) eq "stop" );
      }
      # write UTR output (if there is one)
      &fwd5UTR_introns(\$gene, \@extremes, \%utr5) if ( ($extremes[0] != $genes_span{$gene}[0]) and ( $genes_span{$gene}[2] eq "+")) ;
      &fwd3UTR_introns(\$gene, \@extremes, \%utr3) if ( ($extremes[1] != $genes_span{$gene}[1]) and ( $genes_span{$gene}[2] eq "+")) ;
      &rev5UTR_introns(\$gene, \@extremes, \%utr5) if ( ($extremes[1] != $genes_span{$gene}[1]) and ( $genes_span{$gene}[2] eq "-")) ;
      &rev3UTR_introns(\$gene, \@extremes, \%utr3) if ( ($extremes[0] != $genes_span{$gene}[0]) and ( $genes_span{$gene}[2] eq "-")) ;
    }
  }
}
print "\n\nEND\n" if $track;

close ACE;
close STLBAD;
close CAMBAD;
$db->close;
exit(0);

sub fwd5UTR_introns
  {
    my $gene = shift;
    my $extremes = shift;
    my $utr5 = shift;

    print ACE "\nSequence : $superlink\n";
    print ACE "UTR \"UTR_5:$$gene\" ",$$extremes[0]," ",$genes_span{$$gene}[0]-1," \n";

    # link UTR to CDS
    print ACE "\nSequence : $$gene\nMatching_UTR \"UTR_5:$$gene\"\n";

    print ACE "\nUTR :\"UTR_5:$$gene\"\n";
    print ACE "Method UTR\n";
    print ACE "Source_exons 1 ";

    foreach (sort keys %$utr5) {
      if($$utr5{$_} +1 == $genes_span{$$gene}[0] ) {
	print "$$gene - UTR intron without exon\n";
	next;
      }
      print ACE $_ - $$extremes[0] ,"\nSource_exons ",$$utr5{$_} +1 - $$extremes[0] + 1," ";
    }
    print ACE $genes_span{$$gene}[0] - $$extremes[0] ,"\n";
  }

sub rev5UTR_introns
  {
    my $gene = shift;
    my $extremes = shift;
    my $utr5 = shift;

    # assing UTR to chromosome
    print ACE "\nSequence : $superlink\n";
    print ACE "UTR \"UTR_5:$$gene\" ",$$extremes[1]," ",$genes_span{$$gene}[1]+1," \n";

    # link UTR to CDS
    print ACE "\nSequence : $$gene\nMatching_UTR \"UTR_5:$$gene\"\n";

    # Enter UTR details
    print ACE "\nUTR :\"UTR_5:$$gene\"\n";
    print ACE "Method UTR\n";
    print ACE "Source_exons 1 ";

    foreach (sort {$$utr5{$b}<=>$$utr5{$a}} keys %$utr5) {
      print ACE $$extremes[1] - ($$utr5{$_}),"\nSource_exons ", ($$extremes[1] +1) - ($_ - 1)," ";
    }
    print ACE $$extremes[1] - $genes_span{$$gene}[1] ,"\n";
  }
sub fwd3UTR_introns 
  {
    my $gene = shift;
    my $extremes = shift;
    my $utr3 = shift;

    print ACE "\nSequence : $superlink\n";
    print ACE "UTR \"UTR_3:$$gene\" ",$genes_span{$$gene}[1]+1," ",$$extremes[1]," \n";

    # link UTR to CDS
    print ACE "\nSequence : $$gene\nMatching_UTR \"UTR_3:$$gene\"\n";

    print ACE "\nUTR :\"UTR_3:$$gene\"\n";
    print ACE "Method UTR\n";
    print ACE "Source_exons 1 ";

    foreach (sort keys %$utr3) {
      print ACE ($_)  - ($genes_span{$$gene}[1] + 1) ,"\nSource_exons ",($$utr3{$_} +1) -($genes_span{$$gene}[1])," ";
    }
    print ACE $$extremes[1] - $genes_span{$$gene}[1] ,"\n";
  }

sub rev3UTR_introns 
  {
    my $gene = shift;
    my $extremes = shift;
    my $utr3 = shift;

    print ACE "\nSequence : $superlink\n";
    print ACE "UTR \"UTR_3:$$gene\" ",$genes_span{$$gene}[0]-1," ",$$extremes[0]," \n";

    # link UTR to CDS
    print ACE "\nSequence : $$gene\nMatching_UTR \"UTR_3:$$gene\"\n";

    print ACE "\nUTR :\"UTR_3:$$gene\"\n";
    print ACE "Method UTR\n";
    print ACE "Source_exons 1 ";

    foreach (sort {$$utr3{$b}<=>$$utr3{$a}} keys %$utr3) {
      print ACE ($genes_span{$$gene}[0] - 1) - $$utr3{$_} ,"\nSource_exons ",($genes_span{$$gene}[0] ) - ($_ - 1)," ";
    }
    print ACE $genes_span{$$gene}[0] - $$extremes[0] ,"\n";
  }
sub Check_Introns
  {
    my $strand = shift;
    my $gene = shift;
    my $extremes = shift;
    my $ci = shift;
    my $utr5 = shift;
    my $utr3 = shift;

    if( ( $genes_span{$$gene}[0] > $$ci->[0] ) and ( $$ci->[0] > $$extremes[0] ) ) {
      $$utr5{$$ci->[0]} = $$ci->[1] if $$strand eq "+";
      $$utr3{$$ci->[0]} = $$ci->[1] if $$strand eq "-";
    }
    elsif( ( $genes_span{$$gene}[1] < $$ci->[1] ) and ( $$ci->[1] < $$extremes[1] ) ) {
      $$utr3{$$ci->[0]} = $$ci->[1] if $$strand eq "+";
      $$utr5{$$ci->[0]} = $$ci->[1] if $$strand eq "-";
    }
    # confirmed intron start after the extreme of gene - no point checking any further
    elsif( $$ci->[0] > $$extremes[1] ) {
      return "stop";
    }
  }

sub FindExtremes
 {
   my $extremes = shift;
   my $mcDNA = shift;
   my $gene = shift;

   if( ( ($$extremes[0] - $cDNA{$$mcDNA}[0]) > 20000 ) || 
       ( ($cDNA{$$mcDNA}[1] - $$extremes[1]) > 20000) )  {
     # this determines which file to print to via %BAD_FH based on the gene.
     print {$BAD_FH{ $gene2centre{$$gene} } } "\nSequence : \"$$gene\"\n-D Matching_cDNA \"$$mcDNA\"\n";
   }
   else {
     if( $cDNA{$$mcDNA}[0] < $$extremes[0] ) {
       $$extremes[0] = $cDNA{$$mcDNA}[0];
     }	
     if( $cDNA{$$mcDNA}[1] > $$extremes[1] ) {
       $$extremes[1] = $cDNA{$$mcDNA}[1];
     }
   }
 }


sub reFetchData {

    my ($file,$ref) = @_;
    my $dir = "/nfs/disk100/wormpub/DATABASES/WS98_allele_fix";
    open (FH, "<$dir/$file.dat") or die "can't open $dir/$file.dat";
    undef $/;
    my $VAR1;
    my $data = <FH>;
    eval $data;
    die if $@;
    $/ = "\n";
    close FH;
    %$ref = (%$VAR1);
}

sub reFetchArrayData {

    my ($file,$ref) = @_;
    my $dir = "/nfs/disk100/wormpub/DATABASES/WS98_allele_fix";
    open (FH, "<$dir/$file.dat") or die "can't open $dir/$file.dat";
    undef $/;
    my $VAR1;
    my $data = <FH>;
    eval $data;
    die if $@;
    $/ = "\n";
    close FH;
    @$ref = (@$VAR1);
}
