package check_CDS_tace;
use strict;
use vars qw (@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION $command $gene @genes $CDS_pointer $dna_counter $CDS_seq $bioseq );
use Exporter;
$VERSION=1.00;
@ISA=qw(Exporter);
@EXPORT_OK=('check_CDS');

BEGIN {
  unshift (@INC,"/nfs/disk92/PerlSource/Bioperl/Releases/bioperl-0.05");
}
use Bio::Seq;
$|=1;

sub check_CDS {
$ENV{'ACEDB'} = shift(@_);
my %CDS_len = "";
my %CDS_chksum = "";
my %source = "";
my %protein = "";
my $exec="";
my $name=`uname -sr`;
if    ($name=~/^SunOS/) {($exec)=<~wormpub/acedb/ace4/bin.SUN_4/tace>;}
elsif ($name=~/^IRIX/)  {($exec)=<~wormpub/acedb/ace4/bin.SGI_4/tace>;}
elsif ($name=~/^OSF/)   {($exec)=<~acedb/RELEASE.SUPPORTED/bin.ALPHA_4/giface>;}
elsif ($name=~/^Linux/) {($exec)=<~wormpub/acedb/ace4/bin.LINUX/tace>;}
else {print STDERR "check_CDS_tace: No known binary for $name\n";exit;}
$command=<<EOF;
query find Sequence where Method = \"curated\" 
show -a 
DNA
quit
EOF

open (textace, "echo '$command' | $exec  | ");
while (<textace>) {
   chomp;
   if (/Sequence : \"(\S+)\"/)  {
	$gene = $1;
	push (@genes, $gene);
    }
    if (/Source\s+\"(\S+)\"/)  {
	$source{$gene} = $1;
    }
    if (/Corresponding_protein\s+\"(\S+)\"/)  {
	$protein{$gene} = $1;
    }
    if (/acedb> >(\S+)/) {
	$CDS_pointer = $1;
	$dna_counter = 1;
	$CDS_seq = ">$1\n";
	next;
    }
    if (($dna_counter == 1) && ((/^>/) || (/object dumped/))) {
	$CDS_seq=~tr/a-z/A-Z/;
	$CDS_seq=~s/\>\w+//;
	$CDS_seq=~s/\W+//mg;
	if ($CDS_seq =~ /[^ACGTUMRWSYKVHDBXN]/img) {
	    $CDS_seq=~s/[^ACGTUMRWSYKVHDBXN]//img;
	}
	$CDS_len{$CDS_pointer} = length($CDS_seq);
	$bioseq = Bio::Seq->new(-seq=>$CDS_seq,-id=>$CDS_pointer,-ffmt=>'Fasta',-type=>'Dna',);
	$CDS_chksum{$CDS_pointer}=$bioseq->GCG_checksum;
	undef $bioseq;

	(/>(\S+)/);
	$CDS_seq = ">$1\n";
	$CDS_pointer = $1;
	next;
    }
    if ($dna_counter == 1) {
	$CDS_seq .= $_;
	next;
    }
}
close (textace);

foreach $gene (@genes) {
    if ($protein{$gene} eq "") {
	$protein{$gene} = "n/a     ";
    }
    print "$gene\t$CDS_len{$gene}\t$CDS_chksum{$gene}\t$protein{$gene}\t$source{$gene}\n";
}
}

1;




