package Check_Classes;
use strict;
use vars qw (@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS @TotalClasses $VERSION $class $f_found $s_found $diff);
use Exporter;

# Usage: \&check_classes ("old_db","new_db");

$VERSION=1.00;
@ISA=qw(Exporter);
@EXPORT_OK=('check_classes');

format STDOUT_TOP =
                 Classes Content Comparison

Class Name           First Db     Second Db    First -> Second   
--------------------------------------------------------------
.

format STDOUT =
@<<<<<<<<<<<<<<<<<<<<@<<<<<<<<<<<<@<<<<<<<<<<<<@<<<<<<<<<<<<
$class,              $f_found,    $s_found,    $diff
.          

sub check_classes {
   my ($first,$second,$fh) = @_;
   my ($prog,$name);
   my $name=`/bin/uname -sr`; 
   if    ($name=~/^SunOS/) {($prog)=<~wormpub/acedb/ace4/bin.SUN_4/tace>;}
   elsif ($name=~/^IRIX/)  {($prog)=<~wormpub/acedb/ace4/bin.SGI_4/tace>;}
   elsif ($name=~/^OSF/)   {($prog)=<~acedb/RELEASE.SUPPORTED/bin.ALPHA_4/tace>;}
   elsif ($name=~/^Linux/) {($prog)=<~wormpub/acedb/ace4/bin.LINUX/tace>;}
   else {print "No known binary for $name\n";exit;}
   my $exec_f = "$prog $first";
   my $exec_s = "$prog $second";

# The complete set of classes to dump is between the DATA and END tokens
READARRAY: while (<DATA>) {
    chomp $_;
    last READARRAY if $_ =~ /END/;
    push (@TotalClasses,$_);
}

foreach (@TotalClasses) {
    $class=$_;
    my $command=<<EOF;
query find $class
quit
EOF
    open (FIRSTACE, "echo '$command' | $exec_f  | ");
    while (<FIRSTACE>) {
	if ($_ =~ /Found\s+(\d+)\s+objects/){
            $f_found=$1;
	    chomp $f_found;
	}
    }
    close FIRSTACE;
    open (SECONDACE, "echo '$command' | $exec_s  | ");
    while (<SECONDACE>) {
	if ($_ =~ /Found\s+(\d+)\s+objects/){
            $s_found=$1;
	}
    }
    close SECONDACE;
    $diff = $s_found - $f_found;
    $diff.="\n";
    write ;
}
}

1;

__DATA__
2_point_data
Accession_number
Allele
Author
Balancer
Cell
Cell_group
Class
Clone
Comment
Contig
Database
Display
EMBL_dump_info
EMBL_info
Enzyme
Expr_Pattern
Gene_Class
Jade
Journal
Keyword
Laboratory
Life_stage
Locus
ManyMap
Map
Metabolite
Method
Motif
MultiMap
Multi_counts
Multi_pt_data
Neurodata
Paper
Pathway
PathwayDiagram
Picture
Pos_neg_data
Protein
ReactantInfo
Rearrangement
Reference
RegulatorInfo
Repeat_Info
Sequence
Session
Species
Strain
Table
Table_definition
Tag
Text
Tree
TreeNode
Url
UserSession
View
View_tags
WWW_server
__END__



