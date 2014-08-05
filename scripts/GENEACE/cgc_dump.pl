#!/bin/env perl
# script to create lists for the CGC
# usage: 
#    cgc_dump.pl -database DATABASE_DIRECTORY -labs -geneclass -dir OUTDIR
#    
# it will create a timestamped directory in OUTDIR with dumps and diffs versus
# a last directory, which is a symlink to the last run 

use Ace;
use Getopt::Long;
use IO::File;

my ($db,$geneclass,$labs);
my $originalDir = '.';

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
$mon++;
my $stamp ="$year$mon$mday";

GetOptions(
  '-database=s'  => \$db,
  '-geneclass'   => \$geneclass,
  '-labs'        => \$labs,
  '-dir=s'       => \$originalDir,
)||die(@!);

my $db = Ace->connect(-path => $db,-program => '/software/acedb/bin/giface')||die(@!);

# current directory
$dir = "$originalDir/$stamp";
mkdir($dir)||die(@!);

if($geneclass){
  my @gene_classes = $db->fetch(Gene_class => '*');
  my $outf = IO::File->new(">$dir/geneclass.txt")||die(@!);
  foreach my $gc (sort @gene_classes){
   printf $outf "\"%s\"\t\"%s\"\t\"%s\"\n",$gc,$gc->Designating_laboratory,join(",",$gc->Description);
  }
  $outf->close;
  diff('geneclass');
}

if ($labs){
  my @laboratories = $db->fetch(-query => 'find Laboratory *;CGC');
  my $outf = IO::File->new(">$dir/Laboratory.txt")||die(@!);
  foreach my $l (sort @laboratories){
     # my $rep = $l->Representative?$l->Representative->Standard_name : '';
     printf $outf "\"%s\"\t\"%s\"\t\"%s\"\t\"%s\"\n",$l,$l->Allele_designation,join(',',map {$_->Standard_name} $l->Representative),$l->Mail;
  }
  $outf->close;
  diff('Laboratory');
}

# update last link
system("rm $originalDir/last") && die(@!);
system("cd $originalDir && ln -s $stamp ./last") && die(@!);

sub diff{
   my ($file)=@_;

   # will run a 4k column wide side by side comparison skipping identical lines
   system("diff --side-by-side --suppress-common-lines --width 4000 $originalDir/last/$file.txt $dir/$file.txt > $dir/$file.diff");
   my $inf  = IO::File->new("<$dir/$file.diff")||die(@!);
   my $outf = IO::File->new(">$dir/$file.diff_")||die(@!);
   while (<$inf>){
     s/^\s+//;         # new
     s/.*\s\|\s/\| /; # changed <- will barf on internal pipes
     s/^>\s+/> /;
     print $outf $_;
   }
   $inf->close;
   $outf->close;
   system("mv -f $dir/$file.diff_ $dir/$file.diff") && die(@!);
}

