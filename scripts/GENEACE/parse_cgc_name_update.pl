#!/software/bin/perl -w

#3 WBGene00022629 ZC513.5 algn-12 WBPerson451
#3 WBGene00007043 T24D1.4 algn-10 WBPerson451 tag-179 
#2 WBGene00045189 B0284.6 srpr-1.1

use strict;
use Getopt::Long;
use lib $ENV{'CVS_DIR'};

my $filein;
my $fileout;
my $who;

GetOptions (
    'filein=s'         => \$filein,
    'fileout=s'        => \$fileout,
    'who=s' => \$who,
    );

unless ($who =~ /WBPerson\d+/) {
    die "No/Incorrect WBPerson specified with -who option\n";
    close (IN);
    close (OUT);
    exit(0);
}


my @f;
unless (defined$fileout) {
  $fileout = "${filein}_out";
}

open (OUT, ">$fileout");

print OUT "//Scripted output file.\n";

open (IN, "<$filein") or die("Cannot find file $filein\n");
while (<IN>) {
    chomp;
    @f=split" ";
    unless ($f[1] =~ /WBGene\d{8}/) {next;}
    print OUT "\n\nGene : \"$f[1]\"\nVersion $f[0]\nCGC_name $f[3]";
    if (defined $f[4]){
	#Paper evidence
	if (/WBPaper/) {
	    print OUT " Paper_evidence $f[4]";
	}
	#Person evidence
	elsif (/WBPerson/) {
	    print OUT " Person_evidence $f[4]";
	}
        #laboratory
        elsif ((/\w{2}/) || (/\w{3}/)) {
	    print OUT " Laboratory_evidence $f[4]";
	}
    }
print OUT "\nPublic_name $f[3]\nVersion_change $f[0] now $who Name_change CGC_name $f[3]";

    if (defined $f[5]) {
	if($f[5] =~ /\S+/) {
	    print OUT "\nOther_name $f[5]\nVersion_change $f[0] now $who Name_change Other_name $f[5]";
	} 
    }
    if ($f[3] =~ /(\w+)/) {
	print OUT "\nGene_class $1";
    }
    if (defined $f[5]) {
	if ($f[5] =~ /(\w+)\-\d+/){
	    print OUT "\n\nGene_class $1\nOld_member $f[5]\n"
	}
    }
}

close (IN);
close (OUT);

exit(0);

__END__
