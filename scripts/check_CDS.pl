#!/usr/local/bin/perl
#
# check_CDS.pl v1.0
# dl
# 2000-04-26
#
# Aceperl script to check database gene models. 
#
###########################################################################################
#
# 000426 dl : PP version


########################################
# iniatialise                          #
########################################

use Ace;
BEGIN {
  unshift (@INC,"/nfs/disk92/PerlSource/Bioperl/Releases/bioperl-0.05");
}
use Bio::Seq;

########################################
# command-line parsing                 #
########################################

while ($ARGV[0] =~ /^-/) {
    $_=shift;
    if (/^-c(.*)/) {
       $ace=1;
    }
    elsif (/^-s(.*)/) {
        $ace=2;
    }
    elsif (/^-a(.*)/) {
        $ace=3;
    }
    elsif (/^-d(.*)/) {
        $debug=1;
    }
    elsif (/^-w(.*)/) {
        $html=1;
    }
    else {
        &usage;
    }
}

########################################
# usage subroutine                     #
########################################

sub usage {
    print "Usage: sweeper.pl [-options]/n/n";
    print "Options:\n";
    print "-c   Use Cambridge data  \n";
    print "-s   Use St Louis data   \n";
    print "-a   Use Autoace data  \n";
    print "-d   Debug/Verbose mode\n\n";
    exit;
}

########################################
# Connect with acedb database          #
########################################

my $stlacepath="/nfs/disk100/wormpub/acedb/ace4/stl";
my $camacepath="/nfs/disk100/wormpub/acedb/ace4/cam";
my $autoacepath="/nfs/wormdata1/wormdb2";

$|=1;

if ($ace == 1) {
    if ($debug == 1) {print  "Opening camace ....\n";}
    $db = Ace->connect(-path=>$camacepath) || do { print "Connection failure: ",Ace->error; die();};
}
elsif ($ace == 2) {
    if ($debug == 1) {print  "Opening stlace ....\n";}
    $db = Ace->connect(-path=>$stlacepath) || do { print "Connection failure: ",Ace->error; die();};
}
elsif ($ace == 3) {
    if ($debug == 1) {print  "Opening autoace ....\n";}
    $db = Ace->connect(-path=>$autoacepath) || do { print "Connection failure: ",Ace->error; die();};
}
if ($debug) {print "Connection OK.\n\n";}

########################################
# Main Loop                            #
########################################


#if ($debug) {
#    $count = $db->fetch(-query=> 'find Sequence where Method = curated');
#    print "checking $count curated sequences\n\n";
#}

#print "Gene\tLength (bp)\tCheckSum\tSource\n";


$i = $db->fetch_many(-query=> 'find Sequence where Method = curated');  
while ($obj = $i->next) {
    $gene = $obj;

    $CDS_seq = $obj->asDNA();
    $CDS_seq=~tr/a-z/A-Z/;
    $CDS_seq=~s/\>\w+//;
    $CDS_seq=~s/\W+//mg;
    if ($CDS_seq =~ /[^ACGTUMRWSYKVHDBXN]/img) {
	$CDS_seq=~s/[^ACGTUMRWSYKVHDBXN]//img;
    }
    $CDS_len = length($CDS_seq);
    
    $source  = $obj->Source;
    
    $wormpep = $obj->Corresponding_protein();

    if ($wormpep eq "") {
	$wormpep = "n/a       ";
    }


    $bioseq = Bio::Seq->new(-seq=>$CDS_seq,-id=>$gene,-ffmt=>'Fasta',-type=>'Dna',);
    $chksum=$bioseq->GCG_checksum;
    undef $bioseq;

    print "$obj    \t$CDS_len\t$chksum\t$wormpep\t$source";

#    if (length ($source) < 7) {
#	print "     ";
#    }
#    print "\t";

#    foreach $protein_id ($obj->Protein_id) {
#	($parent,$protein_id_prefix,$protein_id_suffix) = $protein_id->row();
#	print "[$parent -> $protein_id_prefix.$protein_id_suffix] ";
#    }

    print "\n";

    $obj->DESTROY();

} # end object loop    

exit;
