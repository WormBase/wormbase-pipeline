#!/usr/local/bin/perl5.8.0 -w

use Ace;

$campath = "/nfs/disk100/wormpub/acedb/ace4/cam";
$autopath = "/nfs/wormdata1/wormdb2";


%camace = "";
%autoace = "";

# make gene list for camace

print "Connect to database : camace\n";

open (CAMLIST, ">camlist");
$db = Ace->connect(-path=>"$campath") ||
    die "failed to connect to database\n";

print "Connection successful\n";
print "Dump CDS objects from camace\n";

$it = $db->fetch_many(-query=>'find elegans_CDS');
while ($obj=$it->next){

    # discard genewise predictions

    if ($obj =~ /\.gw/) {
	next;
    }

    # extract data

    $source = $obj->Sequence;
    $method = $obj->Method;

    print "$obj\n";
    $camace{$obj}="$method $source";
}
close (CAMLIST);

# make gene list for autoace

print "Connect to database : autoace\n";

open (CAMLIST, ">autolist");
$db = Ace->connect(-path=>"$autopath") ||
    die "failed to connect to database\n";

print "Connection successful\n";
print "Dump CDS objects from autoace\n";

$it = $db->fetch_many(-query=>'find elegans_CDS');
while ($obj=$it->next){


    # discard GSC predictions

    $method = $obj->Method;
    if ($method eq "stl") {
	next;
    }

    $source = $obj->Sequence;

    print "$obj\n";
    $autoace{$obj}="$method $source";
}
close (CAMLIST);

print "Compare the lists\n";

system ("cat camlist | sort > camlist.sorted");
system ("cat autolist | sort > autolist.sorted");

system ("comm -13 camlist.sorted autolist.sorted > intersection_autolist");
system ("comm -23 camlist.sorted autolist.sorted > intersection_camlist");

print "\Genes in camace but not autoace\n";

open (MISSING, "intersection_camlist");
while (<MISSING>) {
    chomp;
    $gene = $_;
    ($source,$method) = split (/\s+/, $camace{$gene});
    print "\nSequence : $gene\nSequence\t$source\nMethod\t$method\n";
}
close (MISSING);

print "\Genes in autoace but not camace\n";

open (MISSING, "intersection_autolist");
while (<MISSING>) {
    chomp;
    $gene = $_;
    ($source,$method) = split (/\s+/, $camace{$gene});
    print "\nSequence : $gene\nSequence\t$source\nMethod\t$method\n";
}
close (MISSING);

# and tidy up




exit;

