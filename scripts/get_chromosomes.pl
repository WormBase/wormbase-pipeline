#!/usr/local/bin/perl
# dl
# 991118
#
# hacked from one of Steve's subroutines


$file = shift;
%chrome=&mapchromes2;

open (FILE, "$file") || die "Can't open sequence file $file\n\n"; 
while (<FILE>) {

    if (/Sequence : \"(\S+)\"/) {
	$gene = $1;
	($project,$num) = split (/\./, $gene);

#	print "Sequence : \"$gene\"\nChromosome \"$chrome{$project}\"\n\n";

	if ($chrome{$project} eq "") {
	    print "Sequence : \"$gene\"\tChromosome \"$chrome{$project}\"\n\n";
	}


    }
}
close (FILE);

exit;



sub mapchromes2 {
    local($command);
    local(%chrome);
    undef %chrome;
        local($map);
    $command=<<EOF;
    find genome_sequence *
    show -a Interpolated_gMap
    quit
EOF

        open(TEXTACE, "echo '$command' | tace /nfs/disk100/wormpub/acedb/ace4/autoace | ");
    while (<TEXTACE>) {
                #print;
                if (/^Sequence\s+:\s+\"(\w+)\"/)  {
                        $sequence=$1;
                }
                if (/^Interpolated_gMap\s+\"(\S+)"\s+/)  {
                        $map=$1;
                        #print "Assigning $sequence to $map\n";
                        $chrome{$sequence}=$map;
        }
    }
    close TEXTACE;     
        #print  "AH6 is on chromosome $chrome{'AH6'}\n";
    return %chrome;
}

