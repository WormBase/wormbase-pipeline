#!/usr/local/bin/perl

$clone = shift;
$mincontig = 1000;
$first = 0;

open (DNA, "echo 'dump2fasta -depad -clone $clone' | /usr/bin/csh - |");	
$dumped="no";
$/=">";
    
while (<DNA>) {
    if (/\S+/) {
	$dumped="yes";
    }
    s/\>//;
    if (/^\d+/) {
	/^(\d+)/;
	$contig=$1;
	s/^\d+.+\n//;
    }
    
    #only want to print contig details if contig => mincontig
    $length=length($_);
    $sequence=$_;
    if ( $length  >= $mincontig){
	$moddate=&time2ace(time);
	print ">$clone Contig ID=$contig; Length=$length; Order=Unknown; Status=Cambridge-Unfinished; Author=$author; LastModification=$moddate; EMBL acc=;\n";
	print $sequence;
    }
}

#Did the dump2fasta fail?  
if ($dumped eq "no") {
    print "dump of $clone unsuccessful\n";
} 

#$/="\n";
close DNA;


##################################
#get the date of a file modifiction
##################################
sub mod {
    local(@date);
    @date=stat $_[0];
    return $date[9];
}

##################################
sub time2ace {
    local($time) = @_;
    local(@lt);

    @lt = localtime($time);

    # yy -> yyyy
    if ($lt[5] < 70) { $lt[5] += 2000; }
    elsif ($lt[5] < 100) { $lt[5] += 1900; }

    sprintf("%4d-%02d-%02d_%02d:%02d:%02d",$lt[5],$lt[4]+1,$lt[3],$lt[2],$lt[1],$lt[0]);

}
##################################

