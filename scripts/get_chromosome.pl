#!/usr/local/bin/perl
# dl
# 991118
#
# hacked from one of Steve's subroutines


$project=shift;
$exec="tace /nfs/disk100/wormpub/acedb/ace4/autoace";

$command=<<EOF;
find genome_sequence $project
show -a Interpolated_gMap
quit
EOF
    
open(TEXTACE, "echo '$command' | $exec -| "); 
while (<TEXTACE>) {
  if (/^Interpolated_gMap\s+\"(\S+)"\s+/)  {
        print "$project is located on chromosome $1\n";
  }
}
close TEXTACE;     

exit;
