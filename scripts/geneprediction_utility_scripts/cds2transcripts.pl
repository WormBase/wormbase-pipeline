#!/usr/bin/env perl

while (<>){
   s/CDS :/Transcript :/;
   next if /^CDS|^Evidence|^Remark|^Predicted|Start_not_found|End_not_found|Genetic_code|^Brief_identification/;
   s/Corresponding_transcript/Corresponding_CDS/;
   s/curated/Coding_transcript/;
   print;
}
