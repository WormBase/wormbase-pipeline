#!/usr/local/bin/perl5.8.0
#
# run_genefinder.pl
#
# Usage: run_genefinder.pl -WSxxx
#
# Last edited by: $Author: dl1 $
# Last edited on: $Date: 2005-07-11 13:16:57 $
 


my $release = shift;                                         # Release number 

my $accession;                                               # Accession for the current genome sequence
our $dna = "/tmp/genefinder.dna";                            # temp file name used 
our $cmid = "./allcmid";                                     # cmid file to utilise

# Copy current allcmid file from /wormsrv2/autoace
system ("scp wormsrv2:/wormsrv2/autoace/allcmid ./allcmid"); # copy allcmid file 

# Touch output file
system ("touch ./${release}_genefinder.ace");

open (CMID, "<$cmid");
while (<CMID>) {

    # parse file name f

    if (/^>(\S+)/) {
	$accession = $1;
	
	# get the DNA
	system ("dbfetch.pl $cmid $accession > $dna");
	
	print "Running genefinder for $accession\n";
	
	# run genefinder
	system("/nfs/disk100/wormpub/BioSW/gf_980506/genefinder/bin/genefinder  -tablenamefile /nfs/disk100/wormpub/BioSW/gf_980506/nemtables/sanger.tablenamefile.cod -intronPenalty /nfs/disk100/wormpub/BioSW/gf_980506/nemtables/intron_penalty.lookup -exonPenalty /nfs/disk100/wormpub/BioSW/gf_980506/nemtables/exon_penalty.lookup $dna > ./OUTPUT/$accession.gf");
	
	# parse to acefile
	system ("genefinder2ace.pl ./OUTPUT/$accession.gf $accession > ./OUTPUT/$accession.ace");
	
	# and push top the full file
	system ("cat ./OUTPUT/$accession.ace >> ./${release}_genefinder.ace");

    }   
    
}
close CMID;

system ("scp ${release}_genefinder.ace wormsrv2:/wormsrv2/wormbase/misc_dynamic/misc_genefinder.ace");

exit(0);
