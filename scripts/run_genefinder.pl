#!/usr/local/bin/perl5.8.0
#
# run_genefinder.pl
#
##

my $release = shift;                                         # Release number 

my $accession;                                               # Accession for the current genome sequence
our $dna = "/tmp/genefinder.dna";                            # temp file name used 
our $cmid = "./allcmid";                                     # cmid file to utilise

# Copy current allcmid file from /wormsrv2/autoace
system ("scp wormsrv2:/wormsrv2/autoace/allcmid ./allcmid"); # copy allcmid file 

# Touch output file
system ("touch ./${release}_genefinder.ace");

open (AGP, "<$cmid");
while (<AGP>) {

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
close AGP;


exit(0);
