#!/usr/local/bin/perl5.8.0 -w

# further_BLAST_dump.pl
#
#
# Author: Chao-Kung CHen
#
# Last updated by: $Author: ar2 $                      
# Last updated on: $Date: 2003-08-01 10:11:51 $        

system("/usr/apps/bin/scp /acari/work2a/wormpipe/dumps/blastp_ensembl.ace wormsrv2:/wormsrv2/wormbase/ensembl_dumps/");
system("/usr/apps/bin/scp /acari/work2a/wormpipe/dumps/blastx_ensembl.ace wormsrv2:/wormsrv2/wormbase/ensembl_dumps/");
system("/usr/apps/bin/scp /acari/work2a/wormpipe/dumps/wormprot_motif_info.ace wormsrv2:/wormsrv2/wormbase/ensembl_dumps/");
system("/usr/apps/bin/scp /acari/work2a/wormpipe/dumps/swissproteins.ace wormsrv2:/wormsrv2/wormbase/ensembl_dumps/");
system("/usr/apps/bin/scp /acari/work2a/wormpipe/dumps/tremblproteins.ace wormsrv2:/wormsrv2/wormbase/ensembl_dumps/");
system("/usr/apps/bin/scp /acari/work2a/wormpipe/dumps/ipi_hits.ace wormsrv2:/wormsrv2/wormbase/ensembl_dumps/");
system("/usr/apps/bin/scp /acari/work2a/wormpipe/dumps/best_blastp_hits wormsrv2:/wormsrv2/wormbase/ensembl_dumps/");
system("/usr/apps/bin/scp /acari/work2a/wormpipe/dumps/brigprot_blastp_ensembl.ace wormsrv2:/wormsrv2/wormbase/ensembl_dumps/");

chdir "/wormsrv2/wormbase/ensembl_dumps";

system("cat ipi_hits.ace fly.ace yeast.ace swissproteins.ace tremblproteins.ace brigpep.ace > ensembl_protein_info.ace");

exit(0);
