package WormBaseConf;

use strict;
use vars qw( %WormBaseConf );


%WormBaseConf = (
		 #location of agp and gff files
		 WB_CHR_INFO => [
				 { 
				  chr_name 	=> 'I',
				  agp_file 	=> 'CHROMOSOME_I.agp',
				  length 	=> '15072418',#lengths can be found at top of *.gff file
				  gff_file 	=> 'CHROMOSOME_I.gff',
				  dna_file 	=> 'CHROMOSOME_I.dna',
				 },
				 { 
				  chr_name 	=> 'II',
				  agp_file 	=> 'CHROMOSOME_II.agp',
				  length 	=> '15279314',
				  gff_file 	=> 'CHROMOSOME_II.gff',
				  dna_file 	=> 'CHROMOSOME_II.dna',
				 },
				 { 
				  chr_name 	=> 'III',
				  agp_file 	=> 'CHROMOSOME_III.agp',
				  length 	=> '13783677',
				  gff_file 	=> 'CHROMOSOME_III.gff',
				  dna_file 	=> 'CHROMOSOME_III.dna',
				 },
				 { 
				  chr_name 	=> 'IV',
				  agp_file 	=> 'CHROMOSOME_IV.agp',
				  length 	=> '17493785',
				  gff_file 	=> 'CHROMOSOME_IV.gff',
				  dna_file 	=> 'CHROMOSOME_IV.dna',
				 },
				 { 
				  chr_name 	=> 'V',
				  agp_file 	=> 'CHROMOSOME_V.agp',
				  length 	=> '20919396',
				  gff_file 	=> 'CHROMOSOME_V.gff',
				  dna_file 	=> 'CHROMOSOME_V.dna',
				 },
				 { 
				  chr_name 	=> 'X',
				  agp_file 	=> 'CHROMOSOME_X.agp',
				  length 	=> '17718851',
				  gff_file 	=> 'CHROMOSOME_X.gff',
				  dna_file 	=> 'CHROMOSOME_X.dna',
				 },
				 { 
				  chr_name 	=> 'MtDNA',
				  agp_file 	=> 'CHROMOSOME_MtDNA.agp',
				  length 	=> '13794',	#sequence-region CHROMOSOME_MtDNA 1 13794
				  gff_file 	=> 'CHROMOSOME_MtDNA.gff',  
				 }
			       ],

		 # agp type which gets written to database meta
		 WB_AGP_TYPE 			=> 'exon',
		 WB_NEW_COORD_SYSTEM_VERSION 	=> 'CEL160',  #CEL150 ???
		 WB_workDIR 			=> '/acari/work2a/wormpipe/ensembl/',  #/ecs/work/username/WormBase150/
		 WB_cvsDIR  			=> '/nfs/acari/wormpipe/ensembl/',  #/nfs/acari/user/cvsDIR/
		 WB_scratchDIR			=> '/acari/scratch1/worm/',  #/ecs/scratch/user/WormBase150/
		
		 
		 # date that this new version was made
		 WB_DDAY    => '0605', #0511  ie.November 2005

		 # database to put sequnece and genes into
		 WB_DBNAME  => 'worm_WB160',
		 WB_DBHOST  => 'ecs1f',
		 WB_DBUSER  => 'wormadmin',
		 WB_DBPASS  => 'XXX',
		 WB_DBPORT  => '3306',

		 # if want the debug statements in wormbase to ensembl scripts printed
		 WB_DEBUG   => 1,

		 # location to write file containing dodgy seq ids
		 WB_SEQ_IDS => '/nfs/acari/wormpipe/ensembl/dodgy_ids',  #ecs/work/user/WormBase150/dodgy_seq_ids.txt

		 # location to write ids of genes which don't translate
		 WB_NON_TRANSLATE => '/nfs/acari/wormpipe/ensembl/dodgy_ids', #/ecs/work/user/WormBase150/non_translateable.txt

		 #coordinate-system / base feature names
		 WB_CLONE_SYSTEM_NAME      => 'clone',
		 WB_CHROMOSOME_SYSTEM_NAME => 'chromosome',

		 # logic name of analysis objectsto be parsed out of the gff
		 WB_LOGIC_NAME             => 'wormbase',
		 WB_OPERON_LOGIC_NAME      => 'Operon',
		 WB_rRNA_LOGIC_NAME	   => 'rRNA',
		 WB_RNAI_LOGIC_NAME        => 'RNAi',
		 WB_EXPR_LOGIC_NAME        => 'Expression_profile',
		 WB_PSEUDO_LOGIC_NAME      => 'Pseudogene',
		 WB_TRNA_LOGIC_NAME        => 'tRNA',
		 #WB_SL1_LOGIC_NAME        => 'SL1',
		 #WB_SL2_LOGIC_NAME        => 'SL2',
		 WB_ce_mrna_LOGIC_NAME     => 'celegans_mrna',
		 WB_ce_est_LOGIC_NAME      => 'celegans_est',
		 WB_other_EST_LOGIC_NAME   => 'other_est',
		 WB_fosmid_LOGIC_NAME      => 'fosmid',
		);


sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else
    # all of GeneConf:
    my @vars = @_ ? @_ : keys( %WormBaseConf );
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;

    foreach (@vars) {
	if ( defined $WormBaseConf{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$WormBaseConf{ $_ };
	} else {
	    die "Error: WormBaseConf: $_ not known\n";
	}
    }
}

1;
