#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  check_genes.pl
#
#        USAGE:  ./check_genes.pl -database DATABASE_NAME -outputfile FILE_NAME [-dna|-transcripts_only]
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:   (), <>
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  11/07/06 17:49:26 BST
#     REVISION:  ---
#===============================================================================

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use IO::File;

my ($database,$dna,$transcript_only,$outfile);
GetOptions( 'database=s'=>\$database,'dna'=>\$dna,'transcript_only'=>\$transcript_only, 'outfile=s' => \$outfile)||die(@!);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host   => 'farmdb1',
        -user   => 'ensro',
        -dbname => $database,
        -pass   => '',
        -port   => 3306,
    );


my $fh = new IO::File "> $outfile" ||die("cannot open $outfile\n");

my $gene_adaptor = $db->get_GeneAdaptor();
my @genes = @{$gene_adaptor->fetch_all()};
foreach my $gene(@genes){
	my $geneId=$gene->stable_id();
	foreach my $trans (@{$gene->get_all_Transcripts()}) {
		my $protein=$trans->translation();
		my $transcriptId=$trans->stable_id();
                my $proteinId=$protein->stable_id();
		my $geneID=$gene->stable_id();
		my $seq;
		if($dna){
			$seq=$trans->translateable_seq()
		}elsif($transcript_only){
			$seq=$trans->seq->seq()
		}else{
			$seq=$protein->seq()
		}
		printf $fh ">$proteinId\t$transcriptId\t$geneId\n%s\n",reformat($seq);
	}
}

sub reformat {
	my ($seq)=@_;
	my @bases = split(//,$seq);
	my $count=1;
	my $reformatted;
	foreach my $b(@bases){
		$reformatted.=$b;
                $reformatted.="\n" if ($count++ % 60 == 0 && $count < scalar(@bases));
	}
	return $reformatted;
}
