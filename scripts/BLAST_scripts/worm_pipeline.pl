#!/usr/local/bin/perl

# ms2@sanger.ac.uk

use strict;
use Getopt::Long;
use DBI;
use Bio::PrimarySeq;
use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::DBSQL::Protein::DBAdaptor;
#use vars qw($opt_f $opt_o $opt_n $opt_m,);

my ($fasta, $database, $brig );

$|=1;

GetOptions (
	    'fasta:s'    => \$fasta,
	     'brig'       => \$brig
	   );

my $usage = "worm_pipeline.pl\n";
$usage .= "-fasta [fasta file of protein sequences to be added to the mysql database]\n";
#$usage .= "-n to match the newly added proteins to swall, using pmatch\n";
#$usage .= "-o to re-match the proteins already in the database to swall, using pmatch (very slow)\n";

unless ($fasta) {
    die "$usage";
}

###############################
# define some variables

my $host = "ecs1f";
my $dbname = $brig ? "worm_brigprot" : "wormprot"; #if -brig opt, use worm_brigprot else use wormprot
my $user = "wormadmin";
my $pass = "worms";
my $species = $brig ? 2 : 1; #assign no based on species

my $pmatch_wrapper = "pmatch.pl";
my $scratchdir = "/tmp";

#################################################################
# datetime subroutine (returns 2001-06-14 11:42:04)

sub now {
    return sprintf ("%04d-%02d-%02d %02d:%02d:%02d",
                    sub {($_[5]+1900, $_[4]+1, $_[3], $_[2], $_[1], $_[0])}->(localtime));
}

#################################################################
# connect to the database, and get some adaptors

my $db = Bio::EnsEMBL::Pipeline::DBSQL::Protein::DBAdaptor->new ( -host   => $host,
                                                                  -dbname => $dbname,
                                                                  -user   => $user,
                                                                  -pass   => $pass,
                                                                );                                                              
my $proteinAdaptor = $db->get_ProteinAdaptor;


#################################################################
# read the sequence fasta file, map the proteins to SWALL,      #
# and populate the protein, peptide and InputIdAnalysis tables  #
#################################################################

print "----------------------------------------------------------\n";
print "starting the WORM PROTEIN analysis pipeline\n";
print "\tprocessing file $fasta\n";
print "\tconnected to $dbname at $host as $user\n";
print "----------------------------------------------------------\n\n";

#################################################################
# old proteins (pmatch + updating swall)

print "check out the proteins already in the DB [".&now."]\n";


my @peptides = $proteinAdaptor->fetch_all_Peptide;

my %exists;
if (@peptides) {
    print "-> there are ".scalar(@peptides)." in the database...\n";
    %exists = map { ($_->id, 1) } @peptides;
}
else {
    print "-> no proteins in the database...\n";
}

##############################################################
# new proteins (pmatch)

my %name2swall;


#################################################################################
# populate the protein, peptide and InputIdAnalysis tables
# (insert with analysisId 1 into the InputIdAnalysis table)
# (corresponds to logic_name 'SubmitProtein')

print "populate the protein, peptide and InputIdAnalysis tables with the new proteins [".&now."]\n";

my $stream = Bio::SeqIO->new ( -file   => $fasta,
                               -format => 'fasta',
                             );
$stream->moltype ('protein');
my $count = 0;
my $count_swall = 0;
while (my $pep = $stream->next_primary_seq) {
    if (exists $exists{$pep->id}) {
        print "\t".$pep->id." already exists in the database...\n";
        next;
    }
    if (exists $name2swall{$pep->id}) {
        $pep->accession_number ($name2swall{$pep->id});
        $count_swall++;
    }
    # !! the organism is currently hard-coded to 1 (= Caenorhabditis elegans in organism table)
    $proteinAdaptor->submit ($pep, $species);
    print "  ".$pep->id."\n";
    $count++;
}

print "-> $count proteins were added to the database ($count_swall map to a swall entry)...\n";

#################################################################################
# run the RuleManager script




#############################

$db->disconnect;


####################################
#   run script to check that whats in the database is whats in wormpep
#$fasta = new_entries.WS$WS_versio
my $version = $fasta =~ /new_entries\.WS(\d+)/;
my $wp_file = "/wormsrv2/WORMPEP/wormpep$version/wp.fasta$version";
system("perl5.8.0 ~/scripts/BLAST_scripts/mysql_vs_db_protein_compare.pl -fasta $wp_file");










