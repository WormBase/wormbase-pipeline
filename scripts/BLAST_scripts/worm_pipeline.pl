#!/usr/local/bin/perl

# ms2@sanger.ac.uk

use strict;
use Getopt::Std;
use DBI;
use Bio::PrimarySeq;
use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::DBSQL::Protein::DBAdaptor;
use vars qw($opt_f $opt_o $opt_n $opt_m);

$|=1;

getopts ("f:onm:");

my $usage = "worm_pipeline.pl\n";
$usage .= "-f [fasta file of protein sequences to be added to the mysql database]\n";
$usage .= "-n to match the newly added proteins to swall, using pmatch\n";
$usage .= "-o to re-match the proteins already in the database to swall, using pmatch (very slow)\n";

unless ($opt_f) {
    die "$usage";
}

###############################
# define some variables

my $host = "ecs1f";
my $dbname = "wormprot";
$dbname = "$opt_m" if $opt_m;
my $user = "wormadmin";
my $pass = "worms";

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
print "\tprocessing file $opt_f\n";
print "\tconnected to $dbname at $host as $user\n";
print "----------------------------------------------------------\n\n";

#################################################################
# old proteins (pmatch + updating swall)

print "check out the proteins already in the DB [".&now."]\n";


my @peptides = $proteinAdaptor->fetch_all_Peptide;

my %exists;
if (@peptides) {
    print "-> there are ".scalar(@peptides)." in the database...\n";
    if ($opt_o) {
        # print the peptides to a fasta file
        open (TMP, ">$scratchdir/$$.fasta") || die "cannot create tmp file\n";
        foreach my $pep (@peptides) {
            $exists{$pep->id} = 1;
            print TMP ">".$pep->id."\n".$pep->seq."\n";
        }
        close TMP;
        # run pmatch
        print "map the proteins already in the DB to swall, updating the previous mapping [".&now."]\n";
        my $cmd = "$pmatch_wrapper -q $scratchdir/$$.fasta -w -o $scratchdir/$$.pmatch -c";
        my $exit_value = system "$cmd";
        unless ($exit_value == 0) {
             die "couldn't make first system call";
        }
        # parse pmatch, updating the swallId's at the same time
        print "make swallId updates [".&now."]\n";
        open (FH1 , "$scratchdir/$$.pmatch") || die "cannot read $scratchdir/$$.pmatch\n";
        while (<FH1>) {
            chomp;
            my @a = split /\t/;
            if ($a[2] eq "match" || $a[2] eq "repeat") {
                $proteinAdaptor->update_swallId ($a[0], $a[1]);
            }
        }
        close FH1;
        unlink "$scratchdir/$$.pmatch";
        unlink "$scratchdir/$$.fasta";
    }
    else {
        foreach my $pep (@peptides) {
            $exists{$pep->id} = 1;
        }
    } 
}
else {
    print "-> no proteins in the database...\n";
}

##############################################################
# new proteins (pmatch)

my %name2swall;

if ($opt_n) {
    print "map the new proteins to swall [".&now."]\n";

    my $cmd = "$pmatch_wrapper -q $opt_f -w -o $scratchdir/$$.pmatch_new -c";
    my $exit_value = system "$cmd";
    unless ($exit_value == 0) {
         die "couldn't make second system call";
    }
    # parse pmatch
    open (FH2 , "$scratchdir/$$.pmatch_new") || die "cannot read $scratchdir/$$.pmatch_new\n";
    while (<FH2>) {
        chomp;
        my @a = split /\t/;
        if ($a[2] eq "match" || $a[2] eq "repeat") {
            $name2swall{$a[0]} = $a[1];
        }
    }
    close FH2;
    unlink "$scratchdir/$$.pmatch_new";
}

#################################################################################
# populate the protein, peptide and InputIdAnalysis tables
# (insert with analysisId 1 into the InputIdAnalysis table)
# (corresponds to logic_name 'SubmitProtein')

print "populate the protein, peptide and InputIdAnalysis tables with the new proteins [".&now."]\n";

my $stream = Bio::SeqIO->new ( -file   => $opt_f,
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
    $proteinAdaptor->submit ($pep, 1);
    print "  ".$pep->id."\n";
    $count++;
}

print "-> $count proteins were added to the database ($count_swall map to a swall entry)...\n";

#################################################################################
# run the RuleManager script




#############################

$db->disconnect;

#################################################################
# run the RuleManager script,
# which creates and run the appropriate jobs
#################################################################













