#!/nfs/team71/worm/mh6/bin/perl
# map_Alleles.pl
#
# by Anthony Rogers
#
# This maps alleles to the genome based on their flanking sequences
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2006-08-08 15:11:43 $
# SubVersion :     $Revision: 1.48 $

use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
use Feature_mapper;
use Modules::GFF_sql;
use Sequence_extract;
use Wormbase;
use Ace;
use Getopt::Long;
use Data::Dumper;
use Memoize;
use threads;
memoize('Sequence_extract::Sub_sequence');

my @opts = @ARGV;

#######################################
# command-line options                #
#######################################

my $debug;       # debug mode, output goes to /tmp
my $limit;       # limit number of alleles to map
my $database;    # specify database to map alleles in, default is autoace
my $release;     # specify release number, defaults to current release number
my $help;        # help mode, show perldoc
my $restart = 'go';    # specify an allele name from which to start mapping
my $no_parse;          # turn off loading of data to $database
my $list;              # read in alleles to map from file rather than from database
my $gff_dir;           # option to use a diferent gff_dir
my $verbose;           # verbose mode, extra output to screen
my $store;             # defines stored configuration~/src/perl/mapper/map_Alleles.pl
my $test;
my $outdir;            # acefile output directory
my $noupdate;          # don't update the mysql database

GetOptions(
    'debug=s'    => \$debug,
    'limit=i'    => \$limit,
    'database=s' => \$database,
    'release=i'  => \$release,
    'help'       => \$help,
    'restart=s'  => \$restart,
    'no_parse'   => \$no_parse,
    'list=s'     => \$list,
    'gff'        => \$gff_dir,
    'verbose'    => \$verbose,
    'store=s'    => \$store,
    'test'       => \$test,
    'outdir=s'   => \$outdir,
    'noupdate'   => \$noupdate
);

# check command line options
if ($help) { print `perldoc $0`; exit }
die "No need for -release option in debug mode" if ( $release && $debug );

###########################
# recreate configuration   #
############################
my $wb;

if ($store) { $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wb = Wormbase->new( -debug => $debug, -test => $test, ) }

if ($gff_dir) { $wb->{'gff_splits'} = $gff_dir }

###########################################
# Variables Part II (depending on $wb)    #
###########################################
$test  = $wb->test  if $wb->test;     # Test mode
$debug = $wb->debug if $wb->debug;    # Debug mode, output only goes to one user

my $maintainers = "All";
my $tace        = $wb->tace;
my $error_count = 0;                  # for tracking alleles that don't map

$release = $wb->get_wormbase_version unless ( defined $release );
my $ace_file = ( $outdir ? $outdir : $wb->acefiles ) . "/allele_mapping.WS$release.ace";

$database = $wb->autoace unless $database;    # ACeDB
my $dbh = $debug ? GFF_sql->new() : GFF_sql->new( { -build => 1 } );    # mySQL

# set variables depending on whether in debug (test) mode or not
if ($debug) {
    $maintainers = "$debug\@sanger.ac.uk";
    print "DEBUG = \"$debug\"\n\n";
}

# check if it needs to submit stuff into lsf

################################################
# Set hashes
################################################

my %worm_gene2class  = $wb->FetchData('worm_gene2class');          # set from Common Data
my %worm_gene2geneID = $wb->FetchData('worm_gene2geneID_name');    # has some issues with double dot names
my %to_map;                                                        # will be set if -list option is specified, stores list of alleles to map
my %original_alleles;                                              # will store details of allele->gene connections already in database
my %allele2gene;                                                   # key is allele name, value is array of CDS/Transcript names
my %allele_data;                                                   # allele data from mapping data, hash of arrays:

# allele => [ (0)type,           (1) 5'flank_seq ,      (2) 3'flank_seq, (3) CDS,
#             (4)end of 5'match, (5) start of 3'match , (6) clone,       (7) chromosome,
#             (8)strains]

my $log = Log_files->make_build_log($wb);                          # create log file, open output file handles

my $extractor = Sequence_extract->invoke( $database, undef, $wb );
###############################################################################################
#
#                                   MAIN part of script
#
###############################################################################################

my $pthread;
$pthread = threads->create( 'prep_db', $debug ) unless $noupdate;    # hmm ... takes quite a while if there are alot of gff files
&process_list if ($list);                                            # read in list of alleles to map if -list is specified
&check_original_mappings;                                            # make hash of allele->gene connections already in database
&map_alleles( $dbh, $pthread );                                      # map alleles, the main routine
&load_alleles_to_database unless ($no_parse);                        # load acefile to autoace (unless -no_parse specified)

#####################################
# tidy up, email log, and exit
#####################################

if ( $error_count > 0 ) {
    print "\nERROR: problems with $error_count alleles\n\n";
    $log->write_to("\nProblems with $error_count alleles\n");
    $log->mail( "$maintainers", "BUILD REPORT: map_Alleles.pl $error_count ERRORS!" );
}
else {
    $log->mail( "$maintainers", "BUILD REPORT: map_Alleles.pl" );
}
exit(0);

##################################################
#                                                #
#            S U B R O U T I N E S               #
#                                                #
##################################################

##################################################
# prepare database                               #
# used to update mySQL while accessing the ACeDB #
# returns 1 on success                           #
##################################################

sub prep_db {
    my ($debug) = @_;
    my $tdbh = $debug ? GFF_sql->new() : GFF_sql->new( { -build => 1 } );    # as database handles cannot be shared across threads ...
    my $gffdir = $wb->gff_splits;
    foreach my $c qw(I II III IV V X) {
        $tdbh->clean("CHROMOSOME_${c}");
        foreach my $file ( glob("$gffdir/CHROMOSOME_${c}_*.gff") ) {
            $tdbh->generate_tags($file);
            $tdbh->load_gff( $file, "CHROMOSOME_${c}" );
        }
    }

    #    $tdbh->DESTROY; # might have to be called explicitely if it doesn't close the connection
    return 1;
}

#########################################################
#             Main allele mapping loop                  #
#########################################################

sub map_alleles {
    my ( $gffdb, $prep_thread ) = @_;
    my $rdy = -1;    # MH ready state for the prepare thread

    # First get allele info from database
    my $db = Ace->connect( -path => $database ) || do { print "$database Connection failure: ", Ace->error; die(); };
    my $mapper = Feature_mapper->new( $database, undef, $wb );
    my @alleles = $db->fetch( -query => 'Find Variation WHERE Flanking_sequences AND species = "Caenorhabditis elegans"' );    #very slow :-(

    open( OUT, ">$ace_file" ) or die "cant open $ace_file\n";

  ALLELE:
    foreach my $allele (@alleles) {
        my $name = $allele->name;
        my ( $sequence, $clone, $left, $right );
        my $go = 0;
        my @affects_CDSs;
        my @affects_Feature;
        my %affects_Splice = ( 'Donor' => [], 'Acceptor' => [] );
        my $orientation    = '.';

        print "\nMapping $name\n" if $verbose;

        # debug facility - this bit is so that it can restart from given allele name
        unless ( "$restart" eq "go" ) {
            if ( "$restart" eq "$name" ) { $restart = 'go' }
            else { print "skipping $name\n" if $verbose; next }
        }

        # debug facility - after the restart means you can specify where to start and how many to do
        if ($limit) { last if $error_count++ >= $limit }
        if ($list) { next unless defined $to_map{$name} }

        # check that allele is attached to a valid sequence object
        $sequence = $allele->Sequence;
        if ( !defined $sequence ) {
            $log->write_to("ERROR: $name has missing Sequence tag\n");
            $error_count++;
            next ALLELE;
        }
        elsif ( !defined $sequence->Source ) {
            $log->write_to("ERROR: $name connects to Sequence $sequence which has no Source tag\n");
            $error_count++;
            next ALLELE;
        }
        unless ( $allele->Flanking_sequences and $allele->Flanking_sequences->right ) {
            $log->write_to( "ERROR: ", $allele->name, " is missing a flanking sequence!\n" );
            next ALLELE;
        }

        #map position on genome
        #(0)type,
        #(1)5'flank_seq ,
        #(2)3'flank_seq
        #(3)CDS
        #(4)end of 5'match
        #(5)start of 3'match
        #(6)clone
        #(7)chromosome
        #(8)strains containing this allele

        # get allele mapping
        # grab both flanking sequences from $database
        $left  = lc $allele->Flanking_sequences->name;
        $right = lc $allele->Flanking_sequences->right->name;

        # warn if flanking sequence is missing
        unless ( $left and $right ) {
            $log->write_to("ERROR: $name does not have two flanking sequences\n");
            $error_count++;
            next ALLELE;
        }
        $allele_data{$name}[0] = $allele->Method;
        $allele_data{$name}[1] = $left;
        $allele_data{$name}[2] = $right;

        # map allele using Feature_mapper.pm, store results of map in @map, warn if mapping failed
        my @map = $mapper->map_feature( $sequence->name, $left, $right );

        if ( "$map[0]" eq "0" ) {
            $log->write_to("ERROR: Couldn't map $name with $left and $right to seq $sequence\n");
            $error_count++;
            next ALLELE;
        }

        unless ( ( $allele->Type_of_mutation and $allele->Type_of_mutation eq "Insertion" )
            or ( $allele->Variation_type and $allele->Variation_type eq "Transposon_insertion" ) )
        {

            # get coords of allele
            if ( $map[2] > $map[1] ) {

                # maps to fwd strand
                $map[1]++;
                $map[2]--;
                $orientation = '+';
            }
            else {
                $map[1]--;
                $map[2]++;
                $orientation = '-';
            }
        }

        $allele_data{$name}[4] = $map[1];
        $allele_data{$name}[5] = $map[2];
        $allele_data{$name}[6] = $map[0];

        # retrieve list of overlapping CDSs/Transcripts

        my ( $chromosome, $start, $stop );
        if ( $map[0] =~ /CHROMOSOME_\d+/ ) {
            ( $chromosome, $start, $stop ) = @map;
        }
        else {
            ( $chromosome, $start ) = $mapper->Coords_2chrom_coords( $map[0], $map[1] );
            ( $chromosome, $stop )  = $mapper->Coords_2chrom_coords( $map[0], $map[2] );
        }

        $rdy = $prep_thread->join() if ( ( $rdy != 1 ) && $prep_thread );
        my @hits = $gffdb->get_chr( $chromosome, { 'start' => $start, 'stop' => $stop, } );    # glorious change introducing GFF database

        # WBGene stuff
        push @affects_CDSs, grep { ( $_->{'source'} eq 'coding_exon'               && $_->{'feature'} eq 'Coding_transcript' ) } @hits;
        push @affects_CDSs, grep { ( $_->{'source'} eq 'miRNA_primary_transcript'  && $_->{'feature'} eq 'miRNA' ) } @hits;
        push @affects_CDSs, grep { ( $_->{'source'} eq 'snoRNA_primary_transcript' && $_->{'feature'} eq 'snoRNA' ) } @hits;
        push @affects_CDSs, grep { ( $_->{'source'} eq 'tRNA_primary_transcript'   && $_->{'feature'} eq 'tRNAscan-SE-1.23' ) } @hits;
        push @affects_CDSs, grep { ( $_->{'source'} eq 'snRNA_primary_transcript'  && $_->{'feature'} eq 'snRNA' ) } @hits;
        push @affects_CDSs, grep { ( $_->{'source'} eq 'rRNA_primary_transcript'   && $_->{'feature'} eq 'rRNA' ) } @hits;
        push @affects_CDSs, grep { ( $_->{'source'} eq 'nc_primary_transcript'     && $_->{'feature'} eq 'Non_coding_transcript' ) } @hits;
        push @affects_CDSs, grep { ( $_->{'source'} eq 'intron'                    && $_->{'feature'} eq 'curated' ) } @hits;
        push @affects_CDSs, grep { ( $_->{'source'} eq 'intron'                    && $_->{'feature'} eq 'Coding_transcript' ) } @hits;
        push @affects_CDSs, grep { ( $_->{'source'} eq 'gene'                      && $_->{'feature'} eq 'gene' ) } @hits;

        # same for Features
        push @affects_Feature, grep { ( $_->{'source'} eq 'SL1_acceptor_site' ) } @hits;
        push @affects_Feature, grep { ( $_->{'source'} eq 'SL2_acceptor_site' ) } @hits;
        push @affects_Feature, grep { ( $_->{'source'} eq 'polyA_signal_sequence' ) } @hits;
        push @affects_Feature, grep { ( $_->{'source'} eq 'polyA_site' ) } @hits;


            # make (unique) array of gene names corresponding to @affects_CDSs
            # first make a hash, then use the keys of the hash to set new array

            my %affects_genes;
            foreach my $hit (@affects_CDSs) {
                next unless ref($hit);
                @hits = ( $hit->{'fluff'} =~ /\"([^\s]+)\"/g );
                map { s/\"//g } @hits;
                foreach my $seq (@hits) {    # connect it to anything in the description
                    if ( $hit->{'source'} eq 'gene' ) { $affects_genes{$seq}->{'gene'} = 1; next }

                    ############## ######## #####

                    if ( $seq =~ /(.*\..*)\.(\d+\w*)$/ ) { $seq = $1 }

                    ############# ######### #####

                    my $gene = $worm_gene2geneID{$seq} || print "WARNING: cannot find $seq  in COMMON_DB\n";

                    $worm_gene2geneID{$seq} || next;
                    $affects_genes{$gene}->{ $hit->{'source'} } = 1;
                }
            }
            my @affects_genes = keys(%affects_genes);    # list of WBGeneIDs

            print Dumper (@affects_CDSs) if $debug;      #bad

            # now compare both arrays (allele->gene connections already in database and
            # allele->gene connections made by script) to look for differences
            my %count;
            foreach my $gene ( @{ $original_alleles{$name} }, @affects_genes ) { $count{$gene}++ }
            foreach my $gene ( keys %count ) {

                # allele->gene connections in both arrays will cause $count{$gene} to be == 2
                # if it only == 1 then there is an error
                if ( $count{$gene} == 1 ) {
                    if ( $affects_genes{$gene} ) {

                        # not so serious, but source database should be updated with connection
                        next                                                      if ( $allele->SNP );
                        print "WARNING: $name->$gene was mapped by script only\n" if ($verbose);
                        $log->write_to("WARNING: $name->$gene was mapped by script only\n");
                    }
                    else {

                        # more serious, source database has a bad allele->gene connection
                        print "ERROR: $name->$gene was in original database only\n" if ($verbose);
                        $log->write_to("ERROR: $name->$gene was in original database only\n");
                        $error_count++;
                    }
                }
            }
            $allele2gene{"$name"} = \@affects_CDSs if ( $affects_CDSs[0] );
            &outputAllele( $name, $chromosome, $start, $stop );
        if ( @affects_Feature >= 1 ) {
            &outputFeature( $name, \@affects_Feature );
            print Dumper (@affects_Feature) if $debug;
        }
    }

    # close database and file handles
    $db->close;
    close OUT;
}

#####################################################################################################

sub outputAllele {
    my ( $allele, $chromosome, $start, $stop ) = @_;
    my %done;

    # only proceed if mapping routine returns the following three pieces of info
    if ( $allele_data{$allele}[6] && $allele_data{$allele}[4] && $allele_data{$allele}[5] ) {

        #allele-sequence
        print OUT "\nSequence : \"$allele_data{$allele}[6]\"\nAllele $allele $allele_data{$allele}[4] $allele_data{$allele}[5]\n\n";

        # process allele if it matches overlapping CDSs/Transcripts
        if ( $allele2gene{$allele} ) {    # ???
            my @affects_CDSs = @{ $allele2gene{$allele} };

            # loop through all CDSs/Transcripts that overlap allele
            foreach my $gff_cds (@affects_CDSs) {
                my @affects = $gff_cds->{'fluff'} =~ /\"[^\s]+\"/g;
                map { s/\"//g } @affects;
                foreach my $cds (@affects) {

                    ############## ######## #####

                    #if ( $cds =~ /([\w\.]*)\.\d+$/ ) { $cds = $1 }

                    ############# ######### #####

                    #allele - WBGene connection

                    # crude hack
                    if ( $gff_cds->{'source'} eq 'gene' ) {
                        next if $done{$cds} && $done{$cds}->{'gene'};    #don't print duplicates
                        $done{$cds}->{'gene'} = 1;

                        print OUT "Variation : $allele\n";
                        print OUT "Gene $cds\n\n";
                        next;
                    }

                    my $WBGene = $worm_gene2geneID{$cds};
                    next if ( !$WBGene );

                    my $type = $gff_cds->{'source'};
                    $type =~ s/^(\w)/\u$1/;

                    # now write relevant acefile output depending on whether allele hits CDS or Transcript
                    if ( !defined( $worm_gene2class{$cds} ) ) {
                        $log->write_to("ERROR: $cds not listed in %worm_gene2class\n");
                        $error_count++;
                        next;
                    }

                    elsif ( $worm_gene2class{$cds} eq "CDS" ) {
                        my ( $splice, $frameshift,$substitution);
                        $splice = &outputSplice( $start, $stop, $gff_cds, $chromosome ) if ( $type eq 'Intron' );
                        $frameshift = &outputFrameshift( $start, $stop, $allele_data{$allele}[0] )
                          if ( $type eq 'Coding_exon'
                            && ( $allele_data{$allele}[0] eq 'Insertion_allele' || $allele_data{$allele}[0] eq 'Deletion_allele' ) && ( abs($start - $stop) < 3 ) );
			$substitution = &outputSubstitution( $start, $stop)
                          if ( $type eq 'Coding_exon' && $allele_data{$allele}[0] eq 'Substitution_allele' && ( $start - $stop <= 3 ) );


                        print OUT "Variation : $allele\n";
                        print OUT "Predicted_CDS $cds $type Inferred_automatically map_Alleles.pl\n";
                        print OUT "Predicted_CDS $cds $frameshift" if $frameshift;
                        print OUT "Predicted_CDS $cds $splice" if $splice;
			print OUT "\n";
                    }
                    elsif ( $worm_gene2class{$cds} eq "Transcript" ) {
                        next if $done{$cds} && $done{$cds}->{ $gff_cds->{'source'} };              #don't print duplicates
                        print OUT "Variation : $allele\n";
                        $type = ' ' if $type =~ /.*_primary_transcript/;                           #crude hack
                        print OUT "Transcript $cds $type\n\n";
                        $done{$cds}->{ $gff_cds->{'source'} } = 1;
                    }
                    else {
                        $log->write_to("ERROR: $cds is not a CDS or Transcript\n");
                        $error_count++;
                    }
                }
            }
        }
    }
}

sub outputSubstitution {
	my ($start,$stop,$from,$top)=@_;
	return;
}


##########################
# to print frameshifts
##########################
sub outputFrameshift {
    my ( $start, $stop, $type ) = @_;
    my $size=abs($start-$stop);
    $type=~s/_allele//;
    return "Frameshift \"($size bp $type)\" Inferred_automatically map_Alleles.pl\n";
}

######################################################################################
# to print the feature stuff
######################################################################################
sub outputFeature {
    my ( $allele, $features ) = @_;
    foreach my $feature ( @{$features} ) {
        $feature->{'fluff'} =~ /^\w+\s"([^\s]+)"/;
        my $featureid = $1;
        print OUT "Feature : $featureid\n";
        print OUT "Associated_with_variation $allele Inferred_automatically map_Alleles.pl\n\n";
        print OUT "Variation : $allele\nAffects Feature $featureid\n\n";
    }
}
################################
# to print splice
###############################
sub outputSplice {
    my ( $type, $sstart, $sstop );
    my ( $start, $stop, $intron, $chromosome ) = @_;

    # donor
    return undef unless ( $start >( $intron->{start} - 11) && $stop < ($intron->{start} + 11) )
      || ( $start > $intron->{stop} - 11 && $stop < $intron->{stop} + 11 );

    if ( $start <= $intron->{start} + 1 && $stop >= $intron->{start} ) {
        $type   = $intron->{orientation} eq '+' ? 'Donor' : 'Acceptor';
        $sstop  = $intron->{start} + 1;
        $sstart = $intron->{start};
    }
    elsif ( $start <= $intron->{stop} && $stop >= $intron->{stop} - 1 ) {
        $type   = ( $intron->{orientation} eq '-' ) ? 'Donor' : 'Acceptor';
        $sstop  = $intron->{stop};
        $sstart = $intron->{stop} - 1;
    }
    else { return undef }
    my $site = &get_seq( $extractor, $chromosome, $sstart, $sstop, $intron->{orientation} );
    return "Splice_site $type $site Inferred_automatically map_Alleles.pl\n";
}

####################################################################################
# small subroutine to load data into %to_map hash from file if -list is specified
####################################################################################

sub process_list {
    open( LIST, "<$list" ) || die "Cannot open file specified by $list\n";
    while (<LIST>) { chomp; $to_map{$_} = 1 }
    close(LIST);
}

#################################################
# get sequence from chromosome through ACeDB
# get_seq($extractor,'I',10,100,'+')
sub get_seq {
    my ( $extractor, $chromosome, $pos1, $pos2, $orientation ) = @_;
    my $seq = $extractor->Sub_sequence($chromosome);
    my $len = length $seq;

    # convert to computer coords
    $pos1--;
    $pos2--;
    my ( $p1, $p2 ) = sort { $a <=> $b } ( $pos1, $pos2 );

    # are we in reverse sense? (i.e. reversed order of positions)
    my $bla = substr( $seq, $p1, $p2 - $p1 + 1 );
    $bla = $extractor->DNA_revcomp($bla) if ( $orientation eq '-' );
    return $bla;
}

##############################################################################################
# get list of original allele->gene connections from source database
#
# this will be a hash of arrays, each array element will be a Gene ID that corresponds
# to an allele (deletion alleles can affect multiple genes).  This hash will be used
# later on to compare the output of the mapping script to what was already in the database.
# some of this info can then be fed back to the source databases
#
##############################################################################################

sub check_original_mappings {

    my $original_db = Ace->connect( -path => "$database" ) || do { print "$database Connection failure: ", Ace->error; die(); };

    my @original_alleles = $original_db->fetch( -query => 'Find Variation WHERE Flanking_sequences' );

    # loop through each allele and then each gene for each allele, adding to hash structure
    foreach my $allele (@original_alleles) {

        # only consider alleles with gene connections
        if ( defined( $allele->Gene ) ) {
            my @gene_IDs = $allele->Gene;
            my $counter  = 0;
            foreach my $gene (@gene_IDs) {
                $original_alleles{ $allele->name }[$counter] = $gene->name;
                $counter++;
            }
        }
    }
    $original_db->close;
}

################################
# read acefiles into autoace   #
#                              #
# use autoace_minder.pl -load  #
################################

sub load_alleles_to_database {

    $log->write_to("\nStart parsing $ace_file in to $database\n\n");
    my $command = "autoace_builder.pl -load $ace_file -tsuser map_Alleles.pl";
    my $status  = system($command);
    if ( ( $status >> 8 ) != 0 ) { die "ERROR: Loading $ace_file file failed \$\? = $status\n" }
    $log->write_to("\nFinished parsing $ace_file in to $database\n\n");
}

__END__

=pod

=head2 NAME - map_Alleles.pl

=head1 USAGE

=over 4

=item map_alleles.pl  [-debug -limit -database=s -release=s -verbose -restart=s -list=s -no_parse -noupdate -test]

=back

This script:

Gets alleles with flanking sequence from designated database and maps them to a clone or superlink using SCAN, then checks which if any CDSs they overlap with.

Also requires that the allele has as a "seed" sequence either a Sequence or CDS (a locus with an associated
genomic_sequence will also work).

Also writes two files for updating allele->Sequence and Allele->CDS in geneace, one to remove the current connections and one to enter the new ones.

Outputs acefiles which are loaded in to the same database.

map_Alleles.pl MANDATORY arguments:

=over 4

=item none

=back

map_Alleles.pl  OPTIONAL arguments:

=over 4

=item -help

this stuff

=back

=over 4

=item -debug 
   
output goes to [user] and uses current_DB

=back

=over 4

=item -noupdate

don't update/touch the mysql database as it might screw up the build

=back

=over 4

=item -limit 

limit the number of alleles mapped (debug tool)

=back

=over 4

=item -list 

list the alleles you want mapped (debug tool) - as filename ie -list "$mapping_dir/to_map"

=back

=over 4

=item -database 

specify which database to read info from and load mapping results in to

=back

=over 4

=item -restart 

choose which allele to start with. all preceding (alphabetically) alleles will be skipped

=back

=over 4

=item -verbose 

greater indication of what procedured are being used to map the allele

=back

=over 4

=item -release

select a version of the database other than that being built, do not use in debug mode

=back




=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk (if run as part
of build.

=item

=item Must be run AFTER gff splitting has produced CHROMOSOME_*.genes.gff

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers ( ar2@sanger.ac.uk)

=back

=cut
