#!/usr/local/bin/perl5.8.0 -w
# map_Alleles.pl
#
# by Anthony Rogers
#
# This maps alleles to the genome based on their flanking sequences
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2006-02-22 10:48:11 $

use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
use Feature_mapper;
use Wormbase;
use Ace;
use Getopt::Long;

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
my $gff;               # option to print output in GFF format as well
my $verbose;           # verbose mode, extra output to screen
my $store;             # to specify stored commandline arguments
my $test;	       # for future testing functionality

GetOptions(
    'debug=s'    => \$debug,
    'limit=i'    => \$limit,
    'database=s' => \$database,
    'release=i'  => \$release,
    'help'       => \$help,
    'restart=s'  => \$restart,
    'no_parse'   => \$no_parse,
    'list=s'     => \$list,
    'gff'        => \$gff,
    'verbose'    => \$verbose,
    'store=s'    => \$store,
    'test'	=> \$test
);

# check command line options
if ($help) { print `perldoc $0`; exit; }
die "No need for -release option in debug mode" if ( $release && $debug );

############################
# recreate configuration   #
############################
my $wb;
if ($store) { $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wb = Wormbase->new( -debug => $debug, -test => $test, ) }

###########################################
# Variables Part II (depending on $wb)    #
###########################################
$test  = $wb->test  if $wb->test;     # Test mode
$debug = $wb->debug if $wb->debug;    # Debug mode, output only goes to one user

###################
# misc variables  #
###################
my $maintainers = "All";
my $tace        = $wb->tace;
my $error_count = 0;           # for tracking alleles that don't map

# set variables depending on whether in debug (test) mode or not
if ($debug) {
    print "DEBUG = \"$debug\"\n\n";
    $maintainers = "$debug\@sanger.ac.uk";
}
$release = $wb->get_wormbase_version unless ( defined $release );
my $mapping_dir = $wb->autoace."/MAPPINGS";
my $ace_file    = "$mapping_dir/allele_mapping.WS$release.ace";
my $gff_file    = "$mapping_dir/allele_mapping.WS$release.gff";
$database = $wb->autoace unless $database;

################################################
# Set hashes
################################################

my %worm_gene2class  = $wb->FetchData('worm_gene2class');         # set from Common Data
my %worm_gene2geneID = $wb->FetchData('worm_gene2geneID_name');
my %to_map;              # will be set if -list option is specified, stores list of alleles to map
my %original_alleles;    # will store details of allele->gene connections already in database
my %allele2gene;         # key is allele name, value is array of CDS/Transcript names
my %allele_data;         # allele data from mapping data, hash of arrays:

# allele => [ (0)type,           (1) 5'flank_seq ,      (2) 3'flank_seq, (3) CDS,
#             (4)end of 5'match, (5) start of 3'match , (6) clone,       (7) chromosome,
#             (8)strains]

###############################################################################################
#
#                                   MAIN part of script
#
###############################################################################################

# read in list of alleles to map if -list is specified
&process_list if ($list);

# make hash of allele->gene connections already in database
&check_original_mappings;

# create log file, open output file handles
my $log = Log_files->make_build_log($wb);


  # map alleles, the main routine
  &map_alleles;

# load acefile to autoace (unless -no_parse specified)
&load_alleles_to_database unless ($no_parse);

#####################################
# tidy up, email log, and exit
#####################################

$log->mail();
exit(0);

##################################################
#                                                #
#            S U B R O U T I N E S               #
#                                                #
##################################################

#########################################################
#             Main allele mapping loop                  #
#########################################################

sub map_alleles {

    # First get allele info from database
    my $db = Ace->connect( -path => $database ) || do { print "$database Connection failure: ", Ace->error; die(); };
    my @alleles = $db->fetch( -query => 'Find Variation WHERE Flanking_sequences' );

    open( OUT, ">$ace_file" ) or die "cant open $ace_file\n";
    open( GFF, ">$gff_file" ) or die "cant open $gff_file\n" if ($gff);

    my $sequence;
    my $clone;
    my $chromosome;
    my $name;
    my $left;
    my $right;
    my $go = 0;
    my @affects_CDSs;

    my $mapper = Feature_mapper->new($database, undef, $wb);

  ALLELE:
    foreach my $allele (@alleles) {
        $name = $allele->name;

        print "\nMapping $name\n" if $verbose;

        # debug facility - this bit is so that it can restart from given allele name
        unless ( "$restart" eq "go" ) {
            if ( "$restart" eq "$name" ) {
                $restart = "go";
            }
            else {
                print "skipping $name\n" if $verbose;
                next;
            }
        }

        # debug facility - after the restart means you can specify where to start and how many to do
        if ($limit) {
            last if $error_count++ >= $limit;
        }

        if ($list) {
            next unless defined $to_map{$name};
        }

        # grab both flanking sequences from $database
        $left  = lc $allele->Flanking_sequences->name;

	# some alleles only have one flank!
	unless ($allele->Flanking_sequences->right) {
	  $log->write_to("ERROR: ",$allele->name," is missing a flanking sequence!\n");
	  next ALLELE;
	}
        $right = lc $allele->Flanking_sequences->right->name;

        # warn if flanking sequence is missing
        unless ( $left and $right ) {
            $log->write_to("ERROR: $name does not have two flanking sequences\n");
            $error_count++;
            next ALLELE;
        }

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

        $allele_data{$name}[1] = $left;
        $allele_data{$name}[2] = $right;

        # map allele using Feature_mapper.pm, store results of map in @map, warn if mapping failed
        my @map = $mapper->map_feature( $sequence->name, $left, $right );

        if ( "$map[0]" eq "0" ) {
            $log->write_to("ERROR: Couldn't map $name with $left and $right to seq $sequence\n");
            $error_count++;
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

        unless ( ( $allele->Type_of_mutation and $allele->Type_of_mutation eq "Insertion" )
            or ( $allele->Variation_type and $allele->Variation_type eq "Transposon_insertion" ) )
        {

            # get coords of allele not 1st / last base of flanks
            if ( $map[2] > $map[1] ) {

                # maps to fwd strand
                $map[1]++;
                $map[2]--;
            }
            else {
                $map[1]--;
                $map[2]++;
            }
        }

        $allele_data{$name}[4] = $map[1];
        $allele_data{$name}[5] = $map[2];
        $allele_data{$name}[6] = $map[0];

        # retrieve list of overlapping CDSs/Transcripts from Feature_mapper.pm
        @affects_CDSs = $mapper->check_overlapping_CDS( $map[0], $map[1], $map[2] );

        # make (unique) array of gene names corresponding to @affects_CDSs
        # first make a hash, then use the keys of the hash to set new array
        my %affects_genes;
        foreach (@affects_CDSs) {
            my $gene = $worm_gene2geneID{$_};
            $affects_genes{$gene} = 1;
        }
        my @affects_genes = keys(%affects_genes);

    # now compare both arrays (allele->gene connections already in database and
    # allele->gene connections made by script) to look for differences
    my %count;
    foreach my $gene (@{$original_alleles{$name}}, @affects_genes){
      $count{$gene}++;
    }
    foreach my $gene (keys %count){
      # allele->gene connections in both arrays will cause $count{$gene} to be == 2
      # if it only == 1 then there is an error
      if ($count{$gene} == 1){
	if($affects_genes{$gene}){
	  # not so serious, but source database should be updated with connection
	  next if ( $allele->SNP or $allele->Transposon_insertion );
	  print "WARNING: $name->$gene was mapped by script. Add connection in Geneace?n" if ($verbose);
	  $log->write_to("WARNING: $name->$gene was mapped by script. Add connection in Geneace?\n");
	  
	}
	else{
	  # more serious, source database has a bad allele->gene connection
	  print "ERROR: $name->$gene was in original database only\n" if ($verbose);
	  $log->write_to("ERROR: $name->$gene is in Geneace. Script indicates this could be wrong\n");
	  $error_count++;       
	}
      }
    }
    
	$allele2gene{"$name"} = \@affects_CDSs if ($affects_CDSs[0]);
	
	&outputAllele($name);
    }

    # close database and file handles
    $db->close;

    close OUT;
    close GFF if $gff;

}

#####################################################################################################

sub outputAllele {
    my $allele = shift;

    # only proceed if mapping routine returns the following three pieces of info
    if ( $allele_data{$allele}[6] and $allele_data{$allele}[4] and $allele_data{$allele}[5] ) {

        print OUT
"\nSequence : \"$allele_data{$allele}[6]\"\nAllele $allele $allele_data{$allele}[4] $allele_data{$allele}[5]\n";
        print GFF
"\n$allele_data{$allele}[6]\tAllele\tTEST\t$allele_data{$allele}[4]\t$allele_data{$allele}[5]\t.\t+\t.\tAllele \"$allele\""
          if $gff;

        # process allele if it matches overlapping CDSs/Transcripts
        if ( $allele2gene{$allele} ) {
            my @affects_CDSs = split( /\s/, "@{$allele2gene{$allele}}" );

            # loop through all CDSs/Transcripts that overlap allele
            foreach my $cds (@affects_CDSs) {

                #allele - WBGene connection
                my $WBGene = $worm_gene2geneID{$cds};
                print OUT "\nGene : $WBGene\n";
                print OUT "Allele $allele\n\n";

                # now write relevant acefile output depending on whether allele hits CDS or Transcript
                if ( !defined( $worm_gene2class{$cds} ) ) {
                    $log->write_to("ERROR: $cds not listed in %worm_gene2class\n");
                    $error_count++;
                    next;
                }
                elsif ( $worm_gene2class{$cds} eq "CDS" ) {
                    print OUT "\nCDS : \"$cds\"\n";
                    print OUT "Variation $allele\n\n";

                    print OUT "Variation : $allele\n";
                    print OUT "Predicted_CDS $cds\n";
                }
                elsif ( $worm_gene2class{$cds} eq "Transcript" ) {
                    print OUT "\nTranscript : \"$cds\"\n";
                    print OUT "Variation $allele\n\n";

                    print OUT "Variation : $allele\n";
                    print OUT "Transcript $cds\n";
                }
                else {
                    $log->write_to("ERROR: $cds is not a CDS or Transcript\n");
                    $error_count++;
                }
		#write Database line for link to KO_consortium
		if( $allele =~ /^[og]k\d/ ) {
		  print "\nVariation : $allele\n";
		  print "Database Gene_Knockout_Consortium GeneID \"<GeneID=$cds>\"\n";
		}
            }
        }
    }
}

####################################################################################
# small subroutine to load data into %to_map hash from file if -list is specified
####################################################################################

sub process_list {
    open( LIST, "<$list" ) || die "Cannot open file specified by $list\n";
    while (<LIST>) {
        chomp;
        $to_map{$_} = 1;
    }
    close(LIST);
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

    my $original_db = Ace->connect( -path => "$database" )
      || do { print "$database Connection failure: ", Ace->error; die(); };

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
# use autoace_builder.pl -load  #
################################

sub load_alleles_to_database {

    $log->write_to("\nStart parsing $ace_file in to $database\n\n");
    $wb->load_to_database($wb->autoace,$ace_file, 'map_alleles');
    $log->write_to("\nFinished parsing $ace_file in to $database\n\n");
}

__END__

=pod

=head2 NAME - map_Alleles.pl

=head1 USAGE

=over 4

=item map_alleles.pl  [-debug -limit -database=s -release=s -verbose -restart=s -list=s -no_parse]

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
   
output goes to ar2 and uses current_DB

=back

=over 4

=item -limit 

limt the number of alleles mapped (debug tool)

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
