#!/usr/local/bin/perl5.8.0 -w                    
#
# map_alleles.pl
#
# by Anthony Rogers
#
# This maps alleles to the genome based on their flanking sequences
#
# Last updated by: $Author: ar2 $                      
# Last updated on: $Date: 2004-08-06 15:42:25 $        

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Feature_mapper;
use Wormbase;
use Ace;
use Getopt::Long;


#######################################
# command-line options                #
#######################################

my $debug;          # debug mode, output goes to /tmp
my $limit;          # limit number of alleles to map
my $database;       # specify database to map alleles in, default is autoace
my $ver;            # specify release number, defaults to current release number
my $help;           # help mode, show perldoc
my $restart = "go"; # specify an allele name from which to start mapping
my $no_parse;       # turn off loading of data to $database
my $list;           # read in alleles to map from file rather than from database
my $gff;            # option to print output in GFF format as well
my $verbose;        # verbose mode, extra output to screen


GetOptions( "debug=s"    => \$debug,
	    "limit=i"    => \$limit,
	    "database=s" => \$database,
	    "WS=i"       => \$ver,
	    "help"       => \$help,
	    "restart=s"  => \$restart,
	    "no_parse"   => \$no_parse,
	    "list=s"     => \$list,
	    "gff"        => \$gff,
	    "verbose"    => \$verbose,
	  );

if ($help) { print `perldoc $0`;exit;}



###################
# misc variables  #
###################

my $maintainers = "All";
my $rundate     = &rundate;
my $runtime     = &runtime;
my $log;
$ver = &get_wormbase_version unless defined $ver;
my $tace = &tace;

my $data_dump_dir;
my $ace_file;
my $gff_file;
my $mapping_dir;

if ($debug)  {
    $data_dump_dir = "/tmp";
    $log           = "$data_dump_dir/map_Alleles.$rundate";
    $database      = "/nfs/disk100/wormpub/DATABASES/current_DB/" unless $database;
    $ver--;
    $ace_file      = "$data_dump_dir/mapped_alleles.ace";
    $gff_file      = "$data_dump_dir/mapped_alleles.gff";
    $mapping_dir   = $data_dump_dir;
    $maintainers   = "$debug\@sanger.ac.uk";
}
else { 
    $log         = "/wormsrv2/logs/map_AllelesWS$ver.$rundate.$$";
    $mapping_dir = "/wormsrv2/autoace/MAPPINGS";
    $ace_file    = "$mapping_dir/allele_mapping.WS$ver.ace";
    $gff_file    = "$mapping_dir/allele_mapping.WS$ver.gff";
    $database    = "/wormsrv2/autoace/" unless $database;
}


# read in list of alleles to map if -list is specified #######################
my %to_map;

if ( $list ) {
  open (LIST, "<$list");
  while (<LIST>) {
    chomp;
    $to_map{$_} = 1;
  }
  close(LIST);
}




########## get list of alleles from geneace ####################

my %geneace_alleles;

# list of alleles will be used later check against for feedback files
my $geneace_db = Ace->connect(-path => "/wormsrv2/geneace") || do { print  "$database Connection failure: ",Ace->error; die();};
my @geneace_alleles = $geneace_db->fetch(-query =>'Find Allele WHERE Flanking_sequences');
foreach ( @geneace_alleles ) {
  my $G_name = $_->name;
  my $wb_gene = $_->Gene;
  $wb_gene = "undef" unless $wb_gene;
  $geneace_alleles{$G_name} = $wb_gene;
}
$geneace_db->close;


##########  File handles etc #############

open (LOG,">$log") or die "cant open $log\n\n";
print LOG "$0 start at $runtime on $rundate\n----------------------------------\n\n";
open (OUT,">$ace_file") or die "cant open $ace_file\n";
open (GFF,">$gff_file") or die "cant open $gff_file\n" if $gff;

my %CDS2gene;
&FetchData('cds2wbgene_id',\%CDS2gene);


#########################################################
#                                                       #
#             Main allele mapping loop                  #
#                                                       #
#########################################################

# First get allele info from database
my $db = Ace->connect(-path  => $database) || do { print  "$database Connection failure: ",Ace->error; die();};
my @alleles = $db->fetch(-query =>'Find Allele WHERE Flanking_sequences');


my %allele_data;   #  allele => [ (0)type, (1)5'flank_seq , (2)3'flank_seq, (3)CDS, (4)end of 5'match, (5)start of 3'match , (6)clone, (7)chromosome, (8)strains]
my $error_count = 0; # for tracking alleles that don't map
my $sequence;
my $clone;
my $chromosome;
my $name;
my $left;
my $right;
my $go = 0;
my %allele2gene;
my @affects_genes;

my $mapper = Feature_mapper->new( $database);

ALLELE:
foreach my $allele (@alleles) {
  $name = $allele->name;

  print "Mapping $name\n" if $verbose;


  # debug facility - this bit is so that it can restart from given allele name
  unless ("$restart" eq "go"){
    if ("$restart" eq "$name") {
      $restart = "go";
    } 
    else { 
      print "skipping $name\n" if $verbose;
      next;
    }
  }

  # debug facility - after the restart means you can specify where to start and how many to do 
  if ( $limit ) {
    last if ++$error_count >= $limit;
  }

  if ( $list ) {
    next unless defined $to_map{$name};
  }

  # grab both flanking sequences from $database
  $left  = lc $allele->Flanking_sequences->name;
  $right = lc $allele->Flanking_sequences->right->name;  

  # warn if flanking sequence is missing
  unless ($left and $right){
    print LOG "ERROR: $name does not have two flanking sequences\n";
    $error_count++;
    next ALLELE;
  }


  # check that allele is attached to a valid sequence object
  $sequence = $allele->Sequence;
  if(!defined $sequence){    
    print LOG "ERROR: $name has missing Sequence tag\n";
    $error_count++;
    next ALLELE;
  }
  elsif(!defined $sequence->Source) {
    print LOG "ERROR: $name connects to Sequence $sequence which has no Source tag\n";
    $error_count++;
    next ALLELE;
  }


  $allele_data{$name}[1] = $left;
  $allele_data{$name}[2] = $right;

  # map allele using Feature_mapper.pm, store results of map in @map, warn if mapping failed
  my @map = $mapper->map_feature($sequence->name,$left, $right);
  if( "$map[0]" eq "0" ) {
    print LOG "ERROR: Couldn't map $name with $left and $right to seq $sequence\n";
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
  

  # get coords of allele not 1st / last base of flanks
  if( $map[2] > $map[1] ) {
    # maps to fwd strand
    $map[1]++; $map[2]--;
  }
  else {
    $map[1]--; $map[2]++;
  }

  $allele_data{$name}[4] = $map[1];
  $allele_data{$name}[5] = $map[2];
  $allele_data{$name}[6] = $map[0];


  @affects_genes = $mapper->check_overlapping_CDS($map[0],$map[1],$map[2]);
  print "$name hits @affects_genes\n" if $verbose;
  $allele2gene{"$name"} = \@affects_genes if ($affects_genes[0]);

  &outputAllele($name);
}

print "\nWARNING: $error_count alleles failed to map\n\n";

$db->close;


close OUT;
close GFF if $gff;

##############################
# read acefiles into autoace #
##############################
unless ( $no_parse ) {
  print LOG "\nStart parsing $ace_file in to $database\n\n";

  my $command =<<END;
pparse $ace_file
save
quit
END

  eval{
    open (TACE,"| $tace -tsuser map_allele $database") || warn "Couldn't open tace connection to $database\n";
    print TACE $command;
    close (TACE);
  };
  if ( $@ ) {
    print LOG "parse failure . . \n$@\n";
  }
  else {
    print LOG "successfully finished parsing\n";
  }
}


# close LOG and send mail

print LOG "ERROR: $error_count alleles failed to map\n" if ($error_count > 0);
print LOG "$0 end at ",&runtime," \n-------------- END --------------------\n\n";
close LOG;
if($error_count > 0){
  &mail_maintainer("BUILD REPORT: map_Alleles.pl $error_count ERRORS!","$maintainers","$log");
}
else{
  &mail_maintainer("BUILD REPORT: map_Alleles.pl","$maintainers","$log");
}
exit(0);

#######################################
#                                     #
#          SUB ROUTINES               #
#                                     #
#######################################

sub outputAllele{
  my $to_dump = shift;

  if( $allele_data{$to_dump}[6] and $allele_data{$to_dump}[4] and  $allele_data{$to_dump}[5]) { 
      
    print OUT "\nSequence : \"$allele_data{$to_dump}[6]\"\nAllele $to_dump $allele_data{$to_dump}[4] $allele_data{$to_dump}[5]\n";
    print GFF "\n$allele_data{$to_dump}[6]\tAllele\tTEST\t$allele_data{$to_dump}[4]\t$allele_data{$to_dump}[5]\t.\t+\t.\tAllele \"$to_dump\"" if $gff;
    if( $allele2gene{$to_dump} ) {
      my @affects_genes = split(/\s/,"@{$allele2gene{$to_dump}}");
      
      # in CDS object
      foreach my $ko (@affects_genes) {
	print OUT "\nCDS : \"$ko\"\nAlleles $to_dump\n";
      }

      # in Allele object
      print OUT "\nAllele : $to_dump\n";
      foreach my $ko (@affects_genes) {
	#allele - CDS connection
	print OUT "Predicted_gene $ko\n";

	#allele - WBGene connection
	my $WBGene = $CDS2gene{$ko};
	print OUT "Gene $WBGene\n";

	if( defined($geneace_alleles{$to_dump}) and $geneace_alleles{$to_dump}->name ne $WBGene ) {
	  print LOG "Geneace Gene error $ko : Geneace $geneace_alleles{$ko}  Mapping $WBGene\n";
	}
      }
    }
    else {
      print "no overlapping gene for $to_dump\n" if $verbose;
    }

    
  }
}



__END__

=pod

=head2 NAME - map_Alleles.pl

=head1 USAGE

=over 4

=item map_alleles.pl  [-debug -limit -database=s -WS=s -verbose -restart=s -list=s -no_parse]

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

=item -ver 

select a version of the database other than that being built

=back




=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=item

=item Must be run AFTER gff splitting has produced CHROMOSOME_*.genes.gff

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers ( ar2@sanger.ac.uk)

=back

=cut
