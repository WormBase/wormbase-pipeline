#!/usr/local/bin/perl5.8.0 -w                    
#
# map_alleles.pl
#
# by Anthony Rogers
#
# This maps alleles to the genome based on their flanking sequence
#
# Last updated by: $Author: krb $                      # These lines will get filled in by cvs and helps us
# Last updated on: $Date: 2003-12-16 16:22:12 $        # quickly see when script was last changed and by whom

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Feature_mapper;
use Wormbase;
use Ace;
use Getopt::Long;

#######################################
# command-line options                #
#######################################

my ($debug, $update, $limit);
my $database;
my $ver;
my $verbose;
my $restart = "go";
my $help;
my $no_geneace;
my $no_parse;
my $list;
my $gff;
my $geneace;

# $debug   -  all output goes to ar/allele_mapping

GetOptions( "debug:s"    => \$debug,
	    "limit=s"    => \$limit,
	    "database=s" => \$database,
	    "WS=s"       => \$ver,
	    "help"       => \$help,
	    "restart=s"  => \$restart,
	    "no_parse"   => \$no_parse,
	    "list=s"     => \$list,
	    "no_geneace" => \$no_geneace,
	    "gff"        => \$gff,
	    "verbose"    => \$verbose,
	    "geneace=s"  => \$geneace
	  );

if ($help) { print `perldoc $0`;exit;}

##############
# variables  #
##############

# Most checking scripts should produce a log file that is a) emailed to us all
#                                                         b) copied to /wormsrv2/logs

my $maintainers = "All";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = &runtime;
my $log;
$ver = &get_wormbase_version unless defined $ver;

my $data_dump_dir;
my $ace_file;
my $strain_file;
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

# geneace selectable for dev
$geneace = "/wormsrv2/geneace" unless $geneace;

##########  File handles etc #############
open (LOG,">$log") or die "cant open $log\n\n";
print LOG "$0 start at $runtime on $rundate\n----------------------------------\n\n";

my $geneace_update        = "$mapping_dir/map_alleles_geneace_update.WS$ver.ace";
my $geneace_update_delete = "$mapping_dir/map_alleles_geneace_update_delete.WS$ver.ace";
my %geneace_alleles;

my $cshl_update           = "$mapping_dir/map_alleles_cshl_update.WS$ver.ace";
my $cshl_update_delete    = "$mapping_dir/map_alleles_cshl_update_delete.WS$ver.ace";

my $KO_overlap_genes      = "$mapping_dir/KO_genes_overlap";

unless ($no_geneace) {
  open (GENEACE,">$geneace_update") or die "can't open $geneace_update: $!\n";
  open (GEN_DEL,">$geneace_update_delete") or die "can't open $geneace_update_delete\n";

  open (CSHLACE,">$cshl_update") or die "can't open $cshl_update: $!\n";
  open (CSHL_DEL,">$cshl_update_delete") or die "can't open $cshl_update_delete\n";

  open (KOC, ">$KO_overlap_genes") or die "can't open $KO_overlap_genes\n";
  
  # get list of alleles from geneace to check against for feedback files
  my $g_db = Ace->connect(-path => $geneace) || do { print  "$database Connection failure: ",Ace->error; die();};
  my @g_alleles = $g_db->fetch(-query =>'Find Allele;flanking_sequences');
  foreach ( @g_alleles ) {
    my $G_name = $_->name;
    $geneace_alleles{$G_name} = 1;
  }
  $g_db->close;
}
open (OUT,">$ace_file") or die "cant open $ace_file\n";
open (GFF,">$gff_file") or die "cant open $gff_file\n" if $gff;
#open (STR,">$strain_file") or die "cant open $strain_file\n";


########### database accesss ####################

#get allele info from database
my $db = Ace->connect(-path  => $database) || do { print  "$database Connection failure: ",Ace->error; die();};
my @alleles = $db->fetch(-query =>'Find Allele;flanking_sequences');
my @KO_alleles = $db->fetch(-query => 'Find Allele; KO_consortium_allele');
my %KO_alleles;
%KO_alleles = map {$_ => 1} @KO_alleles;

# read in list of alleles to map if specified #######################
my %to_map;
if ( $list ) {
  open (LIST, "<$list");
  while (<LIST>) {
    chomp;
    $to_map{$_} = 1;
  }
}



####### allele mapping loop ######################

my %allele_data;   #  allele => [ (0)type, (1)5'flank_seq , (2)3'flank_seq, (3)CDS, (4)end of 5'match, (5)start of 3'match , (6)clone, (7)chromosome, (8)strains]
my $error_count = 0; # for tracking alleles that don't map
my $sequence;
my $clone;
my $chromosome;
my $name;
my $left;
my $right;
my $go = 0;
my $KO_allele;
my %allele2gene;
my @affects_genes;

my $mapper = Feature_mapper->new( $database);

ALLELE:
foreach my $allele (@alleles) {
  undef $KO_allele;
  $name = $allele->name;

  # debug facility - this bit is so that it can restart from given allele name
  unless ("$restart" eq "go"){
    if ("$restart" eq "$name") {
      $restart = "go";
    } else { 
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

  # lets get going . . . 
  $left     = lc $allele->Flanking_sequences->name;
  $right    = lc $allele->Flanking_sequences->right->name;  

  next unless ($left and $right);
  print "mapping $name\n" if $verbose;

  $sequence = $allele->Sequence;
  unless (defined $sequence->Source) {
    print LOG "$name has an invalid sequence\n";
    next ALLELE;
  }

  $allele_data{$name}[1] = $left;
  $allele_data{$name}[2] = $right;

  my @map = $mapper->map_feature($sequence->name,$left, $right);
  if( "$map[0]" eq "0" ) {
    print LOG "ERROR: Couldn't map $name with $left and $right to seq $sequence\n";
    $error_count++;
    next ALLELE;
  }

  $allele_data{$name}[4] = $map[1];
  $allele_data{$name}[5] = $map[2];
  $allele_data{$name}[6] = $map[0];


  @affects_genes = $mapper->check_overlapping_CDS($map[0],$map[1],$map[2]);
  print "$name hits @affects_genes\n" if $verbose;
  $allele2gene{"$name"} = \@affects_genes if ($affects_genes[0]);

  # this identified KO_consortium alleles so they can be fed back to KOAC
  my $method = $allele->Method;
  if ("$method" eq "Knockout_allele") {
    $KO_allele = 1;
  }
  &outputAllele($name);
}

  print "$error_count alleles\n";

print LOG "Update files for geneace allele mappings are available - \n$geneace_update\n$geneace_update_delete\n\n" unless $no_geneace;

$db->close;


close OUT;
close GFF if $gff;
close GEN_DEL unless $no_geneace;
close GENEACE unless $no_geneace;

##############################
# read acefiles into autoace #
##############################
unless ( $no_parse ) {
  print LOG "\nStart parsing $ace_file in to $database\n\n";

  my $command =<<END;
pparse $geneace_update_delete
pparse $ace_file
save
quit
END

  my $tace = &tace;
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
&mail_maintainer("map_alleles","$maintainers","$log");

exit(0);

#######################################
#                                     #
#          SUB ROUTINES               #
#                                     #
#######################################

sub outputAllele
  {
    my $to_dump = shift;
    print KOC "$to_dump @affects_genes\n" if ( defined $KO_alleles{$to_dump} and defined $affects_genes[0]);
    if( $allele_data{$to_dump}[6] and $allele_data{$to_dump}[4] and  $allele_data{$to_dump}[5]) { 
      
      print OUT "\nSequence : \"$allele_data{$to_dump}[6]\"\nAllele $to_dump $allele_data{$to_dump}[4] $allele_data{$to_dump}[5]\n";
      print GFF "\n$allele_data{$to_dump}[6]\tAllele\tTEST\t$allele_data{$to_dump}[4]\t$allele_data{$to_dump}[5]\t.\t+\t.\tAllele \"$to_dump\"" if $gff;
      if( $allele2gene{$to_dump} ) {
	my @affects_genes = split(/\s/,"@{$allele2gene{$to_dump}}");
	
	# in CDS object
	foreach my $ko (@affects_genes) {
	  print OUT "\nCDS : \"$ko\"\nAllele $to_dump\n";
	}

	# in Allele object
	print OUT "\nAllele : $to_dump\n";
	unless ($no_geneace ) {
	  if ( defined $geneace_alleles{$to_dump} ) {
	    print GEN_DEL "\nAllele : $to_dump\n-D Predicted_gene\n-D Sequence\n"; # remove current sequence and predicted genes from Geneace
	    print GENEACE "\nAllele : \"$to_dump\"\nSequence \"$allele_data{$to_dump}[6]\"\n"; # allele -> sequence
	  } else {
	    print CSHL_DEL "\nAllele : $to_dump\n-D Predicted_gene\n-D Sequence\n"; # remove current sequence and predicted genes from Geneace
	    print CSHLACE "\nAllele : \"$to_dump\"\nSequence \"$allele_data{$to_dump}[6]\"\n";
	  }
	}
	foreach my $ko (@affects_genes) {
	  #allele - CDS connection
	  print OUT "Predicted_gene $ko\n";
	  unless ($no_geneace ) {
	    if ( defined $geneace_alleles{$to_dump} ) {
	      print GENEACE "Predicted_gene $ko\n"; # update geneace with allele -> CDS
	    } else {
	      print CSHLACE "Predicted_gene $ko\n";
	    }
	  }
	}
      }
      else {
	print "no overlapping gene for $to_dump\n" if $verbose;
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

    }
  }










# Add perl documentation in POD format
# This should expand on your brief description above and add details of any options
# that can be used with the program.  Such documentation can be viewed using the perldoc
# command.



__END__

=pod

=head2 NAME - map_alleles.pl

=head1 USAGE

=over 4

=item map_alleles.pl  [-debug -limit -database=s -WS=s -verbose -restart=s -list=s -no_parse -no_geneace]

=back

This script:

Gets alleles with flanking sequence from designated database and maps them to a clone or superlink using SCAN, then checks which if any CDSs they overlap with.

Also requires that the allele has as a "seed" sequence either a Sequence or CDS (a locus with an associated
genomic_sequence will also work).

Also writes two files for updating allele->Sequence and Allele->CDS in geneace, one to remove the current connections and one to enter the new ones.

Outputs acefiles which are loaded in to the same database.

map_alleles.pl MANDATORY arguments:

=over 4

=item none

=back

map_alleles.pl  OPTIONAL arguments:

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
