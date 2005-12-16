#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: ar2 $
# Last updated on: $Date: 2005-12-16 11:18:55 $ 

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Getopt::Long;

<<<<<<< update_inferred_multi_pt.pl

##########################################
# command line options
##########################################

my ($help, $database, $debug, $load, $test);

GetOptions ("help"         => \$help,
	    "database=s"   => \$database,   # can specify db for debug purpose (required also when not debug)
	    "debug=s"      => \$debug,
	    "test"         => \$test,       # use test build environment
	    "load"         => \$load,       # load data to autoace
=======
#--------------------
# global variables
#--------------------

my ($help, $database, $debug, $load, $update);

GetOptions ("help"         => \$help,
	    "database=s"  => \$database,   # can specify db for debug purpose (required also when not debug)
	    "debug=s"      => \$debug,
	    "load"         => \$load,       # load CGC approved pseudo map markers
	    "update"       => \$update,
>>>>>>> 1.15.4.1
           );


##############################
# check command line options
###############################

# Display help if required
if ($help){&usage("Help")}


my $user = `whoami`; chomp $user;

# set a default database if not specified
($database = "/wormsrv2/autoace") if (!$database);

if (!$debug && $user ne "wormpub"){print "\nYou need to be wormpub to proceed..exit\n"; exit(0)};

my $maintainers = "All";
$maintainers = $debug if ($debug);


print "\nTarget database for genetic map and multi-pt data uploading is $database...\n\n";



#############################
#
# Misc variables
#
#############################

my $tace = &tace;          # tace executable path

# Set up top level base directory which is different if in test mode
# Make all other directories relative to this
my $basedir   = "/wormsrv2";
$basedir      = glob("~wormpub")."/TEST_BUILD" if ($test);
my $outdir = "$basedir/autoace/acefiles";


# log file information
my $log = Log_files->make_build_log();


# Need input file with genes which have had their Interpolated_map_position promoted to a Map position
# this should have been created by previous build step
# Also need to create two output files, one for new multi point data and one for updated multi point data

my $pseudo_map_positions = "$outdir/pseudo_map_positions.ace";
my $new_output           = "$outdir/new_pseudo_multi_point_data.ace";
my $updated_output       = "$outdir/updated_pseudo_multi_point_data.ace";

# open acedb connection via tace
my $db = Ace->connect(-path    => $database,
                      -program => $tace) || do { print "Connection failure: ",Ace->error; die();};



my %gene2allele; # gene ID is key, first allele attached to that gene is the value
my %locus_order; # key is CGC gene name, value is ordinal position on map
my %order_locus; # reverse of above, key is a integer position, value is CGC gene name



# find out the respective order of genes on genetic map
&get_flanking_loci;


<<<<<<< update_inferred_multi_pt.pl
# now make new multi point data objects
# assumes that promoted map position file has been made earlier in the build
&make_new_inferred_multi_pt_obj; 


# now update existing multi point data objects where necessary
&update_existing_inferred_multi_pt;
=======
#--------------------------------------------------------------------
#   make inferred multi-pt obj for promoted loci during the build
#   update flanking loci of existing old inferred multi-pt obj
#--------------------------------------------------------------------
>>>>>>> 1.15.4.1


#load to database if -load specified
if($load){
  $log->write_to("Loading new pseudo multi point data to $database");
  my $command = "autoace_minder.pl -load $new_output -tsuser new_pseudo_multi_pt_data";

  $log->write_to("Loading updated pseudo multi point data to $database");
  $command = "autoace_minder.pl -load $updated_output -tsuser updated_pseudo_multi_pt_data";

}


################################################
#
# tidy up and exit
#
################################################

$db->close;
$log->mail("$maintainers", "BUILD REPORT: $0");
exit(0);




################################################
#
#            S U B R O U T I N E S
#
################################################



<<<<<<< update_inferred_multi_pt.pl
sub get_flanking_loci {

  # get gene order from last cmp_gmap_with_coord_order_WSXXX.yymmdd.pid file
  my @map_file = glob("$basedir/autoace/MAPPINGS/INTERPOLATED_MAP/cmp_gmap_with_coord_order_*");

  print "Getting flanking loci information from $map_file[-1]\n";

  my $count = 0;
  open(MAP, $map_file[-1]) || die "Can't find $map_file[-1]\n";

  while(<MAP>){
    chomp;
    if ($_ =~ /^(I|V|X)/){
      my($chrom, $map_position, $gene, $cds, $chrom_position) = split(/\s+/, $_);
      $count++;
      if(defined($gene)){       
	$locus_order{$gene}  = $count;
	$order_locus{$count} = $gene;
      }
    }
  }
  close(MAP);
}
=======
>>>>>>> 1.15.4.1



#########################################################################################################

sub make_new_inferred_multi_pt_obj {   # run during the build, when approved pseudo markers are available

  $log->write_to("Creating NEW multi-point data objects\n\n");

  open(PSEUDO, "$pseudo_map_positions") || die $!;

  while (<PSEUDO>){
    chomp;

    #---- first create a hash of gene_id (key) to allele (value)
    #     and also verify that each inferred marker is linked to an allele
    if ($_ =~ /Gene : \"(.+)\"/){
      my $gene_id = $1;
      $gene_id = $db->fetch(-class => 'Gene', -name  => $gene_id);

      if (defined $gene_id->Allele(1)){
	my @alleles = $gene_id->Allele(1);
	foreach my $allele (@alleles){
	  
	  # only want to consider alleles which are not Transposon insertions or SNPs
	  next if(defined($allele->Transposon_insertion));
	  
	  my $method = $allele->Method;
	  next if(defined($method) && ($method eq "SNP" || $method eq "Transposon_insertion"));
	  
	  # now add allele to hash
	  $gene2allele{$gene_id} = $allele; 
	  last;
	}	  
      }
      else {
	$log->write_to("ERROR: $gene_id now has NO allele attached . . . corresponding multi-pt obj needs update. . .\n");
      }      
    }
  }
  close(PSEUDO);


  # get last multipt obj
  my @multi_objs = $db->find('Find Multi_pt_data *');
  my $last_multi = $multi_objs[-1];

  my $multi = $last_multi -1; $multi++;  # last number of multi_pt obj

  # write inferred multi_obj acefile
  open(NEW, ">$new_output") || die $!;


  foreach (keys %gene2allele){

    # store left and right flanking loci
    my $L_locus;
    my $R_locus;

    my $gene = $db->fetch(-class => 'Gene',
			  -name  => $_ );
    my $cgc_name = $gene->CGC_name;
    $multi++;

    # find the loci that are to the left and right of the current locus
    if (defined($cgc_name) && defined($locus_order{$cgc_name})){
      $L_locus = $order_locus{($locus_order{$cgc_name} -1) };
      $R_locus = $order_locus{($locus_order{$cgc_name} +1 )};
    }

    print NEW "\n\nGene : \"$_\"\n";
    print NEW "Multi_point $multi\n\n";

    print NEW "Multi_pt_data : $multi\n";
    print NEW "Gene_A \"$_\" \"$gene2allele{$_}\"\n";
    print NEW "Gene   \"$_\" \"$gene2allele{$_}\"\n";


    # now grab gene IDs for flanking markers based on CGC name
    if ($L_locus && $R_locus){
      my ($L_gene_id)= $db->fetch(-query=>"Find Gene_name $L_locus; Follow CGC_name_for");
      my ($R_gene_id) = $db->fetch(-query=>"Find Gene_name $R_locus; Follow CGC_name_for");
      print NEW "Combined Gene \"$L_gene_id\" 1 Gene \"$_\" 1 Gene \"$R_gene_id\"\n";
      print NEW "Remark \"Data inferred from $gene2allele{$_}, sequence of $cgc_name and interpolated map position (which became genetics map)\" Inferred_automatically\n";

    }
    elsif (!$L_locus){
      $log->write_to("ERROR: Multi-point object $multi has missing left flanking locus for $cgc_name.  Is gene at end of chromosome?\n");
    }
    elsif (!$R_locus){
      $log->write_to("ERROR: Multi-point object $multi has missing right flanking locus for $cgc_name.  Is gene at end of chromosome?\n");
    }

    elsif (!$gene2allele{$_}){
      $log->write_to("ERROR: Multi-point object $multi has no allele information\n");
    }
  }
  close NEW;
}


#############################################################################################################

sub update_existing_inferred_multi_pt {

  $log->write_to("Updating EXISTING multi-point data objects\n\n");

  my $query  = "Find Multi_pt_data * WHERE Remark AND NEXT AND NEXT = \"Inferred_automatically\"";
  push( my @inferred_multi_objs, $db->find($query) );

  open(UPDATE, ">$updated_output") || die $!;

  my ($center_gene, $allele, $L_locus, $R_locus);

  foreach (@inferred_multi_objs){
    $center_gene = $_ -> Combined(5);

    if(!defined($center_gene)) {
      $log->write_to("ERROR: Multi-point object $_ has missing main gene after Combined tag (should be three genes)\n");
    }

    my $gene_id = $db->fetch(-class => 'Gene',
			     -name  => $center_gene );

    my $cgc_name = $gene_id->CGC_name;
    if(!defined($cgc_name)){
      $log->write_to("ERROR: Multi-point object $_ has a gene ($gene_id) which does not have a CGC name\n");
    }

    # get allele associated with center gene
    my $allele = $_ ->Gene(2);
    if(!defined($allele)){
      $log->write_to("ERROR: Multi-point object $_ has a gene ($gene_id) which does not have an allele\n");
    }


    # get flanking gene details from real updated genetic map data 
    # note that this might be different to the existing multi point data if
    # new map markers have been created, e.g. gene order A->B->C is now A->B->D->C
    # where D is a new gene.  Therefore multipoint data for genes B and C have to be
    # updated to reflect new flanking marker

    $L_locus = $order_locus{$locus_order{$cgc_name} -1 };
    $R_locus = $order_locus{$locus_order{$cgc_name} +1 };
    my ($L_gene_id) = $db->fetch(-query=>"Find Gene_name $L_locus; Follow CGC_name_for");
    my ($R_gene_id) = $db->fetch(-query=>"Find Gene_name $R_locus; Follow CGC_name_for");
    
    # Remove old data
    print UPDATE "\nMulti_pt_data : \"$_\"\n";
    print UPDATE "-D Combined\n";
    print UPDATE "-D Remark\n";

    # print new data
    print UPDATE "\nMulti_pt_data : \"$_\"\n";
    if ($L_gene_id && $R_gene_id && $center_gene && $allele){
      print UPDATE "Combined Gene \"$L_gene_id\" 1 Gene \"$center_gene\" 1 Gene \"$R_gene_id\"\n";
      print UPDATE "Remark \"Data inferred from $allele, sequence of $cgc_name and interpolated map position (which became genetics map)\" Inferred_automatically\n";
    }
    elsif (!$L_gene_id) {
      $log->write_to("ERROR: Multi-point object $_ has missing left flanking locus for $gene_id\n");
    }
    elsif (!$R_gene_id) {
      $log->write_to("ERROR: Multi-point object $_ has missing right flanking locus for $gene_id\n");
    }
  }

  close(UPDATE);
}

#####################################################################################################

sub usage {
  my $error = shift;
  if ($error == 0) {
    # Normal help menu
    exec ('perldoc',$0);
  }
}




__END__


=head2 NAME - update_inferred_multi_pt.pl

=head1 USAGE
                                                                                                       
=over 4
                                                                                                       
=item make_pseudo_map_positions.pl -[options]
                                                                                                       
                                                                                                       
=back
                                                                                                       
=head1 DESCRIPTION

An earlier build script (make_pseudo_map_positions.pl) would have 'promoted' certain genes to have a 'Map'
position.  These would have been genes with CGC-names, alleles, and interpolated map positions.  However,
all true map markers, should have associated mapping data stored in Multi_pt_data objects.  These objects
minimally store the name of the gene, an allele used to define it as a mapped marker, and the two markers
that lie either side of the gene.

This script makes two files (in /wormsrv2/autoace/acefiles/).  The first (new_pseudo_multi_point_data.ace)
will contain new multi point data objects to go with the genes with newly promoted markers.  The second 
file (updated_pseudo_multi_point_data.ace) contains modifications to existing multi point objects.  This 
step is necessary because new markers will now lie between markers that were previously flanking markers.

E.g.  assume that gene B is an existing marker which has an associated multi point data object which 
describes two flanking markers (gene A and gene C):

gene A - gene B - gene C

now assume that a new marker (gene D) is added to the map between the positions of gene B and gene D:

gene B - gene D - gene C

the existing multi point data for gene B is now incorrect in stating that gene C is the right hand 
marker.  This second file makes these corrections to multi point data objects.


=back
                                                                                                       
=head1 MANDATORY arguments: <none>
                                                                                                       
=back

=head1 OPTIONAL arguments: -help, -database, -verbose, -test, -load

=over 4

=item -help

Displays this help

=item -database

Specify a path to a valid acedb database to query (defaults to /wormsrv2/autoace)

=item -test

Use the test build environment

=item -load

Load the resulting acefile to autoace (or database specified by -database)

=back

=head1 AUTHOR Keith Bradnam (krb@sanger.ac.uk) - heavily rewritten from a Chao-Kung original


=cut
