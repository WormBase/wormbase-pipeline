#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-05-18 10:01:51 $ 

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'}; 
use Wormbase;
use Ace;
use GENEACE::Geneace;
use Getopt::Long;

#--------------------
# global variables
#--------------------

my ($help, $database, $debug, $version, $update);

GetOptions ("h|help"         => \$help,
	    "db|database=s"  => \$database,
	    "d|debug=s"      => \$debug,
	    "v|version=s"    => \$version,       # the build number, eg 124
	    "u|update"       => \$update,
           );

my $user = `whoami`; chomp $user;
my $rundate = &rundate;
my $tace = glob("~wormpub/ACEDB/bin_ALPHA/tace");          # tace executable path
my $autoace_version = "WS".get_wormbase_version();
my $multi_dir = "/wormsrv1/geneace/JAH_DATA/MULTI_PT_INFERRED";

my $log = "/wormsrv2/logs/new_inferred_multi_pt.$rundate"; `rm -f $log`;
open(LOG, ">>$log") || die $!;
print LOG "# $0 started at ", runtime(), "\n\n";;
print LOG "=============================================================================================\n\n";

if (!$database){print "\nDatabase path not specified..exit.\n"; exit(0)}
if (!$debug && $user ne "wormpub"){print "\nYou need to be wormpub to proceed..exit\n"; exit(0)};

print "\nTarget database for genetic map and multi-pt data uploading is $database...\n\n";


my $db = Ace->connect(-path    => $database,
                      -program => $tace) || do { print "Connection failure: ",Ace->error; die();};

my $ga = init Geneace();

my %Alleles = $ga->get_non_Transposon_alleles($db); # all Geneace alleles which have no Transposon_insertion tag

my (%Gene_info, %gene_id_allele, %locus_order, %order_locus);

# CGC approved promoted loci

my $new_multi_file = `ls $multi_dir/loci_become_genetic_marker_for_$autoace_version`;
chomp $new_multi_file;

#--------------------------------------------------------------------------------------------------------------
#   load promoted marker Gene (loci) to geneace for the next build
#   ENSURES that this gets into the build
#   NOTE: the version number of the loci_become_genetic_marker_for_WSxxx file should be that of the next build,
#         and loci in this file should be approved by Jonathan first
#--------------------------------------------------------------------------------------------------------------

&convert_int_map_to_map_for_promoted_marker_Gene if $version;

#--------------------------------------------------------------------
#   make inferred multi-pt obj for promoted loci during the build
#   update flanking loci of all inferred multi-pt obj
#--------------------------------------------------------------------

if ($update){

  %Gene_info = $ga -> gene_info($database);   ## spcified db

  &get_flanking_loci;

  &make_inferred_multi_pt_obj if $new_multi_file; # do this only if there is new inferred_multi pt obj
  &update_inferred_multi_pt;
}

print LOG "\n$0 finished at ", runtime(), "\n\n";
mail_maintainer("Update inferred multi-pt objects", "ALL", $log) if !$debug;
mail_maintainer("Update inferred multi-pt objects", "$debug\@sanger.ac.uk", $log) if $debug;


#-------------------------
# s u b r o u t i n e s
#-------------------------

 sub convert_int_map_to_map_for_promoted_marker_Gene {

  # load file with CGC approved promoted Gene(loci) to Geneace
  my $upload = "pparse $new_multi_file\nsave\nquit";
  $ga->upload_database($ga->geneace(), $upload, "pseudo_marker_loci", $log);  # upload to Geneace before the start of a WS build
 # $ga->upload_database($database, $upload, "pseudo_marker_loci", $log); # only for testing
}

sub make_inferred_multi_pt_obj {   # run during the build

  open(F, "$new_multi_file") || die $!;

  my $locus;
  while (<F>){
    chomp;
    #---- first create a hash of gene_id (key) to allele (value)
    #     and also verify that each inferred marker is linked to an allele

    if ($_ =~ /Gene : (.+)/){
      my $gene_id = $1;
      $gene_id = $db->fetch(-class => 'Gene',
			    -name  => $gene_id);
      if (defined $gene_id -> CGC_name(1)){
	if (defined $gene_id -> Allele(1)){
	  my @alleles = $gene_id -> Allele(1);
	  foreach my $e (@alleles){
	    if (exists $Alleles{$e} ){
	      $gene_id_allele{$gene_id} = $e; # grep only the allele which has no Transposon_insertion tag
	      #print "$gene_id => $e (1)\n";
	      last;
	    }
	  }
	}
	else {
	  print LOG "ERROR: $gene_id has now NO allele attached . . . corresponding multi-pt obj needs update. . .\n";
	}
      }
    }
  }
  close F;

  # get last multipt obj
  my $multipt = "find Multi_pt_data *";
  my @multi_objs = $db->find($multipt);
  my $last_multi = $multi_objs[-1];
  my $multi = $last_multi -1; $multi++;  # last number of multi_pt obj
  my $error = 0;

  # write inferred multi_obj acefile
  open(NEW, ">$multi_dir/inferred_multi_pt_obj_$autoace_version") || die $!;

  foreach (keys %gene_id_allele){

    $multi++;
    my $L_locus = $order_locus{ $locus_order{ $Gene_info{$_}{'CGC_name'} } -1 } if exists $Gene_info{$_}{'CGC_name'};
       $L_locus = $order_locus{ $locus_order{ $Gene_info{$_}{'Other_name'} } -1 } if exists $Gene_info{$_}{'Other_name'};
    print $L_locus, " (L)\n";
    my $R_locus = $order_locus{ $locus_order{ $Gene_info{$_}{'CGC_name'} } +1 } if exists $Gene_info{$_}{'CGC_name'};
       $R_locus = $order_locus{ $locus_order{ $Gene_info{$_}{'Other_name'} } +1 } if exists $Gene_info{$_}{'Other_name'};
    print $R_locus, " (R)\n";

    print NEW "\n\nGene : \"$_\"\n";
    print NEW "Multi_point $multi\n";
    print NEW "\n\nMulti_pt_data : $multi\n";
    print NEW "Gene_A \"$_\" \"$gene_id_allele{$_}\"\n";
    print NEW "Gene \"$_\" \"$gene_id_allele{$_}\"\n";
    if ($L_locus && $R_locus){
      print NEW "Combined Gene \"$Gene_info{$L_locus}{'Gene'}\" 1 Gene \"$_\" 1 Gene \"$Gene_info{$R_locus}{'Gene'}\"\n";
      print NEW "Remark \"Data inferred from $gene_id_allele{$_}, sequence of $Gene_info{$_}{'CGC_name'} and interpolated map position (which became genetics map)\" Inferred_automatically\n";
    }
    elsif (!$L_locus | !$R_locus){
      $error = 1;
      $L_locus = "NA" if !$L_locus;
      $R_locus = "NA" if !$R_locus;
      print LOG "ERROR: Multi-pt obj $multi has incomplete flanking loci [(L) $L_locus (R) $R_locus] information for $Gene_info{$_}{'CGC_name'}($_)\n";
    }
    elsif (!$gene_id_allele{$_}){
      print LOG "ERROR: Multi-pt obj $multi has NO allele [NA] information\n";
    }
  }
  close NEW;

  print LOG "\nInfo: for incomplete flanking loci: check that they should be the end marker of a chromosome: back to JAH.\n"; 

  # load $multi_dir/inferred_multi_pt_obj_$autoace_version to autoace
  my $upload = "pparse $multi_dir/inferred_multi_pt_obj_$autoace_version\nsave\nquit\n";
  $ga->upload_database($database, $upload, "pseudo_mapping_data", $log) if !$debug;
}

sub update_inferred_multi_pt {

  my $query  = "find Multi_pt_data * where remark AND NEXT AND NEXT = \"inferred_automatically\"";
  push( my @inferred_multi_objs, $db->find($query) );

  open(UPDATE, ">/tmp/updated_multi_pt_flanking_loci_$autoace_version") || die $!;
  `chmod 777 /tmp/updated_multi_pt_flanking_loci_$autoace_version`;

  my ($center_Gene, $allele, $L_locus, $R_locus);

  foreach (@inferred_multi_objs){
    $center_Gene = $_ -> Combined(5);

    my $gene_id = $db->fetch(-class => 'Gene',
			     -name  => $center_Gene );

    $allele =(); my @alleles =();
    if (defined $gene_id -> Allele(1)){
      @alleles = $gene_id -> Allele(1);
      foreach my $e (@alleles){
	if (exists $Alleles{$e} ){
	  $allele = $e; # grep only the allele which has no Transposon_insertion tag
	  last;
	}
      }
    }

    $L_locus = $order_locus{ $locus_order{ $Gene_info{$gene_id}{'Public_name'}} -1 };
    $R_locus = $order_locus{ $locus_order{ $Gene_info{$gene_id}{'Public_name'}} +1 };

    print UPDATE "\nMulti_pt_data : \"$_\"\n";
    print UPDATE "-D Combined\n";
    print UPDATE "-D Remark\n";
    print UPDATE "\nMulti_pt_data : \"$_\"\n";
    if ($L_locus && $R_locus && $allele){
      print UPDATE "Combined Gene \"$Gene_info{$L_locus}{'Gene'}\" 1 Gene \"$center_Gene\" 1 Gene \"$Gene_info{$R_locus}{'Gene'}\"\n";
      print UPDATE "Remark \"Data inferred from $allele, sequence of $Gene_info{$center_Gene}{'CGC_name'} and interpolated map position (which became genetics map)\" Inferred_automatically\n" if exists $Gene_info{$center_Gene}{'CGC_name'};
      print UPDATE "Remark \"Data inferred from $allele, sequence of $Gene_info{$center_Gene}{'Other_name'} and interpolated map position (which became genetics map)\" Inferred_automatically\n" if !exists $Gene_info{$center_Gene}{'CGC_name'};
    }
    elsif (!$L_locus | !$R_locus) {
      $L_locus = "NA" if !$L_locus;
      $R_locus = "NA" if !$R_locus;
      print LOG "ERROR: Multi-pt obj $_ has incomplete flanking loci [(L) $L_locus (R) $R_locus] information for $center_Gene\n";
    }
    elsif (!$allele) {
      $allele  = "NA" if !$allele;
      print LOG "ERROR: Multi-pt obj $_ has NO allele [$allele] information\n";
    }
  }

  print LOG "\n\n";

  my $cmd= "pparse /tmp/updated_multi_pt_flanking_loci_$autoace_version\nsave\nquit\n";
  $ga->upload_database($database, $cmd, "Inferred_multi_pt_data", $log) if !$debug;
}

sub get_flanking_loci {

  # get loci order from last cmp_gmap_with_coord_order_WSXXX.yymmdd.pid file
  my @map_file = glob("/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/cmp_gmap_with_coord_order_$autoace_version*");

  my $count = 0;
  open(MAP, $map_file[-1]) || die $!;

  while(<MAP>){
    chomp;

    if ($_ =~ /^(I|V|X)/){
      my($a, $b, $c, $d, $e) = split(/\s+/, $_);
      my $locus = $c;
      $count++;
      $locus_order{$locus} = $count if $locus;
      $order_locus{$count} = $locus if $locus;
    }
  }
  close MAP;
}


__END__


=head2 NAME - update_inferred_multi_pt.pl

=head3 <USAGE>
 
=head2 Options: [h or help] [db or database]

            -h to display this POD
            -db specify db to upload data
            -d(debug): send email to debugger only


=head3 <DESCRITION>

B< Flowchart of creating inferred multi_pt_data [all these steps are taken care of by script, except steps 2].>

B< There are "before-the-build" and "during-the-build" procedures.>

B<   ################### BEFORE THE BUILD ###################>

B<1.> identify loci which
        (a) have Interpolated_map_position (ie, no Map data)
        (b) have allele info (excluding non specified Tc1 insertion alleles, ie, no mutation info)
        (c) have no mapping_data
        (d) belongs to C. elegans
        (e) are linked to sequence (CDS or Transcript or Pseudogene)
   This should have been generated by geneace_check.pl as a file
   as /wormsrv1/geneace/JAH_DATA/MULTI_PT_INFERRED/loci_become_genetic_marker_for_WSXXX,
   where XXX is the release number of the build about to kick off.
   (this file has already upgraded interpolated_map_positino to Map)

B<2.> Send the generated file in step 1 to JAH for approval of promoted loci.

B<3.> Modified the file of step 1 to include only info for the approved loci from step 2.
   Load this file to Geneace by
   update_inferred_multi_pt.pl -db /wormsrv1/geneace/ -n 123 (this is the release number of
   the build to be started. if it is not available, then there is no promoted loci for that build)
   ----------------------------------------------------------------------------------------------
   Note: this file must already be in Geneace before the build begins, otherwise the map position
   will not be up-to-date.
   ----------------------------------------------------------------------------------------------
B< ###################### DURING THE BUILD #####################>

B<4.> Create inferred multi_pt object for promoted loci from step 3 with minimal "combined results"
   information with the immediate left and right flanking cloned loci as multi_pt A and B
   (see eg, multi_pt object 4134).
   Done by: update_inferred_multi_pt.pl -db /wormsrv1/geneace/ -u

   Step 4 depends on cmp_gmap_with_coord_order_WSXXX.yymmdd.pid file generated by get_interpolated_map.pl.
B<5.> Update flanking loci of all inferred multi-pt objects.





