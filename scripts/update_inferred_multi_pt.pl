#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-06-10 16:43:27 $ 

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use GENEACE::Geneace;
use Getopt::Long;

#--------------------
# global variables
#--------------------

my ($help, $database, $debug, $load, $update);

GetOptions ("h|help"         => \$help,
	    "db|database=s"  => \$database,   # can specify db for debug purpose (required also when not debug)
	    "d|debug=s"      => \$debug,
	    "l|load"         => \$load,       # load CGC approved pseudo map markers
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
                                                    # as many of these are silent and convey no info of gene function

my (%Gene_info, %gene_id_allele, %locus_order, %order_locus);

# file with CGC approved promoted loci
my $new_multi_file = `ls $multi_dir/loci_become_genetic_marker_for_$autoace_version`;
chomp $new_multi_file;

#---------------------------------------------------------------------------------------------------------------------
#   load promoted marker Gene (loci) to geneace for the current build
#   NOTE: the pseudo marker in loci_become_genetic_marker_for_WSxxx file should already be approved by Jonathan first
#---------------------------------------------------------------------------------------------------------------------

&load_pseudo_markers_to_geneace if $load;

#--------------------------------------------------------------------
#   make inferred multi-pt obj for promoted loci during the build
#   update flanking loci of existing old inferred multi-pt obj
#--------------------------------------------------------------------

if ($update){

  %Gene_info = $ga -> gene_info($database);   # spcified db

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

 sub load_pseudo_markers_to_geneace {

  # load file with CGC approved promoted pseudo marker loci to Geneace
  my $upload = "pparse $new_multi_file\nsave\nquit";
  $ga->upload_database($database, $upload, "pseudo_marker_loci", $log);  # upload to specifed database (in non-debug mode, this should be geneace)

}

sub make_inferred_multi_pt_obj {   # run during the build, when approved pseudo markers are available

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

  open(UPDATE, ">$multi_dir/updated_multi_pt_flanking_loci_$autoace_version") || die $!;

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

  my $cmd= "pparse $multi_dir/updated_multi_pt_flanking_loci_$autoace_version\nsave\nquit\n";
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
            -d(debug) user: send email to debugger only
            -l upload approved pseudo markers 


=head3 <DESCRITION>

Flowchart of creating inferred multi_pt_data [all these steps are taken care of by script, except steps 2].
All files are created in the directory /wormsrv1/geneace/JAH_DATA/MULTI_PT_INFERRED/

There are "before-the-build" and "during-the-build" procedures.

   ################### BEFORE THE BUILD ###################

   run the script geneace_check.pl -c ps, which

B<1.>
 identify loci which
        (a) have Interpolated_map_position (ie, no Map data)
        (b) have allele info (excluding non specified Tc1 insertion alleles, ie, no mutation info)
        (c) have no mapping_data
        (d) belongs to C. elegans
        (e) are linked to sequence (CDS or Transcript or Pseudogene)
   And gives CGC roughly 2 weeks time to approve the pseudo markers on the list
   /wormsrv1/geneace/JAH_DATA/MULTI_PT_INFERRED/loci_become_genetic_marker_for_WSxxx (so the release number is
   for the NEXT build; this acefile has already upgraded interpolated_map_positino to Map).

B<2.>
 Send the geneace_check output in step 1 (loci_become_genetic_marker_for_WSxxx: next release number) to JAH
   for approval of promoted loci to be used for the next build.

B<3.>
 Modified the file of step 1 (loci_become_genetic_marker_for_WSxxx: current release number -
   created at start of last build) to include only info for the approved loci.
   Load this file (loci_become_genetic_marker_for_WSxxx; current release number) to Geneace by
   update_inferred_multi_pt.pl -db /wormsrv1/geneace/ -l
   -------------------------------------------------------------------------------------------------------
   Note: this file of step 1 (with only approved ones, loci_become_genetic_marker_for_WSxxx: current release number)
         must already be loaded into Geneace before the build begins, otherwise autoace will not have new pseudo markers
         with map positions (appear in /wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/cmp_gmap_with_coord_order_WS120.yymmdd.psid),
         and later during the build, loci in this approved file will still be picked up and new multip-pt obj generated,
         since the script looks for loci_become_genetic_marker_for_WSxxx (current release number) which is available.

         Consequence(1): an inferred multi-pt obj will be created accordingly, but the pseudo marker, say abc-1 will have
                         no left and right flanking loci that Jonathan (CGC) wants, because abc-1 is not yet on the list
                         of markers without rev. physicals
                         (/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/cmp_gmap_with_coord_order_WSxxx.yymmdd.psid).
         Consequence(2): you will need to remove those inferred multi-pt obj. with incomplete information.
         TO AVOID this: if CGC has not yet approved any marker from the list, then you need to DELETE the file
                        (loci_become_genetic_marker_for_WSxxx; current release number) created in step 1 at start of
                        the last build.
   -------------------------------------------------------------------------------------------------------

   ###################### DURING THE BUILD #####################

B<4.>
 Create inferred multi_pt object for promoted loci based on loci_become_genetic_marker_for_WSxxx (current release)
   from step 3 with minimal "combined results" information with the immediate left and right flanking cloned loci
   as multi_pt A and B (see eg, multi_pt object 4134).
   This is done by: "update_inferred_multi_pt.pl -db /wormsrv2/autoace/ -u" (if approved pseudo markers are available),
   which creates the file inferred_multi_pt_obj_WSxxx, otherwise, it is absent.

   Step 4 depends on cmp_gmap_with_coord_order_WSxxx.yymmdd.pid file generated by get_interpolated_map.pl during
   the build

B<5.>
 Update flanking loci of all inferred multi-pt objects (also done by "update_inferred_multi_pt.pl -db /wormsrv2/autoace/ -u",
   which generates updated_multi_pt_flanking_loci_WSxxx file). This is run even no new approved pseudo markers are there for
   current build.

B<6.>
 Both inferred_multi_pt_obj_WSxxx and updated_multi_pt_flanking_loci_WSxxx will also be uploaded back to geneace by
   running /nfs/team71/worm/ck1/WORMBASE_CVS/scripts/GENEACE/load_related_data_from_Build_to_geneace.pl
   after step "update_inferred_multi_pt.pl -db /wormsrv2/autoace/ -u" of the build guide

