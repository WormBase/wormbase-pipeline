#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-02-06 12:15:04 $ 

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'}; 
use Wormbase;
use Ace;
use lib "/nfs/team71/worm/ck1/WORMBASE_CVS/scripts/";
use Geneace;
use Getopt::Long;

#- ------------------
# global variables
#--------------------

my ($help, $database, $debug);

GetOptions ("h|help"         => \$help,
	    "db|database=s"  => \$database,
	    "d|debug"        => \$debug,
           );

my $user = `whoami`; chomp $user;
my $rundate = &rundate;
my $autoace = "/wormsrv2/autoace";
my $tace = &tace;   
my $multi_dir = "/wormsrv1/geneace/JAH_DATA/MULTI_PT_INFERRED";
my $log = "/wormsrv2/logs/update_inferred_multi_pt.$rundate";
open(LOG, ">$log") || die $!;
print LOG "# $0 started at ", runtime(), "\n\n";;
print LOG "=============================================================================================\n\n";

#-------------------------------------------------
# which database to use for gmap and data upload
#-------------------------------------------------
my $ckdb ="/nfs/disk100/wormpub/DATABASES/TEST_DBs/CK1TEST";

if (!$database && $user ne "wormpub" && !$debug){
  print "\nTarget database for genetic map and data uploading is $autoace.\n\n";
  print "You need to be wormpub to do this!\n";
  exit(0);
}
elsif (!$database && $debug){
  print "\nUsing genetic maps in $autoace\n\n";
}
elsif (!$database && $user eq "wormpub" && !$debug){print "\nUsing genetic maps in $autoace.\n\nTarget database for data uploading is $autoace\n\n"}
elsif ($database && $database ne $autoace){print "\nUsing genetic maps in $database.\n\n"};

print "Target database for data uploading is $ckdb\n\n" if $debug or $database;

if (!$database){$database = $autoace}

my $db = Ace->connect(-path  => $database,
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

#-----------------
# start working
#-----------------

my (%locus_allele, %locus_order, %order_locus);

&int_map_to_map_loci;
&make_inferred_multi_pt_obj; 
&update_inferred_multi_pt;

#-------------------------
# s u b r o u t i n e s
#-------------------------

sub int_map_to_map_loci {

  my $autoace_version = get_wormbase_version();
  # get a list of "promoted" loci from geneace_check output to here
  my @int_loci = `cat $multi_dir/loci_become_genetic_marker_for_WS$autoace_version*`;
 
  # test file
  #my @int_loci = `cat $multi_dir/inferred_multi_pt_obj_WS117`;

  my $locus;

  foreach (@int_loci){
    chomp;
    if ($_ =~ /^Locus : \"(.+)\"/){
      $locus = $1;
      $locus = $db->fetch(-class => 'Locus',
			  -name  => $locus);
      if ($locus){
        $locus_allele{$locus} = $locus -> Allele(1) if defined $locus -> Allele(1); # grep only the first allele
        print LOG "ERROR: $locus has now NO allele attached . . . corresponding multi-pt obj needs update. . .\n" if !defined $locus -> Allele(1); # 
      }
      else {
	$locus = $db->fetch(-class => 'Gene_name',
                            -name  => $1);
	if (defined $locus -> Other_name_for(1)){
	   print LOG "ERROR: $locus has become an Other_name for ", $locus -> Other_name_for(1), " . . . corresponding multi-pt obj needs update. . . \n";  
        }
      }
    }
  }
}

sub make_inferred_multi_pt_obj {

  # get last multipt obj
  my $multipt = "find Multi_pt_data *";
  my @multi_objs = $db->find($multipt);
  my $last_multi = $multi_objs[-1];
  my $multi = $last_multi -1; $multi++;
  
  # check autoace version
  my $autoace_version = get_wormbase_version_name(); # return WSXX

  # get loci order from cmp_gmap_with_coord_order_WSXXX.yymmdd.pid file
  my @map_file = glob("/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/cmp_gmap_with_coord_order_$autoace_version*");
  my $count = 0;
  open(IN, $map_file[-1]) || die $!;
  while(<IN>){
    chomp;
    if ($_ =~ /^(I|V|X)/){
      my($a, $b, $c, $d, $e) = split(/\s+/, $_); 
      my $locus = $c;
      $count++;
  
      $locus_order{$locus} = $count if $locus; # print "$locus -> $count  (1)\n";
      $order_locus{$count} = $locus if $locus; # print "$count->$locus    (2)\n";
    }
  }

  # write inferred multi_obj acefile
  open(NEW, ">/tmp/inferred_multi_pt_obj_to_make") || die $!;

  foreach (keys %locus_allele){ 
    $multi++;
    my $L_locus = $order_locus{$locus_order{$_}-1};
    my $R_locus = $order_locus{$locus_order{$_}+1};

    print NEW "\n\nLocus : \"$_\"\n";
    print NEW "Multi_point $multi\n";
    print NEW "\n\nMulti_pt_data : $multi\n";
    print NEW "Locus_A \"$_\" \"$locus_allele{$_}\"\n";
    print NEW "Locus \"$_\" \"$locus_allele{$_}\"\n";
    print NEW "Combined Locus \"$L_locus\" 1 Locus \"$_\" 1 Locus \"$R_locus\"\n";
    print NEW "Remark \"Data inferred from $locus_allele{$_}, sequence of $_ and interpolated map position (which became genetics map)\" Inferred_automatically\n";
  }
}


sub update_inferred_multi_pt {
  
  my $query  = "find Multi_pt_data * where remark AND NEXT AND NEXT = \"inferred_automatically\"";
  push( my @inferred_multi_objs, $db->find($query) );

  open(UPDATE, ">/tmp/updated_multi_pt_flanking_loci") || die $!;
  my (@center_locus, $locus);

  foreach (@inferred_multi_objs){
    my $center_locus = $_ -> Combined(5);

    $locus = $db->fetch(-class => 'Locus',
		       -name  => $center_locus );

    my $allele = $locus -> Allele(1); # grep only the first allele
    my $L_locus = $order_locus{($locus_order{$center_locus})-1};
    my $R_locus = $order_locus{$locus_order{$center_locus}+1};

    print UPDATE "\nMulti_pt_data : \"$_\"\n";
    print UPDATE "-D Combined\n";
    print UPDATE "-D Remark\n";
    print UPDATE "\nMulti_pt_data : \"$_\"\n";
    if ($L_locus && $R_locus){
      print UPDATE "Combined Locus \"$L_locus\" 1 Locus \"$center_locus\" 1 Locus \"$R_locus\"\n";
      print UPDATE "Remark \"Data inferred from $allele, sequence of $center_locus and interpolated map position (which became genetics map)\" Inferred_automatically\n";	
    }
    else {
      $L_locus = "NA" if !$L_locus;
      $R_locus = "NA" if !$R_locus;
      print LOG "ERROR: Obj $_ has incomplete flanking loci information: (L) $L_locus (R) $R_locus\n";
    }
  }

  print LOG "\n\n";

  # output a updated multi-pt temp file and upload it to database specified
  
  my $command=<<END;
pparse /tmp/updated_multi_pt_flanking_loci
pparse /tmp/inferred_multi_pt_obj_to_make
save
quit
END

  my $ga = init Geneace();
  $ga->upload_database($ckdb, $command, "Inferred_multi_pt_data", $log) if $debug;
  $ga->upload_database($database, $command, "Inferred_multi_pt_data", $log) if !$debug;
}

print LOG "Make sure that all file parsings are OK . . . . .\n\n";
mail_maintainer("Update inferred multi-pt objects", "ALL", $log) if !$debug;
mail_maintainer("Update inferred multi-pt objects", "ck1\@sanger.ac.uk", $log) if $debug;

print LOG "\n$0 finished at ", runtime(), "\n\n";

__END__


=head2 NAME - update_inferred_multi_pt.pl  

=head3 <USAGE> 
 
=head2 Options: [h or help] [db or database]

            -h to display this POD
            -db specify db to upload data, if not specify, the default db is autoace
            - d(debug): use CK1 testdb for data upload, sends email to ck1 only


=head3 <DESCRIPTION> 

B<Flowchart of creating inferred multi_pt_data [all these steps are taken care of by script, except steps 2 and 3(b)]>

B<1.>  Identify loci which
           (a) have Interpolated_map_position (ie, no Map data)
           (b) have allele info
           (c) have no mapping_data
           (d) belongs to C. elegans
           (e) are linked to sequence (CDS or Transcript or Pseudogene)

B<2.>  Send loci in step 1 to JAH for approval

B<3.>  Check geneace_check.pl generated file(s) "loci_become_genetic_marker_for_WSXXX.yymmdd" at 
       /wormsrv1/geneace/JAH_DATA/MULTI_PT_INFERRED/ for CGC approved ones from step 2.
         
   
---------------------------------------------------------------------------
Note: loci created in step 3 need to be already in the build, otherwise 
the map position will not be up-to-date.
---------------------------------------------------------------------------

B<4.>  (a) This script creates inferred multi_pt object for loci from step 3 with minimal 
        "combined results" information, ie, uses the immediate left and right 
        flanking cloned loci as multi_pt A and B loci * (see eg, multi_pt object
        4134).
    (b) The script also updates all previous inferred multi pt obj for flanking loci 
    (c) Hand check to make sure these loci become yellow-highlighted on the right 
        of the scale bar in Gmap.

    * This info is obtained from cmp_gmap_with_coord_order_WSXXX.yymmdd.pid 
    file generated by get_interpolated_map.pl.

    This is best done: 
         (a) after get_interpolated_gmap.pl during the build is finished AND
         (b) geneace gets updated subsequently.


