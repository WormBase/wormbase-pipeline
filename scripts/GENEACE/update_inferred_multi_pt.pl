#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-03-19 15:27:03 $ 

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'}; 
use Wormbase;
use Ace;
use GENEACE::Geneace;
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
my $autoace_version = "WS".get_wormbase_version();
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

#-----------------
# start working
#-----------------

my $db = Ace->connect(-path  => $database,
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

my $ga = init Geneace();
my %Alleles = $ga->get_non_Transposon_alleles($db); # all Geneace alleles which have no Transposon_insertion tag

my (%locus_allele, %locus_order, %order_locus);

my $update = &int_map_to_map_loci;
print $update, "\n";

&make_inferred_multi_pt_obj if $update ne "NA";
#&update_inferred_multi_pt;

#print LOG "Make sure that all file parsings are OK . . . . .\n\n";
##mail_maintainer("Update inferred multi-pt objects", "ALL", $log) if !$debug;
##mail_maintainer("Update inferred multi-pt objects", "ck1\@sanger.ac.uk", $log) if $debug;

#print LOG "\n$0 finished at ", runtime(), "\n\n";

#-------------------------
# s u b r o u t i n e s
#-------------------------

sub int_map_to_map_loci {

  # get a list of "promoted" loci from geneace_check output to here

  # test file
  #my $upt_file = glob -e ("$multi_dir/inferred_multi_pt_obj_WS117");

  #my $upt_file = glob("$multi_dir/loci_become_genetic_marker_for_$autoace_version*");
  my $upt_file = glob("$multi_dir/loci_become_genetic_marker_for_WS123*");
  my ($locus, $allele);

  if (!$upt_file){
    print "\nNO new inferred multi-pt object(s) for $autoace_version . . .\n\n";
    print "Proceed to update flanking loci of previous inferred multi-pt objects . . .\n\n";
    &get_flanking_loci;
    return "NA";
  }
  else {
    print $upt_file, "--------\n";

    open(IN, "$upt_file");
    while(<IN>){  
      chomp;
      if ($_ =~ /^Locus : \"(.+)\"/){
	$locus = $1;
	$locus = $db->fetch(-class => 'Locus',
			    -name  => $locus);

	if ($locus){ 
	  if (defined $locus -> Allele(1)){
	    my @alleles = $locus -> Allele(1); 
	    foreach my $e (@alleles){
	      if (exists $Alleles{$e} ){
		$locus_allele{$locus} = $e; # grep only the allele which has no Transposon_insertion tag
		print "$locus => $e (1)\n";
		last;
	      }
	    }
	  }
	  else {
	    print LOG "ERROR: $locus has now NO allele attached . . . corresponding multi-pt obj needs update. . .\n";
	  }
	}
	else {
	  $locus = $db->fetch(-class => 'Gene_name',
			      -name  => $1);
	  if (defined $locus -> Other_name_for(1)){
	    print LOG "ERROR: $locus has become an Other_name for ", $locus -> Other_name_for(1), " . . . corresponding multi-pt obj needs update. . . \n";  
	  }
	  else {
	    print LOG "ERROR: $locus is NOT in database. . .corresponding multi-pt obj needs update. . . \n";  
	  }
	}
      }
    }
    close IN;
    &get_flanking_loci;
  }
}

sub get_flanking_loci {
  # get loci order from latest cmp_gmap_with_coord_order_WSXXX.yymmdd.pid file
  # my @map_file = glob("/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/cmp_gmap_with_coord_order_$autoace_version*");
  my @map_file = glob("/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/cmp_gmap_with_coord_order_WS121*");
  
  my $count = 0;
  open(MAP, $map_file[-1]) || die $!;
  while(<MAP>){
    chomp;
    if ($_ =~ /^(I|V|X)/){
      my($a, $b, $c, $d, $e) = split(/\s+/, $_); 
      my $locus = $c;
      $count++;
      
      $locus_order{$locus} = $count if $locus; # print "$locus -> $count  (1)\n";
      $order_locus{$count} = $locus if $locus; # print "$count->$locus    (2)\n";
    }
  }
  close MAP;
}

sub make_inferred_multi_pt_obj {
  
  # get last multipt obj
  my $multipt = "find Multi_pt_data *";
  my @multi_objs = $db->find($multipt);
  my $last_multi = $multi_objs[-1];
  my $multi = $last_multi -1; $multi++;
  
  # write inferred multi_obj acefile
  open(NEW, ">/tmp/inferred_multi_pt_obj_$autoace_version") || die $!;
  `chmod 777 /tmp/inferred_multi_pt_obj_$autoace_version`;
  foreach (keys %locus_allele){ 
    $multi++;
    my $L_locus = $order_locus{$locus_order{$_}-1};
    my $R_locus = $order_locus{$locus_order{$_}+1};
    print "\n\nMulti_pt_data : $multi\n";
    print NEW "\n\nLocus : \"$_\"\n";
    print NEW "Multi_point $multi\n";
    print NEW "\n\nMulti_pt_data : $multi\n";
    print NEW "Locus_A \"$_\" \"$locus_allele{$_}\"\n";
    print NEW "Locus \"$_\" \"$locus_allele{$_}\"\n";
    if ($L_locus && $R_locus){
      print NEW "Combined Locus \"$L_locus\" 1 Locus \"$_\" 1 Locus \"$R_locus\"\n";
      print NEW "Remark \"Data inferred from $locus_allele{$_}, sequence of $_ and interpolated map position (which became genetics map)\" Inferred_automatically\n";
    }
    elsif (!$L_locus | !$R_locus){
      print LOG "ERROR: Multi-pt obj $multi has incomplete flanking loci [(L) $L_locus (R) $R_locus] information\n";
    }
    elsif (!$locus_allele{$_}){
      print LOG "ERROR: Multi-pt obj $multi has NO allele [NA] information\n";
    }
  }
}

sub update_inferred_multi_pt {

  my $query  = "find Multi_pt_data * where remark AND NEXT AND NEXT = \"inferred_automatically\"";
  push( my @inferred_multi_objs, $db->find($query) );

  #open(UPDATE, ">/tmp/updated_multi_pt_flanking_loci_$autoace_version") || die $!;
  open(UPDATE, ">/tmp/updated_multi_pt_flanking_loci_WS121") || die $!;
  #`chmod 777 /tmp/updated_multi_pt_flanking_loci_$autoace_version`;
  `chmod 777 /tmp/updated_multi_pt_flanking_loci_WS121`;
  my (@center_locus, $locus);

  foreach (@inferred_multi_objs){
  my $center_locus = $_ -> Combined(5);

    $locus = $db->fetch(-class => 'Locus',
		        -name  => $center_locus );

    my $allele =(); my @alleles =();
    if (defined $locus -> Allele(1)){
      @alleles = $locus -> Allele(1); 
      print "$locus =====> @alleles\n";
      foreach my $e (@alleles){
	if (exists $Alleles{$e} ){
	  $allele = $e; # grep only the allele which has no Transposon_insertion tag 
	  print "$locus => $allele\n";
	  last;
	}
      }
    }
    else {
      print "$locus has no allele===\n";
    }
    my $L_locus = $order_locus{($locus_order{$center_locus})-1};
    my $R_locus = $order_locus{$locus_order{$center_locus}+1};

    print UPDATE "\nMulti_pt_data : \"$_\"\n";
    print UPDATE "-D Combined\n";
    print UPDATE "-D Remark\n";
    print UPDATE "\nMulti_pt_data : \"$_\"\n";
    if ($L_locus && $R_locus && $allele){
      print UPDATE "Combined Locus \"$L_locus\" 1 Locus \"$center_locus\" 1 Locus \"$R_locus\"\n";
      print UPDATE "Remark \"Data inferred from $allele, sequence of $center_locus and interpolated map position (which became genetics map)\" Inferred_automatically\n";	
    }
    elsif (!$L_locus | !$R_locus) {
      $L_locus = "NA" if !$L_locus;
      $R_locus = "NA" if !$R_locus;
      print LOG "ERROR: Multi-pt obj $_ has incomplete flanking loci [(L) $L_locus (R) $R_locus] information\n";
    }
    elsif (!$allele) {
      $allele  = "NA" if !$allele;
      print LOG "ERROR: Multi-pt obj $_ has NO allele [$allele] information\n";
    }
  }

  print LOG "\n\n";

  # output temp file and upload it to database specified and copy $num_multi to $multi_dir

  #`cp /tmp/inferred_multi_pt_obj_$autoace_version $multi_dir` if $update ne "NA"; # if exists
  `cp /tmp/inferred_multi_pt_obj_WS121 $multi_dir` if $update ne "NA"; # if exists

  my $cmd1=<<END;
pparse /tmp/inferred_multi_pt_obj_$autoace_version
pparse /tmp/updated_multi_pt_flanking_loci_$autoace_version
save
quit
END

  my $cmd2=<<END;
pparse /tmp/updated_multi_pt_flanking_loci_$autoace_version
save
quit
END

  $ga->upload_database($ckdb, $cmd1, "Inferred_multi_pt_data", $log) if ($debug && $update ne "NA");
  $ga->upload_database($ckdb, $cmd2, "Inferred_multi_pt_data", $log) if ($debug && $update eq "NA");
  $ga->upload_database($database, $cmd1, "Inferred_multi_pt_data", $log) if (!$debug && $update ne "NA");
  $ga->upload_database($database, $cmd2, "Inferred_multi_pt_data", $log) if (!$debug && $update eq "NA");
}


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
           (b) have allele info (excluding transposon_insertion alleles)
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


