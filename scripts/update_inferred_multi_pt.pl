#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-01-23 16:15:40 $ 

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

my ($help, $database);

GetOptions ("h|help"         => \$help,
	    "db|database=s"  => \$database,
           );

my $user = `whoami`; chomp $user;
my $rundate = &rundate;
my $autoace = "/wormsrv2/autoace";
my $tace = &tace;   
my $multi_dir = "/wormsrv1/geneace/JAH_DATA/MULTI_PT_INFERRED";
my $log = "/wormsrv2/logs/update_inferred_multi_pt.$rundate";
open(LOG, ">$log") || die $!;

if (!$database && $user ne "wormpub"){
  print "\nTarget database for uploading data is $autoace.\n\n";
  print "You need to be wormpub to do this!\n";
  exit(0);
}

elsif (!$database && $user eq "wormpub"){
  print "\nTarget database for uploading data is $autoace.\n\n";
}
elsif ($database && $database ne $autoace){print "\nTarget database for uploading data is $database.\n\n"};

if (!$database){$database = $autoace}


my $db = Ace->connect(-path  => $database,
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

#-----------------
# start working
#-----------------

my %locus_allele = &int_map_to_map_loci;

my($locus_order, $order_locus) = &make_inferred_multi_pt_obj; # return 2 hash refs
my %locus_order = %$locus_order;
my %order_locus = %$order_locus;

&update_inferred_multi_pt;

#-------------------------
# s u b r o u t i n e s
#-------------------------

sub int_map_to_map_loci {

  # get a list of "promoted" loci from geneace_check output to here
  my @int_loci = `cat $multi_dir/loci_become_genetic_marker*`;
  my ($locus, %Locus_allele);

  foreach (@int_loci){
    if ($_ =~ /^Locus : \"(.+)\"/){$locus = $1}
    if ($_ =~ /^-D.+/){
      $locus = $db->fetch(-class => 'Locus',
			  -name  => $locus);
      my $allele = $locus -> Allele(1); # grep only the first allele
      $Locus_allele{$locus} = $allele;
    }
  }
  return %Locus_allele;
}

sub make_inferred_multi_pt_obj {

  # get last multipt obj
  my $multipt = "find Multi_pt_data *";
  my @multi_objs = $db->find($multipt);
  my $last_multi = $multi_objs[-1];
  my $multi = $last_multi -1; $multi++;

  # check autoace version
  my $version = get_wormbase_version_name(); # return WSXX

  # get loci order from cmp_gmap_with_coord_order_WS117.yymmdd.pid file
  my @map_file = `ls /wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/cmp_gmap_with_coord_order_$version*`;
  my $count = 0;
  my (%locus_order, %order_locus);

  open(IN, $map_file[-1]) || die $!;
  while(<IN>){
    chomp;
    my($a, $b, $c, $d, $e) = split(/\s+/, $_); 
    my $locus = $c;
    $count++;
    $locus_order{$locus} = $count if $locus;
    $order_locus{$count} = $locus if $locus;
  }
 
  # write inferred multi_obj acefile
  open(NEW, ">/tmp/temp1") || die $!;

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
  return \%locus_order, \%order_locus;
}


sub update_inferred_multi_pt {
  
  my $query  = "find Multi_pt_data * where remark AND NEXT AND NEXT = \"inferred_automatically\"";
  push( my @inferred_multi_objs, $db->find($query) );
  
  my @Error_multi;
  open(UPDATE, ">/tmp/temp2") || die $!;
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
      push(@Error_multi, $_);
    }
  }

  # output a updated multi-pt temp file and upload it to database specified
 
  print LOG "\nERROR: Check Multi_pt_data object(s) @Error_multi\nfor incomplete flanking loci information\n\n";
  
  my $command=<<END;
pparse /tmp/temp1
pparse /tmp/temp2
save
quit
END

  my $ga = init Geneace();
  $ga->upload_database("/nfs/disk100/wormpub/DATABASES/TEST_DBs/CK1TEST", $command, "Inferred_multi_pt_data", $log);
  
}


__END__

=head2 NAME - update_inferred_multi_pt.pl  

=head3 <USAGE> 
 
=head2 Options: [h or help] [db or database]

            -h to display this POD
            -db specify db to upload data, if not specify, the default db is autoace


=head3 <DESCRIPTION> 

B<Flowchart of creating inferred multi_pt_data [all these steps are taken care of by script, except steps 2 and 3(b)]>

B<1.>  Identify loci which
           (a) have Interpolated_map_position (ie, no Map data)
           (b) have allele info
           (c) have no mapping_data
           (d) belongs to C. elegans
           (e) are linked to sequence (CDS or Transcript or Pseudogene)

B<2.>  Send loci in step 1 to JAH for approval

B<3.>  Upgrade Interpolated_map_info to Map for approved loci from step 2.
   
---------------------------------------------------------------------------
Note: loci created in step 3 need to be already in the build, otherwise 
the map position will not be up-to-date.
---------------------------------------------------------------------------

B<4.>  (a) Create inferred multi_pt object for loci from step 3 with minimal 
        "combined results" information with the immediate left and right 
        flanking cloned loci as multi_pt A and B (see eg, multi_pt object
        4134). 
    (b) Make sure these loci become yellow-highlighted on the right 
        of the scale bar in Gmap.

    This info is obtained from cmp_gmap_with_coord_order_WSXXX.yymmdd.pid 
    file generated by get_interpolated_map.pl.

    This is best done
         (a) after get_interpolated_gmap.pl during the build is finished AND
         (b) geneace gets updated subsequently.


