#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-03-23 17:33:17 $ 

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'}; 
use Wormbase;
use Ace;
use GENEACE::Geneace;
use Getopt::Long;

#--------------------
# global variables
#--------------------

my ($help, $database, $debug, $new, $update);

GetOptions ("h|help"         => \$help,
	    "db|database=s"  => \$database,
	    "d|debug=s"      => \$debug,
	    "n|new=s"        => \$new,
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


my $db = Ace->connect(-path  => $database,
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

my $ga = init Geneace();
my %Alleles = $ga->get_non_Transposon_alleles($db); # all Geneace alleles which have no Transposon_insertion tag

my (%locus_allele, %locus_order, %order_locus);


# CGC approved promoted loci
my $new_multi = glob("$multi_dir/loci_become_genetic_marker_for_WS$new*");


#-----------------------------------
#   load promoted loci to geneace
#-----------------------------------

&convert_int_map_to_map_for_promoted_loci if $new;

#--------------------------------------------------------------------
#   make inferred multi-pt obj for promoted loci during the build
#   update flanking loci of all inferred multi-pt obj
#--------------------------------------------------------------------

if ($update){ 
  &make_inferred_multi_pt_obj; 
  &update_inferred_multi_pt; 
}

print LOG "Make sure that all file parsings are OK . . . . .\n\n";
##mail_maintainer("Update inferred multi-pt objects", "ALL", $log) if !$debug;

print LOG "\n$0 finished at ", runtime(), "\n\n";
mail_maintainer("Update inferred multi-pt objects", "$debug\@sanger.ac.uk", $log) if $debug;


#-------------------------
# s u b r o u t i n e s
#-------------------------

sub convert_int_map_to_map_for_promoted_loci {

  # load file with CGC approved promoted loci to Geneace
  my $upload = "pparse $new_multi\nsave\nquit";
  $ga->upload_database($ga->geneace(), $upload, "pseudo_marker_loci", $log);
}

sub make_inferred_multi_pt_obj {

  my @new_multi = `cat $new_multi`;
  my ($locus, $allele);
  foreach (@new_multi){
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
	   #   print "$locus => $e (1)\n";
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
  &get_flanking_loci;

  # get last multipt obj
  my $multipt = "find Multi_pt_data *";
  my @multi_objs = $db->find($multipt);
  my $last_multi = $multi_objs[-1];
  my $multi = $last_multi -1; $multi++;
  my $error =0;

  # write inferred multi_obj acefile
  open(NEW, ">$multi_dir/inferred_multi_pt_obj_$autoace_version") || die $!;

  foreach (keys %locus_allele){ 

    $multi++;
    my $L_locus = $order_locus{$locus_order{$_}-1}; 
    my $R_locus = $order_locus{$locus_order{$_}+1};

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
      $error = 1;
      $L_locus = "NA" if !$L_locus;
      $R_locus = "NA" if !$R_locus;
      print LOG "ERROR: Multi-pt obj $multi has incomplete flanking loci [(L) $L_locus (R) $R_locus] information for $_\n";
    }
    elsif (!$locus_allele{$_}){
      print LOG "ERROR: Multi-pt obj $multi has NO allele [NA] information\n";
    }
  }
  print LOG "\nInfo: for incomplete flanking loci: they should be the end marker of a chromosome.\n"; 
 
  # load $multi_dir/inferred_multi_pt_obj_$autoace_version to autoace
  my $upload = "pparse $multi_dir/inferred_multi_pt_obj_$autoace_version\nsave\nquit";
  $ga->upload_database($database, $upload, "pseudo_mapping_data", $log) if !$debug;
}

sub get_flanking_loci {
  # get loci order from last cmp_gmap_with_coord_order_WSXXX.yymmdd.pid file
  my $map = $new - 1;
  my @map_file = glob("/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/cmp_gmap_with_coord_order_WS$map*");
  
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


sub update_inferred_multi_pt {

  my $query  = "find Multi_pt_data * where remark AND NEXT AND NEXT = \"inferred_automatically\"";
  push( my @inferred_multi_objs, $db->find($query) );

  open(UPDATE, ">/tmp/updated_multi_pt_flanking_loci_$autoace_version") || die $!;
  `chmod 777 /tmp/updated_multi_pt_flanking_loci_$autoace_version`;

  my (@center_locus, $locus);

  foreach (@inferred_multi_objs){
  my $center_locus = $_ -> Combined(5);

    $locus = $db->fetch(-class => 'Locus',
		        -name  => $center_locus );

    my $allele =(); my @alleles =();
    if (defined $locus -> Allele(1)){
      @alleles = $locus -> Allele(1); 
      foreach my $e (@alleles){
	if (exists $Alleles{$e} ){
	  $allele = $e; # grep only the allele which has no Transposon_insertion tag 
	  last;
	}
      }
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
      print LOG "ERROR: Multi-pt obj $_ has incomplete flanking loci [(L) $L_locus (R) $R_locus] information for $center_locus\n";
    }
    elsif (!$allele) {
      $allele  = "NA" if !$allele;
      print LOG "ERROR: Multi-pt obj $_ has NO allele [$allele] information\n";
    }
  }

  print LOG "\n\n";

  my $cmd= "pparse /tmp/updated_multi_pt_flanking_loci_$autoace_version\nsave\nquit";
  $ga->upload_database($database, $cmd, "Inferred_multi_pt_data", $log) if !$debug;
}


__END__


=head2 NAME - update_inferred_multi_pt.pl  

=head3 <USAGE> 
 
=head2 Options: [h or help] [db or database]

            -h to display this POD
            -db specify db to upload data
            -d(debug): send email to debugger only


=head3 <DESCRIPTION> 

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





