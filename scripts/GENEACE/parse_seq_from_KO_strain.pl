#!/usr/local/bin/perl5.8.0 -w
# 
#
# Last updated by: $Author: ck1 $
# Last updated on: $Date: 2004-03-19 10:48:58 $

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Geneace;

######################
# global variables   # 
######################
print @INC;
__END__
my $tace = "/nfs/disk100/acedb/RELEASE.DEVELOPMENT/bin.ALPHA_5/tace";
my $ga = init Geneace();
my $geneace = $ga->geneace();

# open a connection to database

my $db = Ace->connect(-path  => $geneace,
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

my $query = "Find strain VC*; Remark";
push(my @VC, $db->find($query));

foreach (@VC){
  my $remark = $_ -> Remark(1);
  my @allele =(); my $locus =();
  @allele = $_ -> Allele(1) if defined $_ -> Allele(1);
  $locus = $_ -> Gene if defined $_ -> Gene; 
  foreach my $e (@allele){
    if ($e =~ /^gk\d+|^ok\d+/){
      if ($remark =~ /^([A-Z0-9]+\.\d+|[A-Z0-9]+\.\d+\w)\.\s+.+/){
	my $cds = $1;
	if ($locus =~ /tag-\d+/){
	  print "\nLocus : \"$locus\"\n" if $locus =~ /tag-\d+/;
	  print "Allele  \"$e\"\n";
	  #print "Gene \"$locus\"\n" if $locus =~ /tag-\d+/;
	  print "CDS \"$cds\" Inferred_automatically \"from C. elegans Knockout consortium strain info\"\n";
	  print "Strain \"$_\"\n";
	}	
      }
    }
  }
}
  
