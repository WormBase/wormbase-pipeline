#!/usr/local/bin/perl5.6.1 -w
#
# Check_DNA.pl v1.0
# dl
# 2000-04-26
#

########################################
# iniatialise                          #
########################################

$|=1;
#use IO::Handle;
#use Getopt::Std;
use strict;
use lib '/wormsrv2/scripts';
use Wormbase;
use Ace;

my $gffdir    = "/wormsrv2/autoace/CHROMOSOMES";
my $agpdir    = "/wormsrv2/autoace/yellow_brick_road";
my $scriptdir = "/wormsrv2/scripts";

# prepare array of file names and sort names
my @files = (
	  'CHROMOSOME_I.gff',
	  'CHROMOSOME_II.gff',
	  'CHROMOSOME_III.gff',
	  'CHROMOSOME_IV.gff',
	  'CHROMOSOME_V.gff',
	  'CHROMOSOME_X.gff',
	  );

my @gff_files = sort @files; 
undef @files; 

############################################################
# loop through each GFF file                               #
############################################################

my $debug = 1;

foreach my $file (@gff_files) {
    
    next if ($file eq "");
    my ($chromosome) = $file =~ (/CHROMOSOME\_(\S+)\./);

    # CHROMOSOME_X    Genomic_canonical       Sequence        1196305 1224112 .       +       .       Sequence "C46H3"

    open (OUT, ">$agpdir/CHROMOSOME_$chromosome.clone_path.gff") || die "can't open output file '$agpdir/CHROMOSOME_$chromosome.clone_path.gff'\n";
    open (GFF, "<$gffdir/$file") || die "can't open gff file '$gffdir/$file'\n";
    while (<GFF>) {
	chomp;
	next if ($_ =~ /^\#/);
	my ($name,$method,$feature,$start,$stop,$score,$strand,$other) = split (/\t/,$_);
	
	if (($method eq "Genomic_canonical") && ($feature eq "Sequence")) {
	    print OUT "$_\n";
	}
    }
    close GFF;
    close OUT;
    
    # modify to make the clone_acc lists
    # system ("$scriptdir/GFF_with_accessions $agpdir/CHROMOSOME_$chromosome.clone_path.gff > $agpdir/CHROMOSOME_$chromosome.clone_acc.gff");
    &GFF_with_acc("$agpdir/CHROMOSOME_$chromosome.clone_path.gff", "$agpdir/CHROMOSOME_$chromosome.clone_acc.gff" );

}




exit;

# this was originally a separate script called only by this one, so folded in. Could be improved for greater efficiency :)
sub GFF_with_acc
  {
        my $file   = shift;
	my $output = shift;
	my $wormdb = "/wormsrv2/autoace";
	
	my $db = Ace->connect(-path=>$wormdb) || do { print "Connection failure: ",Ace->error; die();};
	open (OUT, ">$output") or die "cant write output to $output :$!\n";
	open (GFF, "<$file") || die "Can't open GFF file\n\n";
	while (<GFF>) {
	  
	  next if (/^\#/);
	  
	  chomp;
	  
	  my @gff = split (/\t/,$_);
	  
	  my ($gen_can) = $gff[8] =~ /Sequence \"(\S+)\"/; 
	  
	  my $obj = $db->fetch(Sequence=>$gen_can);
	  if (!defined ($obj)) {
	    print "Could not fetch sequence '$gen_can'\n";
	    next;
	  }
	  
	  my @acc = $obj->DB_info->row();
	  
	  print OUT "$_ acc=$acc[3] ver=$acc[4]\n";
	  
	  $obj->DESTROY();
	  
	}
	close(GFF);
	close OUT;

	$db->close;
      }
