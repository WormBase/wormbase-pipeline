#!/usr/local/bin/perl5.8.0 -w
#
# check_DNA.pl
# 
# by Dan Lawson
#
# processes GFF files to make new files which can be used to make agp files 
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2004-03-03 13:52:38 $



use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Getopt::Long;


#################################
# Command-line options          #
#################################

my $test;      # uses test environment
my $quicktest; # same as $test but only runs one chromosome

GetOptions ("test"         => \$test,
	    "quicktest"    => \$quicktest);

# check that -test and -quicktest haven't both been set.  Also if -quicktest is specified, 
# need to make -test true, so that test mode runs for those steps where -quicktest is meaningless
if($test && $quicktest){
  die "both -test and -quicktest specified, only one of these is needed\n";
}
($test = 1) if ($quicktest);


# Set up top level base directory which is different if in test mode
# Make all other directories relative to this
my $basedir   = "/wormsrv2";
$basedir      = glob("~wormpub")."/TEST_BUILD" if ($test); 
my $gffdir    = "$basedir/autoace/CHROMOSOMES";
my $agpdir    = "$basedir/autoace/yellow_brick_road";
my $scriptdir = "$basedir/scripts";

# prepare array of file names and sort names
my @files = (
	  'CHROMOSOME_I.gff',
	  'CHROMOSOME_II.gff',
	  'CHROMOSOME_III.gff',
	  'CHROMOSOME_IV.gff',
	  'CHROMOSOME_V.gff',
	  'CHROMOSOME_X.gff',
	  );

@files = ("CHROMOSOME_III.gff") if ($quicktest);

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
    
    if (($method eq "Genomic_canonical") && ($feature eq "Region")) {
      print OUT "$_\n";
    }
  }
  close GFF;
  close OUT;
  
  # modify to make the clone_acc lists
  &GFF_with_acc("$agpdir/CHROMOSOME_$chromosome.clone_path.gff", "$agpdir/CHROMOSOME_$chromosome.clone_acc.gff" );

}


exit(0);

# this was originally a separate script called only by this one, so folded in. Could be improved for greater efficiency :)
sub GFF_with_acc{
  my $file   = shift;
  my $output = shift;
  my $wormdb = "$basedir/autoace";
  
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

    my $sv;
    if($obj->at('DB_info.Database.EMBL.NDB_SV')){
      ($sv) = $obj->at('DB_info.Database.EMBL.NDB_SV');
    }
    elsif($obj->at('DB_info.Database.GenBank.NDB_SV')){
      ($sv) = $obj->at('DB_info.Database.GenBank.NDB_SV');
    }
    
    # now just want the numerical suffix of sequence version field
    my $suffix = substr($sv,-1,1);  
    substr($sv,-2,2,"");  
    my $acc = $sv;
    print OUT "$_ acc=$acc ver=$suffix\n";
       

    $obj->DESTROY();
    
  }
  close(GFF);
  close OUT;
  
  $db->close;
}


__END__


=pod

=head1 NAME - check_DNA.pl

=head2 USAGE

=over 4

=item check_DNA.pl -[options]

=back

=head1 DESCRIPTION

processes GFF files to make new files which can be used to make agp files.  Creates
two new sets of GFF files which will contain clone path and accession information.

Makes new files in $basedir/autoace/yellow_brick_road/

=back

=head1 MANDATORY arguments:

=over 4

=item none

=back

=head1 OPTIONAL arguments: -test, -quicktest


=over 4

=item -test

Uses test environment in ~wormpub/TEST_BUILD/

=back

=over 4

=item -quicktest

Will only run against one chromosome (for speed) which is CHROMOSOME_III

=back

=head1 AUTHOR Dan Lawson (dl1@sanger.ac.uk) 

=back

=cut

