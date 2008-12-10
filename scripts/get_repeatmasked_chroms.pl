#!/usr/local/bin/perl5.8.0 -w
#
# agp2ensembl
#
# Cared for by Simon Potter
# (C) GRL/EBI 2001
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#
# modified for reading in .agp files for worm ensembl


use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

use lib '/software/worm/lib/bioperl-live';
use lib '/software/worm/ensembl/ensembl/modules';
use Bio::EnsEMBL::DBSQL::DBAdaptor;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $species, $wormbase,$database);

my $agp;
my $out_dir;


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"      => \$store,
	    "output:s"   => \$out_dir,
	    "database:s" => \$database,
	    "species:s"  => \$species,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species,
			     );
}
# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}
# establish log file.
my $log = Log_files->make_build_log($wormbase);

$species = $wormbase->species;
unless ($database =~ /$species/) {
	$log->log_and_die("are you sure you have the right species / database combo\n Species :$species\nDatabase : $database\n");
}

#################################
# Set up some useful paths      #
#################################



$out_dir = $wormbase->chromosomes unless $out_dir;
die "cant write to $out_dir\t$!\n" unless (-w $out_dir );

# open connection to EnsEMBL DB
my $dbobj;

$log->write_to("Connecting to worm_dna\n");

$dbobj = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
						       '-host'   => 'ia64d',
						       '-user'   => 'wormro',
						       '-dbname' => $database
						      )
  or die "Can't connect to Database $database";



$log->write_to("Building chromosomes\n");

print STDERR "Outputting     ";

my $sa=$dbobj->get_SliceAdaptor();

my $chr_assembly;
if ($wormbase->assembly_type eq 'contig') {
	$wormbase->run_command(" rm -f $out_dir/${species}_masked.dna", $log);
	$wormbase->run_command(" rm -f $out_dir/${species}_softmasked.dna", $log);
}else {
	$chr_assembly = 1;
}

foreach my $seq ( @{$sa->fetch_all('toplevel')}) {
	my $name=$seq->seq_region_name();

	my($masked, $soft);
	if( $chr_assembly ) {
		my $outfile = "$out_dir"."/${name}_masked.dna";
		my $outfile2 = "$out_dir"."/${name}_softmasked.dna";
		open($masked,">$outfile") or $log->log_and_die("cant write $outfile :$!\n");
		open($soft,">$outfile2") or $log->log_and_die("cant write $outfile2 :$!\n");
	}
	else {
		my $outfile = "$out_dir"."/${species}_masked.dna";
		my $outfile2 = "$out_dir"."/${species}_softmasked.dna";
		open($masked,">>$outfile") or $log->log_and_die("cant write $outfile :$!\n");
		open($soft,">>$outfile2") or $log->log_and_die("cant write $outfile2 :$!\n");

	}
	print_seq($masked,$name,$seq);
	print_seq($soft,$name,$seq,1);
}

unless ($chr_assembly){
	$wormbase->run_command("gzip -9 $out_dir/${species}_masked.dna", $log);
	$wormbase->run_command("gzip -9 $out_dir/${species}_softmasked.dna", $log);
}

sub print_seq {
  my ($file,$name,$seq,$softmasked) = @_;
 
  $log->write_to("\twriting chromosome $name\n");
  print $file ">$name 1 ",$seq->seq_region_length,"\n";

  my $width = 50;
  my $start_point = 0;
  my $sequence=$softmasked?
         $seq->get_repeatmasked_seq(undef,1)->seq()
	 :$seq->get_repeatmasked_seq->seq();
  
  while ( $start_point + $width < length( $sequence ) ) {
    print $file substr($sequence, $start_point, $width ),"\n";
    $start_point += $width;
  }

  print $file substr($sequence, $start_point),"\n";

  # crude linux hack (will not work on Alphas)
  my $dno=fileno($file);
  my $fname = (`lsof -Fn -a -d $dno -p $$`)[1];
  $fname = substr($fname,1);
  $fname =~ s/\s\(.*\)//;
  ###  lifted from perlmonks ###

  close $file;
  system("gzip -9 $fname") if ($chr_assembly);
}

$log->write_to("Done\n");


$log->mail();
print "Finished.\n" if ($verbose);
exit(0);






##############################################################
#
# Subroutines
#
##############################################################



##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################




# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head1 NAME

agp2ensembl.pl

=head1 SYNOPSIS

get_repeatmasked_chroms.pl -agp ~wormpipe/Elegans/WSXXX.agp

=head1 DESCRIPTION

extracts RepeatMasked sequence from DNA database using specified AGP file

=head1 OPTIONS

    -agp     agp file

=head1 CONTACT

ar2@sanger.ac.uk

=cut
