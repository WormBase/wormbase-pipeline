#!/nfs/disk100/wormpub/bin/perl -w                 
#
# This is to add Confirmed / Predicted Status and RFLP to SNP gff lines as requested by Todd
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2006-11-30 13:38:56 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use Ace;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my $species;

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	   		 "test"       => \$test,
			    "verbose"    => \$verbose,
			    "store:s"      => \$store,
			    "species:s"  => \$species
	    );

$species = 'elegans' unless $species;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
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

##########################
# MAIN BODY OF SCRIPT
##########################

my %SNP;

#load SNP details from table maker query
my $table = $wormbase->table_maker_query($wormbase->autoace, &write_def_file);
while(<$table>) {
	s/\"//g; #"
	next if (/acedb/ or /\/\//);
	my ($snp, $conf, $pred, $rflp, $from_species) = split(/\t/,$_);
	next if (! defined $from_species);
	next unless ($from_species =~ /$species/);
	$SNP{$snp}->{'confirm'} = ($conf or $pred);
	$SNP{$snp}->{'RFLP'} = 1 if ($rflp =~ /\w/);
}

my @chroms = qw(I II III IV V X MtDNA);
my $dir = $wormbase->chromosomes;
my $stat = 0;
foreach my $chrom (@chroms) {
	open(GFF,"<$dir/CHROMOSOME_${chrom}.gff") or $log->log_and_die("cant open $dir/CHROMOSOME_${chrom}.gff");
	open(NEW,">$dir/CHROMOSOME_${chrom}.gff.tmp") or $log->log_and_die("cant open $chrom tmp file\n");
	while( <GFF> ) {
		chomp;	
		print NEW "$_";
		#CHROMOSOME_V    Allele  SNP     155950  155951  .       +       .       Variation "uCE5-508"
		#I       Allele  SNP     126950  126950  .       +       .       Variation "pkP1003"  ;  Status "Confirmed_SNP" ; RFLP "Yes"
		if( /SNP/ and /Allele/) {
			my ($allele) = /Variation \"(.+)\"$/;
			print NEW " ; Status \"",$SNP{$allele}->{'confirm'},"\"" if $SNP{$allele}->{'confirm'};
			print NEW " ; RFLP ", (defined $SNP{$allele}->{'RFLP'}? '"Yes"' : '"No"');
			$stat++;
		}
		print NEW "\n";
	}
	$wormbase->run_command("mv -f $dir/CHROMOSOME_${chrom}.gff.tmp $dir/CHROMOSOME_${chrom}.gff", $log);
}

# Close log files and exit
$log->write_to("\n\nChanged $stat lines\n");
$log->write_to("----------\n\n");

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

sub write_def_file {
	my $def = '/tmp/overload_SNP_GFF.def';
	open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
	 my $txt = <<END;
Sortcolumn 1

Colonne 1
Width 12
Optional
Visible
Class
Class Variation
From 1

Colonne 2
Width 12
Mandatory
Hidden
Show_Tag
From 1
Tag SNP

Colonne 3
Width 12
Optional
Visible
Show_Tag
From 1
Tag Confirmed_SNP

Colonne 4
Width 12
Optional
Visible
Show_Tag
From 1
Tag Predicted_SNP

Colonne 5
Width 12
Optional
Visible
Show_Tag
From 1
Tag RFLP

Colonne 6
Width 12
Optional
Visible
Class
Class Species
From 1
Tag Species
END

	print TMP $txt;
	close TMP;
	return $def;
}

__END__
# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


=pod

=head2 NAME - over_load_SNP_gff.pl

=head1 USAGE

=over 4

=item over_load_SNP_gff.pl  [-options]

=back

This script adds the status (confirmed / predicited) of SNPs and whether are RFLPs to GFF lines.  This is so that the web
page can easily distinguish between them

=over 4

=item None at present.

=back

$0.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
