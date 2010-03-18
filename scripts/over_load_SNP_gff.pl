#!/software/bin/perl -w                 
#
# This is to add Confirmed / Predicted Status and RFLP to SNP gff lines as requested by Todd
#
# Last updated by: $Author: mh6 $     
# Last updated on: $Date: 2010-03-18 17:12:51 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use Ace;
use strict;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($species, $gff_file);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "file:s"     => \$gff_file
	    );

$species = 'elegans' unless $species;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( '-debug'   => $debug,
                             '-test'    => $test,
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

# need to add mol_change details to GFF lines
&get_mol_changes;

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

my @chroms = $wormbase->get_chromosome_names(-mito => 1);
my $dir = $wormbase->chromosomes;
my $stat = 0;
foreach my $chrom (@chroms) {
    my $file = ($gff_file or "$dir/${\$wormbase->chromosome_prefix}${chrom}.gff");
    open(GFF,"<$file") or $log->log_and_die("cant open $file");
    open(NEW,">$file.tmp") or $log->log_and_die("cant open $file tmp file\n");
    while( <GFF> ) {
	chomp;	
	print NEW "$_";
	#CHROMOSOME_V    Allele  SNP     155950  155951  .       +       .       Variation "uCE5-508"
	#I       Allele  SNP     126950  126950  .       +       .       Variation "pkP1003"  ;  Status "Confirmed_SNP" ; RFLP "Yes"
	if( /SNP/ and /Allele/) {
	    my ($allele) = /Variation \"(.+)\"$/;
	    print NEW " ; Status \"",$SNP{$allele}->{'confirm'},"\"" if $SNP{$allele}->{'confirm'};
	    print NEW " ; RFLP ", (defined $SNP{$allele}->{'RFLP'}? '"Yes"' : '"No"');
	    print NEW " ; Mutation_type \"".$SNP{$allele}->{'mol_change'}."\"" if $SNP{$allele}->{'mol_change'};
	    $stat++;
	}
	print NEW "\n";
    }
    $wormbase->run_command("mv -f $file.tmp $file", $log);
    last if $gff_file;
}

##################
# Check the files
##################

if($gff_file){
    $log->write_to("Not checking ad hoc file\n");
}
else { 
    foreach my $chrom (@chroms) {
	$wormbase->check_file("$dir/${\$wormbase->chromosome_prefix}${chrom}.gff", $log,
			      minsize => 1500000,
			      lines => ['^##',
					"^${\$wormbase->chromosome_prefix}${chrom}\\s+\\S+\\s+\\S+\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+"],
			      );
    }
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


sub get_mol_changes{
    $log->write_to("getting molecular change info\n");
    my $table = $wormbase->table_maker_query($wormbase->autoace,&write_mol_change_def );

    my %interested = ('Nonsense'    => 1,
		      'Frameshift'  => 2,
		      'Splice_site' => 3,
		      'Missense'    => 4,
		      );

    while(<$table>) {
	chomp;
	s/\"//g; #"
	next if (/acedb/ or /\/\//);
	my @data = split(/\s+/,$_);
	if($interested{$data[1]}){
	    if( !(defined $SNP{$data[0]}->{'mol_change'}) or  
		 ($interested{$data[1]} < $interested{ $SNP{$data[0]}->{'mol_change'} }) ){
		$SNP{$data[0]}->{'mol_change'} = $data[1];
	    }
	}
    }
    close $table;
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


sub write_mol_change_def {
	my $def = '/tmp/overload_SNP_GFF_mol_chng.def';
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
Mandatory 
Hidden 
Class 
Class CDS 
From 1 
Tag Predicted_CDS 
 
Colonne 4 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # Missense 
 
Colonne 5 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # Nonsense 
 
Colonne 6 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # Splice_site 
 
Colonne  7
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # Frameshift 
 
Colonne 8 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # Intron 
 
Colonne 9 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # Coding_exon 
 
Colonne 10 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # Promoter 
 
Colonne 11 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # UTR_3 
 
Colonne 12 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # UTR_5 
 
Colonne 13 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3 
Tag  HERE  # Regulatory_feature 
 
Colonne 14 
Width 12 
Optional 
Visible 
Show_Tag 
Right_of 3  
Tag  HERE  # Genomic_neighbourhood 
 
 

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
