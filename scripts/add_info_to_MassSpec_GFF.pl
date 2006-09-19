#!/usr/local/bin/perl5.8.0 -w
#
# add_info_to_MassSpec_GFF.pl                           
# 
# by Gary Williams                         
#
# This add some extraneous information to the MassSpec peptides lines in the GFF file
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2006-09-19 08:42:11 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
#use Ace;
#use Sequence_extract;
#use Coords_converter;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"      => \$store,
	    );

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

#################################
# Set up some useful paths      #
#################################


# Set up top level base directories (these are different if in test mode)
my $gff_dir         = $wormbase->gff;         # AUTOACE GFF
###my $gff_dir = glob("~wormpub/DATABASES/current_DB/CHROMOSOMES/");

my $ace_dir = $wormbase->autoace;
###my $ace_dir = glob("~wormpub/DATABASES/current_DB/");

# other paths
my $tace            = $wormbase->tace;        # TACE PATH

# open an ACE connection to parse details
print "Connecting to Ace\n";
my $db = Ace->connect (-path => $ace_dir,
                       -program => $tace) || die "cannot connect to database at $wormbase->database('current')\n";

###################################
# get the protein matches of the sequences
###################################

my $cmd1 = "Query Find Mass_spec_peptide\nshow -a\nquit";

my %matches;
my $id;

print "Finding Mass_spec data\n";
open (TACE, "echo '$cmd1' | $tace $ace_dir |");
while (<TACE>) {
  chomp;
  next if (/acedb\>/);
  next if (/\/\//);
  if (/Mass_spec_peptide\s+:\s+\"(\S+)\"/) {
    $id = $1;
  } elsif (/\"(.+)\"\s+Protein\s+\"(.+)\"/) {
    $matches{$id}{$2}{$1} = 1;	# count the number of experiments this peptide has been seen in
  }
}
close TACE;




##########################
# MAIN BODY OF SCRIPT
##########################
my $count;

# get the set of CDS, history versions and resulting protein IDs
my $protein_history_aref = &get_protein_history;   

# loop through the chromosomes
  my @chromosomes = qw( I II III IV V X );                            # chromosomes
  foreach my $chromosome (@chromosomes) {
    print "Reading chromosome $chromosome\n" if ($verbose);

# loop through the GFF file
    my @f;
    open (GFF, "<$gff_dir/CHROMOSOME_${chromosome}.gff") || die "Failed to open gff file $gff_dir/CHROMOSOME_${chromosome}.gff\n";
    open (OUT, ">$gff_dir/CHROMOSOME_${chromosome}.gff.new") || die "Failed to open gff file $gff_dir/CHROMOSOME_${chromosome}.gff.new\n";
###    open (OUT, ">./CHROMOSOME_${chromosome}.gff.new") || die "Failed to open gff file ./CHROMOSOME_${chromosome}.gff.new\n";
    while (my $line = <GFF>) {
      chomp $line;
      if ($line =~ /^#/ || $line !~ /\S/) {
	print OUT "$line\n";
        next;
      }
      @f = split /\t/, $line;
      my $id;

# is this a MassSpec peptide line?
      if ($f[1] eq 'mass_spec_genome') {
	# get the ID name
	($id) = ($f[8] =~ /Target \"Mass_spec_peptide:(\S+)\"/);

	if (exists $matches{$id}) { # for this peptide id
	  my %times_observed;
	  my $proteins = "";
	  my $cdss = "";

	  foreach my $prot (keys %{ $matches{$id} }) {
	    if ($proteins ne "") {$proteins .= " ";}
	    $proteins .= $prot;

	    # get the CDS name for this protein
	    my $prot_obj = $db->fetch(Protein => $prot);
	    if (! defined $prot_obj) {
	      die "Can't fetch Protein object for $prot\n";
	    }
	    my $cds = $prot_obj->Corresponding_CDS;
# if this is an old CDS not in autoace any more then we have to look in the historical records file
	    if (! defined $cds) { 
#	      print STDERR "can't find CDS for protein $prot\n";
	      $cds = get_previous_wormpep_ids($prot, $protein_history_aref);
	    }
	    $prot_obj->DESTROY();
	    if ($cdss ne "") {$cdss .= " ";}
	    $cdss .= $cds;

	    foreach my $experiment (keys %{ $matches{$id}{$prot} }) {
	      $times_observed{$experiment} = 1;	# count the number of unique experiments that this peptide has been seen in
	    }
	  }
	  $line .= " ; Note \"$id\"";
	  $line .= " ; Protein_matches \"$proteins\"";
	  $line .= " ; CDS_matches \"$cdss\"";
	  my $times_observed = (keys %times_observed); # count the unique experiments
	  $line .= " ; Times_observed \"$times_observed\"";
	  $count++;		# count the number of lines changed for statistics at the end
	  print "$line\n" if ($verbose);
	}
      }

# write out the line
      print OUT "$line\n";

# end of GFF loop
    }

# close files
    close (GFF);
    close (OUT);

# end of chromosome loop
  }

# copy new GFF files over
  foreach my $chromosome (@chromosomes) {
    system("mv -f $gff_dir/CHROMOSOME_${chromosome}.gff.new $gff_dir/CHROMOSOME_${chromosome}.gff");
  }

# Close log files and exit
$log->write_to("\n\nStatistics\n");
$log->write_to("----------\n\n");
$log->write_to("Changed $count lines\n");


$log->mail();
print "Finished.\n" if ($verbose);
exit(0);






##############################################################
#
# Subroutines
#
##############################################################

##########################################
# get the set of CDS, history versions and resulting protein IDs
# my %protein_history = &get_protein_history;
                                                                                                                                          
sub get_protein_history {
  my @protein_history;
                            
  my $release = $wormbase->get_wormbase_version;
                                                                                                              
  my $data_file = "/nfs/disk100/wormpub/BUILD/WORMPEP/wormpep$release/wormpep.history$release";
 
  open (HIST, "< $data_file") || die "Can't open $data_file\n";
  while (my $line = <HIST>) {
    chomp $line;
    my @f = split /\s+/, $line;
    # $cds_id, $wormpep_id, $version1, $version2 <- ($version2 is undef if current version)
    push @protein_history, [@f];
  }
  close(HIST);
 
  return \@protein_history;
}

##########################################
# get the CDS name for a wormpep protein ID
# $cds = get_previous_wormpep_ids($protein_name, $protein_history_aref);
 
sub get_previous_wormpep_ids {
 
  my ($protein_name, $protein_history_aref) = @_;
 
  my @protein_history = @{$protein_history_aref};
 
  my @wormpep_ids;
  my @versions;
 
  my ($protein) = ($protein_name =~ /WP:(\S+)/);

  foreach my $line (@protein_history) {
    my ($cds_id, $wormpep_id, $version1, $version2) = @{$line};
    if ($wormpep_id eq $protein) {
      return $cds_id;
    }
  }
 
  print STDERR "Still can't find the CDS for protein $protein_name\n";
  return "unknown";
}
 

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

=head2 NAME - add_info_to_MassSpec_GFF.pl

=head1 USAGE

=over 4

=item  add_info_to_MassSpec_GFF [-options]

=back

Todd wanted to have some informative data added to the GFF results for mass-spec peptides.
This script gets the information from the autoace database and adds this information to the GFF records.

add_info_to_MassSpec_GFF.pl MANDATORY arguments:

=over 4

=item None at present.

=back

add_info_to_MassSpec_GFF.pl  OPTIONAL arguments:

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

=item Gary Williams (gw3@sanger.ac.uk)

=back

=cut
