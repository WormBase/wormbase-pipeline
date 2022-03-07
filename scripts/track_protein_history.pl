#!/software/bin/perl -w
#
# track_protein_history.pl
# 
# by Gary Williams                        
#
# This gets information from the wormpep history data file
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2011-11-24 12:11:53 $
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
	    "store:s"    => \$store,
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

# establish log file.
my $log = Log_files->make_build_log($wormbase);

# in test mode?                                                                                                                                                                                                          
if ($test) {
    print "In test mode\n";
    $log->write_to("In test mode\n");
}


##########################
# MAIN BODY OF SCRIPT
##########################
my $count = "0";
my $wcount = "0";
print "Finding changes in previous and current wormpep files\n\n";
$log->write_to("Finding changes in previous and current wormpep files\n\n");
my $database_version = $wormbase->get_wormbase_version;

my %hash;

# read in current wormpep CDS and protein IDs
get_wormpep($database_version, \%hash, +1);

# read in previous wormpep CDS and protein IDs
get_wormpep($database_version-1, \%hash, -1);

# get the history file hashref
my $wormbase_history = get_wormpep_history($database_version);

# get the changed CDS/protein pairs from the 'wormpep_current' sequence files and look for problems
foreach my $CDS (keys %hash) {
  foreach my $pep (keys %{$hash{$CDS}}) {
    if ($hash{$CDS}{$pep} == 0)  {
	print "$CDS $pep no change\n" if ($verbose);
	$count++;
      my $status = has_current_history($wormbase_history, $CDS);
      if ($status == 0) {print "$CDS/$pep is in 'wormpep${database_version}' as an unchanged CDS and protein, but the history file thinks $CDS is not live\n";
			 $log->write_to("$CDS/$pep is in 'wormpep${database_version}' as an unchanged CDS and protein, but the history file thinks $CDS is not live\n");
			 $count--;
			 $wcount++;
      }
	
    }
    elsif ($hash{$CDS}{$pep} == 1)  {
      print "$CDS $pep appeared\n" if ($verbose);
      my @cds_history = @{$wormbase_history->{$CDS}};
      foreach my $version (@cds_history) {
	my ($wormpep_id, $version1, $version2) = @{$version};
	my $found = 0;
	if ($wormpep_id eq $pep) { # found when it appeared
	  $found = 1;
	  if ($version1 == $database_version) {
	    $found = 2; # appeared in this version
	  }
	  if ($found == 0) {print "$CDS/$pep is in 'wormpep${database_version}' as newly appeared, but the history file has no record of it\n";
			    $log->write_to("$CDS/$pep is in 'wormpep${database_version}' as newly appeared, but the history file has no record of it\n");
			    $wcount++;
	  }
	  if ($found == 1) {print "$CDS/$pep is in 'wormpep${database_version}' as newly appeared, but the history file says it is not new in this Build\n";
			    $log->write_to("$CDS/$pep is in 'wormpep${database_version}' as newly appeared, but the history file says it is not new in this Build\n");
			    $wcount++;
	  }
	}
      }
    }
    elsif ($hash{$CDS}{$pep} == -1) {
      print "$CDS $pep lost\n" if ($verbose);
      my @cds_history = @{$wormbase_history->{$CDS}};
      foreach my $version (@cds_history) {
	my ($wormpep_id, $version1, $version2) = @{$version};
	my $found = 0;
	if ($wormpep_id eq $pep) { # found it
	  $found = 1;
	  if ($version2 == $database_version) {
	    $found = 2; # lost in this version
	  }
	  if ($found == 0) {print "$CDS/$pep is in 'wormpep${database_version}' as lost, but the history file has no record of it\n";
			    $log->write_to("$CDS/$pep is in 'wormpep${database_version}' as lost, but the history file has no record of it\n");
			    $wcount++;
	  }
	  if ($found == 1) {print "$CDS/$pep is in 'wormpep${database_version}' as lost, but the history file says it is still current\n";
			    $log->write_to("$CDS/$pep is in 'wormpep${database_version}' as lost, but the history file says it is still current\n");
			    $wcount++;
	  }
	}
      }
    }
    else {
	print "$CDS $pep not found at all??!"
    }
  }
}


# get the changed CDS/protein pairs from the history file and look for problems
foreach my $CDS (keys %{$wormbase_history}) {
  my @cds_history = @{$wormbase_history->{$CDS}};
  foreach my $version (@cds_history) {
      my ($pep, $version1, $version2) = @{$version};
      my $pep_tag = "wormpep=$pep";
    if ($version1 == $database_version) { # new in this version
	if (!exists $hash{$CDS}{$pep_tag}) {print "$CDS/$pep is in the history file as new but 'wormpep${database_version}' has no record of it even in the previous version\n";
					$log->write_to("$CDS/$pep is in the history file as new but 'wormpep${database_version}' has no record of it even in the previous version\n");
					$wcount++;
	}
	if ($hash{$CDS}{$pep_tag} == -1)   {print "$CDS/$pep is in the history file as new but 'wormpep${database_version}' thinks it is lost\n";
					$log->write_to("$CDS/$pep is in the history file as new but 'wormpep${database_version}' thinks it is lost\n");
					$wcount++;
	}
	if ($hash{$CDS}{$pep_tag} == 0)    {print "$CDS/$pep is in the history file as new but 'wormpep${database_version}' thinks it is unchanged\n";
					$log->write_to("$CDS/$pep is in the history file as new but 'wormpep${database_version}' thinks it is unchanged\n");
					$wcount++;
	}

    } elsif (!defined $version2) { # still live
	unless (defined $hash{$CDS}{$pep_tag}) {print "$CDS/$pep is in the history file as unchanged but 'wormpep${database_version}' has no record of it even in the previous version\n";
					$log->write_to("$CDS/$pep is in the history file as unchanged but 'wormpep${database_version}' has no record of it even in the previous version\n");
					$wcount++;
	}
	if ($hash{$CDS}{$pep_tag} == -1)   {print "$CDS/$pep is in the history file as unchanged but 'wormpep${database_version}' thinks it is lost\n";
					$log->write_to("$CDS/$pep is in the history file as unchanged but 'wormpep${database_version}' thinks it is lost\n");
					$wcount++;
	}
	if ($hash{$CDS}{$pep_tag} == 1)    {print "$CDS/$pep is in the history file as unchanged but 'wormpep${database_version}' thinks it has just appeared\n";
					$log->write_to("$CDS/$pep is in the history file as unchanged but 'wormpep${database_version}' thinks it has just appeared\n");
					$wcount++;
	}

    } elsif ($version2 == $database_version) { # lost in this version
	if (!exists $hash{$CDS}{$pep_tag}) {print "$CDS/$pep is in the history file as lost but 'wormpep${database_version}' has no record of it even in the previous version\n";
					$log->write_to("$CDS/$pep is in the history file as lost but 'wormpep${database_version}' has no record of it even in the previous version\n");
					$wcount++;
	}
	if ($hash{$CDS}{$pep_tag} == 0)    {print "$CDS/$pep is in the history file as lost but 'wormpep${database_version}' thinks it is unchanged\n";
					$log->write_to("$CDS/$pep is in the history file as lost but 'wormpep${database_version}' thinks it is unchanged\n");
					$wcount++;
	}
	if ($hash{$CDS}{$pep_tag} == 1)    {print "$CDS/$pep is in the history file as lost but 'wormpep${database_version}' thinks it has just appeared\n";
					$log->write_to("$CDS/$pep is in the history file as lost but 'wormpep${database_version}' thinks it has just appeared\n");
					$wcount++;
	}
    }
  }
}
$log->write_to("$count entries unchanged\n");
$log->write_to("$wcount entries have thrown a warning\n");
$log->write_to("Finished.\n");
$log->mail();
print "Finished.\n" if ($verbose);
exit(0);






##############################################################
#
# Subroutines
#
##############################################################

sub get_wormpep {

  my ($version, $hashref, $diff) = @_;

  #get the  CDS and protein IDs in the wormpep sequence file
  open( PEP, $wormbase->basedir."/WORMPEP/wormpep${version}_noah/wormpep$version" );
  my $CDS;
  my $pep;
  while (my $line = <PEP>) {
    #>2L52.1 CE32090
    ($CDS, $pep) = ($line =~ />(\S+)\s+(\S+)/);
    if (defined $CDS) {
      #print "$CDS $pep added\n";
      $hashref->{$CDS}{$pep} += $diff;
    }
  }
  close PEP;

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
# get the wormpep protein IDs used by previous version of the gene

sub get_wormpep_history {

  my ($database_version) = @_;


  # get database version file
  
  my $data_file = $wormbase->basedir . "/WORMPEP/wormpep${database_version}_noah/wormpep.history$database_version";

  my %wormpep_history;

  open (HIST, "< $data_file") || die "Can't open $data_file\n";
  while (my $line = <HIST>) {
    chomp $line;
    my @f = split /\s+/, $line;
    # $cds_id, $wormpep_id, $version1, $version2 <- ($version2 is undef if current version)
    my $cds_id = shift @f;
    push @{$wormpep_history{$cds_id}}, [@f];
  }
  close(HIST);


  return \%wormpep_history;

}
##########################################
# returns true if the CDS has a current history
# returns false if the CDS has been made into an isoform (or pseudogene, etc.)
sub has_current_history {
  my ($history, $CDS_name) = @_;

  # if the CDS_name is not a valid key then this CDS_name has been
  # made into a isoform (or pseudogene, etc.) and no details are known
  # of previous wormpep IDs
  if (!exists $history->{$CDS_name}) {return 0;}

  # if the last line in the history file for this CDS_name has only
  # one value, then there is still a current CDS of this name
  my @cds_history = @{$history->{$CDS_name}};
  my $last_result = pop @cds_history;
  #print "@{$last_result}\n" if ($verbose);
  if (scalar @{$last_result} == 2) {return 1;} 

  # so there isn't a current CDS of this name
  return 0;
}

__END__

