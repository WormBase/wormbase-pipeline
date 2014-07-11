#!/software/bin/perl -w

# script to report stats on the number of curation events for a given species

# types of change reported:
# removed - either a CDS has been converted to a Pseudogene or the Gene has been killed (maybe as part of a Gene Merge).
#         - both the protein ID and its gene sequence name from the previous release have gone.
# changed - a change to A CDS structure has resulted in a new protein or it has been converted to a non-coding isoform.
#         - a new protein ID that was not in the previous release, but its gene sequence name is in the previous and current release and it is not a new isoform.
# new - a new gene has been created or there has been a Gene Split.
#     - a new protein ID was not in the previous release and its gene sequence name was not in the previous release either.
# new isoform - a new isoform has been created.
#             - this is a new isoform sequence name (with a letter) and the protein ID is also new.
# renamed - an existing gene has been renamed but the old protein is still attached to it.
#         - the protein ID still exists, but it is attached to a different gene seqeunce name from the previous release.

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
use DBI;

######################################
# variables and command-line options # 
######################################

my $input;
my ($help, $debug, $test, $verbose, $store, $wormbase, $USER, $release1, $release2, $species);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store"      => \$store,
	    "species=s"  => \$species,
	    "user:s"     => \$USER,
	    "release1=s"	 => \$release1, # WS number to start at
	    "release2=s"	 => \$release2, # WS number to end at
	    );

# always in debug mode
$debug = $USER;


if (!defined $release1) {die "-release1 not specified\n";}

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species,
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

# checks
if (!defined $USER) {$log->log_and_die("-user not specified\n")}
if (!defined $release1) {$log->log_and_die("-release1 not specified\n")}
if (!defined $release2) {$log->log_and_die("-release2 not specified\n")}

# Display help if required
&usage("Help") if ($help);

if (!defined $species) {$species = 'elegans'}

my $prefix = $wormbase->pepdir_prefix . "pep";;

# mysql database parameters
my $dbsn = "DBI:mysql:dbname=wbgene_id;host=shap;port=3303";
my $dbuser = $USER;
my $dbpass = $USER;

my $mysql = DBI -> connect($dbsn, $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to database, $DBI::errstr";



# open an ACE connection
#print "Connecting to Ace\n";
my $database = "/nfs/wormpub/DATABASES/geneace";
my $ace = Ace->connect (-path => $database) || die "cannot connect to database at $database\n";


# release start and end dates 
# NB these are the dates between data-freezes when curation that goes into the next Build is done
my %dates = (
	     200 => ['2009-01-01', '2009-02-01'], # Jan 2009
	     201 => ['2009-02-01', '2009-03-01'], # Feb
	     202 => ['2009-03-01', '2009-04-01'], # Mar
	     203 => ['2009-04-01', '2009-05-01'], # Apr
	     204 => ['2009-05-01', '2009-06-01'], # May
	     205 => ['2009-06-01', '2009-07-01'], # Jun
	     206 => ['2009-07-01', '2009-08-01'], # Jul
	     207 => ['2009-08-01', '2009-09-01'], # Aug
	     208 => ['2009-09-01', '2009-10-01'], # Sep
	     209 => ['2009-10-01', '2009-11-01'], # Oct
	     210 => ['2009-11-01', '2009-12-01'], # Nov
	     211 => ['2009-12-01', '2010-01-01'], # Dec

	     212 => ['2010-01-01', '2010-02-01'], # Jan
	     213 => ['2010-02-01', '2010-03-01'], # Feb
	     214 => ['2010-03-01', '2010-04-01'], # Mar
	     215 => ['2010-04-01', '2010-05-01'], # Apr
	     216 => ['2010-05-01', '2010-06-01'], # May
	     217 => ['2010-06-01', '2010-07-01'], # Jun
	     218 => ['2010-07-01', '2010-08-01'], # Jul
	     219 => ['2010-08-01', '2010-09-01'], # Aug
	     220 => ['2010-09-01', '2010-10-01'], # Sep
	     221 => ['2010-10-01', '2010-11-01'], # Oct
	     222 => ['2010-11-01', '2010-12-01'], # Nov
	     223 => ['2010-12-01', '2011-01-01'], # Dec

	     224 => ['2011-01-01', '2011-02-01'], # Jan
	     225 => ['2011-02-01', '2011-03-01'], # Feb
# bimonthly releases start
	     226 => ['2011-03-01', '2011-05-01'], # Mar-Apr
	     227 => ['2011-05-01', '2011-07-01'], # May-Jun
	     228 => ['2011-07-01', '2011-09-01'], # Jul-Aug
	     229 => ['2011-09-01', '2011-11-01'], # Sep-Oct
	     230 => ['2011-11-01', '2012-01-01'], # Nov-Dec

	     231 => ['2012-01-01', '2012-03-01'], # Jan-Feb
	     232 => ['2012-03-01', '2012-05-01'], # Mar-Apr
	     233 => ['2012-05-01', '2012-07-01'], # May-Jun
	     234 => ['2012-07-01', '2012-09-01'], # Jul-Aug
	     235 => ['2012-09-01', '2012-11-01'], # Sep-Oct
	     236 => ['2012-11-01', '2013-01-01'], # Nov-Dec

	     237 => ['2013-01-01', '2013-03-01'], # Jan-Feb
	     238 => ['2013-03-01', '2013-05-01'], # Mar-Apr
	     239 => ['2013-05-01', '2013-08-01'], # May-Jun-Jul (IWM delayed things)
	     240 => ['2013-08-01', '2013-10-01'], # Aug-Sep
	     241 => ['2013-10-01', '2013-12-01'], # Oct-Nov
	     242 => ['2013-12-01', '2014-02-01'], # Dec-Jan

	     243 => ['2014-02-01', '2014-04-01'], # Feb-Mar
	     244 => ['2014-04-01', '2014-06-01'], # Apr-May (probable dates in the future ...)
	     245 => ['2014-06-01', '2014-08-01'], # Jun-Jul
	     246 => ['2014-08-01', '2014-10-01'], # Aug-Sep
	     247 => ['2014-10-01', '2014-12-01'], # Oct-Nov
	     248 => ['2014-12-01', '2015-02-01'], # Dec-Jan

	     249 => ['2015-02-01', '2015-04-01'], # Feb-Mar
	     250 => ['2015-04-01', '2015-06-01'], # Apr-May
	    );


##########################
# MAIN BODY OF SCRIPT
##########################

my %prot=();
my %seq=();
my %iso=();
my %isoprot=();
my %prev_prot=();
my %prev_seq=();
my %prev_iso=();
my %prev_isoprot=();


my $total_removed=0;
my $total_changed=0;
my $total_new=0;
my $total_new_isoform=0;
my $total_renamed=0;


foreach my $release (200 .. $release2) {
  
  my $count_removed=0;
  my $count_changed=0;
  my $count_new=0;
  my $count_new_isoform=0;
  my $count_renamed=0;
  my $count_changed_to_this_release=0;
  my $count_changed_from_last_release=0;
  my $tmp=0;

  my $file;
  $file = "/nfs/users/nfs_w/wormpub/BUILD/WORMPEP/${prefix}${release}/${prefix}.accession${release}";

  if (open (IN, "<$file")) {
  
    while (my $line = <IN>) {
      if ($line =~ /^\s*$/) {next}
      my @f = split /\s+/, $line;
      my $protein = shift @f;
      foreach my $sequence (@f) {
	# get the sequence name of the gene by removing any isoform letters
	if ($sequence =~ s/(\S+?)([a-z])$/$1/) {
	  $iso{$sequence.$2} = $protein;
	  $isoprot{$protein} = $sequence.$2;
	}
	push @{$prot{$protein}}, $sequence;
	push @{$seq{$sequence}}, $protein;
      }      
    }

    if (%prev_prot) {
      # see what has changed
      # look at proteins in previous release
      foreach my $protein (keys %prev_prot) {
	if (exists $prot{$protein}) {
	  # no change or renamed
	  foreach my $sequence (@{$prev_prot{$protein}}) {
	    #print "checking to see if old sequence $sequence still exists for the current protein $protein in ", @{$prot{$protein}}, "\n" if $verbose && $release >= $release1;
	    if (grep $sequence, @{$prot{$protein}} ) {
	      #print "$release: $protein , $sequence to $protein , $sequence : no change\n" if $verbose && $release >= $release1;
	    } else {
	      print "$release: $protein , $sequence to $protein , - : renamed\n" if $verbose && $release >= $release1;
	      $count_renamed++;
	    }
	  }
	  
	} else {
	  # removed or changed
	  foreach my $sequence (@{$prev_prot{$protein}}) {
	    if (exists $seq{$sequence}) { # this includes CDS converted to non-coding transcript as well as structure changes
	      print "$release: $protein , $sequence to - , $sequence : changed from last release\n" if $verbose && $release >= $release1;
	      $count_changed_from_last_release++;
	    } else {
	      print "$release: $protein , $sequence to - , - : removed\n" if $verbose && $release >= $release1;
	      $count_removed++;
	    }
	  }
	}
	
      }

      # now look at proteins in current release to catch the new genes
      foreach my $protein (keys %prot) {
	if (!exists $prev_prot{$protein}) {
	  # new gene or changed/new isoform
	  foreach my $sequence (@{$prot{$protein}}) {
	    if (exists $prev_seq{$sequence}) {
	      # exclude new isoforms
	      if (exists $isoprot{$protein} && !exists $prev_iso{$isoprot{$protein}}) {
		#print "$release: - , $sequence to $protein , $sequence : new isoform that looked like a change\n" if $verbose && $release >= $release1; # we count the new isoforms properly later on - don't count this twice
		#$new_isoform++;
	      } else {
		print "$release: - , $sequence to $protein , $sequence : changed to this release\n" if $verbose && $release >= $release1; # this includes existing isoforms that are changed
		$count_changed_to_this_release++; # *** this is the form of 'changed' that we are using ***
	      }
	    } else {
	      print "$release: - , - to $protein , $sequence : new gene\n" if $verbose && $release >= $release1;
	      $count_new++;
	    }
	  }
	}
      }

      # now look at the isoforms to see which of these are new
      foreach my $isoform (keys %iso) {
	my $protein = $iso{$isoform};
	if (exists $prev_iso{$isoform}) {
	  if ($prev_iso{$isoform} eq $iso{$isoform}) {
#	    print "$release: - , $isoform to $protein , $isoform : no change isoform\n" if $verbose && $release >= $release1;	      
	  } else {
#	    print "$release: - , $isoform to $protein , $isoform : changed isoform\n" if $verbose && $release >= $release1; # existing isoform changed - just another type of change - don't count this twice
	  }
	} else {
	  if (!exists $prev_prot{$protein}) { # check it isn't an existing CDS renamed to CDS with a 'a'
	    print "$release: - , - to $protein , $isoform : new isoform\n" if $verbose && $release >= $release1;
	    $count_new_isoform++;
	  }
	}
      }
    }
    
    if ($release >= $release1 && $release <= $release2) {
      
      # the two ways of counting changes are nearly the same but one catches changes to non-coding transcripts
      my $count_changed = $count_changed_to_this_release;
      my $new_non_coding_isoforms = $count_changed_from_last_release - $count_changed_to_this_release;

      # the wormpep files only get the removed and new coding genes, so ignore the results we got for those and
      # get the real removed and changed counts from the nameserver - this includes non-coding genes etc.
      if (!exists $dates{$release}) {$log->log_and_die("The release '$release' is not yet set up in the table of release dates.\n")}
      my @dates = @{$dates{$release}};
      my ($startdate, $enddate) = @dates;
      my $full_species = $wormbase->full_name();
      $count_removed = &get_removed($full_species, $startdate, $enddate);
      $count_new     = &get_new    ($full_species, $startdate, $enddate);  

      print "\tWS$release: removed: $count_removed, changed: $count_changed, new: $count_new, new isoform: $count_new_isoform\n";
      
      $total_removed += $count_removed;
      $total_changed += $count_changed;
      $total_new += $count_new;
      $total_new_isoform += $count_new_isoform;
      $total_renamed += $count_renamed;
    }

    %prev_prot = %prot;
    %prev_seq = %seq;
    %prev_iso = %iso;
    %prev_isoprot = %isoprot;

    %prot = ();
    %seq = ();
    %iso = ();
    %isoprot=();

  } else {
    if ($release >= $release1 && $release <= $release2) {
      print "\trelease $release not built\n";

      # the wormpep files only get the removed and new coding genes, so ignore the results we got for those and
      # get the real removed and changed counts from the nameserver - this includes non-coding genes etc.
      if (!exists $dates{$release}) {$log->log_and_die("The release '$release' is not yet set up in the table of release dates.\n")}
      my @dates = @{$dates{$release}};
      my ($startdate, $enddate) = @dates;
      my $full_species = $wormbase->full_name();
      $count_removed = &get_removed($full_species, $startdate, $enddate);
      $count_new     = &get_new    ($full_species, $startdate, $enddate);  

      print "\tWS$release: removed: $count_removed, changed: $count_changed, new: $count_new, new isoform: $count_new_isoform\n";
      
      $total_removed += $count_removed;
      $total_changed += $count_changed;
      $total_new += $count_new;
      $total_new_isoform += $count_new_isoform;
      $total_renamed += $count_renamed;

    }
  }

}


print "Total: removed: $total_removed, changed: $total_changed, new: $total_new, new isoform: $total_new_isoform\n";

# close the ACE connection
$ace->close;

# disconnect from the mysql database
$mysql->disconnect || die "error disconnecting from database", $DBI::errstr;

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);

#################################################################



sub usage {
  print "Script to get stats on curation in a species\n";
  exit(0);
}

#################################################################
# get the number of removed genes in the specified species
sub get_removed {
  my ($full_species, $start_date, $end_date) = @_;

  my %removed;

  my $query="select log.log_what, log.log_when, prim.object_public_id  from identifier_log as log, primary_identifier as prim where log.log_when >= '$start_date' AND log.log_when < '$end_date' AND prim.object_id = log.object_id AND prim.domain_id = 1 AND (log.log_what = 'killed' OR log.log_what = 'mergedTo')";

  my $db_query = $mysql->prepare($query);
  $db_query->execute();
  my $ref_results = $db_query->fetchall_arrayref;
  foreach my $result_row (@$ref_results) {
    my $gene = $result_row->[2];
    my $obj =  $ace->fetch("Gene" => $gene);
    my $gene_species = $obj->Species;
    #print "$gene $gene_species\n";
    if ($gene_species eq $full_species) {
      $removed{$gene} = 1;
    }
  }

  return keys %removed;
}

###################################################################################
# get the number of new genes in the specified species
sub get_new {
  my ($full_species, $start_date, $end_date) = @_;

  my %new;

  my $query="select log.log_what, log.log_when, prim.object_public_id  from identifier_log as log, primary_identifier as prim where log.log_when >= '$start_date' AND log.log_when < '$end_date' AND prim.object_id = log.object_id AND prim.domain_id = 1 AND (log.log_what = 'created' OR log.log_what = 'splitFrom')";

  my $db_query = $mysql->prepare($query);
  $db_query->execute();
  my $ref_results = $db_query->fetchall_arrayref;
  foreach my $result_row (@$ref_results) {
    my $gene = $result_row->[2];
    my $obj =  $ace->fetch("Gene" => $gene);
    if (defined $obj) {   # some very new genes will not yet have been added to geneace - ignore them as they certainly will not be in a recent WS release yet
      my $gene_species = $obj->Species;
      #print "$gene $gene_species\n";
      if ($gene_species eq $full_species) {
	$new{$gene} = 1;
      }
    }
  }

  return keys %new;
}