#!/usr/local/bin/perl5.6.1 -w
#
# whats_new_in_EMBL
#
# dl
#
# checks EMBL for new EST or mRNA entries
#
# Last updated by: $Author: krb $                      
# Last updated on: $Date: 2003-09-05 16:41:16 $        

$|=1;
use strict;
use Getopt::Long;
use lib "/wormsrv2/scripts/";   
use Wormbase;
use Time::Local;

##############################
# Variables                  #
##############################

my ($help, $debug, $days);
my $maintainers = "All";
our $log;

GetOptions ("help"     => \$help,
            "debug=s"  => \$debug,
            "days=s"   => \$days);
	   

# help page	   
&usage("Help") if ($help);
	   
# no debug name
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

&create_log_files;



##########################
# MAIN BODY OF SCRIPT
##########################

my $date     = &get_date;

my $new_elegans_mRNA = 0;
my $new_elegans_EST  = 0;
my $new_nematode_EST = 0;
my $non_elegans_ESTs = 0;
my $new_EMBL_CDS     = 0;



# Elegans mRNA entries
open (NEW_SEQUENCES, "/usr/local/pubseq/bin/getz -c \'([emblnew-Division:inv] & [emblnew-Molecule:rna] & [emblnew-org:caenorhabditis elegans*] & [emblnew-DateCreated#$date:]) | ([emblrelease-Division:inv] & [emblrelease-Molecule:rna] & [emblrelease-org:caenorhabditis elegans*] & [emblrelease-DateCreated#$date:])\' |");

while (<NEW_SEQUENCES>){
  chomp;
  $new_elegans_mRNA = $_;
}
close NEW_SEQUENCES;

print "There are $new_elegans_mRNA new mRNA entries since $date\n" if ($debug);



# Elegans EST entries
open (NEW_SEQUENCES, "/usr/local/pubseq/bin/getz -c \'([emblnew-Division:est] & [emblnew-Molecule:rna] & [emblnew-org:caenorhabditis elegans*] & [emblnew-DateCreated#$date:]) | ([emblrelease-Division:est] & [emblrelease-Molecule:rna] & [emblrelease-org:caenorhabditis elegans*] & [emblrelease-DateCreated#$date:])\' |");
while (<NEW_SEQUENCES>){
  chomp;
  $new_elegans_EST = $_;
}
close NEW_SEQUENCES;

print "There are $new_elegans_EST new EST entries since $date\n" if ($debug);



# All Nematatode EST entries
open (NEW_SEQUENCES, "/usr/local/pubseq/bin/getz -c \'([emblnew-Division:est] & [emblnew-Molecule:rna] & [emblnew-Taxon:Nematoda*] & [emblnew-DateCreated#$date:]) | ([emblrelease-Division:est] & [emblrelease-Molecule:rna] & [emblrelease-Taxon:Nematoda*] & [emblrelease-DateCreated#$date:])\' |");
while (<NEW_SEQUENCES>){
  chomp;
  $new_nematode_EST = $_;
}
close NEW_SEQUENCES;

print "There are $new_nematode_EST new nematode EST entries since $date\n" if ($debug);

$non_elegans_ESTs += ($new_nematode_EST - $new_elegans_EST);
print "There are $non_elegans_ESTs non C. elegans ESTs since $date\n" if ($debug);



# New EMBL CDS entries
open (NEW_SEQUENCES, "/usr/local/pubseq/bin/getz -c \'([emblnew-Division:inv] & [emblnew-Molecule:DNA*] & [emblnew-DateCreated#20030101:] & [emblnew-org:Caenorhabditis elegans*] ! [emblnew-Keywords:HTG*]) | ([emblrelease-Division:inv] & [emblrelease-Molecule:DNA*] & [emblrelease-DateCreated#$date:] & [emblrelease-org:Caenorhabditis elegans*] ! [emblrelease-Keywords:HTG*])\' |");
while(<NEW_SEQUENCES>){
  chomp;
  $new_EMBL_CDS = $_;
}
close NEW_SEQUENCES;
#! [emblrelease-Keywords:HTG*]
print "There are $new_EMBL_CDS new C. elegans non-WormBase CDS entries since $date\n" if ($debug);

##############################
# mail $maintainer report    #
##############################

print LOG "EMBL entries created since $date\n";
print LOG "There are $new_elegans_mRNA new C. elegans mRNA entries\n";
print LOG "There are $new_elegans_EST new C. elegans EST entries\n";
print LOG "There are $non_elegans_ESTs new non-C. elegans nematode ESTs since $date\n\n";
print LOG "There are $new_EMBL_CDS new C. elegans non-WormBase CDS entries\n";


close LOG;

&mail_maintainer("WormBase Report: whats new in EMBL",$maintainers,$log);


##############################
# a extremidade              #
##############################

exit(0);

##############################
# subroutines                #
##############################


sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate     = `date +%y%m%d`; chomp $rundate;
  $log        = "/wormsrv2/logs/$script_name.$rundate.$$";

  open (LOG, ">$log") or die "cant open $log";
  print LOG "$script_name\n";
  print LOG "started at ",`date`,"\n";
  print LOG "=============================================\n";
  print LOG "\n";

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

######################################################

sub get_date {

  my $date;

  # Has -days been specified, if so calculate date for getzc query
  # using that date

  if($days){    
    my %month2num = (
		     'Jan' => '01',
		     'Feb' => '02',
		     'Mar' => '03',
		     'Apr' => '04',
		     'May' => '05',
		     'Jun' => '06',
		     'Jul' => '07',
		     'Aug' => '08',
		     'Sep' => '09',
		     'Oct' => '10',
		     'Nov' => '11',
		     'Dec' => '12'
		    );
    
      
    my $t = time - ($days * 86400);
    my $play = scalar localtime $t;
    my ($month,$day,$year) = $play =~ (/\S+\s+(\S+)\s+(\d+)\s+\S+\s+(\d+)/);
    
    $date = $year . $month2num{$month};
    ($date .= "0") if ($day < 10);
    $date .= $day;
  }
  # Otherwise find out date when last build was started
  else{
    my $last_release_date = &get_wormbase_release_date("short");
    print "Last release date was $last_release_date\n" if $debug;
    print LOG "\nLast release date was $last_release_date\n";
    my ($day, $month, $year) = split(/\//,$last_release_date);
    my $time = timelocal("0","0","0","$day","$month","$year");
    my $start_build_time = $time - (10*86400);
    my ($start_day, $start_month, $start_year) = (localtime($start_build_time))[3,4,5];
    
    $start_day   = sprintf("%02d", $start_day % 100);
    $start_month = sprintf("%02d", $start_month % 100);
#    $start_year  = sprintf("%02d", $start_year % 100);
    $start_year += 1900;
    print "Therefore last build started roughly $start_day/$start_month/$start_year\n" if $debug;
    print LOG "Therefore last build started roughly $start_day/$start_month/$start_year\n\n";
    
    $date  = "$start_year"."$start_month"."$start_day";
    
  }

  return ($date);
}




__END__

=pod

=head1 NAME 

whats_new_in_EMBL.pl

=head2 USAGE

whats_new_in_EMBL.pl [-options]

whats_new_in_EMBL will return the number of created entries in EMBL since
a given number of days ago. If no date is specified, it will attempt to calculate
the difference since the start of the last build (calculated as 10 days before
the date of the creation of the current_release FTP directory).

The script calculates how many new C. elegans EST and mRNA sequences and also how
many new non C. elegans nematode EST sequences.


whats_new_in_EMBL.pl MANDATORY arguments:

=over 4

=item none

=back

whats_new_in_EMBL.pl OPTIONAL arguments:

=item B<-days [int]>, No. of days to check for 

=item B<-debug [username]>, Verbose/debug mode

=item B<-help>, Help

=head2 Dependencies

no dependencies.

=head2 AUTHOR

Dan Lawson (B<dl1@sanger.ac.uk>)

=cut
