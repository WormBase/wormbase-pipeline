#!/usr/local/bin/perl5.6.1 -w
#
# whats_new_in_EMBL
#
# dl
#
# checks EMBL for new EST or mRNA entries
#
# Last updated by: $Author: pad $                      
# Last updated on: $Date: 2006-03-03 11:49:27 $        

use strict;
use Getopt::Long;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Time::Local;
use Modules::Features;

##############################
# Variables                  #
##############################

my ($help, $debug, $days, $wormbase, $test);
my $maintainers = "All";
my $getz   = "/usr/local/pubseq/bin/getzc"; # getz binary

GetOptions ("help"     => \$help,
            "debug=s"  => \$debug,
            "days=s"   => \$days,
	    "test"     => \$test);
	   

# help page	   
&usage("Help") if ($help);

$wormbase = Wormbase->new( -debug => $debug,
                             -test => $test,
                             );
	   

# no debug name
if ($debug) {
    ($maintainers = $debug . '\@sanger.ac.uk');
}

my $log = Log_files->make_build_log($wormbase);

#########################
## MAIN BODY OF SCRIPT ##
#########################

our $buildtime       = 21;                                                # length of build in days
my $date             = &get_date;

my $new_elegans_mRNA = 0;                                                 # Counts for the various classes
my $new_elegans_EST  = 0;
my $new_nematode_EST = 0;
my $non_elegans_ESTs = 0;
my $new_EMBL_CDS     = 0;

my $outdir  = ($debug) ? '/tmp/' : '/nfs/disk100/wormpub/camace_orig';

our $acc;                 # EMBL accession
our $id;                  # EMBL ID
our $sv;                  # EMBL sequence version
our $def = "";            # EMBL description, needs to be initialised to ""
our $protid;              # EMBL Protein_ID
our $protver;             # EMBL Protein_ID version
#our $org;                 # EMBL species (not used mh6)
our $verbose;

our $longtext;

########################
# Elegans mRNA entries #
########################

open (NEW_SEQUENCES, "$getz -c \'([embl-div:inv] & [embl-mol:*rna] & [embl-org:Caenorhabditis elegans] & [embl-DateCreated#$date:])\' |");
while (<NEW_SEQUENCES>){
    chomp;
    $new_elegans_mRNA = $_;
}
close NEW_SEQUENCES;

print "There are $new_elegans_mRNA new mRNA entries since $date\n" if ($debug);

if ($new_elegans_mRNA > 0) {

    open (OUT_ACE,  ">$outdir/new_mRNA.ace");
    open (NEW_SEQUENCES, "$getz -f \"acc\" \'([embl-div:inv] & [embl-mol:*rna] & [embl-org:Caenorhabditis elegans] & [embl-DateCreated#$date:])\' |");
    while (<NEW_SEQUENCES>){
	my $seq='';
        # reset vars
	$def = ''; $id = ''; $acc = ''; $sv = ''; $protid = ''; $protver ='';


	chomp;
	($acc) = (/^AC\s+(\S+)\;/);
        next if ($acc eq "");
        print "Parsing EMBL accession: '$acc'\n" if ($verbose);
                            
	$longtext = 1;
	print OUT_ACE "\nLongText : \"$acc\"\n";

        # pfetch each sequence individually
        open (LOOK, "/usr/local/pubseq/bin/pfetch -F $acc |");
        while (<LOOK>) {
	    chomp;
            print if ($verbose);
                                                                                                                                                    
	    # sequence cleaning <----- SEQ
            if (/^\s/) {
                s/\s+//g;
                s/\d+//g;
                print OUT_ACE  "$_\n";
		$seq.=$_;
            }

	    unless ($longtext == 0) {
		print OUT_ACE "$_\n";
	    }
            # grab various details out of EMBL entry
            if (/^ID\s+(\S+)/)                         {$id  = $1;}  # <-------- ID 
            if (/^SV\s+(\S+\.\d+)/)                    {$sv  = $1;}
            if (/^DE\s+(.+)/)                          {$def = $def . " " . $1;}
            if (/^FT\s+\/protein_id=\"(\S+)\.(\d+)\"/) {$protid = $1; $protver = $2;}
            if (/^SQ/) {


		print OUT_ACE "//\n";
		print OUT_ACE "***LongTextEnd***\n";

		# remove any offending '>' from def line. This is required by transcriptmasker.pl
                $def =~ s/\>//g;
                                                                                                                                                    
		print OUT_ACE "\nSequence : \"$acc\"\n";
		print OUT_ACE "Database EMBL NDB_AC $acc\n";
		print OUT_ACE "Database EMBL NDB_ID $id\n";
		print OUT_ACE "Database EMBL NDB_SV $sv\n";
		print OUT_ACE "Protein_id $acc $protid $protver\n";
		print OUT_ACE "DB_annotation EMBL $acc\n";
		print OUT_ACE "Species \"Caenorhabditis elegans\"\n";
		print OUT_ACE "Title \"$def\"\nMethod NDB\n";
		print OUT_ACE "\nDNA \"$acc\"\n";
		
		# end  LongText object
		$longtext = 0;
	    }


	}
	close LOOK;

	my $feature=Features::annot($seq,$id);
	if ($feature){
	    chomp $feature;
	    print OUT_ACE "\n",$feature;
	}
    }
    # close filehandles
    close SEQUENCES;
    close OUT_ACE;
}

#exit(0);


#######################
# Elegans EST entries #
#######################

open (NEW_SEQUENCES, "$getz -c \'([embl-div:est] & [embl-mol:*rna] & [embl-org:Caenorhabditis elegans] & [embl-DateCreated#$date:])\' |");
while (<NEW_SEQUENCES>){
  chomp;
  $new_elegans_EST = $_;
}
close NEW_SEQUENCES;

print "There are $new_elegans_EST new EST entries since $date\n" if ($debug);


if ($new_elegans_EST > 0) {

    open (OUT_ACE,  ">$outdir/new_EST.ace");
    open (NEW_SEQUENCES, "$getz -f \"acc\" \'([embl-div:est] & [embl-mol:*rna] & [embl-org:Caenorhabditis elegans] & [embl-DateCreated#$date:])\' |");
    while (<NEW_SEQUENCES>){
	my $seq='';
        # reset vars
	$def = ""; $id = ""; $acc = ""; $sv = ""; $protid = ""; $protver ="";


	chomp;
	($acc) = (/^AC\s+(\S+)\;/);
        next if ($acc eq "");
        print "Parsing EMBL accession: '$acc'\n" if ($verbose);
                                                                                                                                                    
        # pfetch each sequence individually
        open (LOOK, "/usr/local/pubseq/bin/pfetch -F $acc |");
        while (<LOOK>) {
            print if ($verbose);
             

            # sequence stuff                                                                                                                                       
            if (/^\s/) {
                s/\s+//g;
                s/\d+//g;
                print OUT_ACE  "$_\n";
		$seq.=$_;
            }
            # grab various details out of EMBL entry
            if (/^ID\s+(\S+)/)                         {$id  = $1;}
            if (/^SV\s+(\S+\.\d+)/)                    {$sv  = $1;}
            if (/^DE\s+(.+)/)                          {$def = $def." ".$1;}
            if (/^FT\s+\/protein_id=\"(\S+)\.(\d+)\"/) {$protid=$1; $protver=$2;}
            if (/^SQ/) {
		# remove any offending '>' from def line. This is required by transcriptmasker.pl
                $def =~ s/\>//g;
                                                                                                                                                    
		print OUT_ACE "\nSequence : \"$acc\"\n";
		print OUT_ACE "Database EMBL NDB_AC $acc\n";
		print OUT_ACE "Database EMBL NDB_ID $id\n";
		print OUT_ACE "Database EMBL NDB_SV $sv\n";
		print OUT_ACE "Protein_id $acc $protid $protver\n";
		print OUT_ACE "Species \"Caenorhabditis elegans\"\n";
		print OUT_ACE "Title \"$def\"\nMethod NDB\n";
		print OUT_ACE "\nDNA \"$acc\"\n";
	    }
	}
	close LOOK;
	my $feature=Features::annot($seq,$id);
	if ($feature){
	    chomp $feature;
	    print OUT_ACE "\n",$feature;
	}

    }
    # close filehandles
    close SEQUENCES;
    close OUT_ACE;
}


############################
# All Nematode EST entries #
############################

open (NEW_SEQUENCES, "$getz -c \'([embl-div:est] & [embl-mol:*rna] & [embl-Taxon:Nematoda] & [embl-DateCreated#$date:])\' |");
while (<NEW_SEQUENCES>){
  chomp;
  $new_nematode_EST = $_;

}
close NEW_SEQUENCES;

print "There are $new_nematode_EST new nematode EST entries since $date\n" if ($debug);

$non_elegans_ESTs += ($new_nematode_EST - $new_elegans_EST);
print "There are $non_elegans_ESTs non C. elegans ESTs since $date\n" if ($debug);



########################
# New EMBL CDS entries #
########################

open (NEW_SEQUENCES, "$getz -c \'([embl-div:inv] & ([embl-mol:genomic dna*] | [embl-mol:unassigned dna]) & [embl-DateCreated#$date:] & [embl-org:Caenorhabditis elegans] ! [embl-Keywords:HTG*]) & ([embl-FtKey:cds] > embl))\' |");
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

$log->write_to("New EMBL sequence entries created since: $date\n");
$log->write_to("$new_elegans_mRNA new C. elegans mRNA entries\n");
$log->write_to("$new_elegans_EST new C. elegans EST entries\n");
$log->write_to("$non_elegans_ESTs new non-C. elegans nematode ESTs\n");
$log->write_to("$new_EMBL_CDS new C. elegans non-WormBase CDS entries\n");


$log->mail();


##############################
# a extremidade              #
##############################

exit(0);

##############################
# subroutines                #
##############################

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
    my $last_release_date = $wormbase->get_wormbase_release_date("short");
    
    print "Last release date was $last_release_date\n" if $debug;
    $log->write_to("\nLast release date was $last_release_date\n");

    my ($day, $month, $year) = split(/\//,$last_release_date);
    $month--;
 
    my $time = timelocal("0","0","0","$day","$month","$year");

    my $start_build_time = $time - ($buildtime  * 86400);
    my ($start_day, $start_month, $start_year) = (localtime($start_build_time))[3,4,5];
    
    $start_day   = sprintf("%02d", $start_day   % 100);
    $start_month = sprintf("%02d", $start_month % 100);

    $start_month++;
    $start_year += 1900;

    print     "Therefore last build started roughly $start_day/$start_month/$start_year\n" if $debug;
    $log->write_to("Therefore last build started roughly $start_day/$start_month/$start_year\n\n");
    
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
