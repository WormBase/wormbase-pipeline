# Wormbase.pm - module for general use by many Wormbase scripts
# adapted from babel.pl
# put together by krb, but mostly using stuff in babel.pl which
# was done by dl1 et al.

package Wormbase;

use lib $ENV{'CVS_DIR'};
use Carp;
use Ace;
use Log_files;
use File::Path;
use Storable;
sub new
  {
    my $class = shift;
    my %params = @_;
    my $self = {};
    bless $self , $class;

    # populate passed parameters ( will overwrite defaults if set )
    foreach ( keys %params ) {
      my $key = $_;
      $key =~ s/-//;
      $self->{$key} = $params{$_};
    }
    $self->establish_paths;
    return $self;
  }

#################################################################################


sub get_wormbase_version {
  my $self = shift;
  unless ( $self->{'version'} ) {
    my $dir = $self->autoace;
    if ( -e ("$dir/wspec/database.wrm") ) {
      my $WS_version = `grep "NAME WS" $dir/wspec/database.wrm`;
      chomp($WS_version);
      $WS_version =~ s/.*WS//;
      $self->version($WS_version);
    }
    else {
	$self->version(666);
    }
  }

  return ( $self->{'version'} );
}

###################################################################################

sub get_wormbase_version_name 
  {
    my $self = shift;
    my $version = $self->get_wormbase_version;
    return("WS$version");
  }


###################################################################################

sub version {
  my $self = shift;
  my $ver  = shift;
  $self->{'version'} = $ver if $ver;
  return $self->get_wormbase_version;
}

sub get_wormbase_release_date {

  my $self   = shift;
  my $format = shift;

  if (!(defined($format))) {
    $format = "long";
  } elsif ($format eq "short") {
    $format = "short";
  } elsif ($format eq "both") {
    $format = "both";
  } else {
    $format = "long";
  }

  
  my $line = `ls -l /nfs/disk69/ftp/pub/wormbase/development_release/md5sum.WS*`;
  my @split = split(/\s+/,$line);

  my $month = $split[5];
  my $month2;

  if ($month eq "Jan") {
    $month = "January";   $month2 = "01";
  } elsif ($month eq "Feb") {
    $month = "February";  $month2 = "02";
  } elsif ($month eq "Mar") {
    $month = "March";     $month2 = "03";
  } elsif ($month eq "Apr") {
    $month = "April";     $month2 = "04";
  } elsif ($month eq "May") {
    $month = "May";       $month2 = "05";
  } elsif ($month eq "Jun") {
    $month = "June";      $month2 = "06";
  } elsif ($month eq "Jul") {
    $month = "July";      $month2 = "07";
  } elsif ($month eq "Aug") {
    $month = "August";    $month2 = "08";
  } elsif ($month eq "Sep") {
    $month = "September"; $month2 = "09";
  } elsif ($month eq "Oct") {
    $month = "October";   $month2 = "10";
  } elsif ($month eq "Nov") {
    $month = "November";  $month2 = "11";
  } elsif ($month eq "Dec") {
    $month = "December";  $month2 = "12";
  }

  my $day = $split[6];

  my $day2;
  if (length($day) == 1) {
    $day2 = "0".$day;
  } else {
    $day2 = $day;
  }

  if ($day eq "1") {
    $day .= "st";
  } elsif ($day eq "2") {
    $day .= "nd";
  } elsif ($day eq "3") {
    $day .= "rd";
  } elsif ($day eq "21") {
    $day .= "st";
  } elsif ($day eq "22") {
    $day .= "nd";
  } elsif ($day eq "23") {
    $day .= "rd";
  } elsif ($day eq "31") {
    $day .= "st";
  } else {
    $day .= "th";
  }

  my $year = `date`;
  $year = substr($year,-3,2);

  # make a text style date
  my $date = $day." ".$month;

  # make a regular xx/xx/xx date
  my $date2 = $day2."/".$month2."/".$year;

  return($date)        if ($format eq "long");
  return($date2)       if ($format eq "short");
  return($date2,$date) if ($format eq "both"); 
}

###################################################################################

sub FetchData {
  my $self = shift;
  my ( $file, $ref, $dir ) = @_;

  # directory to load from can be passed in so that /acari can load files copied over
  unless ($dir) {
    $dir = $self->common_data;
  }
  print STDERR "using $dir for COMMON_DATA\n";
  open( FH, "<$dir/$file.dat" ) or die "can't open $dir/$file.dat\t:$!";
  undef $/;
  my $VAR1;
  my $data = <FH>;
  eval $data;
  die if $@;
  $/ = "\n";
  close FH;
  my $keycount = scalar keys %$VAR1;
  die "$file retrieval through FetchData failed - dat file is empty\n" if $keycount == 0;
  %$ref = (%$VAR1);
}

###################################################################################

sub get_script_version {
  my $self = shift;
  my $script = shift;
  my $script_dir = $ENV{'CVS_DIR'};
  open (GET_SCRIPT_LIST, "/bin/ls -l $script_dir/$script |");
  while (<GET_SCRIPT_LIST>) {
    chomp;
    my $stringlen = length ($_);
    $version = substr ($_,$stringlen-3,3);
    last;
  }
  close GET_SCRIPT_LIST;
  return ($version);
} 


#################################################################################

sub copy_check {
  my $self = shift;
  my ($file1,$file2) = @_;
  my $match = "";
  my $O_SIZE = (stat("$file1"))[7];
  my $N_SIZE = (stat("$file2"))[7];
    
  if ($O_SIZE != $N_SIZE) {
    $match = 0;
  } else {
    $match = 1;
  }
  return ($match);
} 


#################################################################################

sub mail_maintainer {
  my $self = shift;
  my ( $name, $maintainer, $logfile ) = @_;
  $maintainer = "ar2\@sanger.ac.uk, pad\@sanger.ac.uk, mt3\@sanger.ac.uk, gw3\@sanger.ac.uk, mh6\@sanger.ac.uk"
    if ( $maintainer =~ m/All/i );

  croak
    "trying email a log to a file - this will overwrite the existing file -STOPPING\nAre you passing a file name to Log object? \n"
      if ( -e $maintainer );
  open( OUTLOG, "|mailx -s \"$name\" $maintainer " );
  if ($logfile) {
    open( READLOG, "<$logfile" );
    while (<READLOG>) {
      print OUTLOG "$_";
    }
    close READLOG;
  } else {
    print OUTLOG "$name";
  }
  close OUTLOG or die "didn't close mail properly\n\n";
  return;
}


#################################################################################

sub celeaccession 
  {
    my $self = shift;
    local (*text_ace);
    my $seq = shift;
    local ($exec);
    local ($command);
    local ($accession);
    $ENV{'ACEDB'} = $self->autoace;
    $command = <<EOF;
    find sequence $seq
    show DB_info
    quit
EOF

    open( text_ace, "echo '$command' | $tace  | " );
    while (<text_ace>) {
      if (/\s+Database\s+EMBL\s+NDB_AC\s+(\S+)/) {
	$accession=$1;
      }
    }
    close text_ace;
    return $accession;
  }

#################################################################################

sub DNA_string_reverse {
  my $self = shift;
  my $revseq = reverse shift;
  $revseq =~ tr/a/x/;
  $revseq =~ tr/t/a/;
  $revseq =~ tr/x/t/;
  $revseq =~ tr/g/x/;
  $revseq =~ tr/c/g/;
  $revseq =~ tr/x/c/;
  return ($revseq);
}

#################################################################################

sub DNA_string_composition {
  my $self = shift;
  my $seq  = shift;
  $seq =~ tr/[A-Z]/[a-z]/;
  my $A = $seq =~ tr/a/a/;
  my $C = $seq =~ tr/c/c/;
  my $G = $seq =~ tr/g/g/;
  my $T = $seq =~ tr/t/t/;
  my $N = $seq =~ tr/n/n/;
  my $P = $seq =~ tr/-/-/;
  return ( $A, $C, $G, $T, $N, $P );
}

#################################################################################

sub gff_sort {
  my $self = shift;
  while (<>) {
    s/#.*//;
    next unless /\S/;
    @f = split /\t/;
    push @a, $_;
    push @n, $f[0];
    push @s, $f[3];
    push @e, $f[4];
  }

  foreach $i ( sort { $n[$a] cmp $n[$b] or $s[$a] <=> $s[$b] or $e[$a] <=> $e[$b] } 0 .. $#a ) {
    print $a[$i];
  }
}

#################################################################################

sub getseqEMBL {
  my $self = shift;
  my $acc = shift;
  my $absent;
    
  my $querycontent1 = "-e+[embl-acc:'$acc']";
  my $querycontent2 = "-e+[emblnew-acc:'$acc']";
    
  my $request1 = "/srs6bin/cgi-bin/wgetz?$querycontent1";
  my $request2 = "/srs6bin/cgi-bin/wgetz?$querycontent2";
    
  my $server = "srs.ebi.ac.uk";

  if ( !defined( open_TCP( *F, $server, 80 ) ) ) {
    print "Error connecting to server \n";
    exit(-1);
  }

  # get from emblnew
  print F "GET $request2 HTTP/1.0\n\n";
  print F "Accept: */*\n";
  print F "User-Agent: socketsrs/1.0\n\n";

  undef($absent);
  while ( my $return_line = <F> ) {
    if ( $return_line =~ /error/ ) {

      #	    print "Query returned no entry : '$request2'\n";
      $absent = 1;
      last;
    }
  }

  close F;

  # else get from embl
  if ($absent) {
    if ( !defined( open_TCP( *G, $server, 80 ) ) ) {
      print "Error connecting to server \n";
      exit(-1);
    }

    #	print "Query [$j]: '$request1'\n";

    print G "GET $request1 HTTP/1.0\n\n";
    print G "Accept: */*\n";
    print G "User-Agent: socketsrs/1.0\n\n";
  } else {
    if ( !defined( open_TCP( *G, $server, 80 ) ) ) {
      print "Error connecting to server \n";
      exit(-1);
    }

    #	print "Query [$j]: '$request2'\n";

    print G "GET $request2 HTTP/1.0\n\n";
    print G "Accept: */*\n";
    print G "User-Agent: socketsrs/1.0\n\n";
  }

  # Parsing annotation
    
  while (my $return_line=<G>) {
    #	print "$return_line";
	
    # Sequence version
    # SV   AL032679.2
    if ($return_line =~ /^SV\s+(\S+)/) {
      #	    print "parsing SV line : $return_line";
	    
      ($sv_acc,$sv_ver) = split (/\./, $1);
    }
    
    # Sequence
    # SQ   Sequence 2679 BP; 882 A; 510 C; 415 G; 872 T; 0 other;
    if ($return_line =~ /^SQ\s+Sequence\s(\d+)\sBP/) {
      #	    print "parsing SQ line : $return_line";
      $seq_len = $1;
      return ($sv_acc,$sv_ver,$seq_len);
      last;

    }
  }
  close G;
}

#################################################################################

#  Subroutines for generating release letter

##################################################################################
# Databases used in build                                                        #
##################################################################################

sub release_databases
  {
    my $self = shift;
    #GET THE DATABASES USED IN THIS BUILD

    #get the non_cambridge DB info from file
    open( DBS, $self->autoace."/logs/Primary_databases_used_in_build" );
    my ( $stlace, $citace, $brigdb, $cshace );
    my %dates;
    while (<DBS>) {
      chomp;
      my @info = split(/ : /,$_);
      $dates{$info[0]} = $info[1];
    }
    foreach my $key ( keys %dates) {
      my $oldstyle = $dates{$key};
      my $newstyle = "20".substr($oldstyle,0,2)."-".substr($oldstyle,2,2)."-".substr($oldstyle,4);
      $dates{$key} = $newstyle;
    }


    #get the Cambridge dates directly from block file 1
    my @date = $self->find_file_last_modified($self->database("camace")."/database/block1.wrm");
    $dates{camace} = $date[0];
    @date = $self->find_file_last_modified($self->database("geneace")."/database/block1.wrm");
    $dates{genace} = $date[0];

    #PARSE THE RELEASE LETTER FOR LAST BUILD INFO
    my $old_ver = $self->get_wormbase_version - 1;
    my $ver     = $old_ver + 1;

    my %old_dates;
    my $located = 0;
    my $count   = 0;		#determines how many lines to read = no databases
    open( OLD_LETTER, $self->reports."letter.WS$old_ver" );
    while (<OLD_LETTER>) {
      if ( ( $located == 1 ) && ( $count <= 6 ) ) {
	chomp;
	my @info = split( / : | - / , $_ );

	#this will put some crap in the hash from the 1st two line of the section but it no matter
	$old_dates{ $info[0] } = $info[1];
	$count++;
      } elsif ( $_ =~ m/Primary/ ) {
	$located = 1;
      }
    }
    my $dbaseFile = $self->reports."/dbases";
    open( WRITE, ">$dbaseFile" );

    print WRITE "Primary databases used in build WS$ver\n------------------------------------\n";
    foreach my $key ( sort keys %dates ) {
      print WRITE "$key : $dates{$key}";
      if ( "$dates{$key}" gt "$old_dates{$key}" ) {
	print WRITE " - updated\n";
      } elsif ( "$dates{$key}" lt "$old_dates{$key}" ) {
	print WRITE "you're using a older version of $key than for WS$old_ver ! ! \n";
      } else {
	print WRITE "\n";
      }
    }
    close WRITE;

    my $name       = "Database update report";
    my $maintainer = "All";
    $self->mail_maintainer( $name, $maintainer, $dbaseFile );

    return 1;
  }

##################################################################################
# Returns the date yyyy-mm-dd and time hh:mm:ss file was last modified           #
##################################################################################
sub find_file_last_modified {
  my $self     = shift;
  my $filename = shift;
  open( FILE, "<$filename" ) || die "cant open file $filename\n";
  my @fileinfo = stat FILE;
  my @date     = localtime( $fileinfo[9] );
  close(FILE);
  my $year = sprintf( "%d-%02d-%02d\n", $date[5] + 1900, $date[4] + 1, $date[3] );
  my $time = "$date[2]:$date[1]:$date[0]";
  chomp $year;

  my @last_modified = ( $year, $time );
  return @last_modified;
}

##################################################################################
# DNA Sequence composition                                                       #
##################################################################################

sub release_composition
  {
    my $self = shift;
    #get the old info from current_DB
    my $ver     = $self->get_wormbase_version;
    my $old_ver = $ver - 1;
    my %old_data;
    $old_data{"-"} = 0;		# initialise to avoid problems later if no gaps
    my $old_letter = $self->database("WS$old_ver")."/CHROMOSOMES/composition.all";
    open (OLD, "<$old_letter") || die "cant open data file - $old_letter";
    while (<OLD>) {
      chomp;
      if ( $_ =~ m/(\d+)\s+total$/ ) {

	#my $tot = $1;$tot =~ s/,//g; # get rid of commas
	$old_data{Total} = $1;
      } elsif ( $_ =~ m/^\s+([\w-]{1})\s+(\d+)/ ) {
	$old_data{"$1"} = $2;
      }
    }
    close(OLD);

    #now get the new stuff to compare
    my $new_letter = $self->chromosomes . "/composition.all";
    my %new_data;

    $new_data{"-"} = 0;		# initialise to avoid problems later if no gaps
    open( NEW, "<$new_letter" ) || die "cant open data file - $new_letter";
    while (<NEW>) {
      chomp;
      if ( $_ =~ m/(\d+)\s+total$/ ) {
	$new_data{Total} = $1;
      } elsif ( $_ =~ m/^\s+([\w-]{1})\s+(\d+)/ ) {
	$new_data{"$1"} = $2;
      }
    }
    close NEW;

    # Now check the differences
    my %change_data;
    my $compositionFile = $self->autoace . "/REPORTS/composition";
    open( COMP_ANALYSIS, ">$compositionFile" ) || die "cant open $compositionFile";
    print COMP_ANALYSIS "Genome sequence composition:\n----------------------------\n\n";
    print COMP_ANALYSIS "       \tWS$ver       \tWS$old_ver      \tchange\n";
    print COMP_ANALYSIS "----------------------------------------------\n";
    foreach my $key ( keys %old_data ) {
      $change_data{$key} = $new_data{$key} - $old_data{$key};
    }

    my @order = ( "a", "c", "g", "t", "n", "Total" ); # gaps removed
    foreach (@order) {
      if ( "$_" eq "Total" ) {
	print COMP_ANALYSIS "\n";
      }
      printf COMP_ANALYSIS ( "%-5s\t%-8d\t%-8d\t%+4d\n", $_, $new_data{$_}, $old_data{$_}, $change_data{$_} );
    }

    # Report file/email

    my $name       = "BUILD REPORT: Sequence composition";
    my $maintainer = "All";

    if ( $change_data{"-"} > 0 ) {
      print COMP_ANALYSIS "Number of gaps has increased - please investigate ! \n";
      $name = $name . " : Introduced a gap";

    }
    if ( $change_data{"Total"} < 0 ) {

      print COMP_ANALYSIS "Total number of bases has decreased - please investigate ! \n";
      $name = $name . " : Lost sequence";
    }
    if ( $change_data{"Total"} > 0 ) {
      print COMP_ANALYSIS "Total number of bases has increased - please investigate ! \n";
      $name = $name . " : Gained sequence";
    }
    close COMP_ANALYSIS;

    $self->mail_maintainer( $name, $maintainer, $compositionFile );

    return 1;
  }

##################################################################################
#  Wormpep                                                                       #
##################################################################################

sub release_wormpep		#($number_cds $number_total $number_alternate )
  {  
    my $self = shift;
    my ($number_cds, $number_total, $number_alternate) = @_;
    my $ver = $self->get_wormbase_version;
    my $old_ver = $ver -1;

    #extract data from new wormpep files
    my $wormpep = $self->wormpep;
    my $lost     = `more $wormpep/wormpep.diff$ver | grep 'lost' | wc -l`;
    my $new      = `more $wormpep/wormpep.diff$ver | grep 'new' | wc -l`;
    my $changed  = `more $wormpep/wormpep.diff$ver | grep 'changed' | wc -l`;
    my $appeared = `more $wormpep/wormpep.diff$ver | grep 'appear' | wc -l`;
    my $entries  = `cat $wormpep/wormpep.diff$ver | wc -l`;
    my $net      = $new + $appeared - $lost;
    my $codingDNA;

    #get no of coding bases from log file
    open( THIS_LOG, "$wormpep/wormpep_current.log" );
    while (<THIS_LOG>) {
      if ( $_ =~ /No\. of sequences \(letters\) written:\s+\d+\,\d+\s+\((.*)\)/ ) {
	$codingDNA = $1;
      }
    }

    #write new letter
    my $wormpepFile = $self->reports."/wormpep";
    open( LETTER, ">$wormpepFile" ) || die "cant open $wormpepFile\n";

    print LETTER "\n\nWormpep data set:\n----------------------------\n";
    print LETTER
      "\nThere are $number_cds CDS in autoace, $number_total when counting $number_alternate alternate splice forms.\n
The $number_total sequences contain $codingDNA base pairs in total.\n\n";

    print LETTER "Modified entries      $changed";
    print LETTER "Deleted entries       $lost";
    print LETTER "New entries           $new";
    print LETTER "Reappeared entries    $appeared\n";
    printf LETTER "Net change  %+d", $net;

    #get the number of CDS's in the previous build
    open( OLD_LOG, $self->basedir."/WORMPEP/wormpep$old_ver/wormpep_current.log" );
    my $oldCDS;
    while (<OLD_LOG>) {
      if ( $_ =~ /No\. of sequences \(letters\) written:\s+(\d+\,\d+)\s+\(.*\)/ ) {
	$oldCDS = $1;
	$oldCDS =~ s/,//g;
      }
    }
    close OLD_LOG;

    #check
    my $mail;
    if ( $lost + $new + $changed + $appeared != $entries ) {
      print LETTER "cat of wormpep.diff$ver does not add up to the changes (from $0)";
    }
    if ( $oldCDS + $net != $number_total ) {
      print LETTER
	"The differnce between the total CDS's of this ($number_total) and the last build ($oldCDS) does not equal the net change $net\nPlease investigate! ! \n";
    }

    close LETTER;

    my $name       = "Wormpep release stats";
    my $maintainer = "All";
    $self->mail_maintainer( $name, $maintainer, $wormpepFile );

    return 1;
  }

#end of release letter generating subs
#############################################

sub test_user_wormpub {
  my $self = shift;
  my $name = `whoami`;
  chomp $name;
  if ( "$name" eq "wormpub" ) {
    print "running scripts as user wormpub . . . \n\n";
    return;
  } else {
    print
      "You are doing this as $name NOT wormpub ! \n\n If you are going to alter autoace in any way it will break.\nDo you want to continue? (y/n). . ";
    my $response = <STDIN>;
    chomp $response;
    if ( "$response" eq "n" ) {
      die "probably for the best !\n";
    } else {
      print
	"OK - on your head be it !\nBack to the script . .\n#########################################################\n\n\n";
      return;
    }
  }
}

#############################################
sub runtime {
  my $self    = shift;
  my $runtime = `date +%H:%M:%S`;
  chomp $runtime;
  return $runtime;
}
###############################################
sub rundate {
  my $self    = shift;
  my $rundate = `date +%y%m%d`;
  chomp $rundate;
  return $rundate;
}

###################################################
# subs to get the correct version of ACEDB binaries
sub tace {
  my $self = shift;
  return $self->{'tace'};
}

sub giface {
  my $self = shift;
  return $self->{'giface'};
}

####################################
# Check for database write access
####################################

sub check_write_access {

    my $self         = shift;
    my $database     = shift;
    my $write_access = "yes";

    $write_access = "no" if ( -e "${database}/database/lock.wrm" );
  return ($write_access);

}

####################################
# Delete files from directory
####################################
sub delete_files_from {
  my $self = shift;
  my ( $directory, $pattern, $folder ) = @_;
  my $file;
  my $delete_count = 0;
  my $fail_warn    = 1;

  return undef unless ( -e $directory );

  if ( $folder eq "+" ) {
    print "Removing entire dir and subdirs of $directory\n";
    $delete_count = rmtree($directory);
  } else {
    opendir( TO_GO, $directory ) or die "cant get listing of $directory:\t$!\n";

    #   $pattern = "." if $pattern eq "*";

    $pattern = "." unless $pattern;
    $pattern =~ s/\*/\./g;

    while ( $file = readdir(TO_GO) ) {
      next if ( $file eq "." or $file eq ".." );
      if ( $file =~ /$pattern/ ) {
	if ( unlink("$directory/$file") ) {
	  $delete_count++;
	} else {
	  warn "couldn't unlink $directory/$file :\t$!\n";
	  undef $fail_warn;
	}
      }
    }
  }
  return $fail_warn ? $delete_count : $fail_warn; # undef if failed else no. files removed;
}

sub load_to_database {

  my $self     = shift;
  my $database = shift;
  my $file     = shift;

  # tsuser is optional but if set, should replace any dots with underscores just in case
  # if not set im using the filename with dots replaced by '_'
  my $tsuser = shift;

  unless ($tsuser) {

    # remove trailing path of filename
    $tsuser = $file;
    $tsuser =~ s/.*\///;
  }

  $tsuser =~ s/\./_/g;

  unless ( -e "$file" ) {
    die " Couldn't find file named: $file\n";
  }
  my $command = "pparse $file\nsave\nquit\n";

  open( WRITEDB, "| $tace -tsuser $tsuser $database " ) || die "Couldn't open pipe to database\n";
  print WRITEDB $command;
  close(WRITEDB);
}

sub wormpep_files {
  my $self = shift;
  return ( "wormpep", "wormpep.accession", "wormpep.dna", "wormpep.history", "wp.fasta", "wormpep.table",
	   "wormpep.diff" );
}

sub test        { $self = shift; return $self->{'test'}; }
sub debug       { $self = shift; return $self->{'debug'}; }
sub wormpub     { $self = shift; return $self->{'wormpub'}; }
sub basedir     { $self = shift; return $self->{'basedir'}; }
sub autoace     { $self = shift; return $self->{'autoace'}; }
sub wormpep     { $self = shift; return $self->{'wormpep'}; }
sub wormrna     { $self = shift; return $self->{'wormrna'}; }
sub gff         { $self = shift; return $self->{'gff'}; }
sub gff_splits  { $self = shift; return $self->{'gff_splits'}; }
sub chromosomes { $self = shift; return $self->{'chromosomes'}; }
sub logs        { $self = shift; return $self->{'logs'}; }
sub ftp_upload  { $self = shift; return $self->{'ftp_upload'}; }
sub reports     { $self = shift; return $self->{'reports'}; }
sub misc_static { $self = shift; return $self->{'misc_static'}; }
sub misc_dynamic { $self = shift; return $self->{'misc_dynamic'}; }
sub primaries   { $self = shift; return $self->{'primaries'}; }

# this can be modified by calling script
sub common_data {
  $self = shift;
  my $path = shift;
  if ($path) {
    if ( -e $path ) {
      $self->{'common_data'} = $path;
    } else {
      die "$path does not exist\n";
    }
  }
  return $self->{'common_data'};
}

sub database {
  my $self     = shift;
  my $database = shift;
  if ( $self->{'databases'}->{"$database"} ) {
    return $self->{'databases'}->{"$database"};
  } else {

    # try under the usual database path
    my $poss_path = $self->wormpub . "/DATABASES/$database";
    return $poss_path if ( -e $poss_path );

    #build related database
    $poss_path = $self->basedir . "/$database";
    return $poss_path if ( -e $poss_path );
    print STDERR "no such database $database\n";
    return undef;
  }
}

sub primary {
  my $self = shift;
  my $database = shift;
  my $path  = $self->primaries . "/$database";
  print STDERR "no such primary database\n" unless (-e $path);
  return $path;
}

sub establish_paths {
  my $self = shift;
  ( $self->{'wormpub'} ) = glob("~wormpub");
  my $basedir = $self->wormpub . "/BUILD";
  $basedir = $self->wormpub . "/TEST_BUILD" if $self->test;

  $self->{'basedir'}    = $basedir;
  $self->{'autoace'}    = "$basedir/autoace";
  $self->{'ftp_upload'} = "/nfs/ftp_uploads/wormbase";
  $self->{'wormpep'}    = $basedir . "/WORMPEP/wormpep" . $self->get_wormbase_version;
  $self->{'wormrna'}    = $basedir . "/WORMRNA/wormrna" . $self->get_wormbase_version;

  $self->{'logs'}        = $self->autoace . "/logs";
  $self->{'common_data'} = $self->autoace . "/COMMON_DATA";
  $self->{'chromosomes'} = $self->autoace . "/CHROMOSOMES";
  $self->{'reports'}     = $self->autoace . "/REPORTS";
  $self->{'gff'}         = $self->chromosomes; #to maintain backwards compatibility 
  $self->{'gff_splits'}  = $self->autoace . "/GFF_SPLITS";
  $self->{'primaries'}   = $self->basedir . "/PRIMARIES";

  $self->{'tace'}   = glob("~wormpub/ACEDB/bin_ALPHA/tace");
  $self->{'giface'} = glob("~wormpub/ACEDB/bin_ALPHA/giface");

  $self->{'databases'}->{'geneace'} = $self->wormpub . "/DATABASES/geneace";
  $self->{'databases'}->{'camace'}  = $self->wormpub . "/DATABASES/camace";
  $self->{'databases'}->{'current'} = $self->wormpub . "/DATABASES/current_DB";
  $self->{'databases'}->{'autoace'} = $self->autoace;

  $self->{'primary'}->{'camace'}  = $self->primaries ."/camace";
  $self->{'primary'}->{'geneace'} = $self->primaries ."/geneace";
  $self->{'primary'}->{'stlace'}  = $self->primaries ."/stlace";
  $self->{'primary'}->{'citace'}  = $self->primaries ."/citace";
  $self->{'primary'}->{'cshace'}  = $self->primaries ."/cshace";
  $self->{'primary'}->{'brigace'} = $self->primaries ."/brigace";

  $self->{'build_data'} = $self->{'wormpub'} . "/BUILD_DATA";
  $self->{'misc_static'} = $self->{'build_data'} . "/MISC_STATIC";
  $self->{'misc_dynamic'} = $self->{'build_data'} . "/MISC_DYNAMIC";

  # create dirs if missing
  mkpath( $self->logs )        unless ( -e $self->logs );
  mkpath( $self->common_data ) unless ( -e $self->common_data );
  mkpath( $self->wormpep )     unless ( -e $self->wormpep );  system("chmod -R g+w ".$self->wormpep);
  mkpath( $self->wormrna )     unless ( -e $self->wormrna );  system("chmod -R g+w ".$self->wormrna);
  mkpath( $self->chromosomes ) unless ( -e $self->chromosomes );
  mkpath( $self->reports )     unless ( -e $self->reports );
  mkpath( $self->gff )         unless ( -e $self->gff );
  mkpath( $self->gff_splits )  unless ( -e $self->gff_splits );
  mkpath( $self->primaries )   unless ( -e $self->primaries );

  system("chmod -R g+w ".$self->autoace);

}

sub run_script {
  my $self   = shift;
  my $script = shift;
  my $log    = shift;

  my $store = $self->autoace . "/wormbase.store";
  store( $self, $store );
  my $command = "perl $ENV{'CVS_DIR'}/$script -store $store";
  print "$command\n" if $self->test;
  return $self->run_command( "$command", $log );
}


sub run_command {
  my $self    = shift;
  my $command = shift;
  my $log     = shift;
  $log->write_to("running $command\n");
  my $return_status = system("$command");
  if ( ( $return_status >> 8 ) != 0 ) {
    $log->write_to(" WARNING: $script returned non-zero ($return_status)\n") if $log;
    return 1;
  } else {
    $log->write_to("$script exited cleanly\n") if $log;
    return 0;
  }
}

sub wait_for_LSF {
  my $self = shift;
  sleep 10;
  while ( &jobs_left != 0 ) {
    sleep 10;
  }

  print "all jobs finished\n";
  return;

  sub jobs_left {
    my $self  = shift;
    my $count = 0;
    open( JOBS, "bjobs |" );
    while (<JOBS>) {
      print $_;
      next if /JOBID/;		#title line
      $count++;
    }
    close JOBS;
    print "$count jobs left , , \n";
    return $count;
  }
}

################################################################################
#Return a true value
################################################################################

1;

__END__

=pod

=head1 NAME - Wormbase.pm

=head2 DESCRIPTION

The Wormbase.pm module replaces babel.pl which was previouly used
to access some common subroutines for general Wormbase development
work.  

This module provides access to the following subroutines:

=over 4

=item *

get_wormbase_version

This subroutine returns the current WormBase release version.  This is read from
the file: /wormsrv2/autoace/wspec/database.wrm file.  The function returns
the number.

=back

=over 4

=item *

get_wormbase_version_name

As above, but returns the full name, e.g. 'WS47' rather than just '47'

=back

=over 4

=item *

get_wormbase_release_date

Gets the date of the release from the date stamp of the letter.WSxx file in the wormbase ftp 
directory.  Creation of the letter.WSxx file occurs pretty much at the end of the rebuild
so it is really an approximate date.

If no argument is passed to the function it will return the date in 'long' format. E.g.
"21st September".  It will also return this format if the string 'long' is passed to the 
function.

If the string 'short' is passed to the function it will return a six figure date format,
e.g. dd/mm/yy.

If the string 'both' is passed to the function it will return the long and the short versions.

=back


=over 4

=item *

get_wormpep_version

Takes the wormbase version number and adds 10 to it.

=over 4

=item *

get_script_version

This subroutine grabs the version number of the file.  No longer used
and is not exported by default from the module.  Replaced by the
get_cvs_version subroutine.

=back

=over 4

=item *

copy_check

Pass the names of two files to this subroutine and it will return '1' if they
are the same size or '0' if otherwise.

=back

=over 4

=item *

mail_maintainer

Mails the logfile from certain script to desired recipients.

Usage:                                                                    
&mail_maintainer(<title>,<maintainer e-mail list>,<logfile>);                                                                                
No return value.

=back

=over 4

=item *
celeaccession

Pass this subroutine the name of a clone and it will return the 
corresponding accession number

=item *
find_database

Pass a list of clones and 'cam' or 'stl' and it will return the clones 
in camace / stlace

=back

=over 4

=item *

release_databases

writes a file detailing the databases used in the current build. 
The contents of this are included in the release letter.
Checks that databases used are not older than those used in previous build.
Emails output and warns if problems

=back

=over 4

=item *

find_file_last_modified

Passed a filename, this returns the date yyyy-mm-dd and time

=back

=over 4

=item *

release_composition

writes a file detailing the sequence composition in the current build.
This is used in the release letter.
Checks the sequence composition of the current genome compared to the last one.
Does various checks on data integrity and flags any problems in the mail sent.

=back

=over 4

=item *

release_wormpep

Compiles release stats of the current Wormpep and writes to a file, later used in release letter.

=back 

=over 4

=item *

check_write_access

Takes a path to an acedb database and returns "yes" if no database/lock.wrm file is present
(i.e. yes, you have write access) and returns "no" if such a file is present.

=back 

=over 4

<<<<<<< Wormbase.pm
=item *

gff_sort


gff_sort reads a GFF file from STDIN , ignores comment lines and blank
lines and prints the remaining lines to STDOUT sorted by the following
keys:

name (column 0)
start (column 3)
end (column 4)

=back 

=over 4

=item *
=======
=item *

gff_sort


gff_sort reads a GFF file from STDIN , ignores comment lines and blank
lines and prints the remaining lines to STDOUT sorted by the following
keys:

name (column 0)
start (column 3)
end (column 4)

=back 

=over 4
>>>>>>> 1.88.4.7

<<<<<<< Wormbase.pm
tace

tace returns the path for tace that is being used
=======
=item *

tace

tace returns the path for tace that is being used
>>>>>>> 1.88.4.7

=back 

=over 4

=item *

dbfetch

dbfetch takes arguments:

   name of sequence to find
   name of file containing one or more sequences in Fasta format

It returns the first fasta format sequence from the file whose name
matches the input seqeunce name.

=back 

=over 4


=head2 load_to_datase

=head2 SYNOPSIS

=over4

&load_to_database($database, $file, "tsuser" );

=back

=head2 USAGE

=over4

Loads specified file in to specified acedb database with tsuser as specified :)
If tsuser not set then the file name will be used (no path, and '_' replacing and '.'s )

=back

=cut

