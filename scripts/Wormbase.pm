#Wormbase.pm - module for general use by many Wormbase scripts
# adapted from babel.pl
# put together by krb, but mostly using stuff in babel.pl which
# was done by dl1 et al.
# SPECIES BRANCH

package Wormbase;
use strict;
# use lib $ENV{'CVS_DIR'};

$ENV{ACEDB_NO_BANNER}=1;

use Carp;
use Ace;
use Log_files;
use File::Path;
use File::stat;
use Storable;
use Digest::MD5 qw(md5_hex);
use Species;

our @core_organisms=qw(Elegans Briggsae Remanei Brenneri Japonica Pristionchus Brugia Ovolvulus Sratti Tmuris);
#our @tier3_organisms=qw(Mhapla Mincognita Heterorhabditis Hcontortus Hcontortus_gasser Cangaria Tspiralis Ctropicalis Asuum Bxylophilus Csinica Loaloa Asuum_davis Panagrellus Dimmitis Namericanus Acey Tsuis_male Tsuis_female Pexspectatus);
our @tier3_organisms=qw(Cangaria Ctropicalis Csinica Panagrellus Elegans_hawaii Elegans_vc2010 Remanei_px356 Cnigoni Clatens Cinopinata Otipulae Remanei_px439);

our @provisional_organisms = qw();

our @allowed_organisms=(@core_organisms, @tier3_organisms,@provisional_organisms); #class data

sub initialize {
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
  $self->{'version'} = 666 if( $self->test);
  print STDERR "Using ".$self->{'autoace'}."\n" if( $self->{'autoace'} );
  #$self->establish_paths;
  return $self;
}

#################################################################################

# DEMO constructor
#
# should become: initialize(-organism => Elegans,......) 
sub new {
  my $class = shift;
  my %params=@_;

  my $self;
  my $ORGANISM=($params{'-organism'}||'Elegans'); # fixed at the moment
  $ORGANISM = "\u$ORGANISM";	#make sure 1st letter caps
  die "invalid organism: $ORGANISM" if ! grep {$_ eq $ORGANISM} @allowed_organisms;

  $params{'-species'} = lc $ORGANISM;
  $self=$ORGANISM->_new(\%params);
  $self->establish_paths;
  return $self;
}

#######################################################################

sub get_wormbase_version {
  my $self = shift;
  unless ( $self->{'version'} ) {
    my $dir = $self->autoace;
    if ( -e ("$dir/wspec/database.wrm") ) {
      my $WS_version = `grep "NAME WS" $dir/wspec/database.wrm`;
      chomp($WS_version);
      $WS_version =~ s/.*WS//;
      $self->version($WS_version);
    } else {
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

sub get_dev_version {
  my $self = shift;
  unless ( $self->{'version'} ) {
    my $dir = $self->database('current');
    if ( -e ("$dir/wspec/database.wrm") ) {
      my $WS_version = `grep "NAME WS" $dir/wspec/database.wrm`;
      chomp($WS_version);
      $WS_version =~ s/.*WS//;
      $self->version($WS_version);
    } else {
      $self->version(666);
    }
  }

  return ( $self->{'version'} );
}
###################################################################################

sub version { 	 
    my $self = shift; 	 
    my $ver  = shift; 	 
    $self->{'version'} = $ver if $ver; 	 
    return $self->{'version'}; 	 
}

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
  warn "$file retrieval through FetchData failed - dat file is empty\n" if $keycount == 0;
  %$ref = (%$VAR1);
}

###################################################################################

sub get_script_version {
  my $self = shift;
  my $script = shift;
  my $script_dir = $ENV{'CVS_DIR'};
  my $version;
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
  my $O_SIZE = (-s $file1);
  my $N_SIZE = (-s $file2);
  
  return 0 unless ($O_SIZE && $N_SIZE);

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
  $maintainer = "paul.davis\@wormbase.org, magdalena.zarowiecki\@wormbase.org, stavros.diamantakis\@wormbase.org, mark.quintontulloch\@wormbase.org"    if ( $maintainer =~ m/All/i );
  croak "trying to email a log to a file - this will overwrite the existing file -STOPPING\nAre you passing a file name to Log object? \n"  if ( -e $maintainer );
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
  my(@a, @n, @s, @e, @f);
  while (<>) {
    s/#.*//;
    next unless /\S/;
    @f = split /\t/;
    push @a, $_;
    push @n, $f[0];
    push @s, $f[3];
    push @e, $f[4];
  }

  foreach my $i ( sort { $n[$a] cmp $n[$b] or $s[$a] <=> $s[$b] or $e[$a] <=> $e[$b] } 0 .. $#a ) {
    print $a[$i];
  }
}

#################################################################################

#################################################################################

#  Subroutines for generating release letter


##################################################################################
# Returns the date yyyy-mm-dd and time hh:mm:ss file was last modified           #
##################################################################################
sub find_file_last_modified {
  my $self     = shift;
  my $filename = shift;

  my $fileinfo = stat($filename);
  my @date     = localtime( $fileinfo->mtime );

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
    my $log = shift;
    #get the old info from current_DB
    my $ver     = $self->get_wormbase_version;
    my $old_ver = $ver - 1;
    my %old_data;
    $old_data{"-"} = 0;		# initialise to avoid problems later if no gaps
    my $old_letter = $self->database("WS$old_ver")."/CHROMOSOMES/composition.all";
    open (OLD, "<$old_letter") or $log->log_and_die("cant open data file - $old_letter");
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

    my @order = ( "a", "c", "g", "t", "n", "-", "Total" ); # gaps removed
    foreach (@order) {
      if ( "$_" eq "Total" ) {
	print COMP_ANALYSIS "\n";
      }
      if (! exists $old_data{$_}) {print "no data for $_\n"; next}
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

sub release_wormpep {		
  #($number_cds $number_total $number_alternate )
  my $self = shift;
  my ($number_cds, $number_total, $number_alternate, $log) = @_;
  my $ver = $self->get_wormbase_version;
  my $old_ver = $ver -1;
  
  #extract data from new wormpep files
  my $wormpep = $self->wormpep;
  my $lost     = `grep -c 'lost' $wormpep/wormpep.diff$ver`;
  my $new      = `grep -c 'new' $wormpep/wormpep.diff$ver`;
  my $changed  = `grep -c 'changed' $wormpep/wormpep.diff$ver`;
  my $appeared = `grep -c 'appear' $wormpep/wormpep.diff$ver`;

  my $entries  = `cat $wormpep/wormpep.diff$ver | wc -l`;
  my $net      = $new + $appeared - $lost;
  my $codingDNA;
  
  #get no of coding bases 
  open( DNA, "$wormpep/wormpep.dna$ver" );
  while (<DNA>) {
    /^\>/ and next;
    /^(\S+)/ and $codingDNA += length($1);
  }

  
  #write new letter
  my $wormpepFile = $self->reports."/wormpep";
  open( LETTER, ">$wormpepFile" ) || die "cant open $wormpepFile\n";
  
  print LETTER "\n\nWormpep data set:\n----------------------------\n";
  print LETTER "\nThere are $number_total CDSs, from $number_cds protein-coding loci\n";
  print LETTER "\nThe $number_total sequences contain $codingDNA base pairs in total.\n\n";
  print LETTER "Modified entries      $changed";
  print LETTER "Deleted entries       $lost";
  print LETTER "New entries           $new";
  print LETTER "Reappeared entries    $appeared\n";
  printf LETTER "Net change  %+d", $net;
  
  #get the number of CDS's in the previous build
  my $oldCDS = 0;
  open( OLD_PEP, $self->basedir."/WORMPEP/wormpep$old_ver/wormpep${old_ver}.pep" );
  while (<OLD_PEP>) {
    if (/^\>/) {
      $oldCDS++;
    }
  }
  close OLD_PEP;
  
  #check
  my $mail;
  if ( $lost + $new + $changed + $appeared != $entries ) { # this is not very informative - it says that the sum of types of lines in the diff file equals the total number of lines in the diff file - it doesn't say that we have a correct set of history data.
    print LETTER "cat of wormpep.diff$ver does not add up to the changes (from $0)";
  }
  if ( $oldCDS + $net != $number_total ) { # this is comparing the previous release number of proteins + net change of proteins (new + reappeared - lost) to the total number of proteins in this release. This is screwed by the double counting of reappeared proteins in the diff file because these proteins are in the diff file also as lost or changed.
    my $thediff = $number_total - $oldCDS;
    $log->write_to(
	"\nThe difference ($thediff) between the total CDS's of this ($number_total) and the last build ($oldCDS) does not equal the net change $net\nPlease investigate! ! \n");
    $thediff = $number_total + $appeared - $oldCDS;
    $log->write_to("However there is double counting of the $appeared reappeared proteins. Taking this into account the difference is ($thediff).\n");
    if (!$thediff) {
      $log->write_to("So that's probably OK.");
    }
    $log->write_to("We have noticed that there are also some instances of corruption of the history files.\n");
  }
  close LETTER;
  
  my $name       = "Wormpep release stats";
  my $maintainer = $self->debug ? $self->debug : "All";
  $self->mail_maintainer( $name, $maintainer, $wormpepFile );
  
  return 1;
}

#end of release letter generating subs
#############################################

sub test_user_wormpub {
  my $self = shift;
  my $name = `whoami`;
  chomp $name;
  if ( $name eq 'wormpub' or $name eq 'wormbase' or $name eq 'wormpipe' ) {
    print "running scripts as user $name . . . \n\n";
    return 1;
  } elsif ($self->test) {
    print "running in Test\n";
    return 1;
  } else {
    print "You are doing this as $name NOT wormpub ! \n\n If you are going to alter autoace in any way it will break.\nYou have 30 seconds in which to abort before this script will continue!\n";
    sleep 30;
    return 0;
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

sub giface_server {
  my $self = shift;
  return $self->{'giface_server'};
}

sub giface_client {
  my $self = shift;
  return $self->{'giface_client'};
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

####################################
# do various checks on the files output from a script
# (or part of a script) - the tests are read from the config file
# ~wormpub/BUILD/autoace_config/autoace.config
####################################

sub check_files {
  my ($self, $log, $part) = @_;

  # read the filenames and criteria from the config file
  my $config = $self->basedir . "/autoace_config/check_file.config";
  open(FCCONFIG, "<$config") || $log->log_and_die("Can't open $config\n");

  my $species = $self->species;
  my %criteria;
  my $found_species = '';
  my $found_script = 0;
  my $found_a_file = 0;
  my $file;

  my $script_name = $0;
  $script_name =~ s/.+\/(\S+)/$1/; # remove header of path

  while (my $line = <FCCONFIG>) {
    chomp $line;
    if ($line =~ /^#/ || $line =~ /^\s+$/) {  
      next;
    } elsif ($line =~ /^SPECIES/) {
      $found_species = '';
      $found_script = 0;
      if ($line =~ /^SPECIES\s+$species/) {
	$found_species = 'species';
      }
      if ($line =~ /^SPECIES\s+default/) {
	$found_species = 'default';
      }
    } elsif ($found_script && $line !~ /^SCRIPT/) {
      if ($line =~ /^\s*FILE/) {
	($file) = $line =~ /^\s*FILE\s+(\S+)/;
	$file =~ s/wormbase->/self->/g; # convert filenames with '$wormbase->' in to '$self->'
	if ($file =~ /(\$\S+\-\>[\w_\(\)\']+)(\/\S+)/) {
	  $file = eval($1) . $2;         # and expand to the full path
	}
	while ($file =~ /(\S*)\$\{\\(\$\w+\-\>[\w_\(\)\']+)\}(\S*)/) { # look for: anything ${\$object->method} anything
	  $file = $1 . eval($2) . $3;
	}
	# don't want to do a default test it we already have the tests
	# that are specific for this species
	if ($found_species eq 'default' &&
	    exists $criteria{$file}) {next;}
	$found_a_file = 1;
	$criteria{$file}{exists} = 1;
      } elsif ($line =~ /^\s+(\S+)\s+(\S+)/) {
	my $key = $1;
	my $value = $2;
	# get and store the tests for this file
	if ($key eq 'lines' || $key eq 'requires') {
	  push @{$criteria{$file}{$key}}, $value;
	} else {
	  $criteria{$file}{$key} = $value;
	}
      }
    } elsif ($found_species) {
      if ($line =~ /^SCRIPT/) {
	$found_script = 0;
	if ($line =~ /^SCRIPT\s+$script_name\s*$/ || 
	    (defined $part && $line =~ /^SCRIPT\s+$script_name\s+$part$/)) {
	  $found_script = 1;
	}
      }
    }
  }
  close(FCCONFIG);

  # complain if we didn't find an entry for this script
  if (!$found_a_file) {
    if ( $log) {
      $log->write_to("WARNING: Couldn't find any files to test in $config\n");
    }
    carp "WARNING: Couldn''t find any files to test in $config\n";
    return 1;
  }

  my $errors = 0;
  foreach my $file (keys %criteria) {
    $errors++ if ($self->check_file($file, $log, %{$criteria{$file}}));
  }
  
  return $errors;		# number of files failing tests
}

####################################

sub check_file {

  my ($self, $file, $log, %criteria) = @_;
  my @problems;


  # file must not exist
  if (exists $criteria{must_not_exist}) {
    if (-e $file) {
      push @problems,  "the file exists - it should not exist";
    } else {
      delete $criteria{must_not_exist};
    }

  } else {

    # file must exist (the default)
    if (! -e $file) {
      if ( $log) {
	$log->error;
	$log->write_to("ERROR: Couldn't find file named: $file\n");
      }
      carp "ERROR: Couldn't find file named: $file\n";
      return 1;
    }
    delete $criteria{exists}; # this 'exists' criterion is redundant

    if (!-r $file) {
      push @problems,  "file is not readable";
    }

    if (!exists $criteria{readonly}) {
      if (!-w $file) {
	push @problems,  "file is not writeable";
      }
    } else {
      delete $criteria{readonly};
    }

  }

  my $size;
  my $second_file_size;
  if (exists $criteria{samesize}) {
    $size = (-s $file) unless $size;
    $second_file_size = (-s $criteria{samesize});
    if ($second_file_size != $size) {
      push @problems,  "file size ($size) not equal to that of file '$criteria{samesize}' ($second_file_size)";
    }
    delete $criteria{samesize};
  }
  if (exists $criteria{similarsize}) {
    $size = (-s $file);
    $second_file_size = (-s $criteria{similarsize});
    if ($second_file_size < $size * 0.9 || $second_file_size > $size * 1.1) {
      push @problems,  "file size ($size) not similar to that of file '$criteria{similarsize}' ($second_file_size)";
    }
    delete $criteria{similarsize};
  }
  if (exists $criteria{minsize}) {
    $size = (-s $file) unless $size;
    if ($size < $criteria{minsize}) {
      push @problems,  "file size ($size) less than required minimum ($criteria{minsize})";
    }
    delete $criteria{minsize};
  }
  if (exists $criteria{maxsize}) {
    $size = (-s $file) unless $size;
    if ($size > $criteria{maxsize}) {
      push @problems, "file size ($size) greater than required maximum ($criteria{maxsize})";
    }
    delete $criteria{maxsize};
  }
  my $lines;
  my $second_file_lines;
  if (exists $criteria{samelines}) {
    ($lines) = (`wc -l $file` =~ /(\d+)/);
    ($second_file_lines) = (`wc -l $criteria{samelines}` =~ /(\d+)/);
    if ($second_file_lines != $lines) {
      push @problems,  "number of lines ($lines) not equal to that of file '$criteria{samelines}' ($second_file_lines)";
    }
    delete $criteria{samelines};
  }
  if (exists $criteria{similarlines}) {
    ($lines) = (`wc -l $file` =~ /(\d+)/) unless $lines;
    ($second_file_lines) = (`wc -l $criteria{similarlines}` =~ /(\d+)/) unless $lines;
    if ($second_file_lines < $lines * 0.9 || $second_file_lines > $lines * 1.1) {
      push @problems,  "number of lines ($lines) not similar to that of file '$criteria{similarlines}' ($second_file_lines)";
    }
    delete $criteria{similarlines};
  }
  if (exists $criteria{minlines}) {
    ($lines) = (`wc -l $file` =~ /(\d+)/) unless $lines;
    if ($lines < $criteria{minlines}) {
      push @problems, "number of lines ($lines) less than required minimum ($criteria{minlines})";
    }
    delete $criteria{minlines};
  }
  if (exists $criteria{maxlines}) {
    ($lines) = (`wc -l $file` =~ /(\d+)/) unless $lines;
    if ($lines > $criteria{maxlines}) {
      push @problems, "number of lines ($lines) greater than required maximum ($criteria{maxlines})";
    }
    delete $criteria{maxlines};
  }
  
  if (exists $criteria{requires} || exists $criteria{line1} || exists $criteria{line2} || exists $criteria{lines}) {
    open (CHECK_FILE, "< $file") || die "Can't open $file\n";
    my $line_count = 0;
    while ($line_count++, my $line = <CHECK_FILE>) {

      if (exists $criteria{requires}) {
	my $re_count = 0;
	foreach my $regex (@{$criteria{requires}}) { # we need to find at least one of each of these regexps in the file
	  if ($line =~ /$regex/) {
	    splice @{$criteria{requires}}, $re_count, 1; # remove the successful regexp frmo the array
	  }
	  $re_count++;
	}
	if (@{$criteria{requires}} == 0 && !(exists $criteria{line1} || exists $criteria{line2} || exists $criteria{lines})) {last;}
      }

      if ($line_count == 1 && exists $criteria{line1}) {
	if ($line !~ /$criteria{line1}/) {
	  push @problems, "line $line_count:\n$line\ndoesn't match criterion 'line1 => /$criteria{line1}/'";
	  last;
	}
	next;			# don't do 'lines' check on line 1 if 'line1' check exists
      }
      if ($line_count == 2 && exists $criteria{line2}) {
	if ($line !~ /$criteria{line2}/) {
	  push @problems, "line $line_count:\n$line\ndoesn't match criterion 'line2 => /$criteria{line2}/'";
	  last;
	}
	next;			# don't do 'lines' check on line 2 if 'line2' check exists
      }
      if ($line_count > 2 && !exists $criteria{lines} && !@{$criteria{requires}}) {last;}
      if (exists $criteria{lines}) {
	my $line_ok = 0;
	foreach my $regex (@{$criteria{lines}}) { # each line in the file must match one of these regexps
	  if ($line =~ /$regex/) {
	    $line_ok = 1;
	    last;
	  }
	}
	if (!$line_ok) {
	  my $regexp = '/'.join('/ /',@{$criteria{lines}}).'/'; # for easier readability
	  push @problems, "line $line_count:\n$line\ndoesn't match criterion 'lines => [$regexp]]'";
	  last;
	}
      }
    }
    close (CHECK_FILE);

    # check if there are any 'requires' regexps which didn't match
    if (exists $criteria{requires} && @{$criteria{requires}}) {
      push @problems, "the criterion 'requires => [@{$criteria{requires}}]' did not match any line";
    }

    delete $criteria{requires};
    delete $criteria{line1};
    delete $criteria{line2};
    delete $criteria{lines};
  }
  if (exists $criteria{gff}) {
    my ($sequence_name, $sequence_start, $sequence_end);
    my $MAX_FEATURE_LENGTH = 100000;
    open (CHECK_FILE, "< $file") || die "Can't open $file\n";
    my $line_count = 0;
    while ($line_count++, my $line = <CHECK_FILE>) {
      ##sequence-region CHROMOSOME_X 1 17718851
      if ($line =~ /^\#\#sequence-region\s+(\S+)\s+(\d+)\s+(\d+)/) {
	($sequence_name, $sequence_start, $sequence_end) = ($1, $2, $3);
	next;
      }
      #CHROMOSOME_X    Link    region  1       17718851        .       +       .       Sequence "CHROMOSOME_X"
      if (defined $sequence_name && $line =~ /^${sequence_name}\s+Link\s+region\s+(\d+)\s+(\d+)/) {
	next;
      }

      if ($line =~ /^\#/) {next;}

      if (my ($gff_source, $gff_meth, $gff_start, $gff_end) = ($line =~ /^\S+\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+\S+\s+[-+\.]\s+[012\.]/)) {
	if ($gff_end < $gff_start) {
	  push @problems, "line $line_count:\n$line\nGFF feature end is before the feature start";
	  last;
	}
	# there are 'Genomic_canonical' features longer than 100 Kb
	# there are 'BLAT_NEMATODE' features longer than 100 Kb
	# there are 'Vancouver_fosmid' features longer than 100 Kb
	# F47F6.1c.1 is 107835 bases
	# WBGene00018572 is 107835 bases 
	# Locus lin-42 is 107835 bases
	# Oligo_set Aff_Y116F11.ZZ33 is 105 Kb
	# F16H9.2 is 102695 bases
	# WBGene00008901 is 102695 bases
	my $feature_length = $gff_end - $gff_start;
	if ($feature_length > $MAX_FEATURE_LENGTH && 
	    $gff_source ne 'Link' &&
	    $gff_source ne 'Genomic_canonical' &&
	    $gff_source ne 'BLAT_NEMATODE' &&
	    $gff_source ne 'Vancouver_fosmid' &&
            $gff_meth ne 'sequence_alteration' &&
	    !($line =~ /F47F6.1c.1/) &&
	    !($line =~ /WBGene00018572/) &&
	    !($line =~ /Locus\s+lin-42/) &&
	    !($line =~ /Aff_Y116F11.ZZ33/) &&
	    !($line =~ /F16H9.2/) &&
	    !($line =~ /WBGene00008901/) &&
	    !($line =~ /niDf/) &&
	    !($line =~ /PTC01264/) && # BLAT_NEMBASE sequence in briggsae
	    !($line =~ /PT01856/) && # BLAT_WASHU sequence in briggsae
            !($line =~ /CBG16975:bp208/) && # history CDS in briggsae
	    !($line =~ /PT01856/) && # BLAT_WASHU sequence in briggsae
	    !($line =~ /PTC01264/) && # BLAT_NEMBASE sequence in briggsae
	    !($line =~ /CBN23456.2/) && # coding transcript in briggsae
	    !($line =~ /RST3_519777/) && # RST in elegans
	    !($line =~ /FF095578/) && # BLAT_Caen_EST_OTHER in elegans
	    !($line =~ /"WBsf041392"/) && # segmental_duplication
	    !($line =~ /"WBsf041396"/)    # segmental_duplication
	   ) { 
	  push @problems, "line $line_count:\n$line\nGFF feature is longer than $MAX_FEATURE_LENGTH bases ($feature_length bases)";
# report all length errors
#	  last;
	}
	if (defined $sequence_end && $gff_end > $sequence_end) {
	  push @problems, "line $line_count:\n$line\nfeature is off the end of the sequence";
	  last;
	}
	if ((defined $sequence_start && $gff_start < $sequence_start) || $gff_start < 1) {
	  push @problems, "line $line_count:\n$line\nfeature is before the start of the sequence";
	  last;
	}
      } else {
	push @problems, "line $line_count:\n$line\nis a malformed GFF line";
	last;
      }
    }
    close (CHECK_FILE);
    delete $criteria{gff};
  }
  


  foreach my $c (keys %criteria) {
    push @problems, "unknown criterion in check_file() '$c=>$criteria{$c}'";
  }

  foreach my $problem (@problems) {
      if ($log) {
	$log->error;
	$log->write_to("ERROR: $problem found when checking file '$file'\n");
      }
      carp "ERROR: $problem found when checking file '$file'\n";
  }

  if (!@problems) {$log->write_to("Check file: '$file' OK\n");}
  return @problems;
}



####################################



sub load_to_database {

  my $self     = shift;
  my $database = shift;
  my $file     = shift;
  my $tsuser   = shift;
  my $log      = shift;
  my $bk       = shift;
  my $accept_large_differences = shift;

  my $error=0;
  my $species = $self->species;
  my $version = $self->get_wormbase_version;
  my $prev_version = $version-1;
  my $prev_prev_version = $version-2;
  my $pparse_file = $self->build_data . "/COMPARE/pparse_ace.dat"; # file holding pparse details from previous Builds

  unless ( -e "$file" and -e $database) {
    if ( $log) {
      $log->error;
      $log->write_to("Couldn't find file named: $file or database $database\n");
    }
    print STDERR "Couldn't find file named: $file or database $database\n";
    return 1;
  }

  unless ( -r $file) {
    if ( $log) {
      $log->error;
      $log->write_to("Couldn't read file: $file\n");
    }
    print STDERR "Couldn't read file: $file\n";
    return 1;
  }

  # get the base filename without the path
  my $basename = $file;
  $basename =~ s/.*\///;

  my $st = stat($file);
  if ( $st->size > 50000000 and defined ($bk) and $bk ) {
    $log->write_to("backing up block files before loading $file\n") if $log;
    my $db_dir = $database."/database";
    my $tar_file = "backup.".time.".tgz";
    my $tar_cmd = "tar cvfzP $db_dir/$tar_file $db_dir/block* $db_dir/database.map $db_dir/log.wrm";
    $self->run_command("$tar_cmd", $log);

    # remove old backups keeping the one just made and the previous one.
    my @backups = glob("$db_dir/backup*");
    my %details;
    my @sorted;
    @sorted = sort @backups;
    pop @sorted; pop @sorted;	# remove the newest two
    # . . and delete the rest
    foreach (@sorted) {
      unlink;
    }
  }

  #check whether write access is possible.
  if ( $self->check_write_access($database) eq 'no') {
    print STDERR "cant get write access to $database\n";
    if ($log) {
      $log->log_and_die("cant get write access to $database\n");
    } else {
      die;
    }
  }

  # tsuser is optional but if set, should replace any dots with underscores just in case
  # if not set im using the filename with dots replaced by '_'
  unless ($tsuser) {
    # remove trailing path of filename
    $tsuser = $basename
  }

  $tsuser =~ s/\./_/g;

  my $parsed = 0;
  my $active = 0;		# counts of objects


  # split the ace file if it is large as it loads more efficiently
  my @files_to_load;
  if ($st->size > 5000000) {
    # change input separator to paragraph mode;
    my $oldlinesep = $/;
    $/ = "";
    
    # make temp directory
    my $tmpdir = "$database/tmpace";
    mkpath($tmpdir);
    $self->delete_files_from($tmpdir, "*", "-");

    # split ace file
    my $file_count = 0;
    my $entries = 0;
    my $writing = 0;
    my $split_file;

    open (WBTMPIN, "<$file") || $log->log_and_die ("cant open $file\n"); # open original ace file
    while (my $entry = <WBTMPIN>) {
      $entries++;

      # LongText enties are not formed by contiguous lines of text separated by a blank line
      # they can contain blank lines.
      # So, if we find a LongText class entry, stop splitting the file and simply read in the original file.
      if ($entry =~ /LongText\s+\:\s+/) {
	@files_to_load = ($file);
	last;
      }
      
      if (!$writing) {
	$file_count++;
	$writing = 1;
	$split_file = "$tmpdir/${basename}_split_${file_count}";
	open (WBTMPACE, ">$split_file") || $log->log_and_die ("cant open $split_file\n");
	chmod 0666, $split_file;
	push @files_to_load, $split_file;
      }

      print WBTMPACE "\n\n$entry";
      
      if ($entries > 5000) {
	close(WBTMPACE);
	$writing = 0;
	$entries = 0;
      }
    }

    close(WBTMPIN);
    if ($writing) {
      close (WBTMPACE);
    }

    # reset input line separator
    $/ = $oldlinesep;
  } else {
    push @files_to_load, $file;
  }

  foreach my $file_to_load (@files_to_load) {
    print "Loading $file_to_load ...\n";
    my $tace = $self->tace;
    my $command = <<EOF;
pparse $file_to_load
save
quit
EOF
    if (not open( WRITEDB, "echo '$command'| $tace -tsuser $tsuser $database |" )) {
      if ($log) {
	$log->log_and_die("Could not open write pipe to database\n");
      } else {
	die "Couldn't open pipe to database\n";
      }
    }
    
# expect output like:
#
# acedb> // Parsing file /nfs/disk100/wormpub/BUILD/autoace/acefiles/feature_binding_site.ace
# // objects processed: 154 found, 154 parsed ok, 0 parse failed
# // 51 Active Objects
# acedb> // 51 Active Objects
# acedb>
    
    while (my $line = <WRITEDB>) {
      print "$line";
      if ($line =~ 'ERROR') {
	if ($log) {
	  $log->write_to("ERROR while parsing ACE file $file_to_load\n$line\n");
	  $log->error;
	  $error=1;
	}
      } elsif ($line =~ /objects processed:\s+\d+\s+found,\s+(\d+)\s+parsed ok,/) {
	$parsed += $1;
      } elsif ($line =~ /(\d+)\s+Active Objects/) {
	$active += $1;
      }
    }
    if (not close(WRITEDB)) {
      if ($log) {
	$log->log_and_die("Could not close write pipe to database\n");
      } else {
	die "Could not close write pipe to database\n";
      }
    }
  }
  
  if (! $error) {
    # check against previous loads of this file
    my $last_parsed;		# objects parsed on the previous build
    my $last_active;
    my $last_but_one_parsed=0;		# objects parsed on the previous build but one
    my $last_but_one_active=0;
    # get the number of objects in the pparse of this file in the previous Build
    if (open (PPARSE_ACE, "< $pparse_file")) {
      my $got_last = 0;
      my $got_last_but_one = 0;
      while (my $line = <PPARSE_ACE>) {
	my ($pa_version, $pa_file, $pa_species, $pa_parsed, $pa_active) = split /\s+/, $line;
	if ($pa_version == $prev_version && $pa_file eq $file && $species eq $pa_species) {
	  # store to get the last one in the previous build
	  $last_parsed = $pa_parsed;
	  $last_active = $pa_active;
	  $got_last = 1;
	} 
	if ($pa_version == $prev_prev_version && $pa_file eq $file && $species eq $pa_species) {
	  # store to get the last one in the last but one previous build
	  $last_but_one_parsed = $pa_parsed;
	  $last_but_one_active = $pa_active;
	  $got_last_but_one = 1;
	} 
      }

      # Originally we only stored the basename, then we changed to
      # storing the full pathname because there was confusion with
      # files like the acefiles and the merge_all_databases files
      # having the same basename but different contents.
      #
      # Here we look for the basename during the period of transition
      # which will be at least three Builds, so the $pparse_file will
      # contain both the old entries with $basename in and newer entries
      # with $file in.
      #
      # Eventually the following bit looking for matches to $basename
      # will not be required.

      if (!$got_last || !$got_last_but_one) {
	seek PPARSE_ACE, 0, 0;
	while (my $line = <PPARSE_ACE>) {
	  if (!$got_last) {
	    my ($pa_version, $pa_file, $pa_species, $pa_parsed, $pa_active) = split /\s+/, $line;
	    if ($pa_version == $prev_version && $pa_file eq $basename && $species eq $pa_species) {
	      # store to get the last one in the previous build
	      $last_parsed = $pa_parsed;
	      $last_active = $pa_active;
	    } 
	  }
	
	  if (!$got_last_but_one) {
	    my ($pa_version, $pa_file, $pa_species, $pa_parsed, $pa_active) = split /\s+/, $line;
	  
	    if ($pa_version == $prev_prev_version && $pa_file eq $basename && $species eq $pa_species) {
	      # store to get the last one in the last but one previous build
	      $last_but_one_parsed = $pa_parsed;
	      $last_but_one_active = $pa_active;
	    }
	  } 
	}
      }
    

      close (PPARSE_ACE);
    }

    # check the current Build parse object numbers against the previous one
    if (defined $last_parsed) {
      $log->write_to("File: $file\n") if ($log);
      $log->write_to("Version WS$prev_prev_version parsed $last_but_one_parsed objects OK with $last_but_one_active Active Objects\n") if ($log);
      $log->write_to("Version WS$prev_version parsed $last_parsed objects OK with $last_active Active Objects\n") if ($log);
      $log->write_to("Version WS$version parsed $parsed objects OK with $active Active Objects\n\n") if ($log);
      if (!defined $accept_large_differences) {
	if ($parsed < $last_parsed * 0.9 || $parsed > $last_parsed * 1.1
	    ||  $active < $last_active * 0.9 || $active > $last_active * 1.1) {
	  $log->write_to("*** POSSIBLE ERROR found while parsing ACE file $file\n\n") if ($log);
	  $log->error;
	}
      }
    }

    # now store the details for this pparse
    if (open (PPARSE_ACE, ">> $pparse_file")) {
      if ($version && $basename && $species) {
	print PPARSE_ACE "$version $file $species $parsed $active\n";
      } else {
	$log->write_to("*** POSSIBLE ERROR: Couldn't write to $pparse_file because some of the following is blank\nversion=$version, file=$file, species=$species, parsed=$parsed, active=$active\n\n") if ($log);
	$log->error;
      }
      close (PPARSE_ACE);
    } else {
      $log->write_to("WARNING: Couldn't write to $pparse_file\n\n") if ($log);
    }
  }
}

####################################
# remove any blank lines in a sequence file dumped by acedb

sub remove_blank_lines {
  my ($shift, $file, $log) = @_;

  #$log->write_to("Removing blank lines from $file\n");

  $/ = undef;
  open( CHROM,"< $file") or die("cant open $file to read: $!\n");
  my $chrom = <CHROM>;
  close CHROM;
  $chrom =~ s/^\n//;
  $chrom =~ s/\n\n/\n/g;

  open( CHROM,"> $file") or die("cant open $file to write: $!\n");
  print CHROM $chrom;
  close CHROM;
  $/ = "\n";

  # check it starts with a '>'
  if (substr($chrom, 0, 1) ne '>') {
    die("The first line of $file does not start with a '>' character\n");
  }
  print STDERR substr($chrom, 0, 100) if $ENV{TEST};

}


####################################
sub wormpep_files {
  my $self = shift;
  my $prefix = $self->pepdir_prefix;
  return ( $prefix."pep", $prefix."pep.accession", $prefix."pep.dna", $prefix."pep.history", $prefix."pep.fasta", $prefix."pep.table", $prefix."pep.diff" );
}


sub test         { my $self = shift; return $self->{'test'}; }
sub debug        { my $self = shift; return $self->{'debug'}; }
sub wormpub      { my $self = shift; return $self->{'wormpub'}; }
sub scratch_area { my $self = shift; return $self->{'scratch_area'}; }
sub basedir      { my $self = shift; return $self->{'basedir'}; }
sub autoace      { my $self = shift; return $self->{'autoace'}; }
sub wormpep      { my $self = shift; return $self->{'wormpep'}; }
sub peproot      { my $self = shift; return $self->{'peproot'}; }
sub rnaroot      { my $self = shift; return $self->{'rnaroot'}; }
sub brigpep      { my $self = shift; return $self->{'brigpep'}; }
sub wormrna      { my $self = shift; return $self->{'wormrna'}; }
sub gff_splits   { my $self = shift; return $self->{'gff_splits'}; }
sub chromosomes  { my $self = shift; return $self->{'chromosomes'}; }
sub sequences    { my $self = shift; return $self->{'sequences'}; }
sub spell        { my $self = shift; return $self->{'spell_dir'}; }
sub misc_output  { my $self = shift; return $self->{'misc_output'}; }
sub logs         { my $self = shift; return $self->{'logs'}; }
sub ftp_upload   { my $self = shift; return $self->{'ftp_upload'}; }
sub ftp_staging  { my $self = shift; return $self->{'ftp_staging'}; }
sub ftp_site     { my $self = shift; return $self->{'ftp_site'}; }
sub submit_repos { my $self = shift; return $self->{'submit_repos'}; }
sub reports      { my $self = shift; return $self->{'reports'}; }
sub misc_static  { my $self = shift; return $self->{'misc_static'}; }
sub misc_dynamic { my $self = shift; return $self->{'misc_dynamic'}; }
sub primaries    { my $self = shift; return $self->{'primaries'}; }
sub acefiles     { my $self = shift; return $self->{'acefiles'}; }
sub transcripts  { my $self = shift; return $self->{'transcripts'}; }
sub blat         { my $self = shift; return $self->{'blat'}; }
sub farm_dump    { my $self = shift; return $self->{'farm_dump'}; }
sub compare      { my $self = shift; return $self->{'compare'}; }
sub checks       { my $self = shift; return $self->{'checks'}; }
sub build_data   { my $self = shift; return $self->{'build_data'}; }
sub ontology     { my $self = shift; return $self->{'ontology'}; }
sub orgdb        { my $self = shift; return $self->{'orgdb'}; }
sub cdna_dir     { my $self = shift; return $self->{'cdna_dir'};}
sub cdna_acedir  { my $self = shift; return $self->{'cdna_acedir'};}
sub maskedcdna   { my $self = shift; return $self->{'maskedcdna'} ;}
sub seq_db	 { my $self = shift; return $self->database($self->{'species'});}
#sub ebi          { my $self = shift; return $self->{'ebi'} ;}
sub rnaseq       { my $self = shift; return $self->{'rnaseq'} ;}
sub build_lsfout { my $self = shift; return $self->{'build_lsfout'} ;}
sub genome_seq            { my $self = shift; return $self->{'genome_seq'}; }
sub masked_genome_seq     { my $self = shift; return $self->{'masked_genome_seq'} };
sub softmasked_genome_seq { my $self = shift; return $self->{'smasked_genome_seq'} };
sub genome_diffs { my $self = shift; return $self->{'genome_diff'} };

sub gff {
  my $self = shift;

  if ($self->assembly_type eq 'contig') {
    return $self->sequences;
  } else {
    return $self->chromosomes;
  }
}

# this can be modified by calling script
####################################
sub common_data {
  my $self = shift;
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

####################################
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

####################################
sub primary {
  my $self = shift;
  my $database = shift;
  my $path  = $self->{'primary'}->{"$database"};
  print STDERR "no such primary database:$database\n" unless $path&&(-e $path);
  return $path;
}

####################################
# setter methods
sub set_test { 
  my $self = shift; 
  $self->{'test'} = shift; 
  # adjust the paths to point to the TEST_BUILD versions
  $self->establish_paths;
}
sub set_debug { 
  my $self = shift; 
  $self->{'debug'} = shift; 
}

sub establish_paths {
  my $self = shift;

  # are we runing on an EBI machine?
#  if ($ENV{'HOST'} =~ /ebi/) {$self->{'ebi'} = 1}

  # Some farm code uses Wormbase.pm subs so to maintain this farm code needs to be OO Wormbase.pm compliant but we dont want 
  # multiple paths of the build (main and farm) reading/writing the same wormbase.store file . Store the farm version in ~wormpipe
  # and might as well use path retrieval routines as with main build.

  my ($basedir,
      $ftp_uploads,
      $ftp_site);
  
  $self->{'wormpub'} = "/nfs/panda/ensemblgenomes/wormbase";
  $self->{'scratch_area'} = '/nfs/nobackup/ensemblgenomes/wormbase/scratch';

  # if a specified non-build database is being used
  
  if ( $self->autoace ) {
    ($basedir) = $self->autoace =~ /(.*)\/\w+\/*$/;
    $self->{'orgdb'} = $self->{'autoace'};
  } else {
    $basedir = ($self->test) ? $self->wormpub . "/TEST/BUILD" : $self->wormpub . "/BUILD";
    $self->{'autoace'}    = $self->species eq 'elegans' ? "$basedir/autoace" : "$basedir/".$self->species;
    $self->{'orgdb'}      = $self->{'autoace'}; #."/".$self->{'organism'};
  }
 
  $self->{'basedir'}    = $basedir;

  if ($self->test) {
    $self->{'ftp_upload'} = $self->wormpub . "/TEST/ftp_uploads/wormbase";
    $self->{'ftp_staging'} = $self->wormpub . "/TEST/FTP_STAGING";
    $self->{'ftp_site'}   = $self->wormpub . "/TEST/FTP_site/pub/wormbase";
    $self->{'build_data'} = $self->wormpub . "/TEST/BUILD_DATA";
    $self->{'genome_diff'} = $self->wormpub . "/TEST/CHROMOSOME_DIFFERENCES";
  } else {
    $self->{'ftp_upload'} = "/nfs/ftp/private/worm-ftp/upload";
    $self->{'ftp_site'}   = "/nfs/ftp/pub/databases/wormbase";
    $self->{'ftp_staging'} = $self->wormpub . "/FTP_STAGING";
    $self->{'build_data'} = $self->wormpub . "/BUILD_DATA";
    $self->{'genome_diff'} = $self->wormpub . "/CHROMOSOME_DIFFERENCES";
  }    

  $self->{'submit_repos'}  = $self->wormpub . "/analysis/submissions/" . $self->{'species'};  
  $self->{'peproot'}    = $basedir . "/WORMPEP";
  $self->{'rnaroot'}    = $basedir . "/WORMRNA/";
  $self->{'wormrna'}    = $basedir . "/WORMRNA/".$self->pepdir_prefix."rna" . $self->get_wormbase_version;
  $self->{'wormpep'}    = $basedir . "/WORMPEP/".$self->pepdir_prefix."pep" . $self->get_wormbase_version;

  #species specific paths
  $self->{'logs'}        = $self->orgdb . "/logs";
  $self->{'common_data'} = $self->orgdb . "/COMMON_DATA";
  $self->{'chromosomes'} = $self->orgdb . "/CHROMOSOMES";
  $self->{'sequences'}   = $self->orgdb . "/SEQUENCES";
  $self->{'transcripts'} = $self->orgdb . "/TRANSCRIPTS";
  $self->{'spell_dir'}   = $self->orgdb . "/SPELL";
  $self->{'misc_output'} = $self->orgdb . "/MISC_OUTPUT";
  $self->{'reports'}     = $self->orgdb . "/REPORTS";
  $self->{'acefiles'}    = $self->orgdb . "/acefiles";
  $self->{'gff_splits'}  = $self->orgdb . "/GFF_SPLITS";
  $self->{'primaries'}   = $self->basedir . "/PRIMARIES";
  $self->{'blat'}        = $self->orgdb . "/BLAT";
  $self->{'checks'}      = $self->autoace . "/CHECKS";
  $self->{'ontology'}    = $self->autoace . "/ONTOLOGY";
  $self->{'tace'}   = '/nfs/panda/ensemblgenomes/wormbase/software/packages/acedb/RHEL7/4.9.62/tace';
  $self->{'giface'} = '/nfs/panda/ensemblgenomes/wormbase/software/packages/acedb/RHEL7/4.9.62/giface';
  $self->{'giface_server'} = '/nfs/panda/ensemblgenomes/wormbase/software/packages/acedb/RHEL7/4.9.62/sgifaceserver';
  $self->{'giface_client'} = '/nfs/panda/ensemblgenomes/wormbase/software/packages/acedb/RHEL7/4.9.62/saceclient';
  $self->{'databases'}->{'geneace'} = $self->wormpub . "/DATABASES/geneace";
  $self->{'databases'}->{'camace'}  = $self->wormpub . "/DATABASES/camace";
  $self->{'databases'}->{'current'} = $self->wormpub . "/DATABASES/current_DB";
  $self->{'databases'}->{'autoace'} = $self->autoace;
  
  $self->{'primary'}->{'camace'}  = $self->primaries .'/camace';
  $self->{'primary'}->{'geneace'} = $self->primaries .'/geneace';
  $self->{'primary'}->{'citace'}  = $self->primaries .'/citace';
  $self->{'primary'}->{'caltech'} = $self->primaries .'/citace'; # to handle the various names used
  $self->{'primary'}->{'csh'}     = $self->primaries .'/cshace';
  $self->{'primary'}->{'cshace'}  = $self->primaries .'/cshace';
  $self->{'primary'}->{'briggsae'}= $self->primaries .'/briggsae';
  $self->{'primary'}->{'remanei'}  = $self->primaries .'/remanei';
  $self->{'primary'}->{'japonica'}  = $self->primaries .'/japonica';
  $self->{'primary'}->{'brenneri'} = $self->primaries .'/brenneri';
  $self->{'primary'}->{'brugia'} = $self->primaries .'/brugia';
  $self->{'primary'}->{'ovolvulus'} = $self->primaries .'/ovolvulus';
  $self->{'primary'}->{'sratti'} = $self->primaries .'/sratti';
  $self->{'primary'}->{'pristionchus'} = $self->primaries.'/pristionchus';
  $self->{'primary'}->{'tmuris'} = $self->primaries.'/tmuris';
  
  $self->{'misc_static'} = $self->{'build_data'} . "/MISC_STATIC";
  $self->{'misc_dynamic'} = $self->{'build_data'} . "/MISC_DYNAMIC";
  $self->{'compare'}      = $self->{'build_data'} . "/COMPARE";
  $self->{'cdna_dir'}    = $self->{'build_data'} . "/cDNA/".$self->{'species'};
  $self->{'cdna_acedir'} = $self->{'build_data'} . "/cDNAace/".$self->{'species'};
  $self->{'maskedcdna'}  = $basedir . "/cDNA/".$self->{'species'};
  
  $self->{'farm_dump'}    = $self->scratch_area .'/dumps';
  $self->{'rnaseq'}       = $self->scratch_area . '/RNASeq/'.$self->{species}.'/SRA';
  $self->{'build_lsfout'} = $self->scratch_area . "/BUILD/LSF_OUT/" . $self->{species};
  
  $self->{'genome_seq'} = $self->sequences . "/" . $self->{species} . ".genome.fa";
  $self->{'masked_genome_seq'} = $self->sequences . "/" . $self->{species} . ".genome_masked.fa";
  $self->{'smasked_genome_seq'} = $self->sequences . "/" . $self->{species} . ".genome_softmasked.fa";
  
  # create dirs if missing
  mkpath( $self->logs )        unless ( -e $self->logs );
  mkpath( $self->common_data ) unless ( -e $self->common_data );
  mkpath( $self->wormpep )     unless ( -e $self->wormpep );
  mkpath( $self->wormrna )     unless ( -e $self->wormrna ); 
  mkpath( $self->chromosomes ) unless ( -e $self->chromosomes );
  mkpath( $self->sequences )   unless ( -e $self->sequences );
  mkpath( $self->transcripts ) unless ( -e $self->transcripts ); 
  mkpath( $self->spell )       unless ( -e $self->spell ); 
  mkpath( $self->misc_output ) unless ( -e $self->misc_output ); 
  mkpath( $self->reports )     unless ( -e $self->reports );
  mkpath( $self->ontology )    unless ( -e $self->ontology );
  mkpath( $self->gff_splits )  unless ( -e $self->gff_splits );
  mkpath( $self->primaries )   unless ( -e $self->primaries );
  mkpath( $self->acefiles )    unless ( -e $self->acefiles );
  mkpath( $self->blat )        unless ( -e $self->blat );
  mkpath( $self->checks )      unless ( -e $self->checks );
  mkpath( $self->build_lsfout) unless ( -e $self->build_lsfout );
}

####################################
sub run_script {
  my $self   = shift;
  my $script = shift;
  my $log    = shift;  
  
  my $command = $self->build_cmd($script);
  print "$command\n" if $self->test;
  return $self->run_command( "$command", $log );
}

sub run_script_on_farm {
  my $self   = shift;
  my $script = shift;
  my $log = shift;
  
  my $cmd = 'ssh -t farm-login "'.$self->build_cmd($script).'"';
  return $self->run_command( "$cmd", $log );
}

####################################
# abstracted out of run_script so that scripts using LSF::Manager can 
# still use the Storable etc properly
sub build_cmd {
  my $self   = shift;
  my $script = shift;
  
  my $store_file = $self->build_store;
  my $cmd = $self->build_cmd_line($script, $store_file);
  return $cmd;
}

####################################
# abstracted out of run_script so that scripts using LSF::Manager can 
# still use the Storable etc properly
sub build_store {
  my $self   = shift;
  
  my $species = ref $self;
  my $store = $self->autoace . "/$species.store";
  store( $self, $store );
  
  #if user wormpipe this always gives an ERROR and confuses log msgs
  #$self->run_command( "chmod -f 775 $store") unless ($self->test_user_wormpub == 1);
  return $store;
}

####################################
# abstracted out of run_script so that scripts using LSF::Manager can 
# still use the Storable etc properly
sub build_cmd_line {
  my $self   = shift;
  my $script = shift;
  my $store  = shift;
  
  my $command = "perl $ENV{'CVS_DIR'}/$script -store $store";
  print "$command\n" if $self->test;
  return $command;
}

####################################
sub bsub_script  {
	my $self   = shift;
  	my $script = shift;
  	my $script_sp = shift; #species that called script is to operate on.
  	my $log    = shift;
  	my $species = ref $self;
  	my $store;
  	my $wbobj;
  	if(lc $species eq lc $script_sp) {  	
		$store = $self->autoace . "/$species.store";
		$wbobj = $self;
	}
	else {
		#create a WormBase Species object to retain test / debug status
		my $wb = Wormbase->new ('-test' => $self->test,
								'-debug' => $self->debug,
								'-organism' => lc($script_sp)
								);
		$store = $wb->autoace . "/". (ref $wb) .".store";
		$wbobj=$wb;
	}
  	
	store($wbobj,$store) unless -e $store;
  	my $command = "bsub $ENV{'CVS_DIR'}/$script -store $store";
  	print "$command\n" if $self->test;
  	return $self->run_command( "$command", $log );
}


####################################
sub run_command {
  my $self    = shift;
  my $command = shift;
  my $log     = shift;
  if ($log) {
  	if($log eq 'no_log') {undef $log; }
  }
  else {	
	print STDERR "No log obj passed to run_command by ".(caller)."\n";
  }
  $log->write_to("running $command\n") if $log;
  my $return_status = system("$command");
  if ( ( $return_status >> 8 ) != 0 ) {
    if ( $log ) {
      $log->write_to(" WARNING: $command returned non-zero ($return_status)\n");
      $log->error;
    }
    return 1;
  } else {
    $log->write_to("command exited cleanly\n") if $log;
    return 0;
  }
}

####################################
# THIS IS DEPRECATED.
# Use LSF::JobManager instead
sub wait_for_LSF {
  my $self = shift;
  sleep 10;
  my $jobs = &jobs_left;
  while ( $jobs != 0 ) {
    my $time = 100 * $jobs;
    if ($time > 1800) {$time = 1800;}
    sleep $time;
    $jobs = &jobs_left;
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

####################################
sub checkLSF {
  my ($self, $log) = @_;
  
  if (system("bsub -V")) {
    if ($log) {
      $log->log_and_die("You need to be on an LSF enabled system to run this");
    } else {
      die "You need to be on an LSF enabled system to run this";
    }
  }
}

####################################
sub table_maker_query {
  my($self, $database, $def) = @_;
  my $fh;
  if (!-f $def || !-r $def) {die "Cannot find or read the table_maker def file $def\n"}
  open( $fh, "echo \"Table-maker -p $def\" | ". $self->tace." $database |" ) || die "Couldn't access $database\n";
  return $fh;
}

####################################
# accessor for the 'core' species
sub species_accessors {
	my $self = shift;
	my %accessors;
	foreach my $species (@core_organisms ){
		next if (lc($species) eq $self->species); #$wormbase already exists for own species.
		my $wb = Wormbase->new( -debug   => $self->debug,
			     -test     => $self->test,
			     -organism => $species
			   );
		$accessors{lc $species} = $wb;
	}
	return %accessors;
}

#  if you want to find out if you're a core organism check in this list.
sub core_species {
	return @core_organisms;
}

####################################
# accessor for elegans - required sometime for other species builds
sub build_accessor {
	my $self = shift;
	my $wb = Wormbase->new( -debug   => $self->debug,
			     			-test     => $self->test
			     		);
			     		
	return $wb;
}
####################################
# accessor for the 'tier3' species
sub tier3_species_accessors {
	my $self = shift;
	my %accessors;
	foreach my $species (@tier3_organisms ){
		next if (lc($species) eq $self->species); #$wormbase already exists for own species.
		my $wb = Wormbase->new( -debug   => $self->debug,
			     -test     => $self->test,
			     -organism => $species
			   );
		$accessors{lc $species} = $wb;
	}
	return %accessors;
}

# return Species obj for all tierII and III
sub all_species_accessors {
    my $self = shift;
    my %tier2 = $self->species_accessors;
    my %tier3 = $self->tier3_species_accessors;
    foreach my $t3 (keys %tier3){
	$tier2{$t3} = $tier3{$t3};
    }

    return %tier2;
}


sub species {my $self = shift; return $self->{'species'};}

sub format_sequence
{
	my $self = shift;
	my $seq = shift;
	my $length = shift;
	my $new_seq;
	
	$length = $length ? $length : 60;
	my $left;
	$new_seq = $seq if (length($seq) < $length);
	while ($seq =~ /(\S{$length})/g){
		$new_seq .= $&."\n";
		$left = $';#'
	}
	if ($left) {
		$new_seq .= $left;
	}
	else {
		chomp($new_seq);
	}
	return $new_seq;
}

sub get_binned_chroms {
  my $self = shift;
  my %par = @_;
  my $bin_size;
  $bin_size = $par{'-bin_size'} if (defined $par{'-bin_size'});
  $bin_size = 55 unless (defined $par{'-bin_size'}); # that is actually the number of the bins and not the size of each
  my %opt;
  $opt{-prefix} = 1;
  $opt{-mito} = 1 if $par{'-mito'};
  my @chroms = $self->get_chromosome_names(%opt);
  if (scalar @chroms > 50){
    my @bins;
    my $i=0;
    while ($i<scalar @chroms){
      push (@{$bins[$i % $bin_size]},$chroms[$i]);
      $i++;
    }
    map {$_=join(',',@$_)} @bins;
    @chroms = @bins;
  }
  return \@chroms;
}

###################################
# Method: get_chunked_chroms
#
# Given chunk_total and chunk_id, deterministically chunks chromosomes
# into chunk_total chunks and returns the chunk_id'th chunk. 
###################################

sub get_chunked_chroms {
  my ($self, %opt) = @_;

  my @chroms = sort $self->get_chromosome_names(%opt);
  if (exists $opt{-chunk_id} and exists $opt{-chunk_total} and
      $opt{-chunk_id} =~ /^\d+$/ and $opt{-chunk_total} =~ /^\d+$/ and
      $opt{-chunk_id} <= $opt{-chunk_total}) {
    my ($c_id, $c_tot) = ($opt{-chunk_id}, $opt{-chunk_total});

    if (scalar(@chroms) < $c_tot) {
      $c_tot = scalar(@chroms);
    }

    my @chunks;
    for(my $i=0; $i < @chroms; $i++) {
      my $chunk_idx = $i % $c_tot;
      push @{$chunks[$chunk_idx]}, $chroms[$i];
    }

    @chroms = @{$chunks[$c_id - 1]};
  } else {
    warn("get_chunked_chroms: Could not get sensible values for -chunk_id and -chunk_total " . 
         "  - returning full list\n");
    
  }
  return @chroms;
}


##################################################################################
# Wrapper for handling the opening of GFF files for both chromosome or contig
# based assemblies.
# Parameters :
# Sequence to dump ( CHROMOSOME_X, Crem_Contig1245 )
# Method to dump ( curated, BLAT_EST_BEST or undef for full seq)  
# Log_files object
#
# returns a file handle for reading from
#
# eg my $handle = $wormbase->open_GFF_file('Crem_Contig1245', 'curated', $log)
###################################################################################

sub open_GFF_file {
  my ($self, $seq, $method, $log) = @_;

  my $fh;

  if ($self->assembly_type eq 'contig') {      
    # single file for all seqs for contig-based assemblies
    my $file = $self->GFF_file_name($seq, $method);
    if (not defined $seq) {
      open($fh, $file) or $log->log_and_die("Could not open $file for reading\n");
    } else {
      open($fh, "cat $file | grep \"^$seq\\W\" |") or $log->log_and_die("Could not open seq-restricted stream to $file\n");
    }    
  } else {
    if (defined $seq) {
      my $file = $self->GFF_file_name($seq, $method);
      open($fh, $file) or $log->log_and_die("Could not open $file for reading\n");
    } else {
      my @files;
      foreach my $chr ($self->get_chromosome_names(-mito => 1, -prefix => 1)) {
        push @files, $self->GFF_file_name($chr, $method);
      }
      open($fh, "cat @files |") or $log->log_and_die("Could not open cat stream for files: @files\n");
    }
  }

  return $fh;
}


###################################################################################
# Returns the name of the GFF file to use
# Args: full chromosome name, 
#       [optional] GFF method - this forces the use of the GFF_SPLITS/ file, rather than the CHROMOSOMES/ file

sub GFF_file_name {
  my $self = shift;
  my ($chromosome, $method) = @_;

  my $file;
  if ($self->assembly_type ne 'contig') { 
    $file = defined $method ? $self->gff_splits."/${chromosome}_$method.gff" : $self->gff."/$chromosome.gff";
  } else {			# contig based assembly
    $file = defined $method ? $self->gff_splits."/$method.gff" : $self->gff."/".$self->species.".gff";
  }

  return $file;
}


sub open_GFF3_file {
  my $self = shift;
  my $seq = shift;
  my $method = shift;
  my $log = shift;
  my $handle;

  my $file = $self->GFF3_file_name($seq, $method);
  if ($self->assembly_type ne 'contig' ) { 
    open ($handle,"<$file") or $log->log_and_die("cant open $file :$!\n");
  } else {		# contig based assembly
    open($handle,"grep \"^$seq\\W\" $file |") or $log->log_and_die("cant open $file :$!\n");
  }
  
  return $handle;
}

sub GFF3_file_name {
  my $self = shift;
  my ($chromosome, $method) = @_;

  my $file;
  if ($self->assembly_type ne 'contig') { 
    $file = defined $method ? $self->gff_splits."/${chromosome}_$method.gff3" : $self->gff."/$chromosome.gff3";
  } else {			# contig based assembly
    $file = defined $method ? $self->gff_splits."/$method.gff3" : $self->gff."/".$self->species.".gff3";
  }

  return $file;
}

sub processed_GFF_file {
  my ($self, $gff3) = @_;

  my $fname = sprintf("%s/%s.processed.%s", 
                      $self->sequences,
                      $self->species,
                      ($gff3) ? "gff3" : "gff2");
  
  return $fname;
}


sub submit_repos_hash_from_sequence_name {
  my ($self, $seqname) = @_;

  my $hash = (hex(substr(md5_hex($seqname), -2, 2)) % 64) + 1;

  return $hash;

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

=item *

check_file

Checks the existence of the specified file, checks that it is readable
and writeable and optionally checks other things. Note that this can
be used to check on a directory as well as normal files.

Example:
$wormbase->check_file("file.out", $log, 
		      minsize => 10, 
		      maxsize => 10000,
		      minlines => 10,
		      maxlines => 1000,
		      line1 => '^\S+\s+\S+',
		      line2 => '^#',
		      lines => ['^#', '^\s$', '^[a-z]+'],
		      );


Arguments:
    - filename to check
    - $log
    - optional hash containing one or more of the following:
      must_not_exist => 1 (the file must not exist)
      exists => 1 (the default, the file must exist)
      readonly => 1 (allow the file to be readonly)
      samesize => file_name to compare to
      similarsize => file_name to compare to (within 10% of the size)
      minsize => integer number of bytes
      maxsize => integer number of bytes
      minlines => integer number of lines
      maxlines => integer number of lines
      line1 => regular expression that line 1 of the file must match
      line2 => regular expression that line 2 of the file must match
      lines => [reference of list of regular expressions, all other lines must match at least one of these]
      requires => [reference of list of regular expressions, there must be at least one match of each regular expression in the file]

Returns the number of errors found and sets the error flag

=back 

=over 4

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

tace

tace returns the path for tace that is being used

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


=head2 load_to_database

=head2 SYNOPSIS

=over4

&load_to_database($database, $file, "tsuser", $log );

=back

=head2 USAGE

=over4

Loads specified file in to specified acedb database with tsuser as specified :)
If tsuser not set then the file name will be used (no path, and '_' replacing and '.'s )

=back

=cut

