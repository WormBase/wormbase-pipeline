# Wormbase.pm - module for general use by many Wormbase scripts
# adapted from babel.pl
# put together by krb, but mostly using stuff in babel.pl which 
# was done by dl1 et al.

package Wormbase;

use Exporter;
use Carp;
use Ace;
use Log_files;
use File::Path;
@ISA       = qw(Exporter);

@EXPORT    = qw(get_wormbase_version get_wormbase_version_name get_wormbase_release_date copy_check mail_maintainer celeaccession tace gff_sort dbfetch clones_in_database open_TCP DNA_string_reverse DNA_string_composition release_databases find_file_last_modified FetchData release_composition release_wormpep test_user_wormpub runtime rundate giface check_write_access Map_feature scan MapFeature delete_files_from load_to_database);

 


#################################################################################

sub get_wormbase_version {

    my $WS_version = `grep "NAME WS" /wormsrv2/autoace/wspec/database.wrm`;
    chomp($WS_version);
    $WS_version =~ s/.*WS//;
    return($WS_version);
}

###################################################################################

sub get_wormbase_version_name {

    my $WS_version_name = `grep "NAME WS" /wormsrv2/autoace/wspec/database.wrm`;
    chomp($WS_version_name);
    $WS_version_name =~ s/NAME //;    
    $WS_version_name =~ s/ +$//;
    return($WS_version_name);
}

###################################################################################

sub get_wormbase_release_date{

  my $format = shift;

  if(!(defined($format))){
    $format = "long";
  }
  elsif($format eq "short"){
    $format = "short";
  }
  elsif($format eq "both"){
    $format = "both";
  }
  else{$format = "long";}

  
  my $line = `ls -l /nfs/disk69/ftp/pub/wormbase/development_release/md5sum.WS*`;
  my @split = split(/\s+/,$line);

  my $month = $split[5];
  my $month2;

  if($month eq "Jan"){ $month = "January";   $month2 = "01";}
  if($month eq "Feb"){ $month = "February";  $month2 = "02";}
  if($month eq "Mar"){ $month = "March";     $month2 = "03";}
  if($month eq "Apr"){ $month = "April";     $month2 = "04";}
  if($month eq "May"){ $month = "May";       $month2 = "05";}
  if($month eq "Jun"){ $month = "June";      $month2 = "06";}
  if($month eq "Jul"){ $month = "July";      $month2 = "07";}
  if($month eq "Aug"){ $month = "August";    $month2 = "08";}
  if($month eq "Sep"){ $month = "September"; $month2 = "09";}
  if($month eq "Oct"){ $month = "October";   $month2 = "10";}
  if($month eq "Nov"){ $month = "November";  $month2 = "11";}
  if($month eq "Dec"){ $month = "December";  $month2 = "12";}

  my $day = $split[6];

  my $day2;
  if (length($day) == 1){$day2 = "0".$day;}
  else{$day2 = $day;}

  if    ($day eq "1") {$day .= "st";}
  elsif ($day eq "2") {$day .= "nd";}
  elsif ($day eq "3") {$day .= "rd";}
  elsif ($day eq "21"){$day .= "st";}
  elsif ($day eq "22"){$day .= "nd";}
  elsif ($day eq "23"){$day .= "rd";}
  elsif ($day eq "31"){$day .= "st";}
  else                {$day .= "th";} 



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

    my ($file,$ref, $dir) = @_;

    # directory to load from can be passed in so that /acari can load files copied over
    $dir = "/wormsrv2/autoace/COMMON_DATA" unless $dir;
    open (FH, "<$dir/$file.dat") or die "can't open $dir/$file.dat\t:$!";
    undef $/;
    my $VAR1;
    my $data = <FH>;
    eval $data;
    die if $@;
    $/ = "\n";
    close FH;
    %$ref = (%$VAR1);    
}

###################################################################################


sub get_script_version {
    my $script = shift;
    my $script_dir = "/wormsrv2/scripts";
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
    my ($name,$maintainer,$logfile) = @_;
    $maintainer = "dl1\@sanger.ac.uk, ar2\@sanger.ac.uk, ck1\@sanger.ac.uk, krb\@sanger.ac.uk, pad\@sanger.ac.uk" if ($maintainer =~ m/All/i);
    croak "trying email a log to a file - this will overwrite the existing file -STOPPING\nAre you passing a file name to Log object? \n" if ( -e $maintainer );
    open (OUTLOG,  "|/bin/mailx -s \"$name\" $maintainer ");
    if ( $logfile )
      {
	open (READLOG, "<$logfile");
	while (<READLOG>) 
	  { print OUTLOG "$_";
	  }
	close READLOG;
      }
    else {
      print OUTLOG "$name";
    }
    close OUTLOG;
  } 


#################################################################################

sub celeaccession {
    local (*text_ace);
    my $seq = shift;
    local($exec);
    $exec="/nfs/disk100/acedb/RELEASE.SUPPORTED/bin.ALPHA_5/tace";
    local($command);
    local($accession);
    $ENV{'ACEDB'}="/wormsrv2/autoace";
    $command=<<EOF;
    find sequence $seq
    show DB_info
    quit
EOF

    open(text_ace, "echo '$command' | $exec  | ");
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
    my $seq = shift;
    $seq =~ tr/[A-Z]/[a-z]/;
    my $A = $seq =~ tr/a/a/;
    my $C = $seq =~ tr/c/c/;
    my $G = $seq =~ tr/g/g/;
    my $T = $seq =~ tr/t/t/;
    my $N = $seq =~ tr/n/n/;
    my $P = $seq =~ tr/-/-/;
    return ($A,$C,$G,$T,$N,$P);
}    

#################################################################################

sub gff_sort {
  while (<>) {
    s/#.*//;
    next unless /\S/;
    @f = split /\t/;
    push @a, $_;
    push @n, $f[0];
    push @s, $f[3];
    push @e, $f[4];
  }
  
  foreach $i (sort { $n[$a] cmp $n[$b] or $s[$a] <=> $s[$b] or $e[$a] <=> $e[$b] } 0..$#a) { print $a[$i] }
  
}

#################################################################################

sub getseqEMBL {

    my $acc = shift;
    my $absent;
    
    my $querycontent1 = "-e+[embl-acc:'$acc']";
    my $querycontent2 = "-e+[emblnew-acc:'$acc']";
    
    my $request1 = "/srs6bin/cgi-bin/wgetz?$querycontent1";
    my $request2 = "/srs6bin/cgi-bin/wgetz?$querycontent2";
    
    my $server = "srs.ebi.ac.uk";

    if (!defined(open_TCP(*F,$server,80))) {
        print "Error connecting to server \n";
        exit(-1);
    }

    # get from emblnew
    print F "GET $request2 HTTP/1.0\n\n";
    print F "Accept: */*\n";
    print F "User-Agent: socketsrs/1.0\n\n";
    
    undef ($absent);
    while (my $return_line=<F>) {
	if ($return_line =~ /error/) {
#	    print "Query returned no entry : '$request2'\n";
	    $absent = 1;
	    last;
	}
    }
    close F;
    
    # else get from embl
    if ($absent) {
	if (!defined(open_TCP(*G,$server,80))) {
	    print "Error connecting to server \n";
	    exit(-1);
	} 
	
#	print "Query [$j]: '$request1'\n";
	
	print G "GET $request1 HTTP/1.0\n\n";
	print G "Accept: */*\n";
	print G "User-Agent: socketsrs/1.0\n\n";
    }
    else {
	if (!defined(open_TCP(*G,$server,80))) {
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

sub open_TCP {
    use Socket;
    my ($FS,$dest,$port) = @_;
    my $proto = getprotobyname ('tcp');
    socket ($FS,AF_INET,SOCK_STREAM,$proto);
    my $sin = sockaddr_in($port,inet_aton($dest));
    connect($FS,$sin) || return undef;
    my $old_fh = select($FS);
    $| = 1;
    select($old_fh);
}


#################################################################################

sub dbfetch {
    local ($file,$database) = @_;
    my $seq = "";
    $/=">";
  
    open (SHOTGUN, "<$database");
    while (<SHOTGUN>) {
	chop $_;
	/^(\S+)\s+/;
	if ($file eq $1) {
	    $seq = $_;
	    last;
	}
    }
    close SHOTGUN;
    $/="\n";
    return ($seq);
}

#################################################################################


# pass a list of clones and cam or stl

sub clones_in_database {

    my $handle   = ${$_[0]};
    my %count    = ();
    my @output   = ();
    my ($db,$switch);
    if ($handle eq 'cam') {
        $switch = 1;
        $db       = Ace->connect(-path => '/wormsrv2/camace') || die "Couldn't connect to $name\n", Ace->error;
    }
    elsif ($handle eq 'stl') {
        $db       = Ace->connect(-path => '/wormsrv2/stlace') || die "Couldn't connect to $name\n", Ace->error;
    }
    else {
        $db       = Ace->connect(-path => '/nfs/disk100/wormpub/DATABASES/current_DB') || die "Couldn't connect to $name\n", Ace->error;
    }
    my @dbclones = $db->fetch(-query => 'FIND Genome_Sequence');
    
    foreach my $clone (@dbclones) {
	if ($switch == 1) {
            my $string = $clone->Confidential_remark(1);
	    if ((defined $string) && (($string =~ /Louis/) || ($string =~ /not in Cambridge/))) {
		    next;
	    }
	    else {push @output, $clone;}
        }
        else {push @output, $clone;}    
    }
    $db->close;
    return @output;   

}



##################################################################################

#  Subroutines for generating release letter


##################################################################################
# Databases used in build                                                        #
##################################################################################
sub release_databases
  {
    #GET THE DATABASES USED IN THIS BUILD

    #get the non_cambridge DB info from file
    open (DBS,"/wormsrv2/autoace/logs/Primary_databases_used_in_build");
    my ($stlace,$citace, $brigdb, $cshace);
    my %dates;
    while(<DBS>)
      {
	chomp;
	my @info = split(/ : /,$_);
	$dates{$info[0]} = $info[1];
      }
    foreach my $key( keys %dates)
      {
	my $oldstyle = $dates{$key};
	my $newstyle = "20".substr($oldstyle,0,2)."-".substr($oldstyle,2,2)."-".substr($oldstyle,4);
	$dates{$key} = $newstyle;
      }
    
    #get the Cambridge dates directly from block file 1
    my @date = &find_file_last_modified("/wormsrv2/camace/database/block1.wrm");
    $dates{camace} = $date[0];
    @date = &find_file_last_modified("/wormsrv2/geneace/database/block1.wrm");
    $dates{genace} = $date[0];
    

    #PARSE THE RELEASE LETTER FOR LAST BUILD INFO
    my $old_ver = &get_wormbase_version -1;
    my $ver = $old_ver+1;
    my %old_dates;
    my $located = 0;
    my $count =0;  #determines how many lines to read = no databases
    open (OLD_LETTER,"/wormsrv2/autoace/RELEASE_LETTERS/letter.WS$old_ver");
    while(<OLD_LETTER>)
      {
	if( ($located == 1) && ($count <= 6))
	  {
	    chomp;
	    my @info = split(/ : | - /,$_);
	    #this will put some crap in the hash from the 1st two line of the section but it no matter
	    $old_dates{$info[0]} = $info[1];
	    $count++;
	  }
	elsif ($_ =~ m/Primary/){
	  $located = 1;
	}
      }  
    my $dbaseFile = "/wormsrv2/autoace/RELEASE_LETTERS/dbases";
    open (WRITE, ">$dbaseFile");
    
    print WRITE "Primary databases used in build WS$ver\n------------------------------------\n";
    foreach my $key(sort keys %dates)
      {
	print WRITE "$key : $dates{$key}";
	if ("$dates{$key}" gt "$old_dates{$key}") {
	  print WRITE " - updated\n";
	}
	elsif ("$dates{$key}" lt "$old_dates{$key}"){
	  print WRITE "you're using a older version of $key than for WS$old_ver ! ! \n";
	}
	else{
	  print WRITE "\n";
	}
      }
    close WRITE;

    my $name = "Database update report";
    my $maintainer = "All";
    &mail_maintainer($name,$maintainer,$dbaseFile);

    return 1;
 }

##################################################################################
# Returns the date yyyy-mm-dd and time hh:mm:ss file was last modified           #
##################################################################################
sub find_file_last_modified
  { 
    my $filename = shift;
    open (FILE,"<$filename") || die "cant open file $filename\n";
    my @fileinfo = stat FILE;
    my @date = localtime($fileinfo[9]);
    close(FILE);
    my $year = sprintf("%d-%02d-%02d\n",$date[5]+1900,$date[4]+1,$date[3]);
    my $time = "$date[2]:$date[1]:$date[0]";
    chomp $year;

    my @last_modified = ($year, $time);
    return @last_modified;
  }




##################################################################################
# DNA Sequence composition                                                       #
##################################################################################
sub release_composition
  {
    #get the old info from current_DB
    my $ver = &get_wormbase_version;
    my $old_ver = $ver -1;
    my %old_data;
    $old_data{"-"} = 0; # initialise to avoid problems later if no gaps
    my $old_letter ="/wormsrv2/WS$old_ver/CHROMOSOMES/composition.all";
    open (OLD, "<$old_letter") || die "cant open data file - $old_letter";
    while (<OLD>) {
      chomp;
	if ($_ =~ m/(\d+)\s+total$/){
	  #my $tot = $1;$tot =~ s/,//g; # get rid of commas
	  $old_data{Total} = $1;
	}
	elsif ($_ =~ m/^\s+([\w-]{1})\s+(\d+)/){
	  $old_data{"$1"} = $2;
	}
    }
    close (OLD);
    
    #now get the new stuff to compare  
    my $new_letter ="/wormsrv2/autoace/CHROMOSOMES/composition.all";
    my %new_data;
    $new_data{"-"} = 0; # initialise to avoid problems later if no gaps
    open (NEW, "<$new_letter") || die "cant open data file - $new_letter";
    while (<NEW>) {
      chomp;
      if ($_ =~ m/(\d+)\s+total$/){
	$new_data{Total} = $1;
      }
      elsif ($_ =~ m/^\s+([\w-]{1})\s+(\d+)/){
	$new_data{"$1"} = $2;
      }
    }
    close NEW;

    #now check the differences
    my %change_data;
    my $compositionFile = "/wormsrv2/autoace/RELEASE_LETTERS/composition";
    open (COMP_ANALYSIS, ">$compositionFile") || die "cant open $compositionFile";
    print COMP_ANALYSIS "Genome sequence composition:\n----------------------------\n\n";
    print COMP_ANALYSIS "       \tWS$ver       \tWS$old_ver      \tchange\n";
    print COMP_ANALYSIS "----------------------------------------------\n";
    foreach my $key(keys %old_data) {
	$change_data{$key} = $new_data{$key} - $old_data{$key};
      }
   
    my @order = ("a","c","g","t","n","-","Total");
    foreach(@order){
      if ("$_" eq "Total"){
	print COMP_ANALYSIS "\n";
      }
      printf COMP_ANALYSIS ("%-5s\t%-8d\t%-8d\t%+4d\n", $_, $new_data{$_}, $old_data{$_}, $change_data{$_});
    }

    if ($change_data{"-"} > 0){
      print COMP_ANALYSIS "Number of gaps has increased - please investigate ! \n";
    }
    
    if ($change_data{"Total"} < 0) {
      print COMP_ANALYSIS "Total number of bases has decreased - please investigate ! \n";
    }
    close COMP_ANALYSIS;

    my $name = "Sequence composition report";
    my $maintainer = "All";
    &mail_maintainer($name,$maintainer,$compositionFile);
    
    return 1;
  }

##################################################################################
#  Wormpep                                                                       #
##################################################################################
sub release_wormpep       #($number_cds $number_total $number_alternate )
  {  
    my ($number_cds, $number_total, $number_alternate) = @_;
    my $ver = &get_wormbase_version;
    my $old_ver = $ver -1;
    
    #extract data from new wormpep files
    my $lost = `more /wormsrv2/WORMPEP/wormpep$ver/wormpep.diff$ver | grep 'lost' | wc -l`;
    my $new = `more /wormsrv2/WORMPEP/wormpep$ver/wormpep.diff$ver | grep 'new' | wc -l`;
    my $changed = `more /wormsrv2/WORMPEP/wormpep$ver/wormpep.diff$ver | grep 'changed' | wc -l`;
    my $appeared = `more /wormsrv2/WORMPEP/wormpep$ver/wormpep.diff$ver | grep 'appear' | wc -l`;
    my $entries = `cat /wormsrv2/WORMPEP/wormpep$ver/wormpep.diff$ver | wc -l`;
    my $net = $new + $appeared - $lost;
    my $codingDNA;



    #get no of coding bases from log file
    open (THIS_LOG,"/wormsrv2/WORMPEP/wormpep$ver/wormpep_current.log");
    while(<THIS_LOG>){
	if( $_ =~ /No\. of sequences \(letters\) written:\s+\d+\,\d+\s+\((.*)\)/ ){
	 $codingDNA = $1;}
      }

    #write new letter
    my $wormpepFile = "/wormsrv2/autoace/RELEASE_LETTERS/wormpep";
    open (LETTER, ">$wormpepFile") || die "cant open $wormpepFile\n";


    print LETTER "\n\nWormpep data set:\n----------------------------\n";
    print LETTER"\nThere are $number_cds CDS in autoace, $number_total when counting $number_alternate alternate splice forms.\n
The $number_total sequences contain $codingDNA base pairs in total.\n\n";
   
    print LETTER "Modified entries      $changed";
    print LETTER "Deleted entries       $lost";
    print LETTER "New entries           $new";
    print LETTER "Reappeared entries    $appeared\n";
    printf LETTER "Net change  %+d",$net;

    #get the number of CDS's in the previous build
    open (OLD_LOG,"/wormsrv2/WORMPEP/wormpep$old_ver/wormpep_current.log");
    my $oldCDS;
    while(<OLD_LOG>){
	if( $_ =~ /No\. of sequences \(letters\) written:\s+(\d+\,\d+)\s+\(.*\)/ ){
	$oldCDS = $1;
	$oldCDS =~ s/,//g;
      }
    }
    close OLD_LOG;

    #check
    my $mail;
    if ($lost + $new + $changed + $appeared != $entries) {
	print LETTER "cat of wormpep.diff$ver does not add up to the changes (from $0)";
      }
    if ($oldCDS + $net != $number_total){
      print LETTER"The differnce between the total CDS's of this ($number_total) and the last build ($oldCDS) does not equal the net change $net\nPlease investigate! ! \n";
    }
       
    close LETTER; 

    my $name = "Wormpep release stats";
    my $maintainer = "All";
    &mail_maintainer($name,$maintainer,$wormpepFile);

    return 1;
  }
#end of release letter generating subs
#############################################


sub test_user_wormpub
  {
    my $name = `whoami`;
    chomp $name;
    if( "$name" eq "wormpub" ){
      print "running scripts as user wormpub . . . \n\n";
      return;
    }
    else {
      print "You are doing this as $name NOT wormpub ! \n\n If you are going to alter autoace in any way it will break.\nDo you want to continue? (y/n). . ";
	my $response = <STDIN>;
      chomp $response;
      if( "$response" eq "n" ){
	die "probably for the best !\n";
      }
      else {
	print "OK - on your head be it !\nBack to the script . .\n#########################################################\n\n\n";
	return;
      }
    }
  }

#############################################
sub runtime {
  my $runtime    = `date +%H:%M:%S`; chomp $runtime;
  return $runtime;
}
###############################################
sub rundate{
  my $rundate = `date +%y%m%d`; chomp $rundate;
  return $rundate;
}

###################################################
# subs to get the correct version of ACEDB binaries
sub tace {
  my $tace = glob("~wormpub/ACEDB/bin_ALPHA/tace");
  return $tace;
}

sub giface {
  my $giface = glob("~wormpub/ACEDB/bin_ALPHA/giface");
  return $giface;
}

####################################
# Check for database write access
####################################

sub check_write_access{

  my $database = shift;
  my $write_access = "yes";

  $write_access = "no" if (-e "${database}/database/lock.wrm");
  return($write_access);

}


####################################
# Delete files from directory
####################################
sub delete_files_from {
  my ($directory,$pattern,$folder)  = @_;
  my $file;
  my $delete_count = 0;
  my $fail_warn = 1;

  return undef unless (-e $directory );

  if ( $folder eq "+" ) {
    print "Removing entire dir and subdirs of $directory\n";
    $delete_count = rmtree($directory);
  } 
  else {
    opendir (TO_GO,$directory) or die "cant get listing of $directory:\t$!\n";
#   $pattern = "." if $pattern eq "*";

    $pattern = "." unless $pattern;
    $pattern =~ s/\*/\./g;

    while ( $file = readdir(TO_GO)) {
      next if( $file eq "." or $file eq "..") ;
      if ( $file =~ /$pattern/ ) {
	if( unlink ("$directory/$file") ) {
	  $delete_count++;
	}else {
	   warn "couldn't unlink $directory/$file :\t$!\n";
	   undef $fail_warn;
	 }
      }
    }
  }
  return $fail_warn ? $delete_count : $fail_warn ; # undef if failed else no. files removed;
}



sub load_to_database {

  my $database = shift;
  my $file = shift;

  # tsuser is optional but if set, should replace any dots with underscores just in case
  # if not set im using the filename with dots replaced by '_'
  my $tsuser = shift; 

  unless ( $tsuser ) { 
    # remove trailing path of filename
    $tsuser = $file;
    $tsuser =~ s/.*\///;
  }

  $tsuser =~ s/\./_/g;

  unless (-e "$file"){
    die " Couldn't find file named: $file\n";
  }
  my $exe = &tace;
  my $command = "pparse $file\nsave\nquit\n";
 
  open (WRITEDB, "| $exe -tsuser $tsuser $database ") || die "Couldn't open pipe to database\n";
  print WRITEDB $command;
  close (WRITEDB);
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


=cut

=back

=over 4

=item *
gff_sort, tace, and dbfetch

Dont know what these do, sorry.

=back 

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

