# Wormbase.pm - module for general use by many Wormbase scripts
# adapted from babel.pl
# put together by krb, but mostly using stuff in babel.pl which 
# was done by dl1 et al.

package Wormbase;

use Exporter;
use Carp;
use Ace;
@ISA       = qw(Exporter);
@EXPORT    = qw(get_cvs_version get_wormbase_version get_wormbase_version_name get_wormbase_release_date copy_check mail_maintainer celeaccession tace gff_sort dbfetch clones_in_database open_TCP DNA_string_reverse DNA_string_composition);
@EXPORT_OK = qw(get_script_version); 


#################################################################################

sub get_cvs_version{
  my $script_name = shift;
  my $version = `cvs -d /nfs/ensembl/cvsroot/ status $script_name`;
  $version =~ s/.* +Repository revision:\s+([\d\.]*)\s+.*/$1/s; 
  return($version);
}


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

  
  my $line = `ls -l /nfs/disk69/ftp/pub/wormbase/current_release/md5sum.WS??`;
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
    open (OUTLOG,  "|/usr/bin/mailx -s \"$name\" $maintainer ");
    open (READLOG, "<$logfile");
    while (<READLOG>) {
	print OUTLOG "$_";
    }
    close READLOG;
    close OUTLOG;
} 


#################################################################################

sub celeaccession {
    local (*text_ace);
    my $seq = shift;
    local($exec);
    $exec="/nfs/disk100/acedb/RELEASE.SUPPORTED/bin.ALPHA_4/tace";
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
        if (/\s+Database\s+\S+\s+\S+\s+(\S+)/) {
            $accession=$1;
        }
    }
    close text_ace;
    return $accession;
}


#################################################################################

sub tace {
   local($prog);
   local($name); 
   $name=`uname -sr`;
    if ($name=~/^SunOS/)    {($prog)=<~wormpub/acedb/ace4/bin.SUN_4/tace>;}
    elsif ($name=~/^IRIX/)  {($prog)=<~wormpub/acedb/ace4/bin.SGI_4/tace>;}
    elsif ($name=~/^OSF/)   {($prog)=<~acedb/RELEASE.DEVELOPMENT/bin.ALPHA_4/tace>;}
    elsif ($name=~/^Linux/) {($prog)=<~wormpub/acedb/ace4/bin.LINUX/tace>;}
    else {print STDERR "No known binary for $uname\n";exit;}
    return $prog;
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
    $/=">";

    open (SHOTGUN, "<$database");
    while (<SHOTGUN>) {
	chop $_;
	/^(\S+)\s+/;
	if ($file eq $1) {
	    print $_;
	    last;
	}
    }
    close SHOTGUN;
    $/="\n";
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
        $db       = Ace->connect(-path => '/wormsrv2/current_DB') || die "Couldn't connect to $name\n", Ace->error;
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
    return @output;   

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

get_cvs_version

If you pass the name of a script to this subroutine, it will return the
cvs version number.  It is safest to pass the $0 variable but only if
your script is present in /wormsrv2/scripts.  I.e. the script calling
get_cvs_version must be present in a CVS checked-out directory.
Returns the latest CVS version number.

This subroutine replaces the (deprecated) get_script_version (see below).

=back

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

= item *
find_database

Pass a list of clones and 'cam' or 'stl' and it will return the clones 
in camace / stlace

=back

=over 4


=cut

=back

=over 4

=item *
gff_sort, tace, and dbfetch

Don't know what these do, sorry.

=cut
