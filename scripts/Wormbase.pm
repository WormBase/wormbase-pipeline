# Wormbase.pm - module for general use by many Wormbase scripts
# adapted from babel.pl
# put together by krb, but mostly using stuff in babel.pl which 
# was done by dl1 et al.

package Wormbase;

use Exporter;
use Carp;
@ISA       = qw(Exporter);
@EXPORT    = qw(get_cvs_version get_wormbase_version get_old_wormbase_version copy_check mail_maintainer celeaccession tace gff_sort dbfetch find_database);
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
    $WS_version =~ s/NAME //;    
    $WS_version =~ s/ +$//;
    return($WS_version);
}


#################################################################################

sub get_old_wormbase_version {
    my $WS_version = `grep "NAME WS" /wormsrv2/autoace/wspec/database.wrm`;
    chomp($WS_version);
    $WS_version =~ s/NAME //;    
    $WS_version =~ s/ +$//;
    my ($newnumber) = ($WSversion =~ /(\d+$)/);
    my $oldnumber = $newnumber -1;
    my $oldWSversion = "WS".$oldnumber;
    return($oldWSversion);
}


#################################################################################

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

sub find_database {

    my @clones   = @{$_[0]};
    my $handle   = @{$_[1]};
    my $name     = "";
    my %count    = ();
    my @output   = ();
    carp "find_database not called with references\n" if (ref($_[0]) ne 'ARRAY' || ref($_[1]) ne 'HASH');


    if ($handle eq 'cam') {
        my $name = 'camace';
    }
    elsif ($handle eq 'stl') {
        my $name = 'stlace';
    }           


    my $db       = Ace->connect(-path => '/wormsrv2/$name/') || die "Couldn't connect to $name\n", Ace->error;
    my @dbclones = $db->fetch(-query => 'FIND Genome_Sequence');
    foreach my $clone (@clones) {
	my $string = $clone->Confidential_remark(1);
	if ((defined $string) && (($string =~ /Louis/) || ($string =~ /not in Cambridge/))) {
		next;
	}
	else {$count{$clone} = 1;}
    }
    foreach my $clone (@clones) {
        push (@output, $clone) if (exists ($count{$clone})); 
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
the string, i.e. WS48 and not just the number.

=back

=over 4

=item *

get_old_wormbase_version

This subroutine returns the previous WormBase release version.

=back

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
