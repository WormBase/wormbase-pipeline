#!/usr/local/bin/perl
#
# Perl library for C.elegans scripts
#
# 001020 : dl  : Added 'get_script_version'
#              : returns the active version for the script name called
#
#
#################################################################################
#                                                                               #
# CONTENTS:                                                                     #
#                                                                               #
# get_script_version : returns the active version for the script name called    #
# copy_check         : returns +ve if the two files are of the same size        #
# mail_maintainer    : mails the logfile at the end of a script run             #
# celeaccession      : returns the accession number for the given cosmid        #
#                                                                               #
#################################################################################

$script_dir = "/wormsrv2/scripts";

#################################################################################
# get_script_version                                                            #
#                                                                               #
# Returns the version number of the script which is currently active            #
#                                                                               #
# 001020 : dl  : PP version                                                     #
#                                                                               #
# Requires :                                                                    #
#    $script_dir = path of worm scripts                                         #
#                                                                               #
# Usage    :                                                                    #
#    &get_version(<script_name>);                                               #
#                                                                               #
# Returns  :                                                                    #
#    $version = version number of script being used                             #
#                                                                               #
#################################################################################

sub get_script_version {
    my $script = shift;
    open (GET_SCRIPT_LIST, "/bin/ls -l $script_dir/$script |");
    while (<GET_SCRIPT_LIST>) {
	chomp;
	my $stringlen = length ($_);
	$version = substr ($_,$stringlen-3,3);
	last;
    }
    close GET_SCRIPT_LIST;
    return $version;
} # end of sub 'get_script_version'

#################################################################################




#################################################################################
# copy_check                                                                    #
#                                                                               #
# Returns +ve if the two files are of the same size                             #
#                                                                               #
# 001214 : ag3  : PP version                                                    #
#                                                                               #
# Requires :                                                                    #
#    none                                                                       #
#                                                                               #
# Usage    :                                                                    #
#    &copy_check(<file1>,<file2>);                                              #
#                                                                               #
# Returns  :                                                                    #
#    $match = '1' if files are of the same size                                 #
#                                                                               #
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
} # end of sub 'copy_check'

#################################################################################




#################################################################################
# mail_maintainer                                                               #
#                                                                               #
# Mails the logfile at the end of a script run                                  #
#                                                                               #
# 001214 : dl1  : PP version                                                    #
#                                                                               #
# Requires :                                                                    #
#    none                                                                       #
#                                                                               #
# Usage    :                                                                    #
#    &mail_maintainer(<title>,<maintainer e-mail list>,<logfile>);              #
#                                                                               #
# Returns  :                                                                    #
#    none                                                                       #
#                                                                               #
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
} # end of sub 'mail_maintainer'




#################################################################################
# DbWrite                                                                       #
#                                                                               #
# Mails the logfile at the end of a script run                                  #
#                                                                               #
# 001214 : dl1  : PP version                                                    #
#                                                                               #
# Requires :                                                                    #
#    none                                                                       #
#                                                                               #
# Usage    :                                                                    #
#    &mail_maintainer(<title>,<maintainer e-mail list>,<logfile>);              #
#                                                                               #
# Returns  :                                                                    #
#    none                                                                       #
#                                                                               #
#################################################################################

# end of sub 'mail_maintainer'

#################################################################################



#################################################################################
# celeaccession                                                                 #
#                                                                               #
# Returns the accession number for the given cosmid                             #
#                                                                               #
# 001103 : dl  : PP version (from an ag3  subvroutine)                          #
#                                                                               #
# Requires :                                                                    #
#    $acedb_path = path of acedb database to be queried                         #
#    $tace_path  = path of tace executable to use                               #
#                                                                               #
# Usage    :                                                                    #
#    &celeaccession(<clone_name>);                                              #
#                                                                               #
# Returns  :                                                                    #
#    $accession = EMBL accession number of clone                                #
#                                                                               #
#################################################################################

sub celeaccession {
    local (*textace);
    my $seq = shift;
    local($exec);
    $exec="/nfs/disk100/acedb/RELEASE.SUPPORTED/bin.ALPHA_4/tace";
    local($command);
    local($accession);
    $ENV{'ACEDB'}="/wormsrv2/camace";
    $command=<<EOF;
    find sequence $seq
    show DB_info
    quit
EOF

open(textace, "echo '$command' | $exec  | ");
while (<textace>) {
    if (/\s+Database\s+EMBL\s+\S+\s+(\S+)\n/) {$accession=$1;}
    }
close textace;
return $accession;
}



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



################################################################################
#Return a true value
################################################################################

1;

