#!/usr/local/bin/perl
#
# object_count.pl
# v0.1 2000-08-22
# dl
#
# Counts the number of objects in an ACEDB database for each Class stated in the config file
#
#
#
########
#
# 000822 : dl : PP version
#
#
#
#

##################################################
# Connect with acedb database                    #
##################################################

my $acedb = "/wormsrv2/camace";
my $config = "obcount_config";
my $debug = 1;

$|=1;
$ENV{'ACEDB'}="$acedb";
$exec=&tace;

##################################################
# Main Loop                                      #
##################################################

print "# object_count.pl\n# Db = '$acedb'\n\n" if ($debug);

open (CONFIG, "<$config") || die "Failed to open config file '$config'\n\n";
while (<CONFIG>){
    chomp;
    next if ($_ eq "");
    $query = $_;
    &mk_command;
    &count;
}
close (CONFIG);

exit (0);

sub count {
    open (textace, "echo '$command' | $exec  | ");
    while (<textace>) {
	($found = $1)  if (/^\/\/ Found (\d+) objects/);
	($active = $1) if (/^\/\/ (\d+) Active Objects/);
    }
    close (textace);
    (printf "%6s\t'$query'\n",$active) if ($found == $active);
    
}

sub mk_command {
$command=<<EOF;
query find $query 
quit
EOF
}
sub tace {
   local($prog);
   local($name); 
   $name=`uname -sr`;
    
    if ($name=~/^SunOS/) {($prog)=<~wormpub/acedb/ace4/bin.SUN_4/tace>;}
    elsif ($name=~/^IRIX/) {($prog)=<~wormpub/acedb/ace4/bin.SGI_4/tace>;}
    elsif ($name=~/^OSF/)  {($prog)=<~wormpub/ACEDB/bin_ALPHA/giface>;}
    elsif ($name=~/^Linux/)  {($prog)=<~wormpub/acedb/ace4/bin.LINUX/tace>;}

    else {print STDERR "No known binary for $uname\n";exit;}

    return $prog;
}
