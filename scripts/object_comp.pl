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

my $acedb2 = "/wormsrv2/autoace";
my $acedb  = "/wormsrv2/WS16";
my $config = "/wormsrv2/camace/bin/obcomp2_config";
my $error = "obcomp2_error";
my $debug = 0;
my $oblist = 1;
$|=1;
$exec=&tace;

##################################################
# Main Loop                                      #
##################################################

print "# object_comp.pl\n# Db-1 = '$acedb'\n# Db-2 = '$acedb2'\n\n" if ($debug);

my $db1 = $acedb;
my $db2 = $acedb2;

$db1 =~ s/\/wormsrv2\///;
$db2 =~ s/\/wormsrv2\///;


print  " +------------------------------------------+\n";
print  " | Class              |    ACEDB database   |\n";
print  " |                    +----------+----------+--------+\n";
printf " |                    | %8s | %8s | change |\n", $db1,$db2;
print  " +--------------------+----------+----------+--------+\n";


open (ERROR,  ">$error")  || die "Failed to open error file '$config'\n\n";
open (CONFIG, "<$config") || die "Failed to open config file '$config'\n\n";
while (<CONFIG>){
    chomp;
    next if ($_ eq "");
    $query = $_;
    &mk_command;

    printf " | %18s |", $query;

    ########
    # DB-1 #
    ########

    $ENV{'ACEDB'}="$acedb";
    $out = "DB-1_${oblist}";
    &count($out);
    $count_1 = $active;

    print " | ";

    ########
    # DB-2 #
    ########

    $ENV{'ACEDB'}="$acedb2";
    $out = "DB-2_${oblist}";
    &count($out);
    $count_2 = $active;

#    print " | ";

# printf "%6s",$active) if ($found == $active);
   
    my $diff = $count_2 - $count_1;
    if ($diff != 0) {
	printf "| %6s |\n",$diff;
	print ERROR "--> $query\n";
	&diff($oblist);
    }
    else {
	print "|      0 |\n";
    }
    $oblist++;
}
close (CONFIG);
print  " +--------------------+----------+----------+--------+\n";


###########
# TIDY UP #
###########

system ("rm -f DB-1*");
system ("rm -f DB-2*");


exit (0);

sub diff {
    my ($oblist) = @_;
    
    system ("cat DB-1_${oblist} | sort > look-1");
    system ("cat DB-2_${oblist} | sort > look-2");

    open (COMM, "comm -3 look-1 look-2 |");
    while (<COMM>) {
	print ERROR " -> DB-1 $1\n" if (/^(\S+.+)/);
	print ERROR " -> DB-2 $1\n" if (/^\s+(\S+.+)/);
    }
    close (COMM);

    system ("rm -f look-1");
    system ("rm -f look-2");

}



sub count {
    my ($out) = @_;
    open (out, ">$out");
    open (textace, "echo '$command' | $exec  | ");
    while (<textace>) {
	($found = $1)  if (/^\/\/ Found (\d+) objects/);
	($active = $1) if (/^\/\/ (\d+) Active Objects/);
	(print out "$_") if (/\:/);
    }
    close (textace);
    (printf "%9s",$active) if ($found == $active);
    close (out);
    return ($active);
}

sub mk_command {
$command=<<EOF;
query find $query 
list -a 
quit
EOF
}
sub tace {
   local($prog);
   local($name); 
   $name=`uname -sr`;
    
    if ($name=~/^SunOS/)    {($prog)=<~wormpub/acedb/ace4/bin.SUN_4/tace>;}
    elsif ($name=~/^IRIX/)  {($prog)=<~wormpub/acedb/ace4/bin.SGI_4/tace>;}
    elsif ($name=~/^OSF/)   {($prog)=<~wormpub/acedb/ace4/bin.ALPHA_4/giface>;}
    elsif ($name=~/^Linux/) {($prog)=<~wormpub/acedb/ace4/bin.LINUX/tace>;}

    else {print STDERR "No known binary for $uname\n";exit;}

    return $prog;
}
