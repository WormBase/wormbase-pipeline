#!/usr/local/bin/perl


$script_dir = "/nfs/disk100/wormpub/analysis/scripts";

&get_version(make_autoace);
print "$version\n";



###############################################
# get_script_version                          #
###############################################
#
# 001020 : dl  : PP version
#
# Requires :
#    $script_dir = path of worm scripts
#
# Usage    :
#    &get_version(<script_name>);
#
# Output   :
#    $version = version number of script being used
#

sub get_script_version {
    my $script = shift;
    open (LIST, "/bin/ls -l $script_dir/$script |");
    while (<LIST>) {
	chomp;
	my $stringlen = length ($_);
	$version = substr ($_,$stringlen-3,3);
	return $version;
    }
    close (LIST);
}



###############################################
# Prints help and disappears                  #
###############################################

sub PrintHelp {
   exec ('perldoc',$0);
}

__END__

=pod

=head1 NAME - shite.pl

=head2 USAGE

=cut
