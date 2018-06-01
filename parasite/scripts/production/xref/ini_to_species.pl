#!/usr/bin/env perl
use Config::IniFiles;

my @ans;
my $cfg = Config::IniFiles->new( -file => @ARGV ? @ARGV[0] : STDIN );
die unless $cfg;
for my $section ( grep /species/, $cfg->Sections ) {
    for my $name (split "," , $cfg->val( $section, "aliases")){
        if ( grep /wormbase/, $cfg->val( $section, "source" ) ) {
            print "$name\twormbase\n";
        }
        else {
            print "$name\tparasite\n";
        }
    }
}
