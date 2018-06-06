#!/usr/bin/env perl
use JSON;
use Config::IniFiles;

my @ans;
my $cfg = Config::IniFiles->new( 
  -file => @ARGV ? @ARGV[0] : STDIN
) ; 
die unless $cfg;
for my $section ( grep /source/, $cfg -> Sections ) { 
 my @uris = $cfg->val( $section , "data_uri") or ();
 for my $uri ( @uris){ 
   next if $section =~ /brenneri/ or $section =~/remanei/;
   next if $cfg->val($section, "download") eq "N";
   my %ans_here;
   for my $param ($cfg->Parameters($section)){
     next if $param eq "data_uri" or $param eq "download";
     $ans_here{$param} = $cfg->val($section, $param) if $cfg->val($section, $param);
     $ans_here{$param} = int($ans_here{$param}) if $ans_here{$param} =~ /[0-9]+/;
   }
   $uri =~ s/RELEASE/WS$ENV{WORMBASE_VERSION}/g;
   $ans_here{file}= $uri;
   push @ans, \%ans_here;
 }
}
print JSON->new->pretty(1)->encode(\@ans);
