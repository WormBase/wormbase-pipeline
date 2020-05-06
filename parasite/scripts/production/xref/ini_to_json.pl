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
   # this FTP site is available now
   # Tim S. 2020-05-06
   # if ($ENV{WORMBASE_VERSION} > 271) { 
   #   die "$uri Check this case - I think they were missing from FTP but should be back now" if $section =~/remanei/;
   # }
   # I believe this next line also needs to be removed now -- otherwise the line above doesn't make sense
   # Tim S. 2020-05-06
   # next if $section =~/remanei/;
   next if $cfg->val($section, "download") eq "N";
   my %ans_here;
   for my $param ($cfg->Parameters($section)){
     next if $param eq "data_uri" or $param eq "download";
     $param = "release" if $param eq "release_uri";
     $ans_here{$param} = $cfg->val($section, $param) if $cfg->val($section, $param);
     $ans_here{$param} = int($ans_here{$param}) if $ans_here{$param} =~ /^[0-9]+$/;
   }
   $uri =~ s/RELEASE/WS$ENV{WORMBASE_VERSION}/g;
   $ans_here{file}= $uri;
   push @ans, \%ans_here;
 }
}
print JSON->new->pretty(1)->encode(\@ans);
