#!/usr/local/bin/perl -w
# compares clones in different acedb databases                                   
#
# by Kerstin Jekosch
# 10/07/01

use strict;
use Ace;
my (%camace, %stlace, %autoace);


my $camdb     = Ace->connect(-path => '/wormsrv2/camace/') || die "Couldn't connect to  camace\n", Ace->error;
my @camclones = $camdb->fetch(-query => 'FIND Genome_Sequence');
foreach my $camclone (@camclones) {
	my $string = $camclone->Confidential_remark(1);
	if ((defined $string) && (($string =~ /not in Cambridge LINK/) || ($string =~ /Louis/))) {
		next;
	}
	else {$camace{$camclone} = 1;}
}

my $stldb     = Ace->connect(-path => '/wormsrv2/stlace/') || die "Couldn't connect to  stlace\n", Ace->error;
my @stlclones = $stldb->fetch(-query => 'FIND Genome_Sequence');
foreach my $stlclone (@stlclones) {
	$stlace{$stlclone} = 1;
}


my $autodb     = Ace->connect(-path => '/wormsrv2/autoace/') || die "Couldn't connect to  autoace\n", Ace->error;
my @autoclones = $autodb->fetch(-query => 'FIND Genome_Sequence');
foreach my $autoclone (@autoclones) {
	$autoace{$autoclone} = 1;
}


#foreach my $clone (sort keys %autoace) {
#	print "$clone\n" if ((!exists $stlace{$clone}) && (!exists $camace{$clone})) 	
#}
