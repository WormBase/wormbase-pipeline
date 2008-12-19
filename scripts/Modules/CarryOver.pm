#!/software/bin/perl -w
#
# CarryOver.pm                       
# 
# by Anthony Rogers
#
# Carry over data from one build to the next when species not being built
#
# Last updated by: $Author: ar2 $     
# Last updated on: $Date: 2008-12-19 12:32:30 $      

=pod

=head1 NAME

 PWM

=head1 SYNOPSIS



=head1 DESCRIPTION

 This is to move things like wormpep from one release to the next for a species where the database has 
 not been rebuilt.

=head1 CONTACT

Gary ar2@sanger.ac.uk

=head1 APPENDIX

 The rest of the documentation details each of the object methods.
 Internal methods are preceded with a _

=cut
package CarryOver;

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Carp;
use File::Copy;
use File::Path;

sub new {  
	my $class = shift;
	my $self = {};
	bless $self, $class;
	
	$self->{'wormbase'} = shift;
	$self->{'log'} = shift;
	return $self;
}

# copy wormpep to new build version
sub _carry_over {
	my $self = shift;
	my $old_dir = shift;
	my $new_dir = shift;
	my $ver = shift;
	my $oldver = shift;
	
	mkdir("$new_dir",0775) unless( -e $new_dir);
	
	opendir (OLD, $old_dir) or $self->log->log_and_die("cant read $old_dir : $!\n");
	while (my $filename = readdir(OLD)) {
		next if( $filename =~ /^\./);
		my $newname = $filename;
		$newname =~ s/$oldver/$ver/;
		$self->{'log'}->write_to("\tcopying $newname\n");
		copy("$old_dir/$filename","$new_dir/$newname") or $self->log->log_and_die("cant copy $newname : $!\n");
	}
	closedir OLD;
}


sub carry_wormpep {
	my $self = shift;
	my $ver = shift; #version to create
	my $oldver = shift; #version to copy
	
	$ver = $self->{'wormbase'}->get_wormbase_version unless $ver;
	$ver =~ s/WS//;
	$oldver = $ver-1 unless $oldver;
		
	my $newwormpepdir = $self->{'wormbase'}->peproot."/".$self->{'wormbase'}->pepdir_prefix."pep$ver";
	my $oldwormpepdir = $self->{'wormbase'}->peproot."/".$self->{'wormbase'}->pepdir_prefix."pep$oldver";

	$self->_carry_over($oldwormpepdir, $newwormpepdir, $ver, $oldver);
}

sub carry_wormrna {
	my $self = shift;
	my $ver = shift; #version to create
	my $oldver = shift; #version to copy
	
	$ver = $self->{'wormbase'}->get_wormbase_version unless $ver;
	$ver =~ s/WS//;
	$oldver = $ver-1 unless $oldver;
		
	my $newdir = $self->{'wormbase'}->rnaroot."/".$self->{'wormbase'}->pepdir_prefix."rna$ver";
	my $olddir = $self->{'wormbase'}->rnaroot."/".$self->{'wormbase'}->pepdir_prefix."rna$oldver";

	$self->_carry_over($olddir,$newdir, $ver, $oldver);
}


#get/set routine for Log_files object
sub log {
	my $self = shift;
	my $log = shift;
	if( defined $log){
		$self->{'log'} = $log;
	}
	return $self->{'log'};
}

1;
























