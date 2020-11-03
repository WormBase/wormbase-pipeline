#!/usr/bin/env perl

=pod

=head1 NAME

  populate_analysis_descriptions_from_file.pl

=head1 SYNOPSIS

  Put analysis descriptions into a production database

=head1 ARGUMENTS

  perl populate_analysis_descriptions_from_file.pl
         -host
         -port
         -user
         -pass
	 -dbname
         -help
         -descriptions_file

=cut


use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use Pod::Usage qw(pod2usage);

my ($host, $user, $db, $pass, $port, $descs_file, $descs_fh);

my $help = 0;

GetOptions(
	'host=s'			=> \$host,
	'dbname=s'			=> \$db,
	'pass=s'			=> \$pass,
	'port=i'			=> \$port,
	'user=s'			=> \$user,
	'descriptions_file=s'		=> \$descs_file,
	'help|?'			=> \$help,
) or pod2usage(-message => "use -help", -verbose => 1);

pod2usage(-verbose => 2) if $help;

open($descs_fh, '<' ,$descs_file) or die "Could not open $descs_file\n";

my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $host,
    -user   => $user,
    -dbname => $db,
    -pass   => $pass,
    -port   => $port,
      );

# prepare sql


my $add_to_prod = $dba->dbc->prepare("INSERT INTO analysis_description 
			(logic_name, description, display_label, db_version, is_current, web_data_id, displayable)
			VALUES (?,?,?,?,?,?,?)" );


# parse descriptions_file

while(<$descs_fh>){
	my @temp = split('\t',$_);
	my $logic_name 		= 	$temp[0];
	my $description 	= 	$temp[1];
	my $display_label 	=	$temp[2];
	my $db_version 		= 	$temp[3];
	my $is_current		=	$temp[4];
	my $web_data_id 	=	$temp[5];
	if ($web_data_id eq 'NULL'){
		$web_data_id = '103';	# NULL values not allowed - default to coding
	}
	my $displayable  	=	'1';
	$add_to_prod->execute($logic_name, $description, $display_label, $db_version, $is_current, $web_data_id, $displayable);
	$add_to_prod->finish;
	print "Added $logic_name to $db\n";
}


