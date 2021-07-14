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
use Data::Dumper;

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

my $update_prod = $dba->dbc->prepare("UPDATE analysis_description 
			SET description = ?,
			display_label = ?,
			db_version = ?,
			is_current = ?,
			web_data_id = ?,
			displayable = ?
			WHERE logic_name = ?");

my $grab_logic_names = $dba->dbc->prepare("SELECT logic_name FROM analysis_description");


# first retrieve all existing logic names
my %logic_names;
$grab_logic_names->execute();
while (my $ln = $grab_logic_names->fetchrow_array()){
  $logic_names{$ln} = 1;
}
$grab_logic_names->finish();


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

	if (exists $logic_names{$logic_name}){
		print "$logic_name is already in the production database- updating associated data\n";
  		$update_prod->bind_param(1,$description);
		$update_prod->bind_param(2,$display_label);
		$update_prod->bind_param(3,$db_version);
		$update_prod->bind_param(4,$is_current);
		$update_prod->bind_param(5,$web_data_id);
		$update_prod->bind_param(6,$displayable);
		$update_prod->bind_param(7,$logic_name);
 		$update_prod->execute();
		$update_prod->finish();
	}

	else {
	  $add_to_prod->execute($logic_name, $description, $display_label, $db_version, $is_current, $web_data_id, $displayable);
	  $add_to_prod->finish;
	  print "Added $logic_name to $db\n";
	}
		
}


