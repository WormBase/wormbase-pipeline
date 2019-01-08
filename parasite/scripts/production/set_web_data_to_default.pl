#!/usr/bin/env perl
# Copy the tables from previous db to start with a clean slate. Then
# upload an ID mapping, and add kill and creation events for everything else.
# We only add events for genes, while Ensembl also tracks individual exons and transcripts.
use strict;
use warnings;
use ProductionMysql;
use Getopt::Long;
my $db_command;
my %args;
GetOptions (
  'core_db=s'             => \$args{core_db}, # can be a pattern
  'db_command=s'          => \$db_command,
);
my $mysql = ProductionMysql->new($db_command);
$args{core_db} = $mysql->core_db($args{core_db});

my $get_sql = "mysql-pan-prod ensembl_production_parasite -e 'SELECT ad.logic_name, wd.data FROM analysis_description ad LEFT OUTER JOIN web_data wd ON ad.default_web_data_id = wd.web_data_id WHERE ad.is_current = 1 and wd.data is not NULL' ";
my %web_data_for_analyses;
open (my $fh,  "$get_sql | ");
while(<$fh>){
   chomp;
   my ($logic_name, $data_string) = split "\t";
   $web_data_for_analyses{$logic_name} = $data_string;
}

my $insert_web_data_dbh = $mysql->dbc($args{core_db})->prepare("update analysis_description join analysis on (analysis_description.analysis_id = analysis.analysis_id) set web_data=? where logic_name = ?");
my $transcript_adaptor = $mysql->adaptor($args{core_db}, 'Analysis');

my %no_web_data_ok = map {$_ => 1 } qw/dust goa_import interpro2go interpro2pathway ncoils repeatmask repeatmask_customlib seg sfld tmhmm trf xrefchecksum xrefexoneratedna xrefexonerateprotein/;

for my $analysis (grep { not $_->web_data } @{$transcript_adaptor->fetch_all}){
  my $ln = $analysis->logic_name;
  if ($web_data_for_analyses{$ln}){
    $insert_web_data_dbh->execute($web_data_for_analyses{$ln}, $ln);
  } else {
     warn "No web data for: $args{core_db} $ln\n" unless $no_web_data_ok{$ln};
  }
}
