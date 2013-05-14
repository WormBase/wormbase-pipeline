#!/usr/bin/env perl

=head1 SYNOPSIS

Script for writing an analysis table given a configuration file, or the
other way around.

=head1 DESCRIPTION

This script will take a configuration file and write to an analysis table
or take an analysis table and write a configuration file.

=head1 OPTIONS

Database options:

  -dbhost    Host name for database.
  -dbport    Port to connect to (optional).
  -dbname    Database to connect to.
  -dbuser    Username to connect as.
  -dbpass    Password to use.

Other options:

  -read     Read a configuration file and write to the database.
  -write    Read the rule table and write a configuration file.
  -file     File to read from or write too.
  -forceclear By default, fields not specified for existing analyses will retain
              their current values; this option forces them to be cleared

  -help     Displays help text.

This is what a conf entry should like like.  The header for each section
should be the analysis logic name, the rest is key-value pairs:

  [RepeatMask]
  db=repbase
  db_version=020713
  db_file=repbase
  program=RepeatMasker
  program_version=1
  program_file=RepeatMasker
  parameters=-low, -lib, /path/to/file/briggsae.lib
  module=RepeatMasker
  module_version=1
  gff_source=RepeatMasker
  gff_feature=Repeat
  input_id_type=CONTIG

=cut

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);


use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Analysis;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

sub usage {
  exec( 'perldoc', $0 );
  exit;
}

my ($dbhost, $dbuser, $dbpass, $dbport, $dbname);
my ($read, $write, $file, $force_clear, $is_pipeline_db, $help, $do_not_write_new_analyses);

if ( !GetOptions( 'host|dbhost|h:s'  => \$dbhost,
                  'dbname|db|D:s'    => \$dbname,
                  'user|dbuser|u:s'  => \$dbuser,
                  'pass|dbpass|p:s'  => \$dbpass,
                  'port|dbport|P:s'  => \$dbport,
                  'read!'            => \$read,
                  'write!'           => \$write,
                  'forceclear!'      => \$force_clear,
                  'file=s'           => \$file,
                  'nonew'            => \$do_not_write_new_analyses,
                  'help!'            => \$help, ) || $help) {
  usage();
}

if ( !($dbhost) || !($dbuser) || !($dbname) ) {
  print STDERR
    "need to pass in database arguments for script to work\n";
  print STDERR
    "-dbhost $dbhost -dbuser $dbuser -dbpass $dbpass -dbname" .
    " $dbname -dbport $dbport\n";
  usage();
}

if ( ( !($read) && !($write) ) || ( $read && $write ) ) {
  print STDERR "you need to define either read or write on the " .
    "commandline but you shouldn't define both\n";
  usage();
}

if ( !$file ) {
  print STDERR
    "You need to pass a file name which either represents a " .
    "analysis config file to be read and stored in the database or " .
    "written based on the analysis objects in the database\n";
  usage();
}


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dbhost,
                                             -user   => $dbuser,
                                             -pass   => $dbpass,
                                             -dbname => $dbname,
                                             -port   => $dbport, );
if (&check_is_pipeline_db($db)) {
  $db = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor( -host   => $dbhost,
                                                      -user   => $dbuser,
                                                      -pass   => $dbpass,
                                                      -dbname => $dbname,
                                                      -port   => $dbport, );
  $is_pipeline_db = 1;
}

if ($read) {
  my @analyses = @{ parse_files($file) };
  write_into_db( $db, \@analyses);
}
else {
  my $analyses = read_db($db);
  $analyses = [ sort { $a->dbID() <=> $b->dbID() } @{$analyses} ];
  write_file( $file, $analyses );
}



#####################################

=head2 parse_files

  Arg [1]   : array of filenames
  Function  : parse a analysis config file and produce analysis objects
  Returntype: an array ref for an array of analyses objects
  Exceptions: if file doesn't exist
              if config format is in correct
              if key already exists for a particular header'
  Caller    : 
  Example   : my @analyses = @{&parse_files($file)};

=cut

sub parse_files {
  
  my @files = shift;
  
  my %headers;     # will store names of headers and number of keys for each
  
  my $hcounter = 0;
  my %horder; # stores order in which entries were read
  
  
  my $config = {};
  # read each file
  
  foreach my $file (@files) {
    
    if (! -e $file) {
      throw("analysis file $file not found\n");
    }
    my $header = "";
    
    open(FILE, $file) or throw("Could not open file $file");
    while (<FILE>) {
      chomp();
      
      # Comment or blank line
      next if (/^\s$/ || /^\#/);

      # [HEADER]
      if (/^\[(.*)\]\s*$/) {         # $1 will be the header name, without the [] 
	$header = $1;
	$headers{$header} = 0;
        $horder{$header} = $hcounter++;
	#print "Reading stanza $header\n";
      } 

      # key=value
      if (/^([^=\s]+)\s*=\s*(.+?)\s*$/) {   # $1 = key, $2 = value

	my $key = lc($1);           # keys stored as all lowercase, values have case preserved
	my $value = $2;
	if (length($header) == 0) {
	  throw("Found key/value pair $key/$value outside stanza");
	}
	#print "Key: $key Value: $value\n"; 
      	
	# Check if this header/key is already defined
	if (exists($config->{$header}->{$key})) {
	  throw("$key is already defined for [$header]; cannot be redefined");
	} else {
	  # store them in the config hash
	  $config->{$header}->{$key} = $value;
	  #print "$header:$key=$value\n";
	  $headers{$header}++;  # will be used to check for headers with no keys
	}

      }

    } # while <FILE>

    close FILE;
  }

  my @analyses;
  # add a blank key/value for any headers that have no keys
  foreach my $h (sort { $horder{$a} <=> $horder{$b} }  keys (%headers)) {

    my $analysis = Bio::EnsEMBL::Pipeline::Analysis->new(-logic_name => $h);
    my $ext_ana = $db->get_AnalysisAdaptor->fetch_by_logic_name($h);

    foreach my $field ('db_file',
                       'db',
                       'db_version',
                       'program',
                       'program_version',
                       'program_file',
                       'gff_source',
                       'gff_feature',
                       'module',
                       'module_version',
                       'parameters',
                       'input_id_type') {

      if (exists $config->{$h}->{$field}) {
        if (lc($config->{$h}->{$field}) ne 'null') {
          $analysis->$field($config->{$h}->{$field});
        }
      } elsif (defined $ext_ana and not $force_clear) {
        $analysis->$field( $ext_ana->$field );
      }
    }

    push(@analyses, $analysis);
  }

  return \@analyses;

}


=head2 write_into_db

  Arg [1]   : Bio::EnsEMBL::DBSQL::DBAdaptor
      [2]   : Ref to an array of analysis objects
  object already exists
  Function  : Write the analysis objects into the database
  Returntype: N/A
  Exceptions: if dbadaptor is the wrong type of object
  Caller    : 
  Example   : &write_into_db($db, \@analyses);

=cut

sub write_into_db {
  my $db = shift;
  my $analyses = shift;

  #print "have analysis adaptor ".$analysis_adaptor."\n";
  if(!($db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'))){
    throw("need a Pipeline::DBAdaptor not ".$db);
  }

  my $analysis_adaptor = $db->get_AnalysisAdaptor;
  my $sql = "select analysis_id from analysis where logic_name = ?";
  my $sth = $db->prepare($sql);
  ANALYSIS:foreach my $a(@$analyses){ 
    $sth->execute($a->logic_name);
    my ($analysis_id)= $sth->fetchrow;

    if($analysis_id){
      #
      # update of an existing analysis
      #

      if($a->input_id_type){
        if ($is_pipeline_db) {
          my $clean_sql = "DELETE from input_id_type_analysis WHERE analysis_id = $analysis_id";
          $db->dbc->do($clean_sql);
          my $stored_sql = "INSERT into input_id_type_analysis ".
              "(analysis_id, input_id_type) values(?, ?)";
          my $stored_sth = $db->dbc->prepare($stored_sql);
          $stored_sth->execute($analysis_id, $a->input_id_type);
          $stored_sth->finish;
        }
      }
      $a->dbID($analysis_id);
      $a->adaptor($analysis_adaptor);

      my $created_sql = "SELECT NOW();";
      my $created_sth = $db->prepare($created_sql);
      $created_sth->execute();
      my ($current_time) = $created_sth->fetchrow;
      $a->created($current_time);
      $analysis_adaptor->update($a);
      next ANALYSIS;
    }else {
      if (not $do_not_write_new_analyses) {
        $analysis_adaptor->store($a);
      } else {
        warn("Analysis " . $a->logic_name . " not found in database - not writing\n");
      }
    }
  }
}


=head2 read_db

  Arg [1]   : Bio::EnsEMBL::DBAdaptor
  Function  : Read the analysis objects from the database
  Returntype: array ref of analysis objects
  Exceptions: if db isn't the correct type'
  Caller    : 
  Example   : my $analyses = &read_db($db);

=cut

sub read_db{
  my $db = shift;

  if(!($db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'))){
    throw("need a DBAdaptor not ".$db);
  }
 
  my $analysis_adaptor = $db->get_AnalysisAdaptor;

  return $analysis_adaptor->fetch_all;
}


=head2 write_file

  Arg [1]   : filename
      [2]   : arrayref of analysis objects
  Function  : write a config file for the objects given
  Returntype: N/A
  Exceptions: if file doesnt exist
  Caller    : 
  Example   : &write_file($file, $analyses);

=cut

sub write_file{
  my $file = shift;
  my $analyses = shift;
  #print STDERR "opening ".$file."\n";
  open (FH, '>'.$file) or throw ("couldn't open $file to write to");
  foreach my $a(@$analyses){
    print FH "[".$a->logic_name."]\n";
    print FH "db=".$a->db."\n" if($a->db);
    print FH "db_version=".$a->db_version."\n" if($a->db_version);
    print FH "db_file=".$a->db_file."\n" if($a->db_file);
    print FH "program=".$a->program."\n" if($a->program);
    print FH "program_version=".$a->program_version."\n" if($a->program_version);
    print FH "program_file=".$a->program_file."\n" if($a->program_file);
    print FH "parameters=".$a->parameters."\n" if($a->parameters);
    print FH "module=".$a->module."\n" if($a->module);
    print FH "gff_source=".$a->gff_source."\n" if($a->gff_source);
    print FH "gff_feature=".$a->gff_feature."\n" if($a->gff_feature);
    if($a->can("input_id_type")){
      print FH "input_id_type=".$a->input_id_type."\n" if($a->input_id_type);   
    }
    print FH "description ".$a->description."\n" if($a->description);
    print FH "display_name".$a->display_label."\n" if($a->display_label);
    print FH "\n\n";
  }
  
}


############################
sub check_is_pipeline_db {
  my $dbh = shift;

  eval {
    $dbh->dbc->do("SELECT 1 from input_id_type_analysis LIMIT 1");
  };
  if ($@) {
    return 0;
  } else {
    return 1;
  }
}


1;
