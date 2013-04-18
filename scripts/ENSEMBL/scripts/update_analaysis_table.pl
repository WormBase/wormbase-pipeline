#!/usr/bin/env perl
# parse an ENSEMBL analysis config file and wite SQL to update a core database
# tested on e70, if the schema changes the script might also need to be changed

my $file = shift;

my $c = parse_file($file);

while(my($k,$v)=each %$c){
   my @values;
   my @columns;
   my @pairs;
   while (my($k2,$v2)=each %$v){
      $v2="'$v2'" unless /^\d+$/;
      next if $k2 eq 'input_id_type';
      push (@values,$v2);
      push (@columns,$k2);
      push (@pairs,"$k2=$v2")
   }
   # that is for the analysis table
   print 'INSERT INTO analysis (',join(',',('logic_name',@columns)),') VALUES (',join(',',("\"$k\"",@values)),")";
   print ' ON DUPLICATE KEY UPDATE ',join(',',("logic_name=\"$k\"",@pairs)),";\n"; # there is an unique key on logic_name

   # input_id_type_analysis
   print "DELETE FROM input_id_type_analysis WHERE analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name=\"$k\");\n";
   print "INSERT INTO input_id_type_analysis SET analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name=\"$k\"), input_id_type=\"${\$c->{$k}->{input_id_type}}\";\n";
}

# returns a hash with the config data
sub parse_file{
 my ($file)=@_;
 my ($config,$header);
 
 open (FILE, $file) or die("Couldn't open file $file\n");
 while (<FILE>) {

      next if /^\#/;# Comment or blank line
      chomp();

      # [HEADER]
      if (/^\[(.*)\]\s*$/) {
        $header = $1;
      }elsif (/^([^=\s]+)\s*=\s*(.+?)\s*$/) {
        my $key = $1;
        my $value = $2;
       
        # store them in the config hash
        $config->{$header}->{$key} = $value;
      }
  }
  close FILE;
  return $config;
}
