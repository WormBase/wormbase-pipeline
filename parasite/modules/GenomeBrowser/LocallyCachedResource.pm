
package GenomeBrowser::LocallyCachedResource;
use LWP;
use YAML;
use File::Path qw(make_path);
use JSON;
use XML::Simple;
use Text::CSV qw(csv);
sub new {
    my ($class,$root_dir,$species, @other_args) = @_;
    my $dir = "$root_dir/$species";
    my $path_to_local_copy = "$dir/$class";
    make_path $dir;

    YAML::DumpFile($path_to_local_copy, 
       $class->_fetch($species, @other_args)
    ) unless -f $path_to_local_copy;

    return bless YAML::LoadFile($path_to_local_copy), $class;
}

sub get_text {
  my ($class,@urls) = @_;
  my $errors;
  for my $url (@urls){
   my $response = LWP::UserAgent->new->get($url);
   if($response->is_success){
     return $response->decoded_content;
   } else {
     $errors.="$url error:".$response->status_line."\n";
   }
  }
  die $errors;
}
sub get_csv {
  my $class = shift;
  my $text = $class->get_text(@_);
  return csv ({ allow_whitespace => 1, in=>\$text});
}
sub get_json { 
  my $class = shift;
  return from_json($class->get_text(@_));
}

sub get_xml { 
  my $class= shift;
  return XMLin($class->get_text(@_));
}
1;
