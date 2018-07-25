
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
  my ($class,$url) = @_;
  my $response = LWP::UserAgent->new->get($url);
  die "$url error:".$response->status_line unless $response->is_success;
  return $response->decoded_content;
}
sub get_csv {
  my ($class,$url) = @_;
  my $text = $class->get_text($url);
  return csv ({ allow_whitespace => 1, in=>\$text});
}
sub get_json { 
  my ($class,$url) = @_;
  return from_json($class->get_text($url));
}

sub get_xml { 
  my ($class,$url) = @_;
  return XMLin($class->get_text($url));
}
1;
