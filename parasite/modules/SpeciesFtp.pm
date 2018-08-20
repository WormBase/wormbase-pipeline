package SpeciesFtp;

sub new {
  my ($class, $directory_root, $parasite_version) = @_;
  return bless {root => $directory_root, parasite_version => $parasite_version }, $class;
}

sub current_staging {
  return new(shift, join( "/", $ENV{PARASITE_SCRATCH},"dumps", "WBPS$ENV{PARASITE_VERSION}", "FTP"), $ENV{PARASITE_VERSION}); 
}

sub path_to {
  my ($self, $core_db, $extension, $zipped) = @_;
  $core_db = lc ($core_db);
  my ($spe, $cies, $bp ) = split "_", $core_db;
  $bp = uc($bp);
  $bp =~ s/ISS([0-9]+)PRJNA([0-9]+)/ISS$1_PRJNA$2/;
  $zipped //= 1;

  my $dir = join("/", $self->{root}, "${spe}_${cies}", $bp);

  if ($extension) {
     return join "/", $dir, join( ".",
       "${spe}_${cies}",$bp, "WBPS".$self->{parasite_version}, $extension, ($zipped ? "gz" : ())
     );
  } else {
    return $dir;
  } 
}
1;
