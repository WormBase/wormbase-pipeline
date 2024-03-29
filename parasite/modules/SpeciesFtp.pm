package SpeciesFtp;

sub new {
  my ($class, $directory_root, $parasite_version) = @_;
  return bless {root => $directory_root, parasite_version => $parasite_version }, $class;
}

sub current_staging {
  return new(shift, join( "/", $ENV{PARASITE_SCRATCH},"dumps", "WBPS$ENV{PARASITE_VERSION}", "FTP"), $ENV{PARASITE_VERSION});
}
sub previous_staging {
  return new(shift, join( "/", $ENV{PARASITE_SCRATCH},"dumps", "WBPS$ENV{PREVIOUS_PARASITE_VERSION}", "FTP"), $ENV{PREVIOUS_PARASITE_VERSION}); 
}
sub release {
  my ($class, $parasite_version, $release_folder) = @_;
  $release_folder //= "WBPS$parasite_version";
  return new($class, join( "/", "/nfs/ftp/public/databases/wormbase/parasite/releases", $release_folder, "species"), $parasite_version);
}
sub dot_next {
  return release(shift, $ENV{PARASITE_VERSION}, ".next");
}
sub current_release {
  return release(shift, $ENV{PARASITE_VERSION});
}
sub path_to {
  my ($self, $core_db, $extension, $zipped) = @_;
  $core_db = lc ($core_db);
  my ($spe, $cies, $bp ) = split "_", $core_db;
  $bp = uc($bp);
  # fix for special bioproject names assigned to old T. pseudospiralis genome versions
  $bp =~ s/ISS([0-9]+)PRJNA([0-9]+)/ISS$1_PRJNA$2/;
  # fix for special bioproject names assigned to old S. carpocapsae genome versions
  $bp =~ s/V([0-9]+)PRJNA([0-9]+)/V$1_PRJNA$2/;
  # fix for special bioproject names assigned to globodera rhostochiensis genome versions l22 and l19
  $bp =~ s/L([0-9]+)PRJNA([0-9]+)/L$1_PRJNA$2/;
  # fix for special bioproject names assigned to different trematodes
  $bp =~ s/TD([0-9]+)PRJEB([0-9]+)/TD$1_PRJEB$2/;
  # fix for special bioproject names assigned to schmidtea
  $bp =~ s/S2F19H([0-9]+)PRJNA([0-9]+)/S2F19H$1_PRJNA$2/;
  # fix for special bioproject names assigned to aphelenchoides
  $bp =~ s/A([A-Z]+)PRJNA(834627)/A$1_PRJNA$2/;
  # fix for special bioproject names assigned to briggsae
  $bp =~ s/QX(1410)PRJNA(784955)/QX$1_PRJNA$2/;
  $bp =~ s/VX(34)PRJNA(784955)/VX$1_PRJNA$2/;
  # fix for special bioproject names assigned to m chitwoodi
  $bp =~ s/(RACE1)PRJNA(666745)/$1_PRJNA$2/;
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
sub root_exists {
  my ($self) = @_;
  return -d $self->{root};
}
1;
