
package GenomeBrowser::Descriptions;

my $curated_descriptions_canned = {
  schistosoma_mansoni => {
    ERP123000 => {
       ERR001000 => ["short desc", "full desc"],
       ERR001001 => "desc" 
    }
  }
};

sub new {
  my ($class, $species, $curated_descriptions) = @_;
  return bless {
    species => $species,
    curated => $curated_descriptions || $curated_descriptions_canned->{$species},
  }, $class;
}

sub run_description_from_sample_name {
  my ($species, $sample_name) = @_;
  (my $cies = $species) =~ s/.*_//;
  return "" unless $sample_name;
  return "" if ($cies and $sample_name =~ /$cies/i and $sample_name =~ /sample from/i);
  return "" if (scalar(split /\W+/, $sample_name) == 1);
  return "" if $sample_name =~ /private/;
  return "" if length ($sample_name) < 10;
  return $sample_name;
}
# if some of these types can be made nicer in labels,
# e.g. if they have _, add some code that prettifies them
my @types_that_help_explain_values=qw/strain/;
sub run_description_from_factors {
  my ($factors, $attributes) = @_;
  my @result_short;
  my @result_full;

  for my $t (@$factors){
     (my $v = $attributes->{$t} ) or next;
     push @result_short, $v;
     if (grep {$_ eq $t} @types_that_help_explain_values) {
        push @result_full, "$t $v";
     } else {
        push @result_full, $v;
     }
  }
  return "" unless @result_short and @result_full;
  return [join (", ", @result_short), join (", ", @result_full)];
}
sub _get_run_description {
   my ($self, $study_id, $run_id, $factors, $attributes) = @_;
   return (
      $self->{curated}{$study_id}{$run_id}
      or run_description_from_sample_name($self->{species}, $attributes->{sample_name})
      or run_description_from_factors($factors, $attributes)
      or ""
   );
}
sub run_description {
   my ($self, $study_id, $run_id, $factors, $attributes) = @_;
   my $r = _get_run_description(@_);
   return map {"$run_id: $_"} @$r if ref $r eq 'ARRAY';
   return "$run_id: $r", "$run_id: $r", if $r;
   my $species = $self->{species};
   $species =~ s/_/ /g;
   $species = ucfirst($species);
   return "$run_id", "$run_id: sample from $species";
}
# Do we also want to curate these?
sub study_description {
  my ($self, $study_id, $study_attributes) = @_;
  return ( 
    $study_attributes->{"Study description"} 
    or $study_attributes->{study} 
    or $study_id
    or "$species study"
  );
}
1;
