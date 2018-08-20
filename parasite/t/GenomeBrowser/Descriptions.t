use GenomeBrowser::Descriptions;
use Test::More;

my $species = "schistosoma_mansoni";
my $study_id = "SRP1234000";
my $run_id = "SRR1230321";
my $default_long_desc = "$run_id: sample from Schistosoma mansoni";
sub assert_run_description {
  my ($factors, $attributes, $curated, $expected_description_short, $expected_description_full, $desc) = @_;
   my ($description_short, $description_full) = GenomeBrowser::Descriptions->new(
        $species, {$study_id=>$curated})->run_description($study_id, $run_id, $factors, $attributes);
   is_deeply($description_short, $expected_description_short, "description_short $desc");
   is_deeply($description_full, $expected_description_full, "description_full $desc");
}

assert_run_description([], {},{}, $run_id, $default_long_desc, "Default long description ok in null case");
assert_run_description([], {},{$run_id => "curated value"}, "$run_id: curated value", "$run_id: curated value", "curated value");
assert_run_description([], {},{$run_id => ["curated value", "curated value long"]}, "$run_id: curated value", "$run_id: curated value long", "curated value long");
sub sample_name_ok {
  my ($sample_name, $desc) = @_;
  assert_run_description([], {sample_name => $sample_name}, {}, "$run_id: $sample_name", "$run_id: $sample_name", 
    $desc // "sample_name_ok $sample_name");
}
sub sample_name_rejected {
  my ($sample_name, $desc) = @_;
  assert_run_description([], {sample_name => $sample_name}, {}, $run_id, $default_long_desc,
    $desc // "sample_name_rejected $sample_name");
}
sample_name_ok("words words words");
sample_name_rejected("");
sample_name_rejected("123");
sample_name_rejected("short SN");
sample_name_rejected("onlyOneWordIsBadSampleName");
sample_name_rejected("something about being private in GEO");
sample_name_rejected("sample from S. mansoni");

sub _str {
  my ($factors, $attributes) = @_;
  return (join "/ ", @$factors, (join ", ", keys $attributes), (join ", ", values $attributes ));
}
sub factors_ok{
  my ($factors, $attributes, $d_short, $d_full, $desc) = @_;
  assert_run_description($factors, $attributes, {}, "$run_id: $d_short", "$run_id: $d_full",
  $desc // "factors_ok "._str($factors, $attributes) 
  );
}
sub factors_rejected {
  my ($factors, $attributes, $desc) = @_;
  assert_run_description($factors, $attributes, {}, $run_id, $default_long_desc,
  $desc // "factors_ok "._str($factors, $attributes)
  );
}

factors_ok(["type"], {"type"=>"value"}, "value", "value");
factors_ok(["type", "strain"], {"type"=>"value"}, "value", "value");
factors_ok(["type", "other_type"], {"type"=>"value"}, "value", "value");
factors_ok(["type", "strain"], {"type"=>"value", "strain"=>"LE"}, "value, LE", "value, strain LE", "factors_ok full description has strain");
factors_ok(["type", "other_type"], {"type"=>"value", "other_type"=>"other_value"}, "value, other_value", "value, other_value");
factors_rejected([], {});
factors_rejected(["type"], {});
factors_rejected(["type"], {"other_type"=> "value"});

sub assert_study_description {
  my ($study_attributes, $expected, $desc) = @_;
   $desc //= "assert_study_description $study_id -> $expected";
   is_deeply(GenomeBrowser::Descriptions->new(
        $species, {})->study_description($study_id, $study_attributes),
     $expected, $desc);
}

assert_study_description({},$study_id);
assert_study_description({"Study description" => "Desc", "study" => "Study field"},"Desc");
assert_study_description({"study" => "Study field"},"Study field");

done_testing();
