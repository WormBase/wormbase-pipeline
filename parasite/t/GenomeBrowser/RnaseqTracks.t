use GenomeBrowser::RnaseqTracks;
use Test::More;

is_deeply(GenomeBrowser::RnaseqTracks::_label("run id", "", ()),
  "run id", "no attributes means run id"
);

is_deeply(GenomeBrowser::RnaseqTracks::_label("run id", "", ("factor value", "other value")),
  "run id: factor value, other value", "concatenate factor values"
);
is_deeply(GenomeBrowser::RnaseqTracks::_label("run id", "ERS067030", ()),
  "run id",
  "Ignore crappy titles"
);
is_deeply(GenomeBrowser::RnaseqTracks::_label("run id", "very descriptive title", ()),
  "run id: very descriptive title",
  "use study names when they're good"
);

done_testing();
