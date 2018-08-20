use GenomeBrowser::Resources;
use Test::More;

is_deeply(GenomeBrowser::Resources::_run_description("run id", "", ()),
  "run id", "no attributes means run id"
);

is_deeply(GenomeBrowser::Resources::_run_description("run id", "", ("factor value", "other value")),
  "run id: factor value, other value", "concatenate factor values"
);
is_deeply(GenomeBrowser::Resources::_run_description("run id", "ERS067030", ()),
  "run id",
  "Ignore crappy titles"
);
is_deeply(GenomeBrowser::Resources::_run_description("run id", "very descriptive title", ()),
  "run id: very descriptive title",
  "use study names when they're good"
);

done_testing();
