// ----------------------------------------------------------------------------
// 'RunMatcherComparator.C'
// Derek Anderson
// 01.30.2024
//
// root macro to run the STrackMatcherComparator module.
// ----------------------------------------------------------------------------

// analysis specific utilities
#include "MatcherComparatorOptions.h"
#include "/sphenix/user/danderson/install/include/strackmatchercomparator/STrackMatcherComparatorConfig.h"

// load libraries
R__LOAD_LIBRARY(/sphenix/user/danderson/install/lib/libstrackmatchercomparator.so)



// macro body -----------------------------------------------------------------

void RunMatcherComparator() {

  // get configruation
  STrackMatcherComparatorConfig config = MatcherComparatorOptions::GetConfig();

  // run module
  STrackMatcherComparator* comparator = new STrackMatcherComparator(config);
  comparator -> Init();
  comparator -> Analyze();
  comparator -> End();
  return;

}

// end ------------------------------------------------------------------------
