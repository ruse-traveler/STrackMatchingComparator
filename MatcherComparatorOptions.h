// ----------------------------------------------------------------------------
// 'MatcherComparatorOptions.h'
// Derek Anderson
// 01.30.2024
//
// Options for the STrackMatcherComparator module.
// ----------------------------------------------------------------------------

#ifndef MATCHERCOMPARATOROPTIONS_H
#define MATCHERCOMPARATOROPTIONS_H 

// analysis utilities
#include "/sphenix/user/danderson/install/include/strackmatchercomparator/STrackMatcherComparatorConfig.h"

// make common namespaces implicit
using namespace std;



namespace MatcherComparatorOptions {

  // set up configuration -----------------------------------------------------

  STrackMatcherComparatorConfig GetConfig() {

    STrackMatcherComparatorConfig cfg {
      /* configuring goes here */
    };
    return cfg;

  }  // end 'GetConfig()'

}  // end MatcherComparatorOptions namespace

#endif

// end ------------------------------------------------------------------------
