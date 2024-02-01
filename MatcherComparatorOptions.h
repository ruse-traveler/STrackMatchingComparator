// ----------------------------------------------------------------------------
// 'MatcherComparatorOptions.h'
// Derek Anderson
// 01.30.2024
//
// Options for the STrackMatcherComparator module.
// ----------------------------------------------------------------------------

#ifndef MATCHERCOMPARATOROPTIONS_H
#define MATCHERCOMPARATOROPTIONS_H 

// c++ utilities
#include <string>
#include <vector>
// analysis utilities
#include "/sphenix/user/danderson/install/include/strackmatchercomparator/STrackMatcherComparatorConfig.h"

// make common namespaces implicit
using namespace std;



namespace MatcherComparatorOptions {

  // set up configuration -----------------------------------------------------

  STrackMatcherComparatorConfig GetConfig() {

    const vector<string> cfg_info = {
      "#bf{#it{sPHENIX}} Simulation",
      "1 #pi^{-}/event, p_{T} = 10 GeV/c",
      "#bf{Only #pi^{-}}"
    };

    STrackMatcherComparatorConfig cfg {
      .outFileName     = "test.hists.root",
      .newInFileName   = "input/merged/sPhenixG4_oneMatchPerParticle_forCrossCheck_newMatcher.pt10num1evt500pim.d25m1y2024.root",
      .oldInFileName   = "input/merged/sPhenixG4_oneMatchPerParticle_oldEval_forCrossCheck.pt10num1evt500pim.d25m1y2024.root",
      .oddPtFrac       = {0.5, 1.50},
      .info            = cfg_info,
      .doIntNorm       = false,
      .matchVertScales = true
    };
    return cfg;

  }  // end 'GetConfig()'

}  // end MatcherComparatorOptions namespace

#endif

// end ------------------------------------------------------------------------
