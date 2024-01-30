// ----------------------------------------------------------------------------
// 'STrackMatcherComparator.cc'
// Derek Anderson
// 01.30.2024
//
// Small module to produce plots from the output of the
// SvtxEvaluator and FillMatchedClusTree modules.
// ----------------------------------------------------------------------------

// c++ utilities
#include <array>
#include <string>
#include <vector>
#include <utility>
#include <iostream>
// root libraries
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <TFile.h>
#include <TLine.h>
#include <TError.h>
#include <TString.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TDirectory.h>
// analysis utilities
#include "STrackMatcherComparator.h"
#include "STrackMatcherComparatorConfig.h"

// make common namespaces implicit
using namespace std;



// ctor/dtor ------------------------------------------------------------------

STrackMatcherComparator::STrackMatcherComparator(optional<STrackMatcherComparatorConfig> config = nullopt) {

  if (config.has_value()) {
    m_config = config;
  }

}  // end ctor



STrackMatcherComparator::~STrackMatcherComparator() {

  /* nothing to do */

}  // end dtor



// public methods -------------------------------------------------------------

void STrackMatcherComparator::Init() {

  // announce start
  cout << "\n  Beginning track matcher comparison...\n"
       << "    Initialization:"
       << endl;

  // run initialization routines
  OpenOutput();
  OpenInput();
  InitHists();
  return;

}  // end 'Init()'



void STrackMatcherComparator::Analyze() {

  // announce analyzing
  cout << "    Analysis:" << endl;

  // run analysis routines
  GetNewTreeHists();
  GetNewTupleHists();
  GetOldTupleHists();
  return;

}  // end 'Analysis()'



int STrackMatcherComparator::End() {

  // announce end
  cout << "    End:" << endl;

  // run ending routines
  MakeRatiosAndPlots(m_vecTreeHists1D,  m_vecTreeHists2D,  Src::NewTree,  "VsNewTree");
  MakeRatiosAndPlots(m_vecTupleHists1D, m_vecTupleHists2D, Src::NewTuple, "VsNewTuple");
  SaveHistograms();
  CloseInput();
  CloseOutput();

  // exit
  cout << "  Finished track matcher comparison!\n" << endl;
  return;
}



// internal methods -----------------------------------------------------------

void STrackMatcherComparator::OpenOutput() {

  return;

}  // end 'OpenOutput()'



void STrackMatcherComparator::OpenInput() {

  return;

}  // end 'OpenInput()'



void STrackMatcherComparator::InitHists() {

  // output histogram base names
  const vector<vector<string>> vecHistBase = {
    {"hTruthNumTot",  "hTrackNumTot",  "hWeirdNumTot",  "hNormNumTot"},
    {"hTruthNumIntt", "hTrackNumIntt", "hWeirdNumIntt", "hNormNumIntt"},
    {"hTruthNumMvtx", "hTrackNumMvtx", "hWeirdNumMvtx", "hNormNumMvtx"},
    {"hTruthNumTpc",  "hTrackNumTpc",  "hWeirdNumTpc",  "hNormNumTpc"},
    {"hTruthRatTot",  "hTrackRatTot",  "hWeirdRatTot",  "hNormRatTot"},
    {"hTruthRatIntt", "hTrackRatIntt", "hWeirdRatIntt", "hNormRatIntt"},
    {"hTruthRatMvtx", "hTrackRatMvtx", "hWeirdRatMvtx", "hNormRatMvtx"},
    {"hTruthRatTpc",  "hTrackRatTpc",  "hWeirdRatTpc",  "hNormRatTpc"},
    {"hTruthPhi",     "hTrackPhi",     "hWeirdPhi",     "hNormPhi"},
    {"hTruthEta",     "hTrackEta",     "hWeirdEta",     "hNormEta"},
    {"hTruthPt",      "hTrackPt",      "hWeirdPt",      "hNormPt"},
    {"hTruthFrac",    "hTrackFrac",    "hWeirdFrac",    "hNormFrac"},
    {"hTruthQual",    "hTrackQual",    "hWeirdQual",    "hNormQual"},
    {"hTruthPtErr",   "hTrackPtErr",   "hWeirdPtErr",   "hNormPtErr"},
    {"hTruthEtaErr",  "hTrackEtaErr",  "hWeirdEtaErr",  "hNormEtaErr"},
    {"hTruthPhiErr",  "hTrackPhiErr",  "hWeirdPhiErr",  "hNormPhiErr"},
    {"hTruthPtRes",   "hTrackPtRes",   "hWeirdPtRes",   "hNormPtRes"},
    {"hTruthEtaRes",  "hTrackEtaRes",  "hWeirdEtaRes",  "hNormEtaRes"},
    {"hTruthPhiRes",  "hTrackPhiRes",  "hWeirdPhiRes",  "hNormPhiRes"}
  };
  const size_t nHistRows  = vecHistBase.size();
  const size_t nHistTypes = vecHistBase[0].size();

  // 2D histogram name modifiers
  const vector<string> vecVsHistModifiers = {
    "VsTruthPt",
    "VsNumTpc"
  };
  const size_t nVsHistMods = vecVsHistModifiers.size();

  // axis titles
  const string         m_config.count("counts");
  const vector<string> vecBaseAxisVars = {
    "N^{tot} = N_{hit}^{mvtx} + N_{hit}^{intt} + N_{clust}^{tpc}",
    "N_{hit}^{intt}",
    "N_{hit}^{mvtx}",
    "N_{clust}^{tpc}",
    "N_{reco}^{tot} / N_{true}^{tot}",
    "N_{reco}^{intt} / N_{true}^{intt}",
    "N_{reco}^{mvtx} / N_{true}^{mvtx}",
    "N_{reco}^{tpc} / N_{true}^{tpc}",
    "#varphi",
    "#eta",
    "p_{T} [GeV/c]",
    "p_{T}^{reco} / p_{T}^{true}"
    "#chi^{2} / ndf",
    "#deltap_{T} / p_{T}^{reco}",
    "#delta#eta / #eta^{reco}",
    "#delta#varphi / #varphi^{reco}"
    "#Deltap_{T} / p_{T}^{true}",
    "#Delta#eta / #eta^{true}",
    "#Delta#varphi / #varphi^{true}"
  };
  const vector<string> vecVsAxisVars = {
    "p_{T}^{true} [GeV/c]",
    "N_{clust}^{tpc}"
  };

  // output histogram no. of bins
  const uint32_t nNumBins  = 101;
  const uint32_t nRatBins  = 120;
  const uint32_t nEtaBins  = 80;
  const uint32_t nPhiBins  = 360;
  const uint32_t nPtBins   = 101;
  const uint32_t nFracBins = 220;
  const uint32_t nQualBins = 210;
  const uint32_t nResBins  = 110;

  // output histogram bin ranges
  const pair<float, float> xNumBins  = {-0.5,  100.5};
  const pair<float, float> xRatBins  = {-0.5,  5.5};
  const pair<float, float> xEtaBins  = {-2.,   2.};
  const pair<float, float> xPhiBins  = {-3.15, 3.15};
  const pair<float, float> xPtBins   = {-1.,   100.};
  const pair<float, float> xFracBins = {-0.5,  10.5};
  const pair<float, float> xQualBins = {-0.5,  20.5};
  const pair<float, float> xResBins  = {-5.5,  5.5};

  // output histogram base binning definitions
  vector<tuple<uint32_t, pair<float, float>>> vecBaseHistBins = {
    make_tuple(nNumBins,  xNumBins),
    make_tuple(nNumBins,  xNumBins),
    make_tuple(nNumBins,  xNumBins),
    make_tuple(nNumBins,  xNumBins),
    make_tuple(nRatBins,  xRatBins),
    make_tuple(nRatBins,  xRatBins),
    make_tuple(nRatBins,  xRatBins),
    make_tuple(nRatBins,  xRatBins),
    make_tuple(nPhiBins,  xPhiBins),
    make_tuple(nEtaBins,  xEtaBins),
    make_tuple(nPtBins,   xPtBins),
    make_tuple(nFracBins, xFracBins),
    make_tuple(nQualBins, xQualBins),
    make_tuple(nResBins,  xResBins),
    make_tuple(nResBins,  xResBins),
    make_tuple(nResBins,  xResBins),
    make_tuple(nResBins,  xResBins),
    make_tuple(nResBins,  xResBins),
    make_tuple(nResBins,  xResBins)
  };

  // output 2D histogram x-axis binning
  vector<tuple<uint32_t, pair<float, float>>> vecVsHistBins = {
    make_tuple(nPtBins,  xPtBins),
    make_tuple(nNumBins, xNumBins)
  };

  // make 1D base histograms
  vector<vector<TH1D*>> vecBaseHist1D(nHistRows);
  for (size_t iHistRow = 0; iHistRow < nHistRows; iHistRow++) {
    for (const string sHistBase : vecBaseHistBase[iHistRow]) {
      vecBaseHist1D[iHistRow].push_back(
        new TH1D(
          sHistBase.data(),
          "",
          get<0>(vecBaseHistBins[iHistRow]),
          get<1>(vecBaseHistBins[iHistRow]).first,
          get<1>(vecBaseHistBins[iHistRow]).second
        )
      );
    }  // end type loop
  }  // end row loop

  // make 2D base histograms
  vector<vector<vector<TH2D*>>> vecBaseHist2D(nHistRows, vector<vector<TH2D*>>(nHistTypes));
  for (size_t iHistRow = 0; iHistRow < nHistRows; iHistRow++) {
    for (size_t iHistType = 0; iHistType < nHistTypes; iHistType++) {
      for (size_t iVsHistMod = 0; iVsHistMod < nVsHistMods; iVsHistMod++) {
        const string sHistName2D = vecBaseHistBase[iHistRow][iHistType] + vecVsHistModifiers[iVsHistMod];
        vecBaseHist2D[iHistRow][iHistType].push_back(
          new TH2D(
            sHistName2D.data(),
            "",
            get<0>(vecVsHistBins[iVsHistMod]),
            get<1>(vecVsHistBins[iVsHistMod]).first,
            get<1>(vecVsHistBins[iVsHistMod]).second,
            get<0>(vecBaseHistBins[iHistRow]),
            get<1>(vecBaseHistBins[iHistRow]).first,
            get<1>(vecBaseHistBins[iHistRow]).second
          )
        );
      }  // end modifier loop
    }  // end type loop
  }  // end row loop

  // set errors
  for (auto histRow1D : vecBaseHist1D) {
    for (auto hist1D : histRow1D) {
      hist1D -> Sumw2();
    }
  }

  for (auto histRow2D : vecBaseHist2D) {
    for (auto histType2D : histRow2D) {
      for (auto hist2D: histType2D) {
        hist2D -> Sumw2();
      }
    }
  }

  // set axis titles
  size_t iVar = 0;
  for (auto histRow1D : vecBaseHist1D) {
    for (auto hist1D : histRow1D) {
      hist1D -> GetXaxis() -> SetTitle(vecBaseAxisVars.at(iVar).data());
      hist1D -> GetYaxis() -> SetTitle(m_config.count.data());
    }
    ++iVar;
  }
  iVar = 0;

  size_t iComp = 0;
  for (auto histRow2D : vecBaseHist2D) {
    for (auto histType2D : histRow2D) {
      iComp = 0;
      for (auto hist2D : histType2D) {
        hist2D -> GetXaxis() -> SetTitle(vecVsAxisVars.at(iComp).data());
        hist2D -> GetYaxis() -> SetTitle(vecBaseAxisVars.at(iVar).data());
        hist2D -> GetZaxis() -> SetTitle(m_config.count.data());
        ++iComp;
      }
    }
    ++iVar;
  }
  cout << "      Initialized base histograms." << endl;

  /* set up output histograms here */
  return;

}  // end 'InitHists()'



void STrackMatcherComparator::GetNewTreeHists() {

  // announce start of routine
  cout << "      Grabbing new matcher tree histograms:" << endl;

  // declare input leaves (ignore cluster branches for now)
  int   tru_event;
  int   tru_nphg4_part;
  float tru_centrality;
  int   tru_ntrackmatches;
  int   tru_nphg4;
  int   tru_nsvtx;
  int   tru_trackid;
  bool  tru_is_G4track;
  bool  tru_is_Svtrack;
  bool  tru_is_matched;
  float tru_trkpt;
  float tru_trketa;
  float tru_trkphi;
  int   tru_nclus;
  int   tru_nclustpc;
  int   tru_nclusmvtx;
  int   tru_nclusintt;
  float tru_matchrat;
  float tru_matchrat_intt;
  float tru_matchrat_mvtx;
  float tru_matchrat_tpc;

  int   rec_event;
  int   rec_nphg4_part;
  float rec_centrality;
  int   rec_ntrackmatches;
  int   rec_nphg4;
  int   rec_nsvtx;
  int   rec_trackid;
  bool  rec_is_G4track;
  bool  rec_is_Svtrack;
  bool  rec_is_matched;
  float rec_trkpt;
  float rec_trketa;
  float rec_trkphi;
  int   rec_nclus;
  int   rec_nclustpc;
  int   rec_nclusmvtx;
  int   rec_nclusintt;
  float rec_matchrat;
  float rec_matchrat_intt;
  float rec_matchrat_mvtx;
  float rec_matchrat_tpc;

  // set branch addresses (ignore cluster branches for now)
  tInTrue -> SetBranchAddress("event",         &tru_event);
  tInTrue -> SetBranchAddress("nphg4_part",    &tru_nphg4_part);
  tInTrue -> SetBranchAddress("centrality",    &tru_centrality);
  tInTrue -> SetBranchAddress("ntrackmatches", &tru_ntrackmatches);
  tInTrue -> SetBranchAddress("nphg4",         &tru_nphg4);
  tInTrue -> SetBranchAddress("nsvtx",         &tru_nsvtx);
  tInTrue -> SetBranchAddress("trackid",       &tru_trackid);
  tInTrue -> SetBranchAddress("is_G4track",    &tru_is_G4track);
  tInTrue -> SetBranchAddress("is_Svtrack",    &tru_is_Svtrack);
  tInTrue -> SetBranchAddress("is_matched",    &tru_is_matched);
  tInTrue -> SetBranchAddress("trkpt",         &tru_trkpt);
  tInTrue -> SetBranchAddress("trketa",        &tru_trketa);
  tInTrue -> SetBranchAddress("trkphi",        &tru_trkphi);
  tInTrue -> SetBranchAddress("nclus",         &tru_nclus);
  tInTrue -> SetBranchAddress("nclustpc",      &tru_nclustpc);
  tInTrue -> SetBranchAddress("nclusmvtx",     &tru_nclusmvtx);
  tInTrue -> SetBranchAddress("nclusintt",     &tru_nclusintt);
  tInTrue -> SetBranchAddress("matchrat",      &tru_matchrat);
  tInTrue -> SetBranchAddress("matchrat_intt", &tru_matchrat_intt);
  tInTrue -> SetBranchAddress("matchrat_mvtx", &tru_matchrat_mvtx);
  tInTrue -> SetBranchAddress("matchrat_tpc",  &tru_matchrat_tpc);

  tInReco -> SetBranchAddress("event",         &rec_event);
  tInReco -> SetBranchAddress("nphg4_part",    &rec_nphg4_part);
  tInReco -> SetBranchAddress("centrality",    &rec_centrality);
  tInReco -> SetBranchAddress("ntrackmatches", &rec_ntrackmatches);
  tInReco -> SetBranchAddress("nphg4",         &rec_nphg4);
  tInReco -> SetBranchAddress("nsvtx",         &rec_nsvtx);
  tInReco -> SetBranchAddress("trackid",       &rec_trackid);
  tInReco -> SetBranchAddress("is_G4track",    &rec_is_G4track);
  tInReco -> SetBranchAddress("is_Svtrack",    &rec_is_Svtrack);
  tInReco -> SetBranchAddress("is_matched",    &rec_is_matched);
  tInReco -> SetBranchAddress("trkpt",         &rec_trkpt);
  tInReco -> SetBranchAddress("trketa",        &rec_trketa);
  tInReco -> SetBranchAddress("trkphi",        &rec_trkphi);
  tInReco -> SetBranchAddress("nclus",         &rec_nclus);
  tInReco -> SetBranchAddress("nclustpc",      &rec_nclustpc);
  tInReco -> SetBranchAddress("nclusmvtx",     &rec_nclusmvtx);
  tInReco -> SetBranchAddress("nclusintt",     &rec_nclusintt);
  tInReco -> SetBranchAddress("matchrat",      &rec_matchrat);
  tInReco -> SetBranchAddress("matchrat_intt", &rec_matchrat_intt);
  tInReco -> SetBranchAddress("matchrat_mvtx", &rec_matchrat_mvtx);
  tInReco -> SetBranchAddress("matchrat_tpc",  &rec_matchrat_tpc);
  cout << "        Set input branches." << endl;

  // grab no. of entries
  const int64_t nTrueEntries = tInTrue -> GetEntries();
  const int64_t nRecoEntries = tInReco -> GetEntries(); 
  cout << "        Beginning truth particle loop: " << nTrueEntries << " to process" << endl;

  // loop over truth particles
  int64_t nTrueBytes = 0;
  for (int64_t iTrueEntry = 0; iTrueEntry < nTrueEntries; iTrueEntry++) {

    // grab truth particle entry
    const int64_t trueBytes = tInTrue -> GetEntry(iTrueEntry);
    if (trueBytes < 0) {
      cerr << "PANIC: issue with entry " << iTrueEntry << "! Aborting loop!\n" << endl;
      break;
    } else {
      nTrueBytes += trueBytes;
    }

    const int64_t iTrueProg = iTrueEntry + 1;
    if (iTrueProg == nTrueEntries) {
      cout << "          Processing entry " << iTrueProg << "/" << nTrueEntries << "..." << endl;
    } else {
      cout << "          Processing entry " << iTrueProg << "/" << nTrueEntries << "...\r" << flush;
    }

    // select truth particles
    if (!tru_is_G4track) continue;

    // fill truth 1D histograms
    vecHist1D[Var::NTot][Type::Truth]  -> Fill(tru_nclus);
    vecHist1D[Var::NIntt][Type::Truth] -> Fill(tru_nclusintt);
    vecHist1D[Var::NMvtx][Type::Truth] -> Fill(tru_nclusmvtx);
    vecHist1D[Var::NTpc][Type::Truth]  -> Fill(tru_nclustpc);
    vecHist1D[Var::RTot][Type::Truth]  -> Fill(1.);
    vecHist1D[Var::RIntt][Type::Truth] -> Fill(1.);
    vecHist1D[Var::RMvtx][Type::Truth] -> Fill(1.);
    vecHist1D[Var::RTpc][Type::Truth]  -> Fill(1.);
    vecHist1D[Var::Phi][Type::Truth]   -> Fill(tru_trkphi);
    vecHist1D[Var::Eta][Type::Truth]   -> Fill(tru_trketa);
    vecHist1D[Var::Pt][Type::Truth]    -> Fill(tru_trkpt);
    vecHist1D[Var::Frac][Type::Truth]  -> Fill(1.);

    // fill truth 2D histograms
    vecHist2D[Var::NTot][Type::Truth][Comp::VsTruthPt]  -> Fill(tru_trkpt, tru_nclus);
    vecHist2D[Var::NIntt][Type::Truth][Comp::VsTruthPt] -> Fill(tru_trkpt, tru_nclusintt);
    vecHist2D[Var::NMvtx][Type::Truth][Comp::VsTruthPt] -> Fill(tru_trkpt, tru_nclusmvtx);
    vecHist2D[Var::NTpc][Type::Truth][Comp::VsTruthPt]  -> Fill(tru_trkpt, tru_nclustpc);
    vecHist2D[Var::RTot][Type::Truth][Comp::VsTruthPt]  -> Fill(tru_trkpt, 1.);
    vecHist2D[Var::RIntt][Type::Truth][Comp::VsTruthPt] -> Fill(tru_trkpt, 1.);
    vecHist2D[Var::RMvtx][Type::Truth][Comp::VsTruthPt] -> Fill(tru_trkpt, 1.);
    vecHist2D[Var::RTpc][Type::Truth][Comp::VsTruthPt]  -> Fill(tru_trkpt, 1.);
    vecHist2D[Var::Phi][Type::Truth][Comp::VsTruthPt]   -> Fill(tru_trkpt, tru_trkphi);
    vecHist2D[Var::Eta][Type::Truth][Comp::VsTruthPt]   -> Fill(tru_trkpt, tru_trketa);
    vecHist2D[Var::Pt][Type::Truth][Comp::VsTruthPt]    -> Fill(tru_trkpt, tru_trkpt);
    vecHist2D[Var::Frac][Type::Truth][Comp::VsTruthPt]  -> Fill(tru_trkpt, 1.);

    vecHist2D[Var::NTot][Type::Truth][Comp::VsNumTpc]  -> Fill(tru_nclustpc, tru_nclus);
    vecHist2D[Var::NIntt][Type::Truth][Comp::VsNumTpc] -> Fill(tru_nclustpc, tru_nclusintt);
    vecHist2D[Var::NMvtx][Type::Truth][Comp::VsNumTpc] -> Fill(tru_nclustpc, tru_nclusmvtx);
    vecHist2D[Var::NTpc][Type::Truth][Comp::VsNumTpc]  -> Fill(tru_nclustpc, tru_nclustpc);
    vecHist2D[Var::RTot][Type::Truth][Comp::VsNumTpc]  -> Fill(tru_nclustpc, 1.);
    vecHist2D[Var::RIntt][Type::Truth][Comp::VsNumTpc] -> Fill(tru_nclustpc, 1.);
    vecHist2D[Var::RMvtx][Type::Truth][Comp::VsNumTpc] -> Fill(tru_nclustpc, 1.);
    vecHist2D[Var::RTpc][Type::Truth][Comp::VsNumTpc]  -> Fill(tru_nclustpc, 1.);
    vecHist2D[Var::Phi][Type::Truth][Comp::VsNumTpc]   -> Fill(tru_nclustpc, tru_trkphi);
    vecHist2D[Var::Eta][Type::Truth][Comp::VsNumTpc]   -> Fill(tru_nclustpc, tru_trketa);
    vecHist2D[Var::Pt][Type::Truth][Comp::VsNumTpc]    -> Fill(tru_nclustpc, tru_trkpt);
    vecHist2D[Var::Frac][Type::Truth][Comp::VsNumTpc]  -> Fill(tru_nclustpc, 1.);
  }  // end truth particle loop

  // announce next entry loop
  cout << "        Finished truth particle loop.\n"
       << "        Beginning reconstructed track loop: " << nRecoEntries << " to process"
       << endl;

  // loop over reco tracks
  // TODO identify matched truth particles
  int64_t nRecoBytes = 0;
  for (int64_t iRecoEntry = 0; iRecoEntry < nRecoEntries; iRecoEntry++) {

    // grab reco track entry
    const int64_t recoBytes = tInReco -> GetEntry(iRecoEntry);
    if (recoBytes < 0) {
      cerr << "PANIC: issue with entry " << iRecoEntry << "! Aborting loop!\n" << endl;
      break;
    } else {
      nRecoBytes += recoBytes;
    }

    const int64_t iRecoProg = iRecoEntry + 1;
    if (iRecoProg == nRecoEntries) {
      cout << "          Processing entry " << iRecoProg << "/" << nRecoEntries << "..." << endl;
    } else {
      cout << "          Processing entry " << iRecoProg << "/" << nRecoEntries << "...\r" << flush;
    }

    // select only tracks matched to truth particle
    if (!rec_is_matched || !rec_is_Svtrack) continue;

    // fill all matched reco 1D histograms
    // FIXME actually calculate pt fraction
    vecHist1D[Var::NTot][Type::Track]  -> Fill(rec_nclus);
    vecHist1D[Var::NIntt][Type::Track] -> Fill(rec_nclusintt);
    vecHist1D[Var::NMvtx][Type::Track] -> Fill(rec_nclusmvtx);
    vecHist1D[Var::NTpc][Type::Track]  -> Fill(rec_nclustpc);
    vecHist1D[Var::RTot][Type::Track]  -> Fill(rec_matchrat);
    vecHist1D[Var::RIntt][Type::Track] -> Fill(rec_matchrat_intt);
    vecHist1D[Var::RMvtx][Type::Track] -> Fill(rec_matchrat_mvtx);
    vecHist1D[Var::RTpc][Type::Track]  -> Fill(rec_matchrat_tpc);
    vecHist1D[Var::Phi][Type::Track]   -> Fill(rec_trkphi);
    vecHist1D[Var::Eta][Type::Track]   -> Fill(rec_trketa);
    vecHist1D[Var::Pt][Type::Track]    -> Fill(rec_trkpt);
    vecHist1D[Var::Frac][Type::Track]  -> Fill(1.);

    // fill all matched reco 2D histograms
    // FIXME use actual truth pt & ntpc of matched particle
    vecHist2D[Var::NTot][Type::Track][Comp::VsTruthPt]  -> Fill(rec_trkpt, rec_nclus);
    vecHist2D[Var::NIntt][Type::Track][Comp::VsTruthPt] -> Fill(rec_trkpt, rec_nclusintt);
    vecHist2D[Var::NMvtx][Type::Track][Comp::VsTruthPt] -> Fill(rec_trkpt, rec_nclusmvtx);
    vecHist2D[Var::NTpc][Type::Track][Comp::VsTruthPt]  -> Fill(rec_trkpt, rec_nclustpc);
    vecHist2D[Var::RTot][Type::Track][Comp::VsTruthPt]  -> Fill(rec_trkpt, rec_matchrat);
    vecHist2D[Var::RIntt][Type::Track][Comp::VsTruthPt] -> Fill(rec_trkpt, rec_matchrat_intt);
    vecHist2D[Var::RMvtx][Type::Track][Comp::VsTruthPt] -> Fill(rec_trkpt, rec_matchrat_mvtx);
    vecHist2D[Var::RTpc][Type::Track][Comp::VsTruthPt]  -> Fill(rec_trkpt, rec_matchrat_tpc);
    vecHist2D[Var::Phi][Type::Track][Comp::VsTruthPt]   -> Fill(rec_trkpt, rec_trkphi);
    vecHist2D[Var::Eta][Type::Track][Comp::VsTruthPt]   -> Fill(rec_trkpt, rec_trketa);
    vecHist2D[Var::Pt][Type::Track][Comp::VsTruthPt]    -> Fill(rec_trkpt, rec_trkpt);
    vecHist2D[Var::Frac][Type::Track][Comp::VsTruthPt]  -> Fill(rec_trkpt, 1.);

    vecHist2D[Var::NTot][Type::Track][Comp::VsNumTpc]  -> Fill(rec_nclustpc, rec_nclus);
    vecHist2D[Var::NIntt][Type::Track][Comp::VsNumTpc] -> Fill(rec_nclustpc, rec_nclusintt);
    vecHist2D[Var::NMvtx][Type::Track][Comp::VsNumTpc] -> Fill(rec_nclustpc, rec_nclusmvtx);
    vecHist2D[Var::NTpc][Type::Track][Comp::VsNumTpc]  -> Fill(rec_nclustpc, rec_nclustpc);
    vecHist2D[Var::RTot][Type::Track][Comp::VsNumTpc]  -> Fill(rec_nclustpc, rec_matchrat);
    vecHist2D[Var::RIntt][Type::Track][Comp::VsNumTpc] -> Fill(rec_nclustpc, rec_matchrat_intt);
    vecHist2D[Var::RMvtx][Type::Track][Comp::VsNumTpc] -> Fill(rec_nclustpc, rec_matchrat_mvtx);
    vecHist2D[Var::RTpc][Type::Track][Comp::VsNumTpc]  -> Fill(rec_nclustpc, rec_matchrat_tpc);
    vecHist2D[Var::Phi][Type::Track][Comp::VsNumTpc]   -> Fill(rec_nclustpc, rec_trkphi);
    vecHist2D[Var::Eta][Type::Track][Comp::VsNumTpc]   -> Fill(rec_nclustpc, rec_trketa);
    vecHist2D[Var::Pt][Type::Track][Comp::VsNumTpc]    -> Fill(rec_nclustpc, rec_trkpt);
    vecHist2D[Var::Frac][Type::Track][Comp::VsNumTpc]  -> Fill(rec_nclustpc, 1.);

    // fill weird and normal matched reco 1D histograms
    // FIXME actually check if is a weird track
    const bool isWeirdTrack = true;
    if (isWeirdTrack) {
      vecHist1D[Var::NTot][Type::Weird]  -> Fill(rec_nclus);
      vecHist1D[Var::NIntt][Type::Weird] -> Fill(rec_nclusintt);
      vecHist1D[Var::NMvtx][Type::Weird] -> Fill(rec_nclusmvtx);
      vecHist1D[Var::NTpc][Type::Weird]  -> Fill(rec_nclustpc);
      vecHist1D[Var::RTot][Type::Weird]  -> Fill(rec_matchrat);
      vecHist1D[Var::RIntt][Type::Weird] -> Fill(rec_matchrat_intt);
      vecHist1D[Var::RMvtx][Type::Weird] -> Fill(rec_matchrat_mvtx);
      vecHist1D[Var::RTpc][Type::Weird]  -> Fill(rec_matchrat_tpc);
      vecHist1D[Var::Phi][Type::Weird]   -> Fill(rec_trkphi);
      vecHist1D[Var::Eta][Type::Weird]   -> Fill(rec_trketa);
      vecHist1D[Var::Pt][Type::Weird]    -> Fill(rec_trkpt);
      vecHist1D[Var::Frac][Type::Weird]  -> Fill(1.);

      vecHist2D[Var::NTot][Type::Weird][Comp::VsTruthPt]  -> Fill(rec_trkpt, rec_nclus);
      vecHist2D[Var::NIntt][Type::Weird][Comp::VsTruthPt] -> Fill(rec_trkpt, rec_nclusintt);
      vecHist2D[Var::NMvtx][Type::Weird][Comp::VsTruthPt] -> Fill(rec_trkpt, rec_nclusmvtx);
      vecHist2D[Var::NTpc][Type::Weird][Comp::VsTruthPt]  -> Fill(rec_trkpt, rec_nclustpc);
      vecHist2D[Var::RTot][Type::Weird][Comp::VsTruthPt]  -> Fill(rec_trkpt, rec_matchrat);
      vecHist2D[Var::RIntt][Type::Weird][Comp::VsTruthPt] -> Fill(rec_trkpt, rec_matchrat_intt);
      vecHist2D[Var::RMvtx][Type::Weird][Comp::VsTruthPt] -> Fill(rec_trkpt, rec_matchrat_mvtx);
      vecHist2D[Var::RTpc][Type::Weird][Comp::VsTruthPt]  -> Fill(rec_trkpt, rec_matchrat_tpc);
      vecHist2D[Var::Phi][Type::Weird][Comp::VsTruthPt]   -> Fill(rec_trkpt, rec_trkphi);
      vecHist2D[Var::Eta][Type::Weird][Comp::VsTruthPt]   -> Fill(rec_trkpt, rec_trketa);
      vecHist2D[Var::Pt][Type::Weird][Comp::VsTruthPt]    -> Fill(rec_trkpt, rec_trkpt);
      vecHist2D[Var::Frac][Type::Weird][Comp::VsTruthPt]  -> Fill(rec_trkpt, 1.);

      vecHist2D[Var::NTot][Type::Weird][Comp::VsNumTpc]  -> Fill(rec_nclustpc, rec_nclus);
      vecHist2D[Var::NIntt][Type::Weird][Comp::VsNumTpc] -> Fill(rec_nclustpc, rec_nclusintt);
      vecHist2D[Var::NMvtx][Type::Weird][Comp::VsNumTpc] -> Fill(rec_nclustpc, rec_nclusmvtx);
      vecHist2D[Var::NTpc][Type::Weird][Comp::VsNumTpc]  -> Fill(rec_nclustpc, rec_nclustpc);
      vecHist2D[Var::RTot][Type::Weird][Comp::VsNumTpc]  -> Fill(rec_nclustpc, rec_matchrat);
      vecHist2D[Var::RIntt][Type::Weird][Comp::VsNumTpc] -> Fill(rec_nclustpc, rec_matchrat_intt);
      vecHist2D[Var::RMvtx][Type::Weird][Comp::VsNumTpc] -> Fill(rec_nclustpc, rec_matchrat_mvtx);
      vecHist2D[Var::RTpc][Type::Weird][Comp::VsNumTpc]  -> Fill(rec_nclustpc, rec_matchrat_tpc);
      vecHist2D[Var::Phi][Type::Weird][Comp::VsNumTpc]   -> Fill(rec_nclustpc, rec_trkphi);
      vecHist2D[Var::Eta][Type::Weird][Comp::VsNumTpc]   -> Fill(rec_nclustpc, rec_trketa);
      vecHist2D[Var::Pt][Type::Weird][Comp::VsNumTpc]    -> Fill(rec_nclustpc, rec_trkpt);
      vecHist2D[Var::Frac][Type::Weird][Comp::VsNumTpc]  -> Fill(rec_nclustpc, 1.);
    } else {
      vecHist1D[Var::NTot][Type::Normal]  -> Fill(rec_nclus);
      vecHist1D[Var::NIntt][Type::Normal] -> Fill(rec_nclusintt);
      vecHist1D[Var::NMvtx][Type::Normal] -> Fill(rec_nclusmvtx);
      vecHist1D[Var::NTpc][Type::Normal]  -> Fill(rec_nclustpc);
      vecHist1D[Var::RTot][Type::Normal]  -> Fill(rec_matchrat);
      vecHist1D[Var::RIntt][Type::Normal] -> Fill(rec_matchrat_intt);
      vecHist1D[Var::RMvtx][Type::Normal] -> Fill(rec_matchrat_mvtx);
      vecHist1D[Var::RTpc][Type::Normal]  -> Fill(rec_matchrat_tpc);
      vecHist1D[Var::Phi][Type::Normal]   -> Fill(rec_trkphi);
      vecHist1D[Var::Eta][Type::Normal]   -> Fill(rec_trketa);
      vecHist1D[Var::Pt][Type::Normal]    -> Fill(rec_trkpt);
      vecHist1D[Var::Frac][Type::Normal]  -> Fill(1.);

      vecHist2D[Var::NTot][Type::Normal][Comp::VsTruthPt]  -> Fill(rec_trkpt, rec_nclus);
      vecHist2D[Var::NIntt][Type::Normal][Comp::VsTruthPt] -> Fill(rec_trkpt, rec_nclusintt);
      vecHist2D[Var::NMvtx][Type::Normal][Comp::VsTruthPt] -> Fill(rec_trkpt, rec_nclusmvtx);
      vecHist2D[Var::NTpc][Type::Normal][Comp::VsTruthPt]  -> Fill(rec_trkpt, rec_nclustpc);
      vecHist2D[Var::RTot][Type::Normal][Comp::VsTruthPt]  -> Fill(rec_trkpt, rec_matchrat);
      vecHist2D[Var::RIntt][Type::Normal][Comp::VsTruthPt] -> Fill(rec_trkpt, rec_matchrat_intt);
      vecHist2D[Var::RMvtx][Type::Normal][Comp::VsTruthPt] -> Fill(rec_trkpt, rec_matchrat_mvtx);
      vecHist2D[Var::RTpc][Type::Normal][Comp::VsTruthPt]  -> Fill(rec_trkpt, rec_matchrat_tpc);
      vecHist2D[Var::Phi][Type::Normal][Comp::VsTruthPt]   -> Fill(rec_trkpt, rec_trkphi);
      vecHist2D[Var::Eta][Type::Normal][Comp::VsTruthPt]   -> Fill(rec_trkpt, rec_trketa);
      vecHist2D[Var::Pt][Type::Normal][Comp::VsTruthPt]    -> Fill(rec_trkpt, rec_trkpt);
      vecHist2D[Var::Frac][Type::Normal][Comp::VsTruthPt]  -> Fill(rec_trkpt, 1.);

      vecHist2D[Var::NTot][Type::Normal][Comp::VsNumTpc]  -> Fill(rec_nclustpc, rec_nclus);
      vecHist2D[Var::NIntt][Type::Normal][Comp::VsNumTpc] -> Fill(rec_nclustpc, rec_nclusintt);
      vecHist2D[Var::NMvtx][Type::Normal][Comp::VsNumTpc] -> Fill(rec_nclustpc, rec_nclusmvtx);
      vecHist2D[Var::NTpc][Type::Normal][Comp::VsNumTpc]  -> Fill(rec_nclustpc, rec_nclustpc);
      vecHist2D[Var::RTot][Type::Normal][Comp::VsNumTpc]  -> Fill(rec_nclustpc, rec_matchrat);
      vecHist2D[Var::RIntt][Type::Normal][Comp::VsNumTpc] -> Fill(rec_nclustpc, rec_matchrat_intt);
      vecHist2D[Var::RMvtx][Type::Normal][Comp::VsNumTpc] -> Fill(rec_nclustpc, rec_matchrat_mvtx);
      vecHist2D[Var::RTpc][Type::Normal][Comp::VsNumTpc]  -> Fill(rec_nclustpc, rec_matchrat_tpc);
      vecHist2D[Var::Phi][Type::Normal][Comp::VsNumTpc]   -> Fill(rec_nclustpc, rec_trkphi);
      vecHist2D[Var::Eta][Type::Normal][Comp::VsNumTpc]   -> Fill(rec_nclustpc, rec_trketa);
      vecHist2D[Var::Pt][Type::Normal][Comp::VsNumTpc]    -> Fill(rec_nclustpc, rec_trkpt);
      vecHist2D[Var::Frac][Type::Normal][Comp::VsNumTpc]  -> Fill(rec_nclustpc, 1.);
    }
  }  // end reco track loop

  // announce method end and return
  cout << "        Finished reconstructed track loop.\n"
       << "      Finished getting histograms from new matcher cluster tree."
       << endl;
  return;

}  // end 'GetNewTreeHists()'



void STrackMatcherComparator::GetNewTupleHists() {

  // announce start of method
  cout << "      Grabbing new matcher tuple histograms:" << endl;

  // declare input leaves
  float tru_evtid;
  float tru_trkid;
  float tru_pt;
  float tru_eta;
  float tru_phi;
  float tru_nmvtxclust_trkmatcher;
  float tru_ninttclust_trkmatcher;
  float tru_ntpclust_trkmatcher;
  float tru_nmvtxclust_manual;
  float tru_ninttclust_manual;
  float tru_ntpcclust_manual;
  float tru_nmvtxlayer_;
  float tru_ninttlayer;
  float tru_ntpclayer;
  float tru_deltapt;
  float tru_quality;
  float tru_dcaxy;
  float tru_dcaz;
  float tru_sigmadcaxy;
  float tru_sigmadacz;
  float tru_vx;
  float tru_vy;
  float tru_vz;
  float tru_gtrkid;
  float tru_gpt;
  float tru_geta;
  float tru_gphi;
  float tru_gnmvtxclust_trkmatcher;
  float tru_gninttclust_trkmatcher;
  float tru_gntpclust_trkmatchergnmvtxclust_manual;
  float tru_gninttclust_manual;
  float tru_gntpcclust_manual;
  float tru_gnmvtxlayer;
  float tru_gninttlayer;
  float tru_gntpclayer;
  float tru_gdeltapt;
  float tru_gquality;
  float tru_gdcaxy;
  float tru_gdcaz;
  float tru_gsigmadcaxy;
  float tru_gsigmadacz;
  float tru_gvx;
  float tru_gvy;
  float tru_gvz;
  float tru_fracnmvtxmatched;
  float tru_fracninttmatched;
  float tru_fracntpcmatched;

  float rec_evtid;
  float rec_trkid;
  float rec_pt;
  float rec_eta;
  float rec_phi;
  float rec_nmvtxclust_trkmatcher;
  float rec_ninttclust_trkmatcher;
  float rec_ntpclust_trkmatcher;
  float rec_nmvtxclust_manual;
  float rec_ninttclust_manual;
  float rec_ntpcclust_manual;
  float rec_nmvtxlayer_;
  float rec_ninttlayer;
  float rec_ntpclayer;
  float rec_deltapt;
  float rec_quality;
  float rec_dcaxy;
  float rec_dcaz;
  float rec_sigmadcaxy;
  float rec_sigmadacz;
  float rec_vx;
  float rec_vy;
  float rec_vz;
  float rec_gtrkid;
  float rec_gpt;
  float rec_geta;
  float rec_gphi;
  float rec_gnmvtxclust_trkmatcher;
  float rec_gninttclust_trkmatcher;
  float rec_gntpclust_trkmatchergnmvtxclust_manual;
  float rec_gninttclust_manual;
  float rec_gntpcclust_manual;
  float rec_gnmvtxlayer;
  float rec_gninttlayer;
  float rec_gntpclayer;
  float rec_gdeltapt;
  float rec_gquality;
  float rec_gdcaxy;
  float rec_gdcaz;
  float rec_gsigmadcaxy;
  float rec_gsigmadacz;
  float rec_gvx;
  float rec_gvy;
  float rec_gvz;
  float rec_fracnmvtxmatched;
  float rec_fracninttmatched;
  float rec_fracntpcmatched;
 
  // set branch addresses
  tInTrue -> SetBranchAddress("evtid",                                  &tru_evtid);
  tInTrue -> SetBranchAddress("trkid",                                  &tru_trkid);
  tInTrue -> SetBranchAddress("pt",                                     &tru_pt);
  tInTrue -> SetBranchAddress("eta",                                    &tru_eta);
  tInTrue -> SetBranchAddress("phi",                                    &tru_phi);
  tInTrue -> SetBranchAddress("nmvtxclust_trkmatcher",                  &tru_nmvtxclust_trkmatcher);
  tInTrue -> SetBranchAddress("ninttclust_trkmatcher",                  &tru_ninttclust_trkmatcher);
  tInTrue -> SetBranchAddress("ntpclust_trkmatcher",                    &tru_ntpclust_trkmatcher);
  tInTrue -> SetBranchAddress("nmvtxclust_manual",                      &tru_nmvtxclust_manual);
  tInTrue -> SetBranchAddress("ninttclust_manual",                      &tru_ninttclust_manual);
  tInTrue -> SetBranchAddress("ntpcclust_manual",                       &tru_ntpcclust_manual);
  tInTrue -> SetBranchAddress("nmvtxlayer_",                            &tru_nmvtxlayer_);
  tInTrue -> SetBranchAddress("ninttlayer",                             &tru_ninttlayer);
  tInTrue -> SetBranchAddress("ntpclayer",                              &tru_ntpclayer);
  tInTrue -> SetBranchAddress("deltapt",                                &tru_deltapt);
  tInTrue -> SetBranchAddress("quality",                                &tru_quality);
  tInTrue -> SetBranchAddress("dcaxy",                                  &tru_dcaxy);
  tInTrue -> SetBranchAddress("dcaz",                                   &tru_dcaz);
  tInTrue -> SetBranchAddress("sigmadcaxy",                             &tru_sigmadcaxy);
  tInTrue -> SetBranchAddress("sigmadacz",                              &tru_sigmadacz);
  tInTrue -> SetBranchAddress("vx",                                     &tru_vx);
  tInTrue -> SetBranchAddress("vy",                                     &tru_vy);
  tInTrue -> SetBranchAddress("vz",                                     &tru_vz);
  tInTrue -> SetBranchAddress("gtrkid",                                 &tru_gtrkid);
  tInTrue -> SetBranchAddress("gpt",                                    &tru_gpt);
  tInTrue -> SetBranchAddress("geta",                                   &tru_geta);
  tInTrue -> SetBranchAddress("gphi",                                   &tru_gphi);
  tInTrue -> SetBranchAddress("gnmvtxclust_trkmatcher",                 &tru_gnmvtxclust_trkmatcher);
  tInTrue -> SetBranchAddress("gninttclust_trkmatcher",                 &tru_gninttclust_trkmatcher);
  tInTrue -> SetBranchAddress("gntpclust_trkmatchergnmvtxclust_manual", &tru_gntpclust_trkmatchergnmvtxclust_manual);
  tInTrue -> SetBranchAddress("gninttclust_manual",                     &tru_gninttclust_manual);
  tInTrue -> SetBranchAddress("gntpcclust_manual",                      &tru_gntpcclust_manual);
  tInTrue -> SetBranchAddress("gnmvtxlayer",                            &tru_gnmvtxlayer);
  tInTrue -> SetBranchAddress("gninttlayer",                            &tru_gninttlayer);
  tInTrue -> SetBranchAddress("gntpclayer",                             &tru_gntpclayer);
  tInTrue -> SetBranchAddress("gdeltapt",                               &tru_gdeltapt);
  tInTrue -> SetBranchAddress("gquality",                               &tru_gquality);
  tInTrue -> SetBranchAddress("gdcaxy",                                 &tru_gdcaxy);
  tInTrue -> SetBranchAddress("gdcaz",                                  &tru_gdcaz);
  tInTrue -> SetBranchAddress("gsigmadcaxy",                            &tru_gsigmadcaxy);
  tInTrue -> SetBranchAddress("gsigmadacz",                             &tru_gsigmadacz);
  tInTrue -> SetBranchAddress("gvx",                                    &tru_gvx);
  tInTrue -> SetBranchAddress("gvy",                                    &tru_gvy);
  tInTrue -> SetBranchAddress("gvz",                                    &tru_gvz);
  tInTrue -> SetBranchAddress("fracnmvtxmatched",                       &tru_fracnmvtxmatched);
  tInTrue -> SetBranchAddress("fracninttmatched",                       &tru_fracninttmatched);
  tInTrue -> SetBranchAddress("fracntpcmatched",                        &tru_fracntpcmatched);

  tInReco -> SetBranchAddress("evtid",                                  &rec_evtid);
  tInReco -> SetBranchAddress("trkid",                                  &rec_trkid);
  tInReco -> SetBranchAddress("pt",                                     &rec_pt);
  tInReco -> SetBranchAddress("eta",                                    &rec_eta);
  tInReco -> SetBranchAddress("phi",                                    &rec_phi);
  tInReco -> SetBranchAddress("nmvtxclust_trkmatcher",                  &rec_nmvtxclust_trkmatcher);
  tInReco -> SetBranchAddress("ninttclust_trkmatcher",                  &rec_ninttclust_trkmatcher);
  tInReco -> SetBranchAddress("ntpclust_trkmatcher",                    &rec_ntpclust_trkmatcher);
  tInReco -> SetBranchAddress("nmvtxclust_manual",                      &rec_nmvtxclust_manual);
  tInReco -> SetBranchAddress("ninttclust_manual",                      &rec_ninttclust_manual);
  tInReco -> SetBranchAddress("ntpcclust_manual",                       &rec_ntpcclust_manual);
  tInReco -> SetBranchAddress("nmvtxlayer_",                            &rec_nmvtxlayer_);
  tInReco -> SetBranchAddress("ninttlayer",                             &rec_ninttlayer);
  tInReco -> SetBranchAddress("ntpclayer",                              &rec_ntpclayer);
  tInReco -> SetBranchAddress("deltapt",                                &rec_deltapt);
  tInReco -> SetBranchAddress("quality",                                &rec_quality);
  tInReco -> SetBranchAddress("dcaxy",                                  &rec_dcaxy);
  tInReco -> SetBranchAddress("dcaz",                                   &rec_dcaz);
  tInReco -> SetBranchAddress("sigmadcaxy",                             &rec_sigmadcaxy);
  tInReco -> SetBranchAddress("sigmadacz",                              &rec_sigmadacz);
  tInReco -> SetBranchAddress("vx",                                     &rec_vx);
  tInReco -> SetBranchAddress("vy",                                     &rec_vy);
  tInReco -> SetBranchAddress("vz",                                     &rec_vz);
  tInReco -> SetBranchAddress("gtrkid",                                 &rec_gtrkid);
  tInReco -> SetBranchAddress("gpt",                                    &rec_gpt);
  tInReco -> SetBranchAddress("geta",                                   &rec_geta);
  tInReco -> SetBranchAddress("gphi",                                   &rec_gphi);
  tInReco -> SetBranchAddress("gnmvtxclust_trkmatcher",                 &rec_gnmvtxclust_trkmatcher);
  tInReco -> SetBranchAddress("gninttclust_trkmatcher",                 &rec_gninttclust_trkmatcher);
  tInReco -> SetBranchAddress("gntpclust_trkmatchergnmvtxclust_manual", &rec_gntpclust_trkmatchergnmvtxclust_manual);
  tInReco -> SetBranchAddress("gninttclust_manual",                     &rec_gninttclust_manual);
  tInReco -> SetBranchAddress("gntpcclust_manual",                      &rec_gntpcclust_manual);
  tInReco -> SetBranchAddress("gnmvtxlayer",                            &rec_gnmvtxlayer);
  tInReco -> SetBranchAddress("gninttlayer",                            &rec_gninttlayer);
  tInReco -> SetBranchAddress("gntpclayer",                             &rec_gntpclayer);
  tInReco -> SetBranchAddress("gdeltapt",                               &rec_gdeltapt);
  tInReco -> SetBranchAddress("gquality",                               &rec_gquality);
  tInReco -> SetBranchAddress("gdcaxy",                                 &rec_gdcaxy);
  tInReco -> SetBranchAddress("gdcaz",                                  &rec_gdcaz);
  tInReco -> SetBranchAddress("gsigmadcaxy",                            &rec_gsigmadcaxy);
  tInReco -> SetBranchAddress("gsigmadacz",                             &rec_gsigmadacz);
  tInReco -> SetBranchAddress("gvx",                                    &rec_gvx);
  tInReco -> SetBranchAddress("gvy",                                    &rec_gvy);
  tInReco -> SetBranchAddress("gvz",                                    &rec_gvz);
  tInReco -> SetBranchAddress("fracnmvtxmatched",                       &rec_fracnmvtxmatched);
  tInReco -> SetBranchAddress("fracninttmatched",                       &rec_fracninttmatched);
  tInReco -> SetBranchAddress("fracntpcmatched",                        &rec_fracntpcmatched);
  cout << "        Set input branches." << endl;

  // grab no. of entries
  const int64_t nTrueEntries = tInTrue -> GetEntries();
  const int64_t nRecoEntries = tInReco -> GetEntries(); 
  cout << "        Beginning truth particle loop: " << nTrueEntries << " to process" << endl;

  // loop over truth particles
  int64_t nTrueBytes = 0;
  for (int64_t iTrueEntry = 0; iTrueEntry < nTrueEntries; iTrueEntry++) {

    // grab truth particle entry
    const int64_t trueBytes = tInTrue -> GetEntry(iTrueEntry);
    if (trueBytes < 0) {
      cerr << "PANIC: issue with entry " << iTrueEntry << "! Aborting loop!\n" << endl;
      break;
    } else {
      nTrueBytes += trueBytes;
    }

    const int64_t iTrueProg = iTrueEntry + 1;
    if (iTrueProg == nTrueEntries) {
      cout << "          Processing entry " << iTrueProg << "/" << nTrueEntries << "..." << endl;
    } else {
      cout << "          Processing entry " << iTrueProg << "/" << nTrueEntries << "...\r" << flush;
    }

    // run calculations
    const float tru_gnclust = tru_gnmvtxclust_trkmatcher + tru_gninttclust_trkmatcher + tru_gntpcclust_manual;

    // fill truth 1D histograms
    vecHist1D[Var::NTot][Type::Truth]  -> Fill(tru_gnclust);
    vecHist1D[Var::NIntt][Type::Truth] -> Fill(tru_gninttclust_trkmatcher);
    vecHist1D[Var::NMvtx][Type::Truth] -> Fill(tru_gnmvtxclust_trkmatcher);
    vecHist1D[Var::NTpc][Type::Truth]  -> Fill(tru_gntpcclust_manual);
    vecHist1D[Var::RTot][Type::Truth]  -> Fill(1.);
    vecHist1D[Var::RIntt][Type::Truth] -> Fill(1.);
    vecHist1D[Var::RMvtx][Type::Truth] -> Fill(1.);
    vecHist1D[Var::RTpc][Type::Truth]  -> Fill(1.);
    vecHist1D[Var::Phi][Type::Truth]   -> Fill(tru_gphi);
    vecHist1D[Var::Eta][Type::Truth]   -> Fill(tru_geta);
    vecHist1D[Var::Pt][Type::Truth]    -> Fill(tru_gpt);
    vecHist1D[Var::Frac][Type::Truth]  -> Fill(1.);

    // fill truth 2D histograms
    vecHist2D[Var::NTot][Type::Truth][Comp::VsTruthPt]  -> Fill(tru_gpt, tru_gnclust);
    vecHist2D[Var::NIntt][Type::Truth][Comp::VsTruthPt] -> Fill(tru_gpt, tru_gninttclust_trkmatcher);
    vecHist2D[Var::NMvtx][Type::Truth][Comp::VsTruthPt] -> Fill(tru_gpt, tru_gnmvtxclust_trkmatcher);
    vecHist2D[Var::NTpc][Type::Truth][Comp::VsTruthPt]  -> Fill(tru_gpt, tru_gntpcclust_manual);
    vecHist2D[Var::RTot][Type::Truth][Comp::VsTruthPt]  -> Fill(tru_gpt, 1.);
    vecHist2D[Var::RIntt][Type::Truth][Comp::VsTruthPt] -> Fill(tru_gpt, 1.);
    vecHist2D[Var::RMvtx][Type::Truth][Comp::VsTruthPt] -> Fill(tru_gpt, 1.);
    vecHist2D[Var::RTpc][Type::Truth][Comp::VsTruthPt]  -> Fill(tru_gpt, 1.);
    vecHist2D[Var::Phi][Type::Truth][Comp::VsTruthPt]   -> Fill(tru_gpt, tru_gphi);
    vecHist2D[Var::Eta][Type::Truth][Comp::VsTruthPt]   -> Fill(tru_gpt, tru_geta);
    vecHist2D[Var::Pt][Type::Truth][Comp::VsTruthPt]    -> Fill(tru_gpt, tru_gpt);
    vecHist2D[Var::Frac][Type::Truth][Comp::VsTruthPt]  -> Fill(tru_gpt, 1.);

    vecHist2D[Var::NTot][Type::Truth][Comp::VsNumTpc]  -> Fill(tru_ntpclust_trkmatcher, tru_gnclust);
    vecHist2D[Var::NIntt][Type::Truth][Comp::VsNumTpc] -> Fill(tru_ntpclust_trkmatcher, tru_gninttclust_trkmatcher);
    vecHist2D[Var::NMvtx][Type::Truth][Comp::VsNumTpc] -> Fill(tru_ntpclust_trkmatcher, tru_gnmvtxclust_trkmatcher);
    vecHist2D[Var::NTpc][Type::Truth][Comp::VsNumTpc]  -> Fill(tru_ntpclust_trkmatcher, tru_gntpcclust_manual);
    vecHist2D[Var::RTot][Type::Truth][Comp::VsNumTpc]  -> Fill(tru_ntpclust_trkmatcher, 1.);
    vecHist2D[Var::RIntt][Type::Truth][Comp::VsNumTpc] -> Fill(tru_ntpclust_trkmatcher, 1.);
    vecHist2D[Var::RMvtx][Type::Truth][Comp::VsNumTpc] -> Fill(tru_ntpclust_trkmatcher, 1.);
    vecHist2D[Var::RTpc][Type::Truth][Comp::VsNumTpc]  -> Fill(tru_ntpclust_trkmatcher, 1.);
    vecHist2D[Var::Phi][Type::Truth][Comp::VsNumTpc]   -> Fill(tru_ntpclust_trkmatcher, tru_gphi);
    vecHist2D[Var::Eta][Type::Truth][Comp::VsNumTpc]   -> Fill(tru_ntpclust_trkmatcher, tru_geta);
    vecHist2D[Var::Pt][Type::Truth][Comp::VsNumTpc]    -> Fill(tru_ntpclust_trkmatcher, tru_gpt);
    vecHist2D[Var::Frac][Type::Truth][Comp::VsNumTpc]  -> Fill(tru_ntpclust_trkmatcher, 1.);
  }  // end truth particle loop

  // announce next entry loop
  cout << "        Finished truth particle loop.\n"
       << "        Beginning reconstructed track loop: " << nRecoEntries << " to process"
       << endl;

  // loop over reco tracks
  int64_t nRecoBytes = 0;
  for (int64_t iRecoEntry = 0; iRecoEntry < nRecoEntries; iRecoEntry++) {

    // grab reco track entry
    const int64_t recoBytes = tInReco -> GetEntry(iRecoEntry);
    if (recoBytes < 0) {
      cerr << "PANIC: issue with entry " << iRecoEntry << "! Aborting loop!\n" << endl;
      break;
    } else {
      nRecoBytes += recoBytes;
    }

    const int64_t iRecoProg = iRecoEntry + 1;
    if (iRecoProg == nRecoEntries) {
      cout << "          Processing entry " << iRecoProg << "/" << nRecoEntries << "..." << endl;
    } else {
      cout << "          Processing entry " << iRecoProg << "/" << nRecoEntries << "...\r" << flush;
    }

    // run calculations
    const double rec_nclus  = rec_ninttclust_trkmatcher + rec_nmvtxclust_trkmatcher + rec_ntpclust_trkmatcher;
    const double rec_gnclus = rec_gninttclust_trkmatcher + rec_gnmvtxclust_trkmatcher + rec_gntpcclust_manual;
    const double rec_rnclus = rec_nclus / rec_gnclus;
    const double rec_rintt  = rec_ninttclust_trkmatcher / rec_gninttclust_trkmatcher;
    const double rec_rmaps  = rec_nmvtxclust_trkmatcher / rec_gnmvtxclust_trkmatcher;
    const double rec_rtpc   = rec_ntpclust_trkmatcher / rec_gntpcclust_manual;
    const double rec_ptfrac = rec_pt / rec_gpt;

    // fill all matched reco 1D histograms
    vecHist1D[Var::NTot][Type::Track]  -> Fill(rec_nclus);
    vecHist1D[Var::NIntt][Type::Track] -> Fill(rec_ninttclust_trkmatcher);
    vecHist1D[Var::NMvtx][Type::Track] -> Fill(rec_nmvtxclust_trkmatcher);
    vecHist1D[Var::NTpc][Type::Track]  -> Fill(rec_ntpclust_trkmatcher);
    vecHist1D[Var::RTot][Type::Track]  -> Fill(rec_rnclus);
    vecHist1D[Var::RIntt][Type::Track] -> Fill(rec_rintt);
    vecHist1D[Var::RMvtx][Type::Track] -> Fill(rec_rmaps);
    vecHist1D[Var::RTpc][Type::Track]  -> Fill(rec_rtpc);
    vecHist1D[Var::Phi][Type::Track]   -> Fill(rec_phi);
    vecHist1D[Var::Eta][Type::Track]   -> Fill(rec_eta);
    vecHist1D[Var::Pt][Type::Track]    -> Fill(rec_pt);
    vecHist1D[Var::Frac][Type::Track]  -> Fill(rec_ptfrac);

    // fill all matched reco 2D histograms
    vecHist2D[Var::NTot][Type::Track][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_nclus);
    vecHist2D[Var::NIntt][Type::Track][Comp::VsTruthPt] -> Fill(rec_gpt, rec_ninttclust_trkmatcher);
    vecHist2D[Var::NMvtx][Type::Track][Comp::VsTruthPt] -> Fill(rec_gpt, rec_nmvtxclust_trkmatcher);
    vecHist2D[Var::NTpc][Type::Track][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_ntpclust_trkmatcher);
    vecHist2D[Var::RTot][Type::Track][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_rnclus);
    vecHist2D[Var::RIntt][Type::Track][Comp::VsTruthPt] -> Fill(rec_gpt, rec_rintt);
    vecHist2D[Var::RMvtx][Type::Track][Comp::VsTruthPt] -> Fill(rec_gpt, rec_rmaps);
    vecHist2D[Var::RTpc][Type::Track][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_rtpc);
    vecHist2D[Var::Phi][Type::Track][Comp::VsTruthPt]   -> Fill(rec_gpt, rec_phi);
    vecHist2D[Var::Eta][Type::Track][Comp::VsTruthPt]   -> Fill(rec_gpt, rec_eta);
    vecHist2D[Var::Pt][Type::Track][Comp::VsTruthPt]    -> Fill(rec_gpt, rec_pt);
    vecHist2D[Var::Frac][Type::Track][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_ptfrac);

    vecHist2D[Var::NTot][Type::Track][Comp::VsNumTpc]  -> Fill(rec_ntpclust_trkmatcher, rec_nclus);
    vecHist2D[Var::NIntt][Type::Track][Comp::VsNumTpc] -> Fill(rec_ntpclust_trkmatcher, rec_ninttclust_trkmatcher);
    vecHist2D[Var::NMvtx][Type::Track][Comp::VsNumTpc] -> Fill(rec_ntpclust_trkmatcher, rec_nmvtxclust_trkmatcher);
    vecHist2D[Var::NTpc][Type::Track][Comp::VsNumTpc]  -> Fill(rec_ntpclust_trkmatcher, rec_ntpclust_trkmatcher);
    vecHist2D[Var::RTot][Type::Track][Comp::VsNumTpc]  -> Fill(rec_ntpclust_trkmatcher, rec_rnclus);
    vecHist2D[Var::RIntt][Type::Track][Comp::VsNumTpc] -> Fill(rec_ntpclust_trkmatcher, rec_rintt);
    vecHist2D[Var::RMvtx][Type::Track][Comp::VsNumTpc] -> Fill(rec_ntpclust_trkmatcher, rec_rmaps);
    vecHist2D[Var::RTpc][Type::Track][Comp::VsNumTpc]  -> Fill(rec_ntpclust_trkmatcher, rec_rtpc);
    vecHist2D[Var::Phi][Type::Track][Comp::VsNumTpc]   -> Fill(rec_ntpclust_trkmatcher, rec_phi);
    vecHist2D[Var::Eta][Type::Track][Comp::VsNumTpc]   -> Fill(rec_ntpclust_trkmatcher, rec_eta);
    vecHist2D[Var::Pt][Type::Track][Comp::VsNumTpc]    -> Fill(rec_ntpclust_trkmatcher, rec_pt);
    vecHist2D[Var::Frac][Type::Track][Comp::VsNumTpc]  -> Fill(rec_ntpclust_trkmatcher, rec_ptfrac);

    // fill weird and normal matched reco 1D histograms
    const bool im_config.normalTrack = ((rec_ptfrac >= m_config.oddPtFrac.first) && (rec_ptfrac <= m_config.oddPtFrac.second));
    if (im_config.normalTrack) {
      vecHist1D[Var::NTot][Type::Normal]  -> Fill(rec_nclus);
      vecHist1D[Var::NIntt][Type::Normal] -> Fill(rec_ninttclust_trkmatcher);
      vecHist1D[Var::NMvtx][Type::Normal] -> Fill(rec_nmvtxclust_trkmatcher);
      vecHist1D[Var::NTpc][Type::Normal]  -> Fill(rec_ntpclust_trkmatcher);
      vecHist1D[Var::RTot][Type::Normal]  -> Fill(rec_rnclus);
      vecHist1D[Var::RIntt][Type::Normal] -> Fill(rec_rintt);
      vecHist1D[Var::RMvtx][Type::Normal] -> Fill(rec_rmaps);
      vecHist1D[Var::RTpc][Type::Normal]  -> Fill(rec_rtpc);
      vecHist1D[Var::Phi][Type::Normal]   -> Fill(rec_phi);
      vecHist1D[Var::Eta][Type::Normal]   -> Fill(rec_eta);
      vecHist1D[Var::Pt][Type::Normal]    -> Fill(rec_pt);
      vecHist1D[Var::Frac][Type::Normal]  -> Fill(rec_ptfrac);

      vecHist2D[Var::NTot][Type::Normal][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_nclus);
      vecHist2D[Var::NIntt][Type::Normal][Comp::VsTruthPt] -> Fill(rec_gpt, rec_ninttclust_trkmatcher);
      vecHist2D[Var::NMvtx][Type::Normal][Comp::VsTruthPt] -> Fill(rec_gpt, rec_nmvtxclust_trkmatcher);
      vecHist2D[Var::NTpc][Type::Normal][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_ntpclust_trkmatcher);
      vecHist2D[Var::RTot][Type::Normal][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_rnclus);
      vecHist2D[Var::RIntt][Type::Normal][Comp::VsTruthPt] -> Fill(rec_gpt, rec_rintt);
      vecHist2D[Var::RMvtx][Type::Normal][Comp::VsTruthPt] -> Fill(rec_gpt, rec_rmaps);
      vecHist2D[Var::RTpc][Type::Normal][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_rtpc);
      vecHist2D[Var::Phi][Type::Normal][Comp::VsTruthPt]   -> Fill(rec_gpt, rec_phi);
      vecHist2D[Var::Eta][Type::Normal][Comp::VsTruthPt]   -> Fill(rec_gpt, rec_eta);
      vecHist2D[Var::Pt][Type::Normal][Comp::VsTruthPt]    -> Fill(rec_gpt, rec_pt);
      vecHist2D[Var::Frac][Type::Normal][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_ptfrac);

      vecHist2D[Var::NTot][Type::Normal][Comp::VsNumTpc]  -> Fill(rec_ntpclust_trkmatcher, rec_nclus);
      vecHist2D[Var::NIntt][Type::Normal][Comp::VsNumTpc] -> Fill(rec_ntpclust_trkmatcher, rec_ninttclust_trkmatcher);
      vecHist2D[Var::NMvtx][Type::Normal][Comp::VsNumTpc] -> Fill(rec_ntpclust_trkmatcher, rec_nmvtxclust_trkmatcher);
      vecHist2D[Var::NTpc][Type::Normal][Comp::VsNumTpc]  -> Fill(rec_ntpclust_trkmatcher, rec_ntpclust_trkmatcher);
      vecHist2D[Var::RTot][Type::Normal][Comp::VsNumTpc]  -> Fill(rec_ntpclust_trkmatcher, rec_rnclus);
      vecHist2D[Var::RIntt][Type::Normal][Comp::VsNumTpc] -> Fill(rec_ntpclust_trkmatcher, rec_rintt);
      vecHist2D[Var::RMvtx][Type::Normal][Comp::VsNumTpc] -> Fill(rec_ntpclust_trkmatcher, rec_rmaps);
      vecHist2D[Var::RTpc][Type::Normal][Comp::VsNumTpc]  -> Fill(rec_ntpclust_trkmatcher, rec_rtpc);
      vecHist2D[Var::Phi][Type::Normal][Comp::VsNumTpc]   -> Fill(rec_ntpclust_trkmatcher, rec_phi);
      vecHist2D[Var::Eta][Type::Normal][Comp::VsNumTpc]   -> Fill(rec_ntpclust_trkmatcher, rec_eta);
      vecHist2D[Var::Pt][Type::Normal][Comp::VsNumTpc]    -> Fill(rec_ntpclust_trkmatcher, rec_pt);
      vecHist2D[Var::Frac][Type::Normal][Comp::VsNumTpc]  -> Fill(rec_ntpclust_trkmatcher, rec_ptfrac);
    } else {
      vecHist1D[Var::NTot][Type::Weird]  -> Fill(rec_nclus);
      vecHist1D[Var::NIntt][Type::Weird] -> Fill(rec_ninttclust_trkmatcher);
      vecHist1D[Var::NMvtx][Type::Weird] -> Fill(rec_nmvtxclust_trkmatcher);
      vecHist1D[Var::NTpc][Type::Weird]  -> Fill(rec_ntpclust_trkmatcher);
      vecHist1D[Var::RTot][Type::Weird]  -> Fill(rec_rnclus);
      vecHist1D[Var::RIntt][Type::Weird] -> Fill(rec_rintt);
      vecHist1D[Var::RMvtx][Type::Weird] -> Fill(rec_rmaps);
      vecHist1D[Var::RTpc][Type::Weird]  -> Fill(rec_rtpc);
      vecHist1D[Var::Phi][Type::Weird]   -> Fill(rec_phi);
      vecHist1D[Var::Eta][Type::Weird]   -> Fill(rec_eta);
      vecHist1D[Var::Pt][Type::Weird]    -> Fill(rec_pt);
      vecHist1D[Var::Frac][Type::Weird]  -> Fill(rec_ptfrac);

      vecHist2D[Var::NTot][Type::Weird][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_nclus);
      vecHist2D[Var::NIntt][Type::Weird][Comp::VsTruthPt] -> Fill(rec_gpt, rec_ninttclust_trkmatcher);
      vecHist2D[Var::NMvtx][Type::Weird][Comp::VsTruthPt] -> Fill(rec_gpt, rec_nmvtxclust_trkmatcher);
      vecHist2D[Var::NTpc][Type::Weird][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_ntpclust_trkmatcher);
      vecHist2D[Var::RTot][Type::Weird][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_rnclus);
      vecHist2D[Var::RIntt][Type::Weird][Comp::VsTruthPt] -> Fill(rec_gpt, rec_rintt);
      vecHist2D[Var::RMvtx][Type::Weird][Comp::VsTruthPt] -> Fill(rec_gpt, rec_rmaps);
      vecHist2D[Var::RTpc][Type::Weird][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_rtpc);
      vecHist2D[Var::Phi][Type::Weird][Comp::VsTruthPt]   -> Fill(rec_gpt, rec_phi);
      vecHist2D[Var::Eta][Type::Weird][Comp::VsTruthPt]   -> Fill(rec_gpt, rec_eta);
      vecHist2D[Var::Pt][Type::Weird][Comp::VsTruthPt]    -> Fill(rec_gpt, rec_pt);
      vecHist2D[Var::Frac][Type::Weird][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_ptfrac);

      vecHist2D[Var::NTot][Type::Weird][Comp::VsNumTpc]  -> Fill(rec_ntpclust_trkmatcher, rec_nclus);
      vecHist2D[Var::NIntt][Type::Weird][Comp::VsNumTpc] -> Fill(rec_ntpclust_trkmatcher, rec_ninttclust_trkmatcher);
      vecHist2D[Var::NMvtx][Type::Weird][Comp::VsNumTpc] -> Fill(rec_ntpclust_trkmatcher, rec_nmvtxclust_trkmatcher);
      vecHist2D[Var::NTpc][Type::Weird][Comp::VsNumTpc]  -> Fill(rec_ntpclust_trkmatcher, rec_ntpclust_trkmatcher);
      vecHist2D[Var::RTot][Type::Weird][Comp::VsNumTpc]  -> Fill(rec_ntpclust_trkmatcher, rec_rnclus);
      vecHist2D[Var::RIntt][Type::Weird][Comp::VsNumTpc] -> Fill(rec_ntpclust_trkmatcher, rec_rintt);
      vecHist2D[Var::RMvtx][Type::Weird][Comp::VsNumTpc] -> Fill(rec_ntpclust_trkmatcher, rec_rmaps);
      vecHist2D[Var::RTpc][Type::Weird][Comp::VsNumTpc]  -> Fill(rec_ntpclust_trkmatcher, rec_rtpc);
      vecHist2D[Var::Phi][Type::Weird][Comp::VsNumTpc]   -> Fill(rec_ntpclust_trkmatcher, rec_phi);
      vecHist2D[Var::Eta][Type::Weird][Comp::VsNumTpc]   -> Fill(rec_ntpclust_trkmatcher, rec_eta);
      vecHist2D[Var::Pt][Type::Weird][Comp::VsNumTpc]    -> Fill(rec_ntpclust_trkmatcher, rec_pt);
      vecHist2D[Var::Frac][Type::Weird][Comp::VsNumTpc]  -> Fill(rec_ntpclust_trkmatcher, rec_ptfrac);
    }
  }  // end reco track loop

  // announce method end and return
  cout << "        Finished reconstructed track loop.\n"
       << "      Finished getting histograms from new matcher track tuple."
       << endl;
  return;

}  // end 'GetNewTupleHists()'



void STrackMatcherComparator::GetOldTupleHists() {

  // announce start of method
  cout << "      Grabbing old matcher tuple histograms:" << endl;

  // declare input leaves
  float tru_event;
  float tru_seed;
  float tru_gntracks;
  float tru_gnchghad;
  float tru_gtrackID;
  float tru_gflavor;
  float tru_gnhits;
  float tru_gnmaps;
  float tru_gnintt;
  float tru_gnmms;
  float tru_gnintt1;
  float tru_gnintt2;
  float tru_gnintt3;
  float tru_gnintt4;
  float tru_gnintt5;
  float tru_gnintt6;
  float tru_gnintt7;
  float tru_gnintt8;
  float tru_gntpc;
  float tru_gnlmaps;
  float tru_gnlintt;
  float tru_gnltpc;
  float tru_gnlmms;
  float tru_gpx;
  float tru_gpy;
  float tru_gpz;
  float tru_gpt;
  float tru_geta;
  float tru_gphi;
  float tru_gvx;
  float tru_gvy;
  float tru_gvz;
  float tru_gvt;
  float tru_gfpx;
  float tru_gfpy;
  float tru_gfpz;
  float tru_gfx;
  float tru_gfy;
  float tru_gfz;
  float tru_gembed;
  float tru_gprimary;
  float tru_trackID;
  float tru_px;
  float tru_py;
  float tru_pz;
  float tru_pt;
  float tru_eta;
  float tru_phi;
  float tru_deltapt;
  float tru_deltaeta;
  float tru_deltaphi;
  float tru_siqr;
  float tru_siphi;
  float tru_sithe;
  float tru_six0;
  float tru_siy0;
  float tru_tpqr;
  float tru_tpphi;
  float tru_tpthe;
  float tru_tpx0;
  float tru_tpy0;
  float tru_charge;
  float tru_quality;
  float tru_chisq;
  float tru_ndf;
  float tru_nhits;
  float tru_layers;
  float tru_nmaps;
  float tru_nintt;
  float tru_ntpc;
  float tru_nmms;
  float tru_ntpc1;
  float tru_ntpc11;
  float tru_ntpc2;
  float tru_ntpc3;
  float tru_nlmaps;
  float tru_nlintt;
  float tru_nltpc;
  float tru_nlmms;
  float tru_vertexID;
  float tru_vx;
  float tru_vy;
  float tru_vz;
  float tru_dca2d;
  float tru_dca2dsigma;
  float tru_dca3dxy;
  float tru_dca3dxysigma;
  float tru_dca3dz;
  float tru_dca3dzsigma;
  float tru_pcax;
  float tru_pcay;
  float tru_pcaz;
  float tru_nfromtruth;
  float tru_nwrong;
  float tru_ntrumaps;
  float tru_nwrongmaps;
  float tru_ntruintt;
  float tru_nwrongintt;
  float tru_ntrutpc;
  float tru_nwrongtpc;
  float tru_ntrumms;
  float tru_nwrongmms;
  float tru_ntrutpc1;
  float tru_nwrongtpc1;
  float tru_ntrutpc11;
  float tru_nwrongtpc11;
  float tru_ntrutpc2;
  float tru_nwrongtpc2;
  float tru_ntrutpc3;
  float tru_nwrongtpc3;
  float tru_layersfromtruth;
  float tru_npedge;
  float tru_nredge;
  float tru_nbig;
  float tru_novlp;
  float tru_merr;
  float tru_msize;
  float tru_nhittpcall;
  float tru_nhittpcin;
  float tru_nhittpcmid;
  float tru_nhittpcout;
  float tru_nclusall;
  float tru_nclustpc;
  float tru_nclusintt;
  float tru_nclusmaps;
  float tru_nclusmms;

  float rec_event;
  float rec_seed;
  float rec_trackID;
  float rec_crossing;
  float rec_px;
  float rec_py;
  float rec_pz;
  float rec_pt;
  float rec_eta;
  float rec_phi;
  float rec_deltapt;
  float rec_deltaeta;
  float rec_deltaphi;
  float rec_siqr;
  float rec_siphi;
  float rec_sithe;
  float rec_six0;
  float rec_siy0;
  float rec_tpqr;
  float rec_tpphi;
  float rec_tpthe;
  float rec_tpx0;
  float rec_tpy0;
  float rec_charge;
  float rec_quality;
  float rec_chisq;
  float rec_ndf;
  float rec_nhits;
  float rec_nmaps;
  float rec_nintt;
  float rec_ntpc;
  float rec_nmms;
  float rec_ntpc1;
  float rec_ntpc11;
  float rec_ntpc2;
  float rec_ntpc3;
  float rec_nlmaps;
  float rec_nlintt;
  float rec_nltpc;
  float rec_nlmms;
  float rec_layers;
  float rec_vertexID;
  float rec_vx;
  float rec_vy;
  float rec_vz;
  float rec_dca2d;
  float rec_dca2dsigma;
  float rec_dca3dxy;
  float rec_dca3dxysigma;
  float rec_dca3dz;
  float rec_dca3dzsigma;
  float rec_pcax;
  float rec_pcay;
  float rec_pcaz;
  float rec_gtrackID;
  float rec_gflavor;
  float rec_gnhits;
  float rec_gnmaps;
  float rec_gnintt;
  float rec_gntpc;
  float rec_gnmms;
  float rec_gnlmaps;
  float rec_gnlintt;
  float rec_gnltpc;
  float rec_gnlmms;
  float rec_gpx;
  float rec_gpy;
  float rec_gpz;
  float rec_gpt;
  float rec_geta;
  float rec_gphi;
  float rec_gvx;
  float rec_gvy;
  float rec_gvz;
  float rec_gvt;
  float rec_gfpx;
  float rec_gfpy;
  float rec_gfpz;
  float rec_gfx;
  float rec_gfy;
  float rec_gfz;
  float rec_gembed;
  float rec_gprimary;
  float rec_nfromtruth;
  float rec_nwrong;
  float rec_ntrumaps;
  float rec_nwrongmaps;
  float rec_ntruintt;
  float rec_nwrongintt;
  float rec_ntrutpc;
  float rec_nwrongtpc;
  float rec_ntrumms;
  float rec_nwrongmms;
  float rec_ntrutpc1;
  float rec_nwrongtpc1;
  float rec_ntrutpc11;
  float rec_nwrongtpc11;
  float rec_ntrutpc2;
  float rec_nwrongtpc2;
  float rec_ntrutpc3;
  float rec_nwrongtpc3;
  float rec_layersfromtruth;
  float rec_npedge;
  float rec_nredge;
  float rec_nbig;
  float rec_novlp;
  float rec_merr;
  float rec_msize;
  float rec_nhittpcall;
  float rec_nhittpcin;
  float rec_nhittpcmid;
  float rec_nhittpcout;
  float rec_nclusall;
  float rec_nclustpc;
  float rec_nclusintt;
  float rec_nclusmaps;
  float rec_nclusmms;

  // Set branch addresses.
  tInTrue -> SetBranchAddress("event",           &tru_event);
  tInTrue -> SetBranchAddress("seed",            &tru_seed);
  tInTrue -> SetBranchAddress("gntracks",        &tru_gntracks);
  tInTrue -> SetBranchAddress("gnchghad",        &tru_gnchghad);
  tInTrue -> SetBranchAddress("gtrackID",        &tru_gtrackID);
  tInTrue -> SetBranchAddress("gflavor",         &tru_gflavor);
  tInTrue -> SetBranchAddress("gnhits",          &tru_gnhits);
  tInTrue -> SetBranchAddress("gnmaps",          &tru_gnmaps);
  tInTrue -> SetBranchAddress("gnintt",          &tru_gnintt);
  tInTrue -> SetBranchAddress("gnmms",           &tru_gnmms);
  tInTrue -> SetBranchAddress("gnintt1",         &tru_gnintt1);
  tInTrue -> SetBranchAddress("gnintt2",         &tru_gnintt2);
  tInTrue -> SetBranchAddress("gnintt3",         &tru_gnintt3);
  tInTrue -> SetBranchAddress("gnintt4",         &tru_gnintt4);
  tInTrue -> SetBranchAddress("gnintt5",         &tru_gnintt5);
  tInTrue -> SetBranchAddress("gnintt6",         &tru_gnintt6);
  tInTrue -> SetBranchAddress("gnintt7",         &tru_gnintt7);
  tInTrue -> SetBranchAddress("gnintt8",         &tru_gnintt8);
  tInTrue -> SetBranchAddress("gntpc",           &tru_gntpc);
  tInTrue -> SetBranchAddress("gnlmaps",         &tru_gnlmaps);
  tInTrue -> SetBranchAddress("gnlintt",         &tru_gnlintt);
  tInTrue -> SetBranchAddress("gnltpc",          &tru_gnltpc);
  tInTrue -> SetBranchAddress("gnlmms",          &tru_gnlmms);
  tInTrue -> SetBranchAddress("gpx",             &tru_gpx);
  tInTrue -> SetBranchAddress("gpy",             &tru_gpy);
  tInTrue -> SetBranchAddress("gpz",             &tru_gpz);
  tInTrue -> SetBranchAddress("gpt",             &tru_gpt);
  tInTrue -> SetBranchAddress("geta",            &tru_geta);
  tInTrue -> SetBranchAddress("gphi",            &tru_gphi);
  tInTrue -> SetBranchAddress("gvx",             &tru_gvx);
  tInTrue -> SetBranchAddress("gvy",             &tru_gvy);
  tInTrue -> SetBranchAddress("gvz",             &tru_gvz);
  tInTrue -> SetBranchAddress("gvt",             &tru_gvt);
  tInTrue -> SetBranchAddress("gfpx",            &tru_gfpx);
  tInTrue -> SetBranchAddress("gfpy",            &tru_gfpy);
  tInTrue -> SetBranchAddress("gfpz",            &tru_gfpz);
  tInTrue -> SetBranchAddress("gfx",             &tru_gfx);
  tInTrue -> SetBranchAddress("gfy",             &tru_gfy);
  tInTrue -> SetBranchAddress("gfz",             &tru_gfz);
  tInTrue -> SetBranchAddress("gembed",          &tru_gembed);
  tInTrue -> SetBranchAddress("gprimary",        &tru_gprimary);
  tInTrue -> SetBranchAddress("trackID",         &tru_trackID);
  tInTrue -> SetBranchAddress("px",              &tru_px);
  tInTrue -> SetBranchAddress("py",              &tru_py);
  tInTrue -> SetBranchAddress("pz",              &tru_pz);
  tInTrue -> SetBranchAddress("pt",              &tru_pt);
  tInTrue -> SetBranchAddress("eta",             &tru_eta);
  tInTrue -> SetBranchAddress("phi",             &tru_phi);
  tInTrue -> SetBranchAddress("deltapt",         &tru_deltapt);
  tInTrue -> SetBranchAddress("deltaeta",        &tru_deltaeta);
  tInTrue -> SetBranchAddress("deltaphi",        &tru_deltaphi);
  tInTrue -> SetBranchAddress("siqr",            &tru_siqr);
  tInTrue -> SetBranchAddress("siphi",           &tru_siphi);
  tInTrue -> SetBranchAddress("sithe",           &tru_sithe);
  tInTrue -> SetBranchAddress("six0",            &tru_six0);
  tInTrue -> SetBranchAddress("siy0",            &tru_siy0);
  tInTrue -> SetBranchAddress("tpqr",            &tru_tpqr);
  tInTrue -> SetBranchAddress("tpphi",           &tru_tpphi);
  tInTrue -> SetBranchAddress("tpthe",           &tru_tpthe);
  tInTrue -> SetBranchAddress("tpx0",            &tru_tpx0);
  tInTrue -> SetBranchAddress("tpy0",            &tru_tpy0);
  tInTrue -> SetBranchAddress("charge",          &tru_charge);
  tInTrue -> SetBranchAddress("quality",         &tru_quality);
  tInTrue -> SetBranchAddress("chisq",           &tru_chisq);
  tInTrue -> SetBranchAddress("ndf",             &tru_ndf);
  tInTrue -> SetBranchAddress("nhits",           &tru_nhits);
  tInTrue -> SetBranchAddress("layers",          &tru_layers);
  tInTrue -> SetBranchAddress("nmaps",           &tru_nmaps);
  tInTrue -> SetBranchAddress("nintt",           &tru_nintt);
  tInTrue -> SetBranchAddress("ntpc",            &tru_ntpc);
  tInTrue -> SetBranchAddress("nmms",            &tru_nmms);
  tInTrue -> SetBranchAddress("ntpc1",           &tru_ntpc1);
  tInTrue -> SetBranchAddress("ntpc11",          &tru_ntpc11);
  tInTrue -> SetBranchAddress("ntpc2",           &tru_ntpc2);
  tInTrue -> SetBranchAddress("ntpc3",           &tru_ntpc3);
  tInTrue -> SetBranchAddress("nlmaps",          &tru_nlmaps);
  tInTrue -> SetBranchAddress("nlintt",          &tru_nlintt);
  tInTrue -> SetBranchAddress("nltpc",           &tru_nltpc);
  tInTrue -> SetBranchAddress("nlmms",           &tru_nlmms);
  tInTrue -> SetBranchAddress("vertexID",        &tru_vertexID);
  tInTrue -> SetBranchAddress("vx",              &tru_vx);
  tInTrue -> SetBranchAddress("vy",              &tru_vy);
  tInTrue -> SetBranchAddress("vz",              &tru_vz);
  tInTrue -> SetBranchAddress("dca2d",           &tru_dca2d);
  tInTrue -> SetBranchAddress("dca2dsigma",      &tru_dca2dsigma);
  tInTrue -> SetBranchAddress("dca3dxy",         &tru_dca3dxy);
  tInTrue -> SetBranchAddress("dca3dxysigma",    &tru_dca3dxysigma);
  tInTrue -> SetBranchAddress("dca3dz",          &tru_dca3dz);
  tInTrue -> SetBranchAddress("dca3dzsigma",     &tru_dca3dzsigma);
  tInTrue -> SetBranchAddress("pcax",            &tru_pcax);
  tInTrue -> SetBranchAddress("pcay",            &tru_pcay);
  tInTrue -> SetBranchAddress("pcaz",            &tru_pcaz);
  tInTrue -> SetBranchAddress("nfromtruth",      &tru_nfromtruth);
  tInTrue -> SetBranchAddress("nwrong",          &tru_nwrong);
  tInTrue -> SetBranchAddress("ntrumaps",        &tru_ntrumaps);
  tInTrue -> SetBranchAddress("nwrongmaps",      &tru_nwrongmaps);
  tInTrue -> SetBranchAddress("ntruintt",        &tru_ntruintt);
  tInTrue -> SetBranchAddress("nwrongintt",      &tru_nwrongintt);
  tInTrue -> SetBranchAddress("ntrutpc",         &tru_ntrutpc);
  tInTrue -> SetBranchAddress("nwrongtpc",       &tru_nwrongtpc);
  tInTrue -> SetBranchAddress("ntrumms",         &tru_ntrumms);
  tInTrue -> SetBranchAddress("nwrongmms",       &tru_nwrongmms);
  tInTrue -> SetBranchAddress("ntrutpc1",        &tru_ntrutpc1);
  tInTrue -> SetBranchAddress("nwrongtpc1",      &tru_nwrongtpc1);
  tInTrue -> SetBranchAddress("ntrutpc11",       &tru_ntrutpc11);
  tInTrue -> SetBranchAddress("nwrongtpc11",     &tru_nwrongtpc11);
  tInTrue -> SetBranchAddress("ntrutpc2",        &tru_ntrutpc2);
  tInTrue -> SetBranchAddress("nwrongtpc2",      &tru_nwrongtpc2);
  tInTrue -> SetBranchAddress("ntrutpc3",        &tru_ntrutpc3);
  tInTrue -> SetBranchAddress("nwrongtpc3",      &tru_nwrongtpc3);
  tInTrue -> SetBranchAddress("layersfromtruth", &tru_layersfromtruth);
  tInTrue -> SetBranchAddress("npedge",          &tru_npedge);
  tInTrue -> SetBranchAddress("nredge",          &tru_nredge);
  tInTrue -> SetBranchAddress("nbig",            &tru_nbig);
  tInTrue -> SetBranchAddress("novlp",           &tru_novlp);
  tInTrue -> SetBranchAddress("merr",            &tru_merr);
  tInTrue -> SetBranchAddress("msize",           &tru_msize);
  tInTrue -> SetBranchAddress("nhittpcall",      &tru_nhittpcall);
  tInTrue -> SetBranchAddress("nhittpcin",       &tru_nhittpcin);
  tInTrue -> SetBranchAddress("nhittpcmid",      &tru_nhittpcmid);
  tInTrue -> SetBranchAddress("nhittpcout",      &tru_nhittpcout);
  tInTrue -> SetBranchAddress("nclusall",        &tru_nclusall);
  tInTrue -> SetBranchAddress("nclustpc",        &tru_nclustpc);
  tInTrue -> SetBranchAddress("nclusintt",       &tru_nclusintt);
  tInTrue -> SetBranchAddress("nclusmaps",       &tru_nclusmaps);
  tInTrue -> SetBranchAddress("nclusmms",        &tru_nclusmms);

  tInReco -> SetBranchAddress("event",           &rec_event);
  tInReco -> SetBranchAddress("seed",            &rec_seed);
  tInReco -> SetBranchAddress("trackID",         &rec_trackID);
  tInReco -> SetBranchAddress("crossing",        &rec_crossing);
  tInReco -> SetBranchAddress("px",              &rec_px);
  tInReco -> SetBranchAddress("py",              &rec_py);
  tInReco -> SetBranchAddress("pz",              &rec_pz);
  tInReco -> SetBranchAddress("pt",              &rec_pt);
  tInReco -> SetBranchAddress("eta",             &rec_eta);
  tInReco -> SetBranchAddress("phi",             &rec_phi);
  tInReco -> SetBranchAddress("deltapt",         &rec_deltapt);
  tInReco -> SetBranchAddress("deltaeta",        &rec_deltaeta);
  tInReco -> SetBranchAddress("deltaphi",        &rec_deltaphi);
  tInReco -> SetBranchAddress("siqr",            &rec_siqr);
  tInReco -> SetBranchAddress("siphi",           &rec_siphi);
  tInReco -> SetBranchAddress("sithe",           &rec_sithe);
  tInReco -> SetBranchAddress("six0",            &rec_six0);
  tInReco -> SetBranchAddress("siy0",            &rec_siy0);
  tInReco -> SetBranchAddress("tpqr",            &rec_tpqr);
  tInReco -> SetBranchAddress("tpphi",           &rec_tpphi);
  tInReco -> SetBranchAddress("tpthe",           &rec_tpthe);
  tInReco -> SetBranchAddress("tpx0",            &rec_tpx0);
  tInReco -> SetBranchAddress("tpy0",            &rec_tpy0);
  tInReco -> SetBranchAddress("charge",          &rec_charge);
  tInReco -> SetBranchAddress("quality",         &rec_quality);
  tInReco -> SetBranchAddress("chisq",           &rec_chisq);
  tInReco -> SetBranchAddress("ndf",             &rec_ndf);
  tInReco -> SetBranchAddress("nhits",           &rec_nhits);
  tInReco -> SetBranchAddress("nmaps",           &rec_nmaps);
  tInReco -> SetBranchAddress("nintt",           &rec_nintt);
  tInReco -> SetBranchAddress("ntpc",            &rec_ntpc);
  tInReco -> SetBranchAddress("nmms",            &rec_nmms);
  tInReco -> SetBranchAddress("ntpc1",           &rec_ntpc1);
  tInReco -> SetBranchAddress("ntpc11",          &rec_ntpc11);
  tInReco -> SetBranchAddress("ntpc2",           &rec_ntpc2);
  tInReco -> SetBranchAddress("ntpc3",           &rec_ntpc3);
  tInReco -> SetBranchAddress("nlmaps",          &rec_nlmaps);
  tInReco -> SetBranchAddress("nlintt",          &rec_nlintt);
  tInReco -> SetBranchAddress("nltpc",           &rec_nltpc);
  tInReco -> SetBranchAddress("nlmms",           &rec_nlmms);
  tInReco -> SetBranchAddress("layers",          &rec_layers);
  tInReco -> SetBranchAddress("vertexID",        &rec_vertexID);
  tInReco -> SetBranchAddress("vx",              &rec_vx);
  tInReco -> SetBranchAddress("vy",              &rec_vy);
  tInReco -> SetBranchAddress("vz",              &rec_vz);
  tInReco -> SetBranchAddress("dca2d",           &rec_dca2d);
  tInReco -> SetBranchAddress("dca2dsigma",      &rec_dca2dsigma);
  tInReco -> SetBranchAddress("dca3dxy",         &rec_dca3dxy);
  tInReco -> SetBranchAddress("dca3dxysigma",    &rec_dca3dxysigma);
  tInReco -> SetBranchAddress("dca3dz",          &rec_dca3dz);
  tInReco -> SetBranchAddress("dca3dzsigma",     &rec_dca3dzsigma);
  tInReco -> SetBranchAddress("pcax",            &rec_pcax);
  tInReco -> SetBranchAddress("pcay",            &rec_pcay);
  tInReco -> SetBranchAddress("pcaz",            &rec_pcaz);
  tInReco -> SetBranchAddress("gtrackID",        &rec_gtrackID);
  tInReco -> SetBranchAddress("gflavor",         &rec_gflavor);
  tInReco -> SetBranchAddress("gnhits",          &rec_gnhits);
  tInReco -> SetBranchAddress("gnmaps",          &rec_gnmaps);
  tInReco -> SetBranchAddress("gnintt",          &rec_gnintt);
  tInReco -> SetBranchAddress("gntpc",           &rec_gntpc);
  tInReco -> SetBranchAddress("gnmms",           &rec_gnmms);
  tInReco -> SetBranchAddress("gnlmaps",         &rec_gnlmaps);
  tInReco -> SetBranchAddress("gnlintt",         &rec_gnlintt);
  tInReco -> SetBranchAddress("gnltpc",          &rec_gnltpc);
  tInReco -> SetBranchAddress("gnlmms",          &rec_gnlmms);
  tInReco -> SetBranchAddress("gpx",             &rec_gpx);
  tInReco -> SetBranchAddress("gpy",             &rec_gpy);
  tInReco -> SetBranchAddress("gpz",             &rec_gpz);
  tInReco -> SetBranchAddress("gpt",             &rec_gpt);
  tInReco -> SetBranchAddress("geta",            &rec_geta);
  tInReco -> SetBranchAddress("gphi",            &rec_gphi);
  tInReco -> SetBranchAddress("gvx",             &rec_gvx);
  tInReco -> SetBranchAddress("gvy",             &rec_gvy);
  tInReco -> SetBranchAddress("gvz",             &rec_gvz);
  tInReco -> SetBranchAddress("gvt",             &rec_gvt);
  tInReco -> SetBranchAddress("gfpx",            &rec_gfpx);
  tInReco -> SetBranchAddress("gfpy",            &rec_gfpy);
  tInReco -> SetBranchAddress("gfpz",            &rec_gfpz);
  tInReco -> SetBranchAddress("gfx",             &rec_gfx);
  tInReco -> SetBranchAddress("gfy",             &rec_gfy);
  tInReco -> SetBranchAddress("gfz",             &rec_gfz);
  tInReco -> SetBranchAddress("gembed",          &rec_gembed);
  tInReco -> SetBranchAddress("gprimary",        &rec_gprimary);
  tInReco -> SetBranchAddress("nfromtruth",      &rec_nfromtruth);
  tInReco -> SetBranchAddress("nwrong",          &rec_nwrong);
  tInReco -> SetBranchAddress("ntrumaps",        &rec_ntrumaps);
  tInReco -> SetBranchAddress("nwrongmaps",      &rec_nwrongmaps);
  tInReco -> SetBranchAddress("ntruintt",        &rec_ntruintt);
  tInReco -> SetBranchAddress("nwrongintt",      &rec_nwrongintt);
  tInReco -> SetBranchAddress("ntrutpc",         &rec_ntrutpc);
  tInReco -> SetBranchAddress("nwrongtpc",       &rec_nwrongtpc);
  tInReco -> SetBranchAddress("ntrumms",         &rec_ntrumms);
  tInReco -> SetBranchAddress("nwrongmms",       &rec_nwrongmms);
  tInReco -> SetBranchAddress("ntrutpc1",        &rec_ntrutpc1);
  tInReco -> SetBranchAddress("nwrongtpc1",      &rec_nwrongtpc1);
  tInReco -> SetBranchAddress("ntrutpc11",       &rec_ntrutpc11);
  tInReco -> SetBranchAddress("nwrongtpc11",     &rec_nwrongtpc11);
  tInReco -> SetBranchAddress("ntrutpc2",        &rec_ntrutpc2);
  tInReco -> SetBranchAddress("nwrongtpc2",      &rec_nwrongtpc2);
  tInReco -> SetBranchAddress("ntrutpc3",        &rec_ntrutpc3);
  tInReco -> SetBranchAddress("nwrongtpc3",      &rec_nwrongtpc3);
  tInReco -> SetBranchAddress("layersfromtruth", &rec_layersfromtruth);
  tInReco -> SetBranchAddress("npedge",          &rec_npedge);
  tInReco -> SetBranchAddress("nredge",          &rec_nredge);
  tInReco -> SetBranchAddress("nbig",            &rec_nbig);
  tInReco -> SetBranchAddress("novlp",           &rec_novlp);
  tInReco -> SetBranchAddress("merr",            &rec_merr);
  tInReco -> SetBranchAddress("msize",           &rec_msize);
  tInReco -> SetBranchAddress("nhittpcall",      &rec_nhittpcall);
  tInReco -> SetBranchAddress("nhittpcin",       &rec_nhittpcin);
  tInReco -> SetBranchAddress("nhittpcmid",      &rec_nhittpcmid);
  tInReco -> SetBranchAddress("nhittpcout",      &rec_nhittpcout);
  tInReco -> SetBranchAddress("nclusall",        &rec_nclusall);
  tInReco -> SetBranchAddress("nclustpc",        &rec_nclustpc);
  tInReco -> SetBranchAddress("nclusintt",       &rec_nclusintt);
  tInReco -> SetBranchAddress("nclusmaps",       &rec_nclusmaps);
  tInReco -> SetBranchAddress("nclusmms",        &rec_nclusmms);
  cout << "        Set input branches." << endl;

  // grab no. of entries
  const int64_t nTrueEntries = tInTrue -> GetEntries();
  const int64_t nRecoEntries = tInReco -> GetEntries(); 
  cout << "        Beginning truth particle loop: " << nTrueEntries << " to process" << endl;

  // loop over truth particles
  int64_t nTrueBytes = 0;
  for (int64_t iTrueEntry = 0; iTrueEntry < nTrueEntries; iTrueEntry++) {

    // grab truth particle entry
    const int64_t trueBytes = tInTrue -> GetEntry(iTrueEntry);
    if (trueBytes < 0) {
      cerr << "PANIC: issue with entry " << iTrueEntry << "! Aborting loop!\n" << endl;
      break;
    } else {
      nTrueBytes += trueBytes;
    }

    const int64_t iTrueProg = iTrueEntry + 1;
    if (iTrueProg == nTrueEntries) {
      cout << "          Processing entry " << iTrueProg << "/" << nTrueEntries << "..." << endl;
    } else {
      cout << "          Processing entry " << iTrueProg << "/" << nTrueEntries << "...\r" << flush;
    }

    // select only primary truth particles
    const bool isPrimary = ((tru_gprimary == 1) && !isnan(tru_trackID));
    if (!isPrimary) continue;

    // run calculations
    const double tru_gntot = tru_gnintt + tru_gnmaps + tru_gntpc;

    // fill truth 1D histograms
    vecHist1D[Var::NTot][Type::Truth]  -> Fill(tru_gntot);
    vecHist1D[Var::NIntt][Type::Truth] -> Fill(tru_gnintt);
    vecHist1D[Var::NMvtx][Type::Truth] -> Fill(tru_gnmaps);
    vecHist1D[Var::NTpc][Type::Truth]  -> Fill(tru_gntpc);
    vecHist1D[Var::RTot][Type::Truth]  -> Fill(1.);
    vecHist1D[Var::RIntt][Type::Truth] -> Fill(1.);
    vecHist1D[Var::RMvtx][Type::Truth] -> Fill(1.);
    vecHist1D[Var::RTpc][Type::Truth]  -> Fill(1.);
    vecHist1D[Var::Phi][Type::Truth]   -> Fill(tru_gphi);
    vecHist1D[Var::Eta][Type::Truth]   -> Fill(tru_geta);
    vecHist1D[Var::Pt][Type::Truth]    -> Fill(tru_gpt);
    vecHist1D[Var::Frac][Type::Truth]  -> Fill(1.);

    // fill truth 2D histograms
    vecHist2D[Var::NTot][Type::Truth][Comp::VsTruthPt]  -> Fill(tru_gpt, tru_gntot);
    vecHist2D[Var::NIntt][Type::Truth][Comp::VsTruthPt] -> Fill(tru_gpt, tru_gnintt);
    vecHist2D[Var::NMvtx][Type::Truth][Comp::VsTruthPt] -> Fill(tru_gpt, tru_gnmaps);
    vecHist2D[Var::NTpc][Type::Truth][Comp::VsTruthPt]  -> Fill(tru_gpt, tru_gntpc);
    vecHist2D[Var::RTot][Type::Truth][Comp::VsTruthPt]  -> Fill(tru_gpt, 1.);
    vecHist2D[Var::RIntt][Type::Truth][Comp::VsTruthPt] -> Fill(tru_gpt, 1.);
    vecHist2D[Var::RMvtx][Type::Truth][Comp::VsTruthPt] -> Fill(tru_gpt, 1.);
    vecHist2D[Var::RTpc][Type::Truth][Comp::VsTruthPt]  -> Fill(tru_gpt, 1.);
    vecHist2D[Var::Phi][Type::Truth][Comp::VsTruthPt]   -> Fill(tru_gpt, tru_gphi);
    vecHist2D[Var::Eta][Type::Truth][Comp::VsTruthPt]   -> Fill(tru_gpt, tru_geta);
    vecHist2D[Var::Pt][Type::Truth][Comp::VsTruthPt]    -> Fill(tru_gpt, tru_gpt);
    vecHist2D[Var::Frac][Type::Truth][Comp::VsTruthPt]  -> Fill(tru_gpt, 1.);

    vecHist2D[Var::NTot][Type::Truth][Comp::VsNumTpc]  -> Fill(tru_gntpc, tru_gntot);
    vecHist2D[Var::NIntt][Type::Truth][Comp::VsNumTpc] -> Fill(tru_gntpc, tru_gnintt);
    vecHist2D[Var::NMvtx][Type::Truth][Comp::VsNumTpc] -> Fill(tru_gntpc, tru_gnmaps);
    vecHist2D[Var::NTpc][Type::Truth][Comp::VsNumTpc]  -> Fill(tru_gntpc, tru_gntpc);
    vecHist2D[Var::RTot][Type::Truth][Comp::VsNumTpc]  -> Fill(tru_gntpc, 1.);
    vecHist2D[Var::RIntt][Type::Truth][Comp::VsNumTpc] -> Fill(tru_gntpc, 1.);
    vecHist2D[Var::RMvtx][Type::Truth][Comp::VsNumTpc] -> Fill(tru_gntpc, 1.);
    vecHist2D[Var::RTpc][Type::Truth][Comp::VsNumTpc]  -> Fill(tru_gntpc, 1.);
    vecHist2D[Var::Phi][Type::Truth][Comp::VsNumTpc]   -> Fill(tru_gntpc, tru_gphi);
    vecHist2D[Var::Eta][Type::Truth][Comp::VsNumTpc]   -> Fill(tru_gntpc, tru_geta);
    vecHist2D[Var::Pt][Type::Truth][Comp::VsNumTpc]    -> Fill(tru_gntpc, tru_gpt);
    vecHist2D[Var::Frac][Type::Truth][Comp::VsNumTpc]  -> Fill(tru_gntpc, 1.);
  }  // end truth particle loop

  // announce next entry loop
  cout << "        Finished truth particle loop.\n"
       << "        Beginning reconstructed track loop: " << nRecoEntries << " to process"
       << endl;

  // loop over reco tracks
  int64_t nRecoBytes = 0;
  for (int64_t iRecoEntry = 0; iRecoEntry < nRecoEntries; iRecoEntry++) {

    // grab reco track entry
    const int64_t recoBytes = tInReco -> GetEntry(iRecoEntry);
    if (recoBytes < 0) {
      cerr << "PANIC: issue with entry " << iRecoEntry << "! Aborting loop!\n" << endl;
      break;
    } else {
      nRecoBytes += recoBytes;
    }

    const int64_t iRecoProg = iRecoEntry + 1;
    if (iRecoProg == nRecoEntries) {
      cout << "          Processing entry " << iRecoProg << "/" << nRecoEntries << "..." << endl;
    } else {
      cout << "          Processing entry " << iRecoProg << "/" << nRecoEntries << "...\r" << flush;
    }

    // skip nan's
    //   TODO also add option to filter out non-primary tracks
    if (isnan(rec_gpt)) continue;

    // run calculations
    const double rec_ntot   = rec_nintt + rec_nmaps + rec_ntpc;
    const double rec_gntot  = rec_gnintt + rec_gnmaps + rec_gntpc;
    const double rec_rtot   = rec_ntot / rec_gntot;
    const double rec_rintt  = rec_nintt / rec_gnintt;
    const double rec_rmaps  = rec_nmaps / rec_gnmaps;
    const double rec_rtpc   = rec_ntpc / rec_gntpc;
    const double rec_ptfrac = rec_pt / rec_gpt;

    // fill all matched reco 1D histograms
    vecHist1D[Var::NTot][Type::Track]  -> Fill(rec_ntot);
    vecHist1D[Var::NIntt][Type::Track] -> Fill(rec_nintt);
    vecHist1D[Var::NMvtx][Type::Track] -> Fill(rec_nmaps);
    vecHist1D[Var::NTpc][Type::Track]  -> Fill(rec_ntpc);
    vecHist1D[Var::RTot][Type::Track]  -> Fill(rec_rtot);
    vecHist1D[Var::RIntt][Type::Track] -> Fill(rec_rintt);
    vecHist1D[Var::RMvtx][Type::Track] -> Fill(rec_rmaps);
    vecHist1D[Var::RTpc][Type::Track]  -> Fill(rec_rtpc);
    vecHist1D[Var::Phi][Type::Track]   -> Fill(rec_phi);
    vecHist1D[Var::Eta][Type::Track]   -> Fill(rec_eta);
    vecHist1D[Var::Pt][Type::Track]    -> Fill(rec_pt);
    vecHist1D[Var::Frac][Type::Track]  -> Fill(rec_ptfrac);

    // fill all matched reco 2D histograms
    vecHist2D[Var::NTot][Type::Track][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_ntot);
    vecHist2D[Var::NIntt][Type::Track][Comp::VsTruthPt] -> Fill(rec_gpt, rec_nintt);
    vecHist2D[Var::NMvtx][Type::Track][Comp::VsTruthPt] -> Fill(rec_gpt, rec_nmaps);
    vecHist2D[Var::NTpc][Type::Track][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_ntpc);
    vecHist2D[Var::RTot][Type::Track][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_rtot);
    vecHist2D[Var::RIntt][Type::Track][Comp::VsTruthPt] -> Fill(rec_gpt, rec_rintt);
    vecHist2D[Var::RMvtx][Type::Track][Comp::VsTruthPt] -> Fill(rec_gpt, rec_rmaps);
    vecHist2D[Var::RTpc][Type::Track][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_rtpc);
    vecHist2D[Var::Phi][Type::Track][Comp::VsTruthPt]   -> Fill(rec_gpt, rec_phi);
    vecHist2D[Var::Eta][Type::Track][Comp::VsTruthPt]   -> Fill(rec_gpt, rec_eta);
    vecHist2D[Var::Pt][Type::Track][Comp::VsTruthPt]    -> Fill(rec_gpt, rec_pt);
    vecHist2D[Var::Frac][Type::Track][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_ptfrac);

    vecHist2D[Var::NTot][Type::Track][Comp::VsNumTpc]  -> Fill(rec_gntpc, rec_ntot);
    vecHist2D[Var::NIntt][Type::Track][Comp::VsNumTpc] -> Fill(rec_gntpc, rec_nintt);
    vecHist2D[Var::NMvtx][Type::Track][Comp::VsNumTpc] -> Fill(rec_gntpc, rec_nmaps);
    vecHist2D[Var::NTpc][Type::Track][Comp::VsNumTpc]  -> Fill(rec_gntpc, rec_ntpc);
    vecHist2D[Var::RTot][Type::Track][Comp::VsNumTpc]  -> Fill(rec_gntpc, rec_rtpc);
    vecHist2D[Var::RIntt][Type::Track][Comp::VsNumTpc] -> Fill(rec_gntpc, rec_rintt);
    vecHist2D[Var::RMvtx][Type::Track][Comp::VsNumTpc] -> Fill(rec_gntpc, rec_rmaps);
    vecHist2D[Var::RTpc][Type::Track][Comp::VsNumTpc]  -> Fill(rec_gntpc, rec_rtpc);
    vecHist2D[Var::Phi][Type::Track][Comp::VsNumTpc]   -> Fill(rec_gntpc, rec_phi);
    vecHist2D[Var::Eta][Type::Track][Comp::VsNumTpc]   -> Fill(rec_gntpc, rec_eta);
    vecHist2D[Var::Pt][Type::Track][Comp::VsNumTpc]    -> Fill(rec_gntpc, rec_pt);
    vecHist2D[Var::Frac][Type::Track][Comp::VsNumTpc]  -> Fill(rec_gntpc, rec_ptfrac);

    // fill weird and normal matched reco 1D histograms
    const bool im_config.normalTrack = ((rec_ptfrac >= m_config.oddPtFrac.first) && (rec_ptfrac <= m_config.oddPtFrac.second));
    if (im_config.normalTrack) {
      vecHist1D[Var::NTot][Type::Normal]  -> Fill(rec_ntot);
      vecHist1D[Var::NIntt][Type::Normal] -> Fill(rec_nintt);
      vecHist1D[Var::NMvtx][Type::Normal] -> Fill(rec_nmaps);
      vecHist1D[Var::NTpc][Type::Normal]  -> Fill(rec_ntpc);
      vecHist1D[Var::RTot][Type::Normal]  -> Fill(rec_rtot);
      vecHist1D[Var::RIntt][Type::Normal] -> Fill(rec_rintt);
      vecHist1D[Var::RMvtx][Type::Normal] -> Fill(rec_rmaps);
      vecHist1D[Var::RTpc][Type::Normal]  -> Fill(rec_rtpc);
      vecHist1D[Var::Phi][Type::Normal]   -> Fill(rec_phi);
      vecHist1D[Var::Eta][Type::Normal]   -> Fill(rec_eta);
      vecHist1D[Var::Pt][Type::Normal]    -> Fill(rec_pt);
      vecHist1D[Var::Frac][Type::Normal]  -> Fill(rec_ptfrac);

      vecHist2D[Var::NTot][Type::Normal][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_ntot);
      vecHist2D[Var::NIntt][Type::Normal][Comp::VsTruthPt] -> Fill(rec_gpt, rec_nintt);
      vecHist2D[Var::NMvtx][Type::Normal][Comp::VsTruthPt] -> Fill(rec_gpt, rec_nmaps);
      vecHist2D[Var::NTpc][Type::Normal][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_ntpc);
      vecHist2D[Var::RTot][Type::Normal][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_rtot);
      vecHist2D[Var::RIntt][Type::Normal][Comp::VsTruthPt] -> Fill(rec_gpt, rec_rintt);
      vecHist2D[Var::RMvtx][Type::Normal][Comp::VsTruthPt] -> Fill(rec_gpt, rec_rmaps);
      vecHist2D[Var::RTpc][Type::Normal][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_rtpc);
      vecHist2D[Var::Phi][Type::Normal][Comp::VsTruthPt]   -> Fill(rec_gpt, rec_phi);
      vecHist2D[Var::Eta][Type::Normal][Comp::VsTruthPt]   -> Fill(rec_gpt, rec_eta);
      vecHist2D[Var::Pt][Type::Normal][Comp::VsTruthPt]    -> Fill(rec_gpt, rec_pt);
      vecHist2D[Var::Frac][Type::Normal][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_ptfrac);

      vecHist2D[Var::NTot][Type::Normal][Comp::VsNumTpc]  -> Fill(rec_gntpc, rec_ntot);
      vecHist2D[Var::NIntt][Type::Normal][Comp::VsNumTpc] -> Fill(rec_gntpc, rec_nintt);
      vecHist2D[Var::NMvtx][Type::Normal][Comp::VsNumTpc] -> Fill(rec_gntpc, rec_nmaps);
      vecHist2D[Var::NTpc][Type::Normal][Comp::VsNumTpc]  -> Fill(rec_gntpc, rec_ntpc);
      vecHist2D[Var::RTot][Type::Normal][Comp::VsNumTpc]  -> Fill(rec_gntpc, rec_rtot);
      vecHist2D[Var::RIntt][Type::Normal][Comp::VsNumTpc] -> Fill(rec_gntpc, rec_rintt);
      vecHist2D[Var::RMvtx][Type::Normal][Comp::VsNumTpc] -> Fill(rec_gntpc, rec_rmaps);
      vecHist2D[Var::RTpc][Type::Normal][Comp::VsNumTpc]  -> Fill(rec_gntpc, rec_rtpc);
      vecHist2D[Var::Phi][Type::Normal][Comp::VsNumTpc]   -> Fill(rec_gntpc, rec_phi);
      vecHist2D[Var::Eta][Type::Normal][Comp::VsNumTpc]   -> Fill(rec_gntpc, rec_eta);
      vecHist2D[Var::Pt][Type::Normal][Comp::VsNumTpc]    -> Fill(rec_gntpc, rec_pt);
      vecHist2D[Var::Frac][Type::Normal][Comp::VsNumTpc]  -> Fill(rec_gntpc, rec_ptfrac);
    } else {
      vecHist1D[Var::NTot][Type::Weird]  -> Fill(rec_ntot);
      vecHist1D[Var::NIntt][Type::Weird] -> Fill(rec_nintt);
      vecHist1D[Var::NMvtx][Type::Weird] -> Fill(rec_nmaps);
      vecHist1D[Var::NTpc][Type::Weird]  -> Fill(rec_ntpc);
      vecHist1D[Var::RTot][Type::Weird]  -> Fill(rec_rtot);
      vecHist1D[Var::RIntt][Type::Weird] -> Fill(rec_rintt);
      vecHist1D[Var::RMvtx][Type::Weird] -> Fill(rec_rmaps);
      vecHist1D[Var::RTpc][Type::Weird]  -> Fill(rec_rtpc);
      vecHist1D[Var::Phi][Type::Weird]   -> Fill(rec_phi);
      vecHist1D[Var::Eta][Type::Weird]   -> Fill(rec_eta);
      vecHist1D[Var::Pt][Type::Weird]    -> Fill(rec_pt);
      vecHist1D[Var::Frac][Type::Weird]  -> Fill(rec_ptfrac);

      vecHist2D[Var::NTot][Type::Weird][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_ntot);
      vecHist2D[Var::NIntt][Type::Weird][Comp::VsTruthPt] -> Fill(rec_gpt, rec_nintt);
      vecHist2D[Var::NMvtx][Type::Weird][Comp::VsTruthPt] -> Fill(rec_gpt, rec_nmaps);
      vecHist2D[Var::NTpc][Type::Weird][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_ntpc);
      vecHist2D[Var::RTot][Type::Weird][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_rtot);
      vecHist2D[Var::RIntt][Type::Weird][Comp::VsTruthPt] -> Fill(rec_gpt, rec_rintt);
      vecHist2D[Var::RMvtx][Type::Weird][Comp::VsTruthPt] -> Fill(rec_gpt, rec_rmaps);
      vecHist2D[Var::RTpc][Type::Weird][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_rtpc);
      vecHist2D[Var::Phi][Type::Weird][Comp::VsTruthPt]   -> Fill(rec_gpt, rec_phi);
      vecHist2D[Var::Eta][Type::Weird][Comp::VsTruthPt]   -> Fill(rec_gpt, rec_eta);
      vecHist2D[Var::Pt][Type::Weird][Comp::VsTruthPt]    -> Fill(rec_gpt, rec_pt);
      vecHist2D[Var::Frac][Type::Weird][Comp::VsTruthPt]  -> Fill(rec_gpt, rec_ptfrac);

      vecHist2D[Var::NTot][Type::Weird][Comp::VsNumTpc]  -> Fill(rec_gntpc, rec_ntot);
      vecHist2D[Var::NIntt][Type::Weird][Comp::VsNumTpc] -> Fill(rec_gntpc, rec_nintt);
      vecHist2D[Var::NMvtx][Type::Weird][Comp::VsNumTpc] -> Fill(rec_gntpc, rec_nmaps);
      vecHist2D[Var::NTpc][Type::Weird][Comp::VsNumTpc]  -> Fill(rec_gntpc, rec_ntpc);
      vecHist2D[Var::RTot][Type::Weird][Comp::VsNumTpc]  -> Fill(rec_gntpc, rec_rtpc);
      vecHist2D[Var::RIntt][Type::Weird][Comp::VsNumTpc] -> Fill(rec_gntpc, rec_rintt);
      vecHist2D[Var::RMvtx][Type::Weird][Comp::VsNumTpc] -> Fill(rec_gntpc, rec_rmaps);
      vecHist2D[Var::RTpc][Type::Weird][Comp::VsNumTpc]  -> Fill(rec_gntpc, rec_rtpc);
      vecHist2D[Var::Phi][Type::Weird][Comp::VsNumTpc]   -> Fill(rec_gntpc, rec_phi);
      vecHist2D[Var::Eta][Type::Weird][Comp::VsNumTpc]   -> Fill(rec_gntpc, rec_eta);
      vecHist2D[Var::Pt][Type::Weird][Comp::VsNumTpc]    -> Fill(rec_gntpc, rec_pt);
      vecHist2D[Var::Frac][Type::Weird][Comp::VsNumTpc]  -> Fill(rec_gntpc, rec_ptfrac);
    }
  }  // end reco track loop

  // announce method end and return
  cout << "        Finished reconstructed track loop.\n"
       << "      Finished getting histograms from old matcher track tuple."
       << endl;
  return;

}  // end 'GetOldTupleHists()'



void STrackMatcherComparator::MakeRatiosAndPlots(
  const vector<TH1D*> vecNewHists1D,
  const vector<TH2D*> vecNewHists2D,
  const int iDir,
  const string sLabel
) {

  // announce start of routine
  cout << "      Making ratios and plots:" << endl;

  // histograms to compare
  const vector<string> vecInHists1D = {
    "hTruthNumTot",
    "hTrackNumTot",
    "hWeirdNumTot",
    "hNormNumTot",
    "hTruthNumIntt",
    "hTrackNumIntt",
    "hWeirdNumIntt",
    "hNormNumIntt",
    "hTruthNumMvtx",
    "hTrackNumMvtx",
    "hWeirdNumMvtx",
    "hNormNumMvtx",
    "hTruthNumTpc",
    "hTrackNumTpc",
    "hWeirdNumTpc",
    "hNormNumTpc",
    "hTruthRatTot",
    "hTrackRatTot",
    "hWeirdRatTot",
    "hNormRatTot",
    "hTruthRatIntt",
    "hTrackRatIntt",
    "hWeirdRatIntt",
    "hNormRatIntt",
    "hTruthRatMvtx",
    "hTrackRatMvtx",
    "hWeirdRatMvtx",
    "hNormRatMvtx",
    "hTruthRatTpc",
    "hTrackRatTpc",
    "hWeirdRatTpc",
    "hNormRatTpc",
    "hTruthPhi",
    "hTrackPhi",
    "hWeirdPhi",
    "hNormPhi",
    "hTruthEta",
    "hTrackEta",
    "hWeirdEta",
    "hNormEta",
    "hTruthPt",
    "hTrackPt",
    "hWeirdPt",
    "hNormPt",
    "hTruthFrac",
    "hTrackFrac",
    "hWeirdFrac",
    "hNormFrac"
  };
  const vector<string> vecInHists2D = {
    "hTruthNumTotVsTruthPt",
    "hTruthNumTotVsNumTpc",
    "hTrackNumTotVsTruthPt",
    "hTrackNumTotVsNumTpc",
    "hWeirdNumTotVsTruthPt",
    "hWeirdNumTotVsNumTpc",
    "hNormNumTotVsTruthPt",
    "hNormNumTotVsNumTpc",
    "hTruthNumInttVsTruthPt",
    "hTruthNumInttVsNumTpc",
    "hTrackNumInttVsTruthPt",
    "hTrackNumInttVsNumTpc",
    "hWeirdNumInttVsTruthPt",
    "hWeirdNumInttVsNumTpc",
    "hNormNumInttVsTruthPt",
    "hNormNumInttVsNumTpc",
    "hTruthNumMvtxVsTruthPt",
    "hTruthNumMvtxVsNumTpc",
    "hTrackNumMvtxVsTruthPt",
    "hTrackNumMvtxVsNumTpc",
    "hWeirdNumMvtxVsTruthPt",
    "hWeirdNumMvtxVsNumTpc",
    "hNormNumMvtxVsTruthPt",
    "hNormNumMvtxVsNumTpc",
    "hTruthNumTpcVsTruthPt",
    "hTruthNumTpcVsNumTpc",
    "hTrackNumTpcVsTruthPt",
    "hTrackNumTpcVsNumTpc",
    "hWeirdNumTpcVsTruthPt",
    "hWeirdNumTpcVsNumTpc",
    "hNormNumTpcVsTruthPt",
    "hNormNumTpcVsNumTpc",
    "hTruthRatTotVsTruthPt",
    "hTruthRatTotVsNumTpc",
    "hTrackRatTotVsTruthPt",
    "hTrackRatTotVsNumTpc",
    "hWeirdRatTotVsTruthPt",
    "hWeirdRatTotVsNumTpc",
    "hNormRatTotVsTruthPt",
    "hNormRatTotVsNumTpc",
    "hTruthRatInttVsTruthPt",
    "hTruthRatInttVsNumTpc",
    "hTrackRatInttVsTruthPt",
    "hTrackRatInttVsNumTpc",
    "hWeirdRatInttVsTruthPt",
    "hWeirdRatInttVsNumTpc",
    "hNormRatInttVsTruthPt",
    "hNormRatInttVsNumTpc",
    "hTruthRatMvtxVsTruthPt",
    "hTruthRatMvtxVsNumTpc",
    "hTrackRatMvtxVsTruthPt",
    "hTrackRatMvtxVsNumTpc",
    "hWeirdRatMvtxVsTruthPt",
    "hWeirdRatMvtxVsNumTpc",
    "hNormRatMvtxVsTruthPt",
    "hNormRatMvtxVsNumTpc",
    "hTruthRatTpcVsTruthPt",
    "hTruthRatTpcVsNumTpc",
    "hTrackRatTpcVsTruthPt",
    "hTrackRatTpcVsNumTpc",
    "hWeirdRatTpcVsTruthPt",
    "hWeirdRatTpcVsNumTpc",
    "hNormRatTpcVsTruthPt",
    "hNormRatTpcVsNumTpc",
    "hTruthPhiVsTruthPt",
    "hTruthPhiVsNumTpc",
    "hTrackPhiVsTruthPt",
    "hTrackPhiVsNumTpc",
    "hWeirdPhiVsTruthPt",
    "hWeirdPhiVsNumTpc",
    "hNormPhiVsTruthPt",
    "hNormPhiVsNumTpc",
    "hTruthEtaVsTruthPt",
    "hTruthEtaVsNumTpc",
    "hTrackEtaVsTruthPt",
    "hTrackEtaVsNumTpc",
    "hWeirdEtaVsTruthPt",
    "hWeirdEtaVsNumTpc",
    "hNormEtaVsTruthPt",
    "hNormEtaVsNumTpc",
    "hTruthPtVsTruthPt",
    "hTruthPtVsNumTpc",
    "hTrackPtVsTruthPt",
    "hTrackPtVsNumTpc",
    "hWeirdPtVsTruthPt",
    "hWeirdPtVsNumTpc",
    "hNormPtVsTruthPt",
    "hNormPtVsNumTpc",
    "hTruthFracVsTruthPt",
    "hTruthFracVsNumTpc",
    "hTrackFracVsTruthPt",
    "hTrackFracVsNumTpc",
    "hWeirdFracVsTruthPt",
    "hWeirdFracVsNumTpc",
    "hNormFracVsTruthPt",
    "hNormFracVsNumTpc"
  };

  // canvas names
  vector<string> vecCanvasNames1D = {
    "cTruthNumTot",
    "cTrackNumTot",
    "cWeirdNumTot",
    "cNormNumTot",
    "cTruthNumIntt",
    "cTrackNumIntt",
    "cWeirdNumIntt",
    "cNormNumIntt",
    "cTruthNumMvtx",
    "cTrackNumMvtx",
    "cWeirdNumMvtx",
    "cNormNumMvtx",
    "cTruthNumTpc",
    "cTrackNumTpc",
    "cWeirdNumTpc",
    "cNormNumTpc",
    "cTruthRatTot",
    "cTrackRatTot",
    "cWeirdRatTot",
    "cNormRatTot",
    "cTruthRatIntt",
    "cTrackRatIntt",
    "cWeirdRatIntt",
    "cNormRatIntt",
    "cTruthRatMvtx",
    "cTrackRatMvtx",
    "cWeirdRatMvtx",
    "cNormRatMvtx",
    "cTruthRatTpc",
    "cTrackRatTpc",
    "cWeirdRatTpc",
    "cNormRatTpc",
    "cTruthPhi",
    "cTrackPhi",
    "cWeirdPhi",
    "cNormPhi",
    "cTruthEta",
    "cTrackEta",
    "cWeirdEta",
    "cNormEta",
    "cTruthPt",
    "cTrackPt",
    "cWeirdPt",
    "cNormPt",
    "cTruthFrac",
    "cTrackFrac",
    "cWeirdFrac",
    "cNormFrac"
  };
  vector<string> vecCanvasNames2D = {
    "cTruthNumTotVsTruthPt",
    "cTruthNumTotVsNumTpc",
    "cTrackNumTotVsTruthPt",
    "cTrackNumTotVsNumTpc",
    "cWeirdNumTotVsTruthPt",
    "cWeirdNumTotVsNumTpc",
    "cNormNumTotVsTruthPt",
    "cNormNumTotVsNumTpc",
    "cTruthNumInttVsTruthPt",
    "cTruthNumInttVsNumTpc",
    "cTrackNumInttVsTruthPt",
    "cTrackNumInttVsNumTpc",
    "cWeirdNumInttVsTruthPt",
    "cWeirdNumInttVsNumTpc",
    "cNormNumInttVsTruthPt",
    "cNormNumInttVsNumTpc",
    "cTruthNumMvtxVsTruthPt",
    "cTruthNumMvtxVsNumTpc",
    "cTrackNumMvtxVsTruthPt",
    "cTrackNumMvtxVsNumTpc",
    "cWeirdNumMvtxVsTruthPt",
    "cWeirdNumMvtxVsNumTpc",
    "cNormNumMvtxVsTruthPt",
    "cNormNumMvtxVsNumTpc",
    "cTruthNumTpcVsTruthPt",
    "cTruthNumTpcVsNumTpc",
    "cTrackNumTpcVsTruthPt",
    "cTrackNumTpcVsNumTpc",
    "cWeirdNumTpcVsTruthPt",
    "cWeirdNumTpcVsNumTpc",
    "cNormNumTpcVsTruthPt",
    "cNormNumTpcVsNumTpc",
    "cTruthRatTotVsTruthPt",
    "cTruthRatTotVsNumTpc",
    "cTrackRatTotVsTruthPt",
    "cTrackRatTotVsNumTpc",
    "cWeirdRatTotVsTruthPt",
    "cWeirdRatTotVsNumTpc",
    "cNormRatTotVsTruthPt",
    "cNormRatTotVsNumTpc",
    "cTruthRatInttVsTruthPt",
    "cTruthRatInttVsNumTpc",
    "cTrackRatInttVsTruthPt",
    "cTrackRatInttVsNumTpc",
    "cWeirdRatInttVsTruthPt",
    "cWeirdRatInttVsNumTpc",
    "cNormRatInttVsTruthPt",
    "cNormRatInttVsNumTpc",
    "cTruthRatMvtxVsTruthPt",
    "cTruthRatMvtxVsNumTpc",
    "cTrackRatMvtxVsTruthPt",
    "cTrackRatMvtxVsNumTpc",
    "cWeirdRatMvtxVsTruthPt",
    "cWeirdRatMvtxVsNumTpc",
    "cNormRatMvtxVsTruthPt",
    "cNormRatMvtxVsNumTpc",
    "cTruthRatTpcVsTruthPt",
    "cTruthRatTpcVsNumTpc",
    "cTrackRatTpcVsTruthPt",
    "cTrackRatTpcVsNumTpc",
    "cWeirdRatTpcVsTruthPt",
    "cWeirdRatTpcVsNumTpc",
    "cNormRatTpcVsTruthPt",
    "cNormRatTpcVsNumTpc",
    "cTruthPhiVsTruthPt",
    "cTruthPhiVsNumTpc",
    "cTrackPhiVsTruthPt",
    "cTrackPhiVsNumTpc",
    "cWeirdPhiVsTruthPt",
    "cWeirdPhiVsNumTpc",
    "cNormPhiVsTruthPt",
    "cNormPhiVsNumTpc",
    "cTruthEtaVsTruthPt",
    "cTruthEtaVsNumTpc",
    "cTrackEtaVsTruthPt",
    "cTrackEtaVsNumTpc",
    "cWeirdEtaVsTruthPt",
    "cWeirdEtaVsNumTpc",
    "cNormEtaVsTruthPt",
    "cNormEtaVsNumTpc",
    "cTruthPtVsTruthPt",
    "cTruthPtVsNumTpc",
    "cTrackPtVsTruthPt",
    "cTrackPtVsNumTpc",
    "cWeirdPtVsTruthPt",
    "cWeirdPtVsNumTpc",
    "cNormPtVsTruthPt",
    "cNormPtVsNumTpc",
    "cTruthFracVsTruthPt",
    "cTruthFracVsNumTpc",
    "cTrackFracVsTruthPt",
    "cTrackFracVsNumTpc",
    "cWeirdFracVsTruthPt",
    "cWeirdFracVsNumTpc",
    "cNormFracVsTruthPt",
    "cNormFracVsNumTpc"
  };

  // normalize histograms and set axes' ranges as needed ----------------------

  // normalize by integral if needed
  if (m_config.doIntNorm) {
    for (TH1D* hOldHist1D : vecOldHists1D) {
      const double intOld1D = hOldHist1D -> Integral();
      if (intOld1D > 0.) {
        hOldHist1D -> Scale(1. / intOld1D);
      }
    }
    for (TH1D* hNewHist1D : vecNewHists1D) {
      const double intNew1D = hNewHist1D -> Integral();
      if (intNew1D > 0.) {
        hNewHist1D -> Scale(1. / intNew1D);
      }
    }
    for (TH2D* hOldHist2D : vecOldHists2D) {
      const double intOld2D = hOldHist2D -> Integral();
      if (intOld2D > 0.) {
        hOldHist2D -> Scale(1. / intOld2D);
      }
    }
    for (TH2D* hNewHist2D : vecNewHists2D) {
      const double intNew2D = hNewHist2D -> Integral();
      if (intNew2D > 0.) {
        hNewHist2D -> Scale(1. / intNew2D);
      }
    }
    cout << "    Normalized histograms." << endl;
  }

  // set z-axis ranges if needed
  if (m_config.matchVertScales) {
    size_t nHist2D = vecOldHists2D.size();
    for (size_t iHist2D = 0; iHist2D < nHist2D; iHist2D++) {
      const float oldMin = vecOldHists2D[iHist2D] -> GetMinimum(0.);
      const float oldMax = vecOldHists2D[iHist2D] -> GetMaximum();
      const float newMin = vecNewHists2D[iHist2D] -> GetMinimum(0.);
      const float newMax = vecNewHists2D[iHist2D] -> GetMaximum();
      const float setMin = min(oldMin, newMin);
      const float setMax = max(oldMax, newMax);
      vecOldHists2D[iHist2D] -> GetZaxis() -> SetRangeUser(setMin, setMax);
      vecNewHists2D[iHist2D] -> GetZaxis() -> SetRangeUser(setMin, setMax);
    }
    cout << "    Adjusted z-axis scales to match." << endl;
  }

  // pick relevant count
  string m_config.countUse("");
  if (m_config.doIntNorm) {
    m_config.countUse = m_config.norm;
  } else {
    m_config.countUse = m_config.count;
  }

  // calculate ratios ---------------------------------------------------------

  vector<TH1D*> vecRatios1D(vecOldHists1D.size());
  vector<TH2D*> vecRatios2D(vecOldHists2D.size());
  for (size_t iHist1D = 0; iHist1D < vecOldHists1D.size(); iHist1D++) {

    // make histogram name
    TString sRatioName = vecOldHists1D.at(iHist1D) -> GetName();
    sRatioName.Append("_");
    sRatioName.Append(sLabel.data());

    vecRatios1D.at(iHist1D) = (TH1D*) vecOldHists1D.at(iHist1D) -> Clone();
    vecRatios1D.at(iHist1D) -> SetName(sRatioName.Data());
    vecRatios1D.at(iHist1D) -> Reset("ICES");
    vecRatios1D.at(iHist1D) -> Divide(vecOldHists1D.at(iHist1D), vecNewHists1D.at(iHist1D), 1., 1.);
  }
  for (size_t iHist2D = 0; iHist2D < vecOldHists2D.size(); iHist2D++) {

    // make histogram name
    TString sRatioName = vecOldHists2D.at(iHist2D) -> GetName();
    sRatioName.Append("_");
    sRatioName.Append(sLabel.data());

    vecRatios2D.at(iHist2D) = (TH2D*) vecOldHists2D.at(iHist2D) -> Clone();
    vecRatios2D.at(iHist2D) -> SetName(sRatioName.Data());
    vecRatios2D.at(iHist2D) -> Reset("ICES");
    vecRatios2D.at(iHist2D) -> Divide(vecOldHists2D.at(iHist2D), vecNewHists2D.at(iHist2D), 1., 1.);
  }
  cout << "    Calculated ratios." << endl;

  // set histogram styles -----------------------------------------------------

  // style options
  const uint32_t fTxt(42);
  const uint32_t fAln(12);
  const uint32_t fCnt(1);
  const float    fLabH[m_const.nAxes]   = {0.04,  0.04,  0.03};
  const float    fLabR1[m_const.nAxes]  = {0.074, 0.074, 0.056};
  const float    fLabR2[m_const.nAxes]  = {0.04,  0.04,  0.03};
  const float    fTitH[m_const.nAxes]   = {0.04,  0.04,  0.04};
  const float    fTitR1[m_const.nAxes]  = {0.074, 0.074, 0.056};
  const float    fTitR2[m_const.nAxes]  = {0.04,  0.04,  0.04};
  const float    fOffTH[m_const.nAxes]  = {1.0,   1.3,   1.2};
  const float    fOffTR1[m_const.nAxes] = {0.8,   0.8,   1.0};
  const float    fOffTR2[m_const.nAxes] = {1.0,   1.3,   1.2};
  const float    fOffLH[m_const.nAxes]  = {0.005, 0.005, 0.000};
  const float    fOffLR1[m_const.nAxes] = {0.005, 0.005, 0.000};
  const float    fOffLR2[m_const.nAxes] = {0.005, 0.005, 0.000};

  // set old histogram styles
  for (TH1D* hOldHist1D : vecOldHists1D) {
    hOldHist1D -> SetName(sOldName.Data());
    hOldHist1D -> SetMarkerColor(m_config.fCol.first);
    hOldHist1D -> SetMarkerStyle(m_config.fMar.first);
    hOldHist1D -> SetFillColor(m_config.fCol.first);
    hOldHist1D -> SetFillStyle(m_config.fFil.first);
    hOldHist1D -> SetLineColor(m_config.fCol.first);
    hOldHist1D -> SetLineStyle(m_config.fLin.first);
    hOldHist1D -> SetLineWidth(m_config.fWid.first);
    hOldHist1D -> SetTitle("");
    hOldHist1D -> SetTitleFont(fTxt);
    hOldHist1D -> GetXaxis() -> SetTitleFont(fTxt);
    hOldHist1D -> GetXaxis() -> SetTitleSize(fTitH[0]);
    hOldHist1D -> GetXaxis() -> SetTitleOffset(fOffTH[0]);
    hOldHist1D -> GetXaxis() -> SetLabelFont(fTxt);
    hOldHist1D -> GetXaxis() -> SetLabelSize(fLabH[0]);
    hOldHist1D -> GetXaxis() -> SetLabelOffset(fOffLH[0]);
    hOldHist1D -> GetXaxis() -> CenterTitle(fCnt);
    hOldHist1D -> GetYaxis() -> SetTitle(m_config.countUse.data());
    hOldHist1D -> GetYaxis() -> SetTitleFont(fTxt);
    hOldHist1D -> GetYaxis() -> SetTitleSize(fTitH[1]);
    hOldHist1D -> GetYaxis() -> SetTitleOffset(fOffTH[1]);
    hOldHist1D -> GetYaxis() -> SetLabelFont(fTxt);
    hOldHist1D -> GetYaxis() -> SetLabelSize(fLabH[1]);
    hOldHist1D -> GetYaxis() -> SetLabelOffset(fOffLH[1]);
    hOldHist1D -> GetYaxis() -> CenterTitle(fCnt);
  }
  for (TH2D* hOldHist2D : vecOldHists2D) {
    hOldHist2D -> SetName(sOldName.Data());
    hOldHist2D -> SetMarkerColor(m_config.fCol.first);
    hOldHist2D -> SetMarkerStyle(m_config.fMar.first);
    hOldHist2D -> SetFillColor(m_config.fCol.first);
    hOldHist2D -> SetFillStyle(m_config.fFil.first);
    hOldHist2D -> SetLineColor(m_config.fCol.first);
    hOldHist2D -> SetLineStyle(m_config.fLin.first);
    hOldHist2D -> SetLineWidth(m_config.fWid.first);
    hOldHist2D -> SetTitle(m_config.legOld.data());
    hOldHist2D -> SetTitleFont(fTxt);
    hOldHist2D -> GetXaxis() -> SetTitleFont(fTxt);
    hOldHist2D -> GetXaxis() -> SetTitleSize(fTitH[0]);
    hOldHist2D -> GetXaxis() -> SetTitleOffset(fOffTH[0]);
    hOldHist2D -> GetXaxis() -> SetLabelFont(fTxt);
    hOldHist2D -> GetXaxis() -> SetLabelSize(fLabH[0]);
    hOldHist2D -> GetXaxis() -> SetLabelOffset(fOffLH[0]);
    hOldHist2D -> GetXaxis() -> CenterTitle(fCnt);
    hOldHist2D -> GetYaxis() -> SetTitleFont(fTxt);
    hOldHist2D -> GetYaxis() -> SetTitleSize(fTitH[1]);
    hOldHist2D -> GetYaxis() -> SetTitleOffset(fOffTH[1]);
    hOldHist2D -> GetYaxis() -> SetLabelFont(fTxt);
    hOldHist2D -> GetYaxis() -> SetLabelSize(fLabH[1]);
    hOldHist2D -> GetYaxis() -> SetLabelOffset(fOffLH[1]);
    hOldHist2D -> GetYaxis() -> CenterTitle(fCnt);
    hOldHist2D -> GetZaxis() -> SetTitle(m_config.countUse.data());
    hOldHist2D -> GetZaxis() -> SetTitleFont(fTxt);
    hOldHist2D -> GetZaxis() -> SetTitleSize(fTitH[2]);
    hOldHist2D -> GetZaxis() -> SetTitleOffset(fOffTH[2]);
    hOldHist2D -> GetZaxis() -> SetLabelFont(fTxt);
    hOldHist2D -> GetZaxis() -> SetLabelSize(fLabH[2]);
    hOldHist2D -> GetZaxis() -> SetLabelOffset(fOffLH[2]);
    hOldHist2D -> GetZaxis() -> CenterTitle(fCnt);
  }

  // set new histogram styles
  for (TH1D* hNewHist1D : vecNewHists1D) {
    hNewHist1D -> SetName(sNewName.Data());
    hNewHist1D -> SetMarkerColor(m_config.fCol.second);
    hNewHist1D -> SetMarkerStyle(m_config.fMar.second);
    hNewHist1D -> SetFillColor(m_config.fCol.second);
    hNewHist1D -> SetFillStyle(m_config.fFil.second);
    hNewHist1D -> SetLineColor(m_config.fCol.second);
    hNewHist1D -> SetLineStyle(m_config.fLin.second);
    hNewHist1D -> SetLineWidth(m_config.fWid.second);
    hNewHist1D -> SetTitle("");
    hNewHist1D -> SetTitleFont(fTxt);
    hNewHist1D -> GetXaxis() -> SetTitleFont(fTxt);
    hNewHist1D -> GetXaxis() -> SetTitleSize(fTitH[0]);
    hNewHist1D -> GetXaxis() -> SetTitleOffset(fOffTH[0]);
    hNewHist1D -> GetXaxis() -> SetLabelFont(fTxt);
    hNewHist1D -> GetXaxis() -> SetLabelSize(fLabH[0]);
    hNewHist1D -> GetXaxis() -> SetLabelOffset(fOffLH[0]);
    hNewHist1D -> GetXaxis() -> CenterTitle(fCnt);
    hNewHist1D -> GetYaxis() -> SetTitle(m_config.countUse.data());
    hNewHist1D -> GetYaxis() -> SetTitleFont(fTxt);
    hNewHist1D -> GetYaxis() -> SetTitleSize(fTitH[1]);
    hNewHist1D -> GetYaxis() -> SetTitleOffset(fOffTH[1]);
    hNewHist1D -> GetYaxis() -> SetLabelFont(fTxt);
    hNewHist1D -> GetYaxis() -> SetLabelSize(fLabH[1]);
    hNewHist1D -> GetYaxis() -> SetLabelOffset(fOffLH[1]);
    hNewHist1D -> GetYaxis() -> CenterTitle(fCnt);
  }
  for (TH2D* hNewHist2D : vecNewHists2D) {
    hNewHist2D -> SetName(sNewName.Data());
    hNewHist2D -> SetMarkerColor(m_config.fCol.first);
    hNewHist2D -> SetMarkerStyle(m_config.fMar.first);
    hNewHist2D -> SetFillColor(m_config.fCol.first);
    hNewHist2D -> SetFillStyle(m_config.fFil.first);
    hNewHist2D -> SetLineColor(m_config.fCol.first);
    hNewHist2D -> SetLineStyle(m_config.fLin.first);
    hNewHist2D -> SetLineWidth(m_config.fWid.first);
    hNewHist2D -> SetTitle(m_config.legNew.data());
    hNewHist2D -> SetTitleFont(fTxt);
    hNewHist2D -> GetXaxis() -> SetTitleFont(fTxt);
    hNewHist2D -> GetXaxis() -> SetTitleSize(fTitH[0]);
    hNewHist2D -> GetXaxis() -> SetTitleOffset(fOffTH[0]);
    hNewHist2D -> GetXaxis() -> SetLabelFont(fTxt);
    hNewHist2D -> GetXaxis() -> SetLabelSize(fLabH[0]);
    hNewHist2D -> GetXaxis() -> SetLabelOffset(fOffLH[0]);
    hNewHist2D -> GetXaxis() -> CenterTitle(fCnt);
    hNewHist2D -> GetYaxis() -> SetTitleFont(fTxt);
    hNewHist2D -> GetYaxis() -> SetTitleSize(fTitH[1]);
    hNewHist2D -> GetYaxis() -> SetTitleOffset(fOffTH[1]);
    hNewHist2D -> GetYaxis() -> SetLabelFont(fTxt);
    hNewHist2D -> GetYaxis() -> SetLabelSize(fLabH[1]);
    hNewHist2D -> GetYaxis() -> SetLabelOffset(fOffLH[1]);
    hNewHist2D -> GetYaxis() -> CenterTitle(fCnt);
    hNewHist2D -> GetZaxis() -> SetTitle(m_config.countUse.data());
    hNewHist2D -> GetZaxis() -> SetTitleFont(fTxt);
    hNewHist2D -> GetZaxis() -> SetTitleSize(fTitH[2]);
    hNewHist2D -> GetZaxis() -> SetTitleOffset(fOffTH[2]);
    hNewHist2D -> GetZaxis() -> SetLabelFont(fTxt);
    hNewHist2D -> GetZaxis() -> SetLabelSize(fLabH[2]);
    hNewHist2D -> GetZaxis() -> SetLabelOffset(fOffLH[2]);
    hNewHist2D -> GetZaxis() -> CenterTitle(fCnt);
  }


  // set ratio styles
  for (TH1D* hRatio1D : vecRatios1D) {
    hRatio1D -> SetMarkerColor(m_config.fCol.first);
    hRatio1D -> SetMarkerStyle(m_config.fMar.first);
    hRatio1D -> SetFillColor(m_config.fCol.first);
    hRatio1D -> SetFillStyle(m_config.fFil.first);
    hRatio1D -> SetLineColor(m_config.fCol.first);
    hRatio1D -> SetLineStyle(m_config.fLin.first);
    hRatio1D -> SetLineWidth(m_config.fWid.first);
    hRatio1D -> SetTitle("");
    hRatio1D -> SetTitleFont(fTxt);
    hRatio1D -> GetXaxis() -> SetTitleFont(fTxt);
    hRatio1D -> GetXaxis() -> SetTitleSize(fTitR1[0]);
    hRatio1D -> GetXaxis() -> SetTitleOffset(fOffTR1[0]);
    hRatio1D -> GetXaxis() -> SetLabelFont(fTxt);
    hRatio1D -> GetXaxis() -> SetLabelSize(fLabR1[0]);
    hRatio1D -> GetXaxis() -> SetLabelOffset(fOffLR1[0]);
    hRatio1D -> GetXaxis() -> CenterTitle(fCnt);
    hRatio1D -> GetYaxis() -> SetTitle(sRatio.data());
    hRatio1D -> GetYaxis() -> SetTitleFont(fTxt);
    hRatio1D -> GetYaxis() -> SetTitleSize(fTitR1[1]);
    hRatio1D -> GetYaxis() -> SetTitleOffset(fOffTR1[1]);
    hRatio1D -> GetYaxis() -> SetLabelFont(fTxt);
    hRatio1D -> GetYaxis() -> SetLabelSize(fLabR1[1]);
    hRatio1D -> GetYaxis() -> SetLabelOffset(fOffLR1[1]);
    hRatio1D -> GetYaxis() -> CenterTitle(fCnt);
  }
  for (TH2D* hRatio2D : vecRatios2D) {
    hRatio2D -> SetMarkerColor(m_config.fCol.first);
    hRatio2D -> SetMarkerStyle(m_config.fMar.first);
    hRatio2D -> SetFillColor(m_config.fCol.first);
    hRatio2D -> SetFillStyle(m_config.fFil.first);
    hRatio2D -> SetLineColor(m_config.fCol.first);
    hRatio2D -> SetLineStyle(m_config.fLin.first);
    hRatio2D -> SetLineWidth(m_config.fWid.first);
    hRatio2D -> SetTitle("");
    hRatio2D -> SetTitleFont(fTxt);
    hRatio2D -> GetXaxis() -> SetTitleFont(fTxt);
    hRatio2D -> GetXaxis() -> SetTitleSize(fTitR2[0]);
    hRatio2D -> GetXaxis() -> SetTitleOffset(fOffTR2[0]);
    hRatio2D -> GetXaxis() -> SetLabelFont(fTxt);
    hRatio2D -> GetXaxis() -> SetLabelSize(fLabR2[0]);
    hRatio2D -> GetXaxis() -> SetLabelOffset(fOffLR2[0]);
    hRatio2D -> GetXaxis() -> CenterTitle(fCnt);
    hRatio2D -> GetYaxis() -> SetTitleFont(fTxt);
    hRatio2D -> GetYaxis() -> SetTitleSize(fTitR2[1]);
    hRatio2D -> GetYaxis() -> SetTitleOffset(fOffTR2[1]);
    hRatio2D -> GetYaxis() -> SetLabelFont(fTxt);
    hRatio2D -> GetYaxis() -> SetLabelSize(fLabR2[1]);
    hRatio2D -> GetYaxis() -> SetLabelOffset(fOffLR2[1]);
    hRatio2D -> GetYaxis() -> CenterTitle(fCnt);
    hRatio2D -> GetZaxis() -> SetTitle(sRatio.data());
    hRatio2D -> GetZaxis() -> SetTitleFont(fTxt);
    hRatio2D -> GetZaxis() -> SetTitleSize(fTitR2[2]);
    hRatio2D -> GetZaxis() -> SetTitleOffset(fOffTR2[2]);
    hRatio2D -> GetZaxis() -> SetLabelFont(fTxt);
    hRatio2D -> GetZaxis() -> SetLabelSize(fLabR2[2]);
    hRatio2D -> GetZaxis() -> SetLabelOffset(fOffLR2[2]);
    hRatio2D -> GetZaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

  // make legends and text boxes ----------------------------------------------

  // make legend
  const uint32_t fColLe(0);
  const uint32_t fFilLe(0);
  const uint32_t fLinLe(0);
  const float    fLegXY[m_const.nVtx] = {0.1, 0.1, 0.3, 0.2};

  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], m_config.header.data());
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  leg -> AddEntry(vecOldHists1D.front(), m_config.legOld.data(), "pf");
  leg -> AddEntry(vecNewHists1D.front(), m_config.legNew.data(), "pf");
  cout << "    Made legend." << endl;

  // make text
  const uint32_t fColTx(0);
  const uint32_t fFilTx(0);
  const uint32_t fLinTx(0);
  const float    fTxtXY[m_const.nVtx] = {0.3, 0.1, 0.5, 0.25};

  TPaveText *txt = new TPaveText(fTxtXY[0], fTxtXY[1], fTxtXY[2], fTxtXY[3], "NDC NB");
  txt -> SetFillColor(fColTx);
  txt -> SetFillStyle(fFilTx);
  txt -> SetLineColor(fColTx);
  txt -> SetLineStyle(fLinTx);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAln);
  for (const string txtLine : m_config.info) {
    txt -> AddText(txtLine.data());
  }
  cout << "    Made text." << endl;

  // make plots ---------------------------------------------------------------

  // create output directories
  array<TDirectory*, m_const.nDir> arrOutDirs;
  for (size_t iDir = 0; iDir < m_const.nDir; iDir++) {
    arrOutDirs[iDir] = fOutput -> mkdir(arrDirNames[iDir].data());
  }
  cout << "    Made directories." << endl;

  // canvas parameters
  const uint32_t width1D(750);
  const uint32_t width2D(1500);
  const uint32_t width2DR(2250);
  const uint32_t height(750);
  const uint32_t heightR1(900);
  const uint32_t heightR2(500);
  const uint32_t fMode(0);
  const uint32_t fBord(2);
  const uint32_t fGrid(0);
  const uint32_t fTick(1);
  const uint32_t fLogX(0);
  const uint32_t fLogY(1);
  const uint32_t fLogZ(1);
  const uint32_t fFrame(0);
  const string   sOldPadName("pOld");
  const string   sNewPadName("pNew");
  const string   sHistPadName("pHist");
  const string   sRatPadName("pRatio");
  const float    fMargin1D[m_const.nSide]  = {0.02,  0.02, 0.15,  0.15};
  const float    fMargin1DH[m_const.nSide] = {0.02,  0.02, 0.005, 0.15};
  const float    fMargin1DR[m_const.nSide] = {0.005, 0.02, 0.2,   0.15};
  const float    fMargin2D[m_const.nSide]  = {0.10, 0.15, 0.15, 0.15};
  const float    xyOldPad[m_const.nVtx]    = {0.0,  0.0,  0.5,  1.0};
  const float    xyOldPadR[m_const.nVtx]   = {0.0,  0.0,  0.33, 1.0};
  const float    xyNewPad[m_const.nVtx]    = {0.5,  0.0,  1.0,  1.0};  
  const float    xyNewPadR[m_const.nVtx]   = {0.33, 0.0,  0.66, 1.0};
  const float    xyHistPadR[m_const.nVtx]  = {0.0,  0.33, 1.0,  1.0};
  const float    xyRatPadR1[m_const.nVtx]  = {0.0,  0.0,  1.0,  0.33};
  const float    xyRatPadR2[m_const.nVtx]  = {0.66, 0.0,  1.0,  1.0};

  // make 1D plots
  size_t iDraw1D = 0;
  for (const string sCanvasName : vecCanvasNames1D) {

    // construct canvas
    TCanvas* cPlot1D = new TCanvas(sCanvasName.data(), "", width1D, height);
    cPlot1D -> SetGrid(fGrid, fGrid);
    cPlot1D -> SetTicks(fTick, fTick);
    cPlot1D -> SetBorderMode(fMode);
    cPlot1D -> SetBorderSize(fBord);
    cPlot1D -> SetFrameBorderMode(fFrame);
    cPlot1D -> SetTopMargin(fMargin1D[0]);
    cPlot1D -> SetRightMargin(fMargin1D[1]);
    cPlot1D -> SetBottomMargin(fMargin1D[2]);
    cPlot1D -> SetLeftMargin(fMargin1D[3]);
    cPlot1D -> SetLogx(fLogX);
    cPlot1D -> SetLogy(fLogY);
    cPlot1D -> cd();

    // draw old vs. new histograms
    vecOldHists1D.at(iDraw1D) -> Draw();
    vecNewHists1D.at(iDraw1D) -> Draw("same");

    // draw text and save
    leg                   -> Draw();
    txt                   -> Draw();
    arrOutDirs[Dir::Plot] -> cd();
    cPlot1D               -> Write();
    cPlot1D               -> Close();

    // increment counter
    ++iDraw1D;
  }

  iDraw1D = 0;
  for (const string sCanvasName : vecCanvasNames1D) {

    // make new name
    const string sNameWithRatio = sCanvasName + "_" + sLabel;

    // construct canvas
    TCanvas* cPlot1D = new TCanvas(sNameWithRatio.data(), "", width1D, heightR1);
    TPad*    pPadH1D = new TPad(sHistPadName.data(), "", xyHistPadR[0], xyHistPadR[1], xyHistPadR[2], xyHistPadR[3]);
    TPad*    pPadR1D = new TPad(sRatPadName.data(),  "", xyRatPadR1[0], xyRatPadR1[1], xyRatPadR1[2], xyRatPadR1[3]);
    cPlot1D -> SetGrid(fGrid, fGrid);
    cPlot1D -> SetTicks(fTick, fTick);
    cPlot1D -> SetBorderMode(fMode);
    cPlot1D -> SetBorderSize(fBord);
    pPadH1D -> SetGrid(fGrid, fGrid);
    pPadH1D -> SetTicks(fTick, fTick);
    pPadH1D -> SetBorderMode(fMode);
    pPadH1D -> SetBorderSize(fBord);
    pPadH1D -> SetFrameBorderMode(fFrame);
    pPadH1D -> SetTopMargin(fMargin1DH[0]);
    pPadH1D -> SetRightMargin(fMargin1DH[1]);
    pPadH1D -> SetBottomMargin(fMargin1DH[2]);
    pPadH1D -> SetLeftMargin(fMargin1DH[3]);
    pPadH1D -> SetLogx(fLogX);
    pPadH1D -> SetLogy(fLogY);
    pPadR1D -> SetGrid(fGrid, fGrid);
    pPadR1D -> SetTicks(fTick, fTick);
    pPadR1D -> SetBorderMode(fMode);
    pPadR1D -> SetBorderSize(fBord);
    pPadR1D -> SetFrameBorderMode(fFrame);
    pPadR1D -> SetTopMargin(fMargin1DR[0]);
    pPadR1D -> SetRightMargin(fMargin1DR[1]);
    pPadR1D -> SetBottomMargin(fMargin1DR[2]);
    pPadR1D -> SetLeftMargin(fMargin1DR[3]);
    pPadR1D -> SetLogx(fLogX);
    pPadR1D -> SetLogy(fLogY);
    cPlot1D -> cd();
    pPadH1D -> Draw();
    pPadR1D -> Draw();

    // draw old vs. new histograms
    pPadH1D                   -> cd();
    vecOldHists1D.at(iDraw1D) -> Draw();
    vecNewHists1D.at(iDraw1D) -> Draw("same");
    leg                       -> Draw();
    txt                       -> Draw();

    // draw ratio
    pPadR1D                 -> cd();
    vecRatios1D.at(iDraw1D) -> Draw();

    // draw text and save
    arrOutDirs[Dir::Plot] -> cd();
    cPlot1D               -> Write();
    cPlot1D               -> Close();

    // increment counter
    ++iDraw1D;
  }

  // make 2D plots
  size_t iDraw2D = 0;
  for (const string sCanvasName : vecCanvasNames2D) {

    // construct canvas
    TCanvas* cPlot2D = new TCanvas(sCanvasName.data(), "", width2D, height);
    TPad*    pOld    = new TPad(sOldPadName.data(), "", xyOldPad[0], xyOldPad[1], xyOldPad[2], xyOldPad[3]);
    TPad*    pNew    = new TPad(sNewPadName.data(), "", xyNewPad[0], xyNewPad[1], xyNewPad[2], xyNewPad[3]);
    cPlot2D -> SetGrid(fGrid, fGrid);
    cPlot2D -> SetTicks(fTick, fTick);
    cPlot2D -> SetBorderMode(fMode);
    cPlot2D -> SetBorderSize(fBord);
    pOld    -> SetGrid(fGrid, fGrid);
    pOld    -> SetTicks(fTick, fTick);
    pOld    -> SetBorderMode(fMode);
    pOld    -> SetBorderSize(fBord);
    pOld    -> SetFrameBorderMode(fFrame);
    pOld    -> SetTopMargin(fMargin2D[0]);
    pOld    -> SetRightMargin(fMargin2D[1]);
    pOld    -> SetBottomMargin(fMargin2D[2]);
    pOld    -> SetLeftMargin(fMargin2D[3]);
    pOld    -> SetLogx(fLogX);
    pOld    -> SetLogy(fLogY);
    pOld    -> SetLogz(fLogZ);
    pNew    -> SetGrid(fGrid, fGrid);
    pNew    -> SetTicks(fTick, fTick);
    pNew    -> SetBorderMode(fMode);
    pNew    -> SetBorderSize(fBord);
    pNew    -> SetFrameBorderMode(fFrame);
    pNew    -> SetTopMargin(fMargin2D[0]);
    pNew    -> SetRightMargin(fMargin2D[1]);
    pNew    -> SetBottomMargin(fMargin2D[2]);
    pNew    -> SetLeftMargin(fMargin2D[3]);
    pNew    -> SetLogx(fLogX);
    pNew    -> SetLogy(fLogY);
    pNew    -> SetLogz(fLogZ);
    cPlot2D -> cd();
    pOld    -> Draw();
    pNew    -> Draw();

    // draw old vs. new histograms
    pOld                      -> cd();
    vecOldHists2D.at(iDraw2D) -> Draw("colz");
    pNew                      -> cd();
    vecNewHists2D.at(iDraw2D) -> Draw("colz");

    // draw text and save
    pNew                  -> cd();
    txt                   -> Draw();
    arrOutDirs[Dir::Plot] -> cd();
    cPlot2D               -> Write();
    cPlot2D               -> Close();

    // increment counter
    ++iDraw2D;
  }

  iDraw2D = 0;
  for (const string sCanvasName : vecCanvasNames2D) {

    // make new name
    const string sNameWithRatio = sCanvasName + "_" + sLabel;

    // construct canvas
    TCanvas* cPlot2D = new TCanvas(sNameWithRatio.data(), "", width2D, heightR2);
    TPad*    pOld    = new TPad(sOldPadName.data(), "", xyOldPadR[0],  xyOldPadR[1],  xyOldPadR[2],  xyOldPadR[3]);
    TPad*    pNew    = new TPad(sNewPadName.data(), "", xyNewPadR[0],  xyNewPadR[1],  xyNewPadR[2],  xyNewPadR[3]);
    TPad*    pRat    = new TPad(sRatPadName.data(), "", xyRatPadR2[0], xyRatPadR2[1], xyRatPadR2[2], xyRatPadR2[3]);
    cPlot2D -> SetGrid(fGrid, fGrid);
    cPlot2D -> SetTicks(fTick, fTick);
    cPlot2D -> SetBorderMode(fMode);
    cPlot2D -> SetBorderSize(fBord);
    pOld    -> SetGrid(fGrid, fGrid);
    pOld    -> SetTicks(fTick, fTick);
    pOld    -> SetBorderMode(fMode);
    pOld    -> SetBorderSize(fBord);
    pOld    -> SetFrameBorderMode(fFrame);
    pOld    -> SetTopMargin(fMargin2D[0]);
    pOld    -> SetRightMargin(fMargin2D[1]);
    pOld    -> SetBottomMargin(fMargin2D[2]);
    pOld    -> SetLeftMargin(fMargin2D[3]);
    pOld    -> SetLogx(fLogX);
    pOld    -> SetLogy(fLogY);
    pOld    -> SetLogz(fLogZ);
    pNew    -> SetGrid(fGrid, fGrid);
    pNew    -> SetTicks(fTick, fTick);
    pNew    -> SetBorderMode(fMode);
    pNew    -> SetBorderSize(fBord);
    pNew    -> SetFrameBorderMode(fFrame);
    pNew    -> SetTopMargin(fMargin2D[0]);
    pNew    -> SetRightMargin(fMargin2D[1]);
    pNew    -> SetBottomMargin(fMargin2D[2]);
    pNew    -> SetLeftMargin(fMargin2D[3]);
    pNew    -> SetLogx(fLogX);
    pNew    -> SetLogy(fLogY);
    pNew    -> SetLogz(fLogZ);
    pRat    -> SetGrid(fGrid, fGrid);
    pRat    -> SetTicks(fTick, fTick);
    pRat    -> SetBorderMode(fMode);
    pRat    -> SetBorderSize(fBord);
    pRat    -> SetFrameBorderMode(fFrame);
    pRat    -> SetTopMargin(fMargin2D[0]);
    pRat    -> SetRightMargin(fMargin2D[1]);
    pRat    -> SetBottomMargin(fMargin2D[2]);
    pRat    -> SetLeftMargin(fMargin2D[3]);
    pRat    -> SetLogx(fLogX);
    pRat    -> SetLogy(fLogY);
    pRat    -> SetLogz(fLogZ);
    cPlot2D -> cd();
    pOld    -> Draw();
    pNew    -> Draw();
    pRat    -> Draw();

    // draw old vs. new vs. ratio histograms
    pOld                      -> cd();
    vecOldHists2D.at(iDraw2D) -> Draw("colz");
    txt                       -> Draw();
    pNew                      -> cd();
    vecNewHists2D.at(iDraw2D) -> Draw("colz");
    pRat                      -> cd();
    vecRatios2D.at(iDraw2D)   -> Draw("colz");

    // draw text and save
    arrOutDirs[Dir::Plot] -> cd();
    cPlot2D               -> Write();
    cPlot2D               -> Close();

    // increment counter
    ++iDraw2D;
  }

  // announce end of routine and return
  cout << "    Made plots.\n"
       << endl; 

  // exit
  cout << "  Finished old vs. new comparison script!\n" << endl;
  return;

}  // end 'MakeRatiosAndPlots(vector<TH1D*>, vector<TH2D*>, int, string)'



void STrackMatcherComparator::SaveHistograms() {

  return;

}  // end 'SaveHistograms()'



void STrackMatcherComparator::CloseInput() {

  return;

}  // end 'CloseInput()'



void STrackMatcherComparator::CloseOutput() {

  return;

}  // end 'CloseOutput()'

// end ------------------------------------------------------------------------
