// ----------------------------------------------------------------------------
// 'STrackMatcherComparator.h'
// Derek Anderson
// 01.30.2024
//
// Small module to produce plots from the output of the
// SvtxEvaluator and FillMatchedClusTree modules.
// ----------------------------------------------------------------------------

#ifndef STRACKMATCHERCOMPARATOR_H
#define STRACKMATCHERCOMPARATOR_H

// c++ utilities
#include <string>
#include <optional>
// root libraries
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TDirectory.h>
// analysis utilities
#include "STrackMatcherComparatorConfig.h"

// make common namespaces implicit
using namespace std;

// forward declarations
class PHCompositeNode;



// STrackMatcherComparator definition -----------------------------------------

class STrackMatcherComparator {

  public:

    // ctor/dtor
    STrackMatcherComparator(optional<STrackMatcherComparatorConfig> config = nullopt);
    ~STrackMatcherComparator();

    // public methods
    void Init();
    void Analyze();
    void End();

     // setters
     void SetConfig(STrackMatcherComparatorConfig& config) {m_config = config;}

     // getters
     STrackMatcherComparatorConfig GetConfig() {return m_config;}

   private:

    // internal methods
    void OpenOutput();
    void OpenInput();
    void InitHists();
    void GetNewTreeHists();
    void GetNewTupleHists();
    void GetOldTupleHists();
    void MakeRatiosAndPlots(const vector<TH1D*> vecNewHists1D, const vector<TH2D*> vecNewHists2D, const int iDir, const string sLabel);
    void SaveOutput();
    void CloseInput();
    void CloseOutput();

    // configuration
    STrackMatcherComparatorConfig m_config;

    // root members
    TFile*   m_outFile      = NULL;
    TFile*   m_treeInFile   = NULL;
    TFile*   m_tupleInFile  = NULL;
    TFile*   m_oldInFile    = NULL;
    TTree*   m_tTreeTrue    = NULL;
    TTree*   m_tTreeReco    = NULL;
    TNtuple* m_ntTupleTruth = NULL;
    TNtuple* m_ntTupleReco  = NULL;
    TNtuple* m_ntOldTruth   = NULL;
    TNtuple* m_ntNewTrack   = NULL;
    TNtuple* m_ntOldTrack   = NULL;

    // output directories
    vector<TDirectory*> m_vecHistDirs;
    vector<TDirectory*> m_vecRatioDirs;
    vector<TDirectory*> m_vecPlotDirs;

    // histogram vectors
    vector<TH1D*> m_vecTreeHists1D;
    vector<TH1D*> m_vecTupleHists1D;
    vector<TH1D*> m_vecOldHists1D;
    vector<TH2D*> m_vecTreeHists2D;
    vector<TH2D*> m_vecTupleHists2D;
    vector<TH2D*> m_vecOldHists2D;

    // class-wide constants
    struct Consts {
      size_t nDir;
      size_t nVtx;
      size_t nSide;
      size_t nAxes;
    } m_const = {4, 4, 4, 3};

    // accessors
    enum Src {
      NewTree,
      NewTuple,
      OldTuple
    };
    enum Var {
      NTot,
      NIntt,
      NMvtx,
      NTpc,
      RTot,
      RIntt,
      RMvtx,
      RTpc,
      Phi,
      Eta,
      Pt,
      Frac,
      Qual,
      PtErr,
      EtaErr,
      PhiErr,
      PtRes,
      EtaRes,
      PhiRes
    };
    enum Type {
      Truth,
      Track,
      Weird,
      Normal
    };
    enum Comp {
      VsTruthPt,
      VsNumTpc
    };

};

#endif

// end ------------------------------------------------------------------------
