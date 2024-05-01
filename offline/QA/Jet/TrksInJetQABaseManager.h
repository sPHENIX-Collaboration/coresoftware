// ----------------------------------------------------------------------------
// 'TrksInJetQABaseManager.h'
// Derek Anderson
// 04.03.2024
//
// Base hist manager submodule for the TrksInJetQA module which
// consolidates methods/data common to all of the hist managers
// ----------------------------------------------------------------------------

#ifndef TRKSINJETQABASEMANAGER_H
#define TRKSINJETQABASEMANAGER_H

// c++ utilities
#include <string>
#include <vector>
#include <utility>
#include <iostream>
// root libraries
#include <TH1.h>
#include <TH2.h>
#include <TDirectory.h>
// phool libraries
#include <phool/phool.h>
// module utilities
#include "TrksInJetQAHist.h"
#include "TrksInJetQATypes.h"
#include "TrksInJetQAConfig.h"



// TrksInJetQABaseManager definition ------------------------------------------

class TrksInJetQABaseManager {

  public:

    // ctor/dtor
    TrksInJetQABaseManager(TrksInJetQAConfig& config, TrksInJetQAHist& hist);
    ~TrksInJetQABaseManager();

    // public methods
    void MakeHistograms(std::string label = "");
    void SaveHistograms(TDirectory* outFile, std::string outDirName);
    void GrabHistograms(std::vector<TH1D*>& vecOutHist1D, std::vector<TH2D*>& vecOutHist2D);

  protected:

    // private methods
    void BuildHistograms(std::string label = "");
    void ResetVectors();

    // private helper methods
    bool IsInMvtx(const uint16_t layer);
    bool IsInIntt(const uint16_t layer);
    bool IsInTpc(const uint16_t layer);

    // virtual private methods
    virtual void DefineHistograms() = 0;

    // histograms
    VecHist1D m_vecHist1D;
    VecHist2D m_vecHist2D;

    // histogram definitions
    VecHistTypes m_vecHistTypes;
    VecHistDef1D m_vecHistDef1D;
    VecHistDef2D m_vecHistDef2D;

    // module utilities
    TrksInJetQAConfig m_config;
    TrksInJetQAHist   m_hist;

};  // end TrksInJetQABaseManager


#endif

// end ------------------------------------------------------------------------
