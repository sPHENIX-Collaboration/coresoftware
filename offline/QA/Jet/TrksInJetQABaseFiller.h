// ----------------------------------------------------------------------------
// 'TrksInJetQABaseFiller.h'
// Derek Anderson
// 04.11.2024
//
// A submodule for the TrksInJetQA F4A module to produce
// QA histograms for tracks and more in jets
// ----------------------------------------------------------------------------

#ifndef TRKSINJETQABASEFILLER_H
#define TRKSINJETQABASEFILLER_H

// c++ utilities
#include <string>
// root libraries
#include <TFile.h>
// phool libraries
#include <phool/phool.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
// tracking libraries
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrackMap.h>
// jet libraries
#include <jetbase/JetContainer.h>
// submodule definitions
#include "TrksInJetQAHitManager.h"
#include "TrksInJetQAClustManager.h"
#include "TrksInJetQATrkManager.h"
#include "TrksInJetQAJetManager.h"
// module utilties
#include "TrksInJetQAHist.h"
#include "TrksInJetQAConfig.h"



// TrksInJetQABaseFiller ------------------------------------------------------

class TrksInJetQABaseFiller {

  public:

    // ctor/dtor
    TrksInJetQABaseFiller(TrksInJetQAConfig& config, TrksInJetQAHist& hist);
    ~TrksInJetQABaseFiller();

    // public methods
    void MakeHistograms(std::string label = "");
    void SaveHistograms(TFile* outFile, std::string outDirName);
    void GrabHistograms(std::vector<TH1D*>& vecOutHist1D, std::vector<TH2D*>& vecOutHist2D);

    // virtual public methods
    virtual void Fill(PHCompositeNode* topNode) = 0;

  protected:

    // private methods
    void GetNodes(PHCompositeNode* topNode);

    // necessary dst nodes
    //   - FIXME these should be smart pointers!
    ActsGeometry*         m_actsGeom = NULL;
    TrkrHitSetContainer*  m_hitMap   = NULL;
    TrkrClusterContainer* m_clustMap = NULL;
    SvtxTrackMap*         m_trkMap   = NULL;
    JetContainer*         m_jetMap   = NULL;

    // submodules to use
    std::unique_ptr<TrksInJetQAHitManager>   m_hitManager   = NULL;
    std::unique_ptr<TrksInJetQAClustManager> m_clustManager = NULL;
    std::unique_ptr<TrksInJetQATrkManager>   m_trackManager = NULL;
    std::unique_ptr<TrksInJetQAJetManager>   m_jetManager   = NULL;

    // module utilities
    TrksInJetQAConfig m_config;
    TrksInJetQAHist   m_hist;

};  // end TrksInJetQABaseFiller

#endif

// end ------------------------------------------------------------------------
