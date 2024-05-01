// ----------------------------------------------------------------------------
// 'TrksInJetQAClustManager.h'
// Derek Anderson
// 03.25.2024
//
// A submodule for the TrksInJetQA module to generate
// QA plots for track clusters
// ----------------------------------------------------------------------------

#ifndef TRKSINJETQACLUSTMANAGER_H
#define TRKSINJETQACLUSTMANAGER_H

// c++ utilities
#include <limits>
#include <vector>
#include <utility>
// root libraries
#include <TH1.h>
#include <TH2.h>
// tracking libraries
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/ActsGeometry.h>
// submodule definitions
#include "TrksInJetQABaseManager.h"



// TrksInJetQAClustManager definition -----------------------------------------

class TrksInJetQAClustManager : public TrksInJetQABaseManager {

  public:

    // histogram accessors
    enum Type {Mvtx, Intt, Tpc, All};
    enum H1D  {PosX, PosY, PosZ, PosR};
    enum H2D  {PosYvsX, PosRvsZ};

    // histogram content
    struct ClustQAContent {
      double x = std::numeric_limits<double>::max();
      double y = std::numeric_limits<double>::max();
      double z = std::numeric_limits<double>::max();
      double r = std::numeric_limits<double>::max();
    };

    // ctor/dtor
    using TrksInJetQABaseManager::TrksInJetQABaseManager;
    ~TrksInJetQAClustManager() {};

    // public methods
    void GetInfo(TrkrCluster* cluster, TrkrDefs::cluskey& clustKey, ActsGeometry* actsGeom);

  private:

    // private methods
    void FillHistograms(const int type, ClustQAContent& content);

    // inherited private methods
    void DefineHistograms() override;

};  // end TrksInJetQAClustManager

#endif

// end ------------------------------------------------------------------------
