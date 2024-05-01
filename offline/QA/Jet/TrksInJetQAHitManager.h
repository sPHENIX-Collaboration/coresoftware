// ----------------------------------------------------------------------------
// 'TrksInJetQAHitManager.h'
// Derek Anderson
// 03.25.2024
//
// A submodule for the TrksInJetQA module
// to generate QA plots for track hits
// ----------------------------------------------------------------------------

#ifndef TRKSINJETQAHITMANAGER_H
#define TRKSINJETQAHITMANAGER_H

// c++ utilities
#include <limits>
#include <vector>
#include <utility>
// root libraries
#include <TH1.h>
#include <TH2.h>
// tracking libraries
#include <trackbase/TrkrHit.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrDefs.h>
// submodule definitions
#include "TrksInJetQABaseManager.h"



// TrksInJetQAHitManager definition -------------------------------------------

class TrksInJetQAHitManager : public TrksInJetQABaseManager {

  public:

    // histogram accessors
    enum Type {Mvtx, Intt, Tpc, All};
    enum H1D  {Ene, ADC, Layer, PhiBin, ZBin};
    enum H2D  {EneVsLayer, EneVsADC, PhiVsZBin};

    // histogram content
    struct HitQAContent {
      double   ene    = std::numeric_limits<double>::max();
      uint64_t adc    = std::numeric_limits<uint64_t>::max();
      uint16_t layer  = std::numeric_limits<uint16_t>::max();
      uint16_t phiBin = std::numeric_limits<uint16_t>::max();
      uint16_t zBin   = std::numeric_limits<uint16_t>::max();
    };

    // ctor/dtor
    using TrksInJetQABaseManager::TrksInJetQABaseManager;
    ~TrksInJetQAHitManager() {};

    // public methods
    void GetInfo(TrkrHit* hit, TrkrDefs::hitsetkey& setKey, TrkrDefs::hitkey& hitKey);

  private:

    // private methods
    void FillHistograms(const int type, HitQAContent& content);

    // inherited private methods
    void DefineHistograms() override;

};  // end TrksInJetQAHitManager

#endif

// end ------------------------------------------------------------------------
