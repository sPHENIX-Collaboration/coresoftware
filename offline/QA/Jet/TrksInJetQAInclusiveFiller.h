// ----------------------------------------------------------------------------
// 'TrksInJetQAInclusiveFiller.h'
// Derek Anderson
// 04.03.2024
//
// A submodule for the TrksInJetQA F4A module to produce
// QA histograms for tracks and more in jets
// ----------------------------------------------------------------------------

#ifndef TRACKSINJETSQAMAKER_INCLUSIVEQAHISTFILLER_H
#define TRACKSINJETSQAMAKER_INCLUSIVEQAHISTFILLER_H

// c+ utilities
#include <cassert>
// phool libraries
#include <phool/phool.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
// tracking libraries
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
// jet libraries
#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>
// submodule definitions
#include "TrksInJetQABaseFiller.h"



// TrksInJetQAInclusiveFiller -------------------------------------------------

class TrksInJetQAInclusiveFiller : public TrksInJetQABaseFiller {

  public:

    // ctor/dtor
    using TrksInJetQABaseFiller::TrksInJetQABaseFiller;
    ~TrksInJetQAInclusiveFiller() {};

    // inherited public methods
    void Fill(PHCompositeNode* topNode) override;

  private:

    // private methods
    void FillHitQAHists();
    void FillClustQAHists();
    void FillTrackQAHists();
    void FillJetQAHists();

};  // end TrksInJetQAInclusiveFiller

#endif

// end ------------------------------------------------------------------------
