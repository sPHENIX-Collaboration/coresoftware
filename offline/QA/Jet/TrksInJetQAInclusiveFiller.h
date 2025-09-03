/// ===========================================================================
/*! \file   TrksInJetQAInclusiveFiller.h
 *  \author Derek Anderson
 *  \date   04.03.2024
 *
 *  A submodule for the TrksInJetQA F4A module to produce
 *  QA histograms for tracks and more in jets
 */
/// ===========================================================================

#ifndef TRKSINJETQAINCLUSIVEQAFILLER_H
#define TRKSINJETQAINCLUSIVEQAFILLER_H

// submodule definitions
#include "TrksInJetQABaseFiller.h"

// jet libraries
#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>

// phool libraries
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// tracking libraries
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

// c+ utilities
#include <cassert>

// ============================================================================
//! Inclusive histogram filler for TrksInJetQA module
// ============================================================================
/*! This histogram filler defines how to fill histograms
 *  for inclusive populations (i.e. all objects regardless
 *  of if they're in jets or not).
 */
class TrksInJetQAInclusiveFiller : public TrksInJetQABaseFiller
{
 public:
  // ctor/dtor
  using TrksInJetQABaseFiller::TrksInJetQABaseFiller;
  ~TrksInJetQAInclusiveFiller() override = default;

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

/// end =======================================================================
