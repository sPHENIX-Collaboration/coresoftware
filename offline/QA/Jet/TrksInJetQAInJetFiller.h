// ----------------------------------------------------------------------------
// 'TrksInJetQAInJetFiller.h'
// Derek Anderson
// 04.03.2024
//
// A submodule for the TrksInJetsQA F4A module to produce
// QA histograms for tracks and more in jets
// ----------------------------------------------------------------------------

#ifndef TRKSINJETQAINJETFILLER_H
#define TRKSINJETQAINJETFILLER_H

// c+ utilities
#include <vector>
#include <cassert>
// root libraries
#include <TMath.h>
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
// particle flow libraries
#include <particleflowreco/ParticleFlowElement.h>
#include <particleflowreco/ParticleFlowElementContainer.h>
// jet libraries
#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>
// g4eval libraries
#include <g4eval/ClusKeyIter.h>
// module utilities
#include "TrksInJetQATypes.h"
// submodule definitions
#include "TrksInJetQABaseFiller.h"



// TrksInJetQAInJetFiller -----------------------------------------------------

class TrksInJetQAInJetFiller : public TrksInJetQABaseFiller {

  public:

    // additional nodes to grab
    enum Node {Flow};

    // ctor/dtor
    using TrksInJetQABaseFiller::TrksInJetQABaseFiller;
    ~TrksInJetQAInJetFiller() {};

    // inherited public methods
    void Fill(PHCompositeNode* topNode) override;

  private:

    // private methods
    void       GetNode(const int node, PHCompositeNode* topNode);
    void       FillJetAndTrackQAHists(PHCompositeNode* topNode);
    void       FillClustAndHitQAHists(SvtxTrack* track, PHCompositeNode* topNode);
    void       GetCstTracks(Jet* jet, PHCompositeNode* topNode);
    void       GetNonCstTracks(Jet* jet, PHCompositeNode* topNode);
    bool       IsCstNotRelevant(const uint32_t type);
    bool       IsTrkInList(const uint32_t id);
    double     GetTrackJetDist(SvtxTrack* track, Jet* jet);
    PFObject*  GetPFObject(const uint32_t id, PHCompositeNode* topNode);
    SvtxTrack* GetTrkFromPFO(PFObject* pfo);

    // additional dst nodes needed
    PFObjectStore* m_flowStore = NULL;

    // for tracks in jet
    std::vector<SvtxTrack*> m_trksInJet;

};  // end TrksInJetQAInJetFiller

#endif

// end ------------------------------------------------------------------------
