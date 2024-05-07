// ----------------------------------------------------------------------------
// 'TrksInJetQATrkManager.h'
// Derek Anderson
// 03.25.2024
//
// A submodule for the TrksInJetQA module to generate
// QA plots for tracks
// ----------------------------------------------------------------------------

#ifndef TRKSINJETQATRKMANAGER_H
#define TRKSINJETQATRKMANAGER_H

// submodule definitions
#include "TrksInJetQABaseManager.h"

// tracking includes
#include <trackbase_historic/SvtxTrack.h>

// root includes
#include <TH1.h>
#include <TH2.h>

// c++ utilities
#include <limits>
#include <utility>
#include <vector>

// TrksInJetQATrkManager definition -------------------------------------------

class TrksInJetQATrkManager : public TrksInJetQABaseManager
{
 public:
  // histogram accessors
  //   - TODO split tracks into seed types
  enum Type
  {
    All
  };
  enum H1D
  {
    Eta,
    Phi,
    Pt,
    Qual
  };
  enum H2D
  {
    EtaVsPhi,
    PtVsQual
  };

  // histogram content
  struct TrackQAContent
  {
    double eta = std::numeric_limits<double>::max();
    double phi = std::numeric_limits<double>::max();
    double pt = std::numeric_limits<double>::max();
    double qual = std::numeric_limits<double>::max();
  };

  // ctor/dtor
  using TrksInJetQABaseManager::TrksInJetQABaseManager;
  ~TrksInJetQATrkManager(){};

  // public methods
  void GetInfo(SvtxTrack* track);

 private:
  // private methods
  void FillHistograms(const int type, TrackQAContent& content);

  // inherited private methods
  void DefineHistograms() override;

};  // end TrksInJetQATrkManager

#endif

// end ------------------------------------------------------------------------
