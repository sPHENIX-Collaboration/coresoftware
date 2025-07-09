/// ===========================================================================
/*! \file   TrksInJetQATrkManager.h
 *  \author Derek Anderson
 *  \date   03.25.2024
 *
 *  A submodule for the TrksInJetQA module to generate
 *  QA plots for tracks
 */
/// ===========================================================================

#ifndef TRKSINJETQATRKMANAGER_H
#define TRKSINJETQATRKMANAGER_H

// submodule definitions
#include "TrksInJetQABaseManager.h"

// jet includes
#include <jetbase/Jet.h>

// tracking includes
#include <trackbase_historic/SvtxTrack.h>

// root includes
#include <TH1.h>
#include <TH2.h>

// c++ utilities
#include <limits>
#include <optional>
#include <utility>
#include <vector>

// ============================================================================
//! Track histogram manager for TrksInJetQA module
// ============================================================================
/*! This histogram manager defines what to histogram
 *  from tracks.
 */
class TrksInJetQATrkManager : public TrksInJetQABaseManager
{
 public:
  ///! enumerates types of tracks
  // - TODO split tracks into seed types
  enum Type
  {
    All
  };

  ///! enumerates 1D histograms
  enum H1D
  {
    Eta,
    Phi,
    Pt,
    Z,
    Jt,
    Qual
  };

  ///! enumerates 2D histograms
  enum H2D
  {
    EtaVsPhi,
    PtVsQual
  };

  // --------------------------------------------------------------------------
  // Track histogram content
  // --------------------------------------------------------------------------
  /*! A small struct to consolidate what variables
   *  to histogram for tracks.
   */ 
  struct TrackQAContent
  {
    double eta = std::numeric_limits<double>::max();
    double phi = std::numeric_limits<double>::max();
    double pt = std::numeric_limits<double>::max();
    double z = std::numeric_limits<double>::max();
    double jt = std::numeric_limits<double>::max();
    double qual = std::numeric_limits<double>::max();
  };

  // ctor/dtor
  using TrksInJetQABaseManager::TrksInJetQABaseManager;
  ~TrksInJetQATrkManager(){};

  // public methods
  void GetInfo(SvtxTrack* track, std::optional<Jet*> jet = std::nullopt);

 private:
  // private methods
  void FillHistograms(const int type, TrackQAContent& content);

  // inherited private methods
  void DefineHistograms() override;

};  // end TrksInJetQATrkManager

#endif

// end ========================================================================
