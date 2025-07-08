/// ===========================================================================
/*! \file   TrksInJetQAClustManager.h
 *  \author Derek Anderson
 *  \date   03.25.2024
 *
 *  A submodule for the TrksInJetQA module to generate
 *  QA plots for track clusters
 */
/// ===========================================================================

#ifndef TRKSINJETQACLUSTMANAGER_H
#define TRKSINJETQACLUSTMANAGER_H

// submodule definitions
#include "TrksInJetQABaseManager.h"

// tracking includes
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>

// root includes
#include <TH1.h>
#include <TH2.h>

// c++ utilities
#include <limits>
#include <utility>
#include <vector>

// ============================================================================
//! Tracker cluster histogram manager for TrksInJetQA module
// ============================================================================
/*! This histogram manager defines what to histogram
 *  from tracker clusters.
 */
class TrksInJetQAClustManager : public TrksInJetQABaseManager
{
 public:
  ///! enumerate cluster subsystems
  enum Type
  {
    All,
    Mvtx,
    Intt,
    Tpc
  };

  ///! enumerate 1D histograms
  enum H1D
  {
    PosX,
    PosY,
    PosZ,
    PosR
  };

  ///! enumerate 2D histograms
  enum H2D
  {
    PosYvsX,
    PosRvsZ
  };

  // --------------------------------------------------------------------------
  // Cluster histogram content
  // --------------------------------------------------------------------------
  /*! A small struct to consolidate what variables
   *  to histogram for clusters.
   */ 
  struct ClustQAContent
  {
    double x = std::numeric_limits<double>::max();
    double y = std::numeric_limits<double>::max();
    double z = std::numeric_limits<double>::max();
    double r = std::numeric_limits<double>::max();
  };

  // ctor/dtor
  using TrksInJetQABaseManager::TrksInJetQABaseManager;
  ~TrksInJetQAClustManager(){};

  // public methods
  void GetInfo(TrkrCluster* cluster, TrkrDefs::cluskey& clustKey, ActsGeometry* actsGeom);

 private:
  // private methods
  void FillHistograms(const int type, ClustQAContent& content);

  // inherited private methods
  void DefineHistograms() override;

};  // end TrksInJetQAClustManager

#endif

// end ========================================================================
