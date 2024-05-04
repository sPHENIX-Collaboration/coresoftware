// ----------------------------------------------------------------------------
// 'TrksInJetQAConfig.h'
// Derek Anderson
// 03.25.2024
//
// Configurable parameters for the TrksInJetQA module.
// ----------------------------------------------------------------------------

#ifndef TRKSINJETQACONFIG_H
#define TRKSINJETQACONFIG_H

// c++ utilities
#include <string>



// TrksInJetQAConfig definition -----------------------------------------------

struct TrksInJetQAConfig {

  // system options
  int  outMode     = 0;
  int  verbose     = 0;
  bool doDebug     = false;
  bool doInclusive = true;
  bool doInJet     = true;
  bool doHitQA     = false;
  bool doClustQA   = false;
  bool doTrackQA   = true;
  bool doJetQA     = true;

  // jet options
  double rJet = 0.4;

  // input options
  std::string jetInNode   = "AntiKt_Track_r04";
  std::string trkInNode   = "SvtxTrackMap";
  std::string clustInNode = "TRKR_CLUSTER";
  std::string hitInNode   = "TRKR_HITSET";

  // output options
  std::string inclusiveDir = "Inclusive";
  std::string inJetDir     = "InJet";
  std::string hitOutDir    = "Hit";
  std::string clustOutDir  = "Clust";
  std::string trackOutDir  = "Track";
  std::string jetOutDir    = "Jet";

  // tracker parameters
  uint16_t nMvtxLayer = 3;
  uint16_t nInttLayer = 4;
  uint16_t nTpcLayer  = 48;

};  // end TrksInJetQAConfig

#endif

// end ------------------------------------------------------------------------
