/// ===========================================================================
/*! \file   TrksInJetQABaseManager.h
 *  \author Derek Anderson
 *  \date   04.03.2024
 *
 *  Base hist manager submodule for the TrksInJetQA module which
 *  consolidates methods/data common to all of the hist managers
 */
/// ===========================================================================

#ifndef TRKSINJETQABASEMANAGER_H
#define TRKSINJETQABASEMANAGER_H

// module utilities
#include "TrksInJetQAConfig.h"
#include "TrksInJetQADefs.h"
#include "TrksInJetQAHist.h"

// phool includes
#include <phool/phool.h>

// root includes
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>

// c++ utilities
#include <iostream>
#include <regex>
#include <string>
#include <utility>
#include <vector>

// ============================================================================
//! Base histogram manager for TrksInJetQA module
// ============================================================================
/*! This class consolidates all common functionality
 *  for Histogram Managers deployed in the TrksInJetQA
 *  module. The histogram managers define what variables
 *  to histogram for a given object (e.g. tracks) and
 *  how to extract the relevant information from that
 *  object.
 */
class TrksInJetQABaseManager
{
 public:
  // ctor/dtor
  TrksInJetQABaseManager(TrksInJetQAConfig& config, TrksInJetQAHist& hist);
  virtual ~TrksInJetQABaseManager();

  // public methods
  void MakeHistograms(const std::string& prefix = "", const std::string& suffix = "");
  void SaveHistograms(TDirectory* topDir, const std::string& outDirName);
  void GrabHistograms(std::vector<TH1D*>& vecOutHist1D, std::vector<TH2D*>& vecOutHist2D);

 protected:
  // private methods
  void BuildHistograms(const std::string& prefix = "", const std::string& suffix = "");
  void ResetVectors();

  // private helper methods
  bool IsInMvtx(const uint16_t layer) const;
  bool IsInIntt(const uint16_t layer) const;
  bool IsInTpc(const uint16_t layer) const;

  // virtual private methods
  virtual void DefineHistograms() = 0;

  // histograms
  // FIXME THIS IS WHAT I NEED TO CHANGE
  TrksInJetQADefs::VecHist1D m_vecHist1D;
  TrksInJetQADefs::VecHist2D m_vecHist2D;

  // histogram definitions
  // FIXME THIS IS WHAT I NEED TO CHANGE
  TrksInJetQADefs::VecHistTypes m_vecHistTypes;
  TrksInJetQADefs::VecHistDef1D m_vecHistDef1D;
  TrksInJetQADefs::VecHistDef2D m_vecHistDef2D;

  // module utilities
  TrksInJetQAConfig m_config;
  TrksInJetQAHist m_hist;

};  // end TrksInJetQABaseManager

#endif

// end ========================================================================
