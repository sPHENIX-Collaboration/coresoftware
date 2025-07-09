/// ===========================================================================
/*! \file   TrksInJetQADefs.h
 *  \author Derek Anderson
 *  \date   04.30.2024
 *
 *  Conveninent type definitions for the TrksInJetQA module.
 */
/// ===========================================================================

#ifndef TRKSINJETSQADEFS_H
#define TRKSINJETSQADEFS_H

// particle flow libraries
#include <particleflowreco/ParticleFlowElement.h>
#include <particleflowreco/ParticleFlowElementContainer.h>

// root libraries
#include <TH1.h>
#include <TH2.h>

// c++ utilities
#include <map>
#include <string>
#include <utility>
#include <vector>

// ============================================================================
//! Definitions for the TrksInJetQA module
// ============================================================================
namespace TrksInJetQADefs
{
  typedef std::pair<float, float> BinRange;
  typedef std::pair<uint32_t, BinRange> BinDef;
  typedef std::tuple<std::string, BinDef> HistDef1D;
  typedef std::tuple<std::string, BinDef, BinDef> HistDef2D;
  typedef std::vector<HistDef1D> VecHistDef1D;
  typedef std::vector<HistDef2D> VecHistDef2D;
  typedef std::vector<std::vector<TH1D*>> VecHist1D;
  typedef std::vector<std::vector<TH2D*>> VecHist2D;
  typedef std::vector<std::string> VecHistTypes;
  typedef std::map<int, HistDef1D> MapHistDef1D;
  typedef std::map<int, HistDef2D> MapHistDef2D;
  typedef std::map<int, TH1D*> MapHist1D;
  typedef std::map<int, TH2D*> MapHist2D;
  typedef std::map<int, std::string> MapHistTypes;
  typedef ParticleFlowElement PFObject;
  typedef ParticleFlowElementContainer PFObjectStore;
}

#endif

// end ========================================================================
