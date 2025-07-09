/// ===========================================================================
/*! \file   TrksInJetQAHist.h
 *  \author Derek Anderson
 *  \date   03.25.2024
 *
 *  Configurable parameters for histograms (like binning, etc.)
 *  for the TrksInJetQA module.
 */
/// ===========================================================================

#ifndef TRKSINJETSQAHIST_H
#define TRKSINJETSQAHIST_H

// module utilities
#include "TrksInJetQADefs.h"

// c++ utilities
#include <string>
#include <utility>

// ============================================================================
//! Histogram parameters for the TrksInJetQA module
// ============================================================================
struct TrksInJetQAHist
{
  ///! Enumerates possible histograms
  enum Var
  {
    Num,
    Adc,
    ZBin,
    PhiBin,
    Layer,
    PosXY,
    PosZ,
    PosR,
    Phi,
    Eta,
    Ene,
    Qual,
    Frac
  };

  // no. of bins
  //   - FIXME the phi/z bin binning needs
  //     some more thought...
  //   - FIXME same with the ADC...
  uint32_t nNumBins = 101;
  uint32_t nAdcBins = 101;
  uint32_t nZBinBins = 101;
  uint32_t nPhiBinBins = 101;
  uint32_t nLayerBins = 101;
  uint32_t nPosXYBins = 6000;
  uint32_t nPosZBins = 6000;
  uint32_t nPosRBins = 3000;
  uint32_t nPhiBins = 180;
  uint32_t nEtaBins = 180;
  uint32_t nEneBins = 505;
  uint32_t nQualBins = 22;
  uint32_t nFracBins = 44;

  // bin ranges
  TrksInJetQADefs::BinRange rNumBins = {-0.5, (float) nNumBins + 0.5};
  TrksInJetQADefs::BinRange rAdcBins = {-0.5, (float) nAdcBins + 0.5};
  TrksInJetQADefs::BinRange rZBinBins = {-0.5, (float) nZBinBins + 0.5};
  TrksInJetQADefs::BinRange rPhiBinBins = {-0.5, (float) nPhiBinBins + 0.5};
  TrksInJetQADefs::BinRange rLayerBins = {-0.5, (float) nLayerBins + 0.5};
  TrksInJetQADefs::BinRange rPosXYBins = {-300., 300.};
  TrksInJetQADefs::BinRange rPosZBins = {-300., 300.};
  TrksInJetQADefs::BinRange rPosRBins = {0., 300.};
  TrksInJetQADefs::BinRange rPhiBins = {-3.15, 3.15};
  TrksInJetQADefs::BinRange rEtaBins = {-4.0, 4.0};
  TrksInJetQADefs::BinRange rEneBins = {-0.5, 100.5};
  TrksInJetQADefs::BinRange rQualBins = {-0.5, 10.5};
  TrksInJetQADefs::BinRange rFracBins = {-0.1, 2.1};

  // construct list of binnings
  std::vector<TrksInJetQADefs::BinDef> GetVecHistBins()
  {
    std::vector<TrksInJetQADefs::BinDef> vecHistBins = {
        std::make_pair(nNumBins, rNumBins),
        std::make_pair(nAdcBins, rAdcBins),
        std::make_pair(nZBinBins, rZBinBins),
        std::make_pair(nPhiBinBins, rPhiBinBins),
        std::make_pair(nLayerBins, rLayerBins),
        std::make_pair(nPosXYBins, rPosXYBins),
        std::make_pair(nPosZBins, rPosZBins),
        std::make_pair(nPosRBins, rPosRBins),
        std::make_pair(nPhiBins, rPhiBins),
        std::make_pair(nEtaBins, rEtaBins),
        std::make_pair(nEneBins, rEneBins),
        std::make_pair(nQualBins, rQualBins),
        std::make_pair(nFracBins, rFracBins)};
    return vecHistBins;
  }  // end 'GetVecHistBins()'
};  // end TrksInJetQAHist

#endif

// end ========================================================================
