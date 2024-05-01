// ----------------------------------------------------------------------------
// 'TrksInJetQAHist.h'
// Derek Anderson
// 03.25.2024
//
// Configurable parameters for histograms (like binning, etc.)
// for the TrksInJetQA module.
// ----------------------------------------------------------------------------

#ifndef TRKSINJETSQAHIST_H
#define TRKSINJETSQAHIST_H

// c++ utilities
#include <string>
#include <utility>
// module utilities
#include "TrksInJetQATypes.h"



// TrksInJetQAHist definition -------------------------------------------------

struct TrksInJetQAHist {

  enum Var {
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
    Qual
  };

  // no. of bins
  //   - FIXME the phi/z bin binning needs
  //     some more thought...
  //   - FIXME same with the ADC...
  uint32_t nNumBins    = 101;
  uint32_t nAdcBins    = 101;
  uint32_t nZBinBins   = 101;
  uint32_t nPhiBinBins = 101;
  uint32_t nLayerBins  = 101;
  uint32_t nPosXYBins  = 6000;
  uint32_t nPosZBins   = 6000;
  uint32_t nPosRBins   = 3000;
  uint32_t nPhiBins    = 180;
  uint32_t nEtaBins    = 180;
  uint32_t nEneBins    = 505;
  uint32_t nQualBins   = 22;

  // bin ranges
  BinRange rNumBins    = {-0.5, (float) nNumBins + 0.5};
  BinRange rAdcBins    = {-0.5, (float) nAdcBins + 0.5};
  BinRange rZBinBins   = {-0.5, (float) nZBinBins + 0.5};
  BinRange rPhiBinBins = {-0.5, (float) nPhiBinBins + 0.5};
  BinRange rLayerBins  = {-0.5, (float) nLayerBins + 0.5};
  BinRange rPosXYBins  = {-300., 300.};
  BinRange rPosZBins   = {-300., 300.};
  BinRange rPosRBins   = {0., 300.};
  BinRange rPhiBins    = {-3.15, 3.15};
  BinRange rEtaBins    = {-4.0, 4.0};
  BinRange rEneBins    = {-0.5, 100.5};
  BinRange rQualBins   = {-0.5, 10.5};

  // construct list of binnings
  std::vector<BinDef> GetVecHistBins() {

    std::vector<BinDef> vecHistBins = {
      std::make_pair(nNumBins,    rNumBins),
      std::make_pair(nAdcBins,    rAdcBins),
      std::make_pair(nZBinBins,   rZBinBins),
      std::make_pair(nPhiBinBins, rPhiBinBins),
      std::make_pair(nLayerBins,  rLayerBins),
      std::make_pair(nPosXYBins,  rPosXYBins),
      std::make_pair(nPosZBins,   rPosZBins),
      std::make_pair(nPosRBins,   rPosRBins),
      std::make_pair(nPhiBins,    rPhiBins),
      std::make_pair(nEtaBins,    rEtaBins),
      std::make_pair(nEneBins,    rEneBins),
      std::make_pair(nQualBins,   rQualBins)
    };
    return vecHistBins;

  }  // end 'GetVecHistBins()'

};  // end TrksInJetQAHist

#endif

// end ------------------------------------------------------------------------
