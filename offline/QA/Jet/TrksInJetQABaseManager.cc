/// ===========================================================================
/*! \file   TrksInJetQABaseManager.cc
 *  \author Derek Anderson
 *  \date   04.03.2024
 *
 *  Base hist manager submodule for the TrksInJetQA module which
 *  consolidates methods/data common to all of the hist managers
 */
/// ===========================================================================

#define TRKSINJETQABASEMANAGER_CC

#include "TrksInJetQABaseManager.h"

// ctor/dtor ==================================================================

// ----------------------------------------------------------------------------
//! Default ctor
// ----------------------------------------------------------------------------
TrksInJetQABaseManager::TrksInJetQABaseManager(TrksInJetQAConfig& config,
                                               TrksInJetQAHist& hist)
  : m_config(config)
  , m_hist(hist)
{};

// ----------------------------------------------------------------------------
//! Default dtor
// ----------------------------------------------------------------------------
TrksInJetQABaseManager::~TrksInJetQABaseManager() {};

// public methods =============================================================

// ----------------------------------------------------------------------------
//! Define and generate histograms for manager
// ----------------------------------------------------------------------------
void TrksInJetQABaseManager::MakeHistograms(const std::string& prefix,
                                            const std::string& suffix)
{
  DefineHistograms();
  BuildHistograms(prefix, suffix);
}  // end 'MakeHistograms(std::string)'

// ----------------------------------------------------------------------------
//! Save histograms to output file
// ----------------------------------------------------------------------------
/*! Note that this only relevant if output
 *  mode is OutMode::File.
 */
void TrksInJetQABaseManager::SaveHistograms(TDirectory* topDir, const std::string& outDirName)
{
  TDirectory* outDir = topDir->mkdir(outDirName.c_str());
  if (!outDir)
  {
    std::cerr << PHWHERE << ": PANIC: unable to make output directory!" << std::endl;
    exit(1);
  }

  outDir->cd();
  for (const auto& hist1D : m_mapHist1D)
  {
    hist1D.second->Write();
  }
  for (const auto& hist2D : m_mapHist2D)
  {
    hist2D.second->Write();
  }
}  // end 'SaveHistograms(TDirectory*, std::string)'

// ----------------------------------------------------------------------------
//! Grab histograms from manager
// ----------------------------------------------------------------------------
void TrksInJetQABaseManager::GrabHistograms(std::vector<TH1D*>& vecOutHist1D,
                                            std::vector<TH2D*>& vecOutHist2D)
{
  for (const auto& hist1D : m_mapHist1D)
  {
    vecOutHist1D.push_back(hist1D.second);
  }
  for (const auto& hist2D : m_mapHist2D)
  {
    vecOutHist2D.push_back(hist2D.second);
  }
}  // end 'GrabHistograms(std::vector<TH1D*>&, std::vector<TH2D*>&)'

// private methods ============================================================

// ----------------------------------------------------------------------------
//! Build histograms from definitions
// ----------------------------------------------------------------------------
/*! Note that the specific histogram definitions are
 *  implemented in the derived classes.
 */
void TrksInJetQABaseManager::BuildHistograms(const std::string& prefix,
                                             const std::string& suffix)
{
  // build 1d histograms
  for (const auto& histType : m_mapHistTypes)
  {
    for (const auto& histDef1D : m_mapHistDef1D)
    {
      // make name
      std::string sHistName(prefix + "_");
      sHistName += histType.second;
      sHistName += std::get<0>(histDef1D.second);
      sHistName += "_";
      sHistName += suffix;

      // make sure histogram name is lower case
      std::transform(sHistName.begin(),
                     sHistName.end(),
                     sHistName.begin(),
                     ::tolower);
      std::regex_replace(
          sHistName,
          std::regex("__"),
          "_");

      // create histogram
      m_mapHist1D.emplace(Index(histType.first, histDef1D.first),
                          new TH1D(sHistName.data(),
                                   "",
                                   std::get<1>(histDef1D.second).first,
                                   std::get<1>(histDef1D.second).second.first,
                                   std::get<1>(histDef1D.second).second.second));
    }  // end 1D hist loop

    // build 2d histograms
    for (const auto& histDef2D : m_mapHistDef2D)
    {
      // make name
      std::string sHistName(prefix + "_");
      sHistName += histType.second;
      sHistName += std::get<0>(histDef2D.second);
      sHistName += "_";
      sHistName += suffix;

      // make sure histogram name is lower case
      std::transform(sHistName.begin(),
                     sHistName.end(),
                     sHistName.begin(),
                     ::tolower);
      std::regex_replace(
          sHistName,
          std::regex("__"),
          "_");

      // create histogram
      m_mapHist2D.emplace(Index(histType.first, histDef2D.first),
                          new TH2D(sHistName.data(),
                                   "",
                                   std::get<1>(histDef2D.second).first,
                                   std::get<1>(histDef2D.second).second.first,
                                   std::get<1>(histDef2D.second).second.second,
                                   std::get<2>(histDef2D.second).first,
                                   std::get<2>(histDef2D.second).second.first,
                                   std::get<2>(histDef2D.second).second.second));
    }  // end hist loop
  }  // end type loop
}  // end 'BuildHistograms(std::string)'

// private helper methods =====================================================

// ----------------------------------------------------------------------------
//! Check if a layer is in the MVTX
// ----------------------------------------------------------------------------
bool TrksInJetQABaseManager::IsInMvtx(const uint16_t layer) const
{
  return (layer < m_config.nMvtxLayer);
}  // end 'IsInMvtx(uint16_t)'

// ----------------------------------------------------------------------------
//! Check if a layer is in the INTT
// ----------------------------------------------------------------------------
bool TrksInJetQABaseManager::IsInIntt(const uint16_t layer) const
{
  return ((layer >= m_config.nMvtxLayer) &&
          (layer < m_config.nInttLayer + m_config.nMvtxLayer));
}  // end 'IsInIntt(uint16_t)'

// ----------------------------------------------------------------------------
//! Check if a layer is in the TPC
// ----------------------------------------------------------------------------
bool TrksInJetQABaseManager::IsInTpc(const uint16_t layer) const
{
  return (layer >= m_config.nMvtxLayer + m_config.nInttLayer);
}  // end 'IsInTpc(uint16_t)'

// ----------------------------------------------------------------------------
//! Create an index based on object type and histogram index
// ----------------------------------------------------------------------------
int TrksInJetQABaseManager::Index(const int type, const int hist) const
{
  return (100 * type) + hist;
}  // end 'Index(int, int)'

// end ========================================================================
