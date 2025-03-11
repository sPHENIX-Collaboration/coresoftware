/// ===========================================================================
/*! \file   CaloStatusMapperLinkDef.h
 *  \author Derek Anderson
 *  \date   05.22.2024
 *
 *  A Fun4All QA module to plot no. of towers per event
 *  and vs. eta, phi as a function of status.
 */
/// ===========================================================================

#define CLUSTERSTATUSMAPPER_CC

// module definition
#include "CaloStatusMapper.h"

// calo base
#include <calobase/TowerInfov2.h>

// calo trigger
#include <calotrigger/TriggerAnalyzer.h>

// f4a libraries
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>

// phool libraries
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>

// qa utilities
#include <qautils/QAHistManagerDef.h>

// root libraries
#include <TH1.h>
#include <TH2.h>

// c++ utiilites
#include <algorithm>
#include <cassert>
#include <iostream>

// abbreviate namespace for convenience
namespace CSMD = CaloStatusMapperDefs;



// ctor/dtor ==================================================================

// ----------------------------------------------------------------------------
//! Default module constructor
// ----------------------------------------------------------------------------
CaloStatusMapper::CaloStatusMapper(const std::string& modulename, const bool debug)
  : SubsysReco(modulename)
{

  // print debug message
  if (debug && (Verbosity() > 1))
  {
    std::cout << "CaloStatusMapper::CaloStatusMapper(std::string&, bool) Calling ctor" << std::endl;
  }

  // make sure node vector is empty
  m_inNodes.clear();

}  // end ctor(std::string&, bool)



// ----------------------------------------------------------------------------
//! Module constructor accepting a configuration
// ----------------------------------------------------------------------------
CaloStatusMapper::CaloStatusMapper(const Config& config)
  : SubsysReco(config.moduleName)
  , m_config(config)
{

  // print debug message
  if (m_config.debug && (Verbosity() > 1))
  {
    std::cout << "CaloStatusMapper::CaloStatusMapper(Config&) Calling ctor" << std::endl;
  }

}  // end ctor(Config&)



// ----------------------------------------------------------------------------
//! Module destructor
// ----------------------------------------------------------------------------
CaloStatusMapper::~CaloStatusMapper()
{

  // print debug message
  if (m_config.debug && (Verbosity() > 1))
  {
    std::cout << "CaloStatusMapper::~CaloStatusMapper() Calling dtor" << std::endl;
  }
  delete m_analyzer;

}  // end dtor



// fun4all methods ============================================================

// ----------------------------------------------------------------------------
//! Initialize module
// ----------------------------------------------------------------------------
int CaloStatusMapper::Init(PHCompositeNode* topNode)
{

  if (m_config.debug)
  {
    std::cout << "CaloStatusMapper::Init(PHCompositeNode*) Initializing" << std::endl;
  }

  // initialize trigger analyzer
  delete m_analyzer;
  m_analyzer = new TriggerAnalyzer();

  // initialize manager/histograms
  InitHistManager();
  BuildHistograms();

  // make sure event no. is set to 0
  m_nEvent = 0;
  return Fun4AllReturnCodes::EVENT_OK;

}  // end 'Init(PHCompositeNode*)'



// ----------------------------------------------------------------------------
//! Grab inputs and fills histograms
// ----------------------------------------------------------------------------
int CaloStatusMapper::process_event(PHCompositeNode* topNode)
{

  if (m_config.debug)
  {
    std::cout << "CaloStatusMapper::process_event(PHCompositeNode* topNode) Processing Event" << std::endl;
  }

  // if needed, check if selected trigger fired
  if (m_config.doTrgSelect)
  {
    m_analyzer -> decodeTriggers(topNode);
    bool hasTrigger = JetQADefs::DidTriggerFire(m_config.trgToSelect, m_analyzer);
    if (!hasTrigger)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  // grab input nodes
  GrabNodes(topNode);

  // loop over input nodes
  for (size_t iNode = 0; iNode < m_inNodes.size(); ++iNode)
  {

    // grab node name & make status base
    const std::string nodeName = m_config.inNodeNames[iNode].first;
    const std::string statBase = MakeBaseName("Status", nodeName);

    // loop over towers
    TowerInfoContainer* towers = m_inNodes[iNode];
    for (size_t iTower = 0; iTower < towers->size(); ++iTower)
    {

      // grab eta, phi indices
      const int32_t key  = towers -> encode_key(iTower);
      const int32_t iEta = towers -> getTowerEtaBin(key);
      const int32_t iPhi = towers -> getTowerPhiBin(key);

      // get status
      const auto tower  = towers -> get_tower_at_channel(iTower);
      const auto status = CSMD::GetTowerStatus(tower);
      if (status == CSMD::Stat::Unknown)
      {
        std::cout << PHWHERE << ": Warning! Tower has an unknown status!\n"
                  << "  channel = " << iTower << ", key = " << key << "\n"
                  << "  node = " << m_config.inNodeNames[iNode].first
                  << std::endl; 
        continue;
      } 

      // make base eta/phi hist name
      const std::string statLabel  = m_mapStatLabels[status];
      const std::string perEtaBase = MakeBaseName("NPerEta", nodeName, statLabel);
      const std::string perPhiBase = MakeBaseName("NPerPhi", nodeName, statLabel);
      const std::string phiEtaBase = MakeBaseName("PhiVsEta", nodeName, statLabel);

      // fill histograms accordingly
      m_hists[statBase]   -> Fill(status);
      m_hists[perEtaBase] -> Fill(iEta);
      m_hists[perPhiBase] -> Fill(iPhi);
      m_hists[phiEtaBase] -> Fill(iEta, iPhi);

    }  // end tower loop
  }  // end node loop

  // increment event no. and return
  ++m_nEvent;
  return Fun4AllReturnCodes::EVENT_OK;

}  // end 'process_event(PHCompositeNode*)'



// ----------------------------------------------------------------------------
//! Run final calculations
// ----------------------------------------------------------------------------
int CaloStatusMapper::End(PHCompositeNode *topNode)
{

  if (m_config.debug)
  {
    std::cout << "CaloStatusMapper::End(PHCompositeNode* topNode) This is the End..." << std::endl;
  }

  // normalize avg. status no.s
  for (const auto& nodeName : m_config.inNodeNames)
  {
    const std::string statBase = MakeBaseName("Status", nodeName.first);
    m_hists[statBase] -> Scale(1. / (double) m_nEvent);
  }

  // register hists and exit
  for (const auto& hist : m_hists) {
    m_manager -> registerHisto(hist.second);
  }
  return Fun4AllReturnCodes::EVENT_OK;

}  // end 'End(PHCompositeNode*)'



// private methods ============================================================

// ----------------------------------------------------------------------------
//! Initialize histogram manager
// ----------------------------------------------------------------------------
void CaloStatusMapper::InitHistManager()
{

  // print debug message
  if (m_config.debug && (Verbosity() > 0))
  {
    std::cout << "CaloStatusMapper::InitHistManager() Initializing histogram manager" << std::endl;
  }

  m_manager = QAHistManagerDef::getHistoManager();
  if (!m_manager)
  {
    std::cerr << PHWHERE << ": PANIC! Couldn't grab histogram manager!" << std::endl;
    assert(m_manager);
  }
  return;

}  // end 'InitHistManager()'



// ----------------------------------------------------------------------------
//! Build histograms
// ----------------------------------------------------------------------------
void CaloStatusMapper::BuildHistograms()
{

  // print debug message
  if (m_config.debug && (Verbosity() > 0))
  {
    std::cout << "CaloStatusMapper::BuildHistograms() Creating histograms" << std::endl;
  }

  // instantiate histogram definitions
  const CSMD::EMCalHistDef emHistDef;
  const CSMD::HCalHistDef  hcHistDef;

  // loop over input node names
  for (const auto& nodeName : m_config.inNodeNames)
  {

    // make status hist name
    const std::string statBase = MakeBaseName("Status", nodeName.first);
    const std::string statName = CSMD::MakeQAHistName(statBase, m_config.moduleName, m_config.histTag);

    // create status hist
    //   - n.b. calo type doesn't matter here
    m_hists[statBase] = emHistDef.MakeStatus1D(statName);

    // loop over status labels
    for (const auto& statLabel : m_mapStatLabels)
    {

      // set relevant bin label for status histogram
      m_hists[statBase] -> GetXaxis() -> SetBinLabel(statLabel.first + 1, statLabel.second.data());

      // make base eta/phi hist name
      const std::string perEtaBase = MakeBaseName("NPerEta", nodeName.first, statLabel.second);
      const std::string perPhiBase = MakeBaseName("NPerPhi", nodeName.first, statLabel.second);
      const std::string phiEtaBase = MakeBaseName("PhiVsEta", nodeName.first, statLabel.second);

      // make full eta/phi hist name
      const std::string namePerEta = CSMD::MakeQAHistName(perEtaBase, m_config.moduleName, m_config.histTag);
      const std::string namePerPhi = CSMD::MakeQAHistName(perPhiBase, m_config.moduleName, m_config.histTag);
      const std::string namePhiEta = CSMD::MakeQAHistName(phiEtaBase, m_config.moduleName, m_config.histTag);

      // make eta/phi hists
      switch (nodeName.second)
      {
        case CSMD::Calo::HCal:
          m_hists[perEtaBase] = hcHistDef.MakeEta1D(namePerEta);
          m_hists[perPhiBase] = hcHistDef.MakePhi1D(namePerPhi);
          m_hists[phiEtaBase] = hcHistDef.MakePhiEta2D(namePhiEta);
          break;
        case CSMD::Calo::EMCal:
          [[fallthrough]];
        default:
          m_hists[perEtaBase] = emHistDef.MakeEta1D(namePerEta);
          m_hists[perPhiBase] = emHistDef.MakePhi1D(namePerPhi);
          m_hists[phiEtaBase] = emHistDef.MakePhiEta2D(namePhiEta);
          break;
      }

    }  // end status loop
  }  // end node loop
  return;

}  // end 'BuildHistograms()'



// ----------------------------------------------------------------------------
//! Grab input nodes
// ----------------------------------------------------------------------------
void CaloStatusMapper::GrabNodes(PHCompositeNode* topNode)
{

  // print debug message
  if (m_config.debug && (Verbosity() > 0))
  {
    std::cout << "CaloStatusMapper::GrabNodes(PHCompositeNode*) Grabbing input nodes" << std::endl;
  }

  // make sure node vector is empty
  m_inNodes.clear();

  // loop over nodes to grab
  for (auto inNodeName : m_config.inNodeNames)
  {
    m_inNodes.push_back(
      findNode::getClass<TowerInfoContainer>(topNode, inNodeName.first)
    );
    if (!m_inNodes.back())
    {
      std::cerr << PHWHERE << ":" << " PANIC! Not able to grab node " << inNodeName.first << "! Aborting!" << std::endl;
      assert(m_inNodes.back());
    }
  }  // end input node name loop
  return;

}  // end 'GrabNodes(PHCompositeNode*)'



// ----------------------------------------------------------------------------
//! Make base histogram name
// ----------------------------------------------------------------------------
std::string CaloStatusMapper::MakeBaseName(
  const std::string& base,
  const std::string& node,
  const std::string& stat) const
{

  if (m_config.debug && (Verbosity() > 2))
  {
    std::cout << "CaloStatusMapper::MakeBaseName(std::string& x 3) Making base histogram name" << std::endl;
  }

  std::string name = base + "_" + node;
  if (!stat.empty())
  {
    name.insert(0, stat + "_");
  }
  return name;

}  // end 'MakeBaseName(std::string& x 3)'

// end ========================================================================
