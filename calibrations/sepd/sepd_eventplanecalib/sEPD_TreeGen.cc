#include "sEPD_TreeGen.h"
#include "QVecDefs.h"

// -- c++
#include <format>
#include <iostream>

// -- event
#include <ffaobjects/EventHeader.h>

// -- Fun4All
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

// -- Nodes
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

// -- Calo
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

// -- Vtx
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

// -- MB
#include <calotrigger/MinimumBiasClassifier.h>
#include <calotrigger/MinimumBiasInfo.h>
#include <centrality/CentralityInfo.h>

// -- sEPD
#include <epd/EpdGeom.h>

//____________________________________________________________________________..
sEPD_TreeGen::sEPD_TreeGen(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int sEPD_TreeGen::Init([[maybe_unused]] PHCompositeNode *topNode)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Print("NODETREE");

  unsigned int bins_sepd_totalcharge{100};
  double sepd_totalcharge_low{0};
  double sepd_totalcharge_high{2e4};

  unsigned int bins_centrality{80};
  double centrality_low{-0.5};
  double centrality_high{79.5};

  hSEPD_Charge = std::make_unique<TProfile>("hSEPD_Charge", "|z| < 10 cm and MB; Channel; Avg Charge", QVecShared::sepd_channels, 0, QVecShared::sepd_channels);
  hSEPD_Charge->Sumw2();

  h2SEPD_totalcharge_centrality = std::make_unique<TH2F>("h2SEPD_totalcharge_centrality", "|z| < 10 cm and MB; sEPD Total Charge; Centrality [%]", bins_sepd_totalcharge, sepd_totalcharge_low, sepd_totalcharge_high, bins_centrality, centrality_low, centrality_high);

  m_output = std::make_unique<TFile>(m_outtree_name.c_str(), "recreate");
  m_output->cd();

  // TTree
  m_tree = std::make_unique<TTree>("T", "T");
  m_tree->SetDirectory(m_output.get());
  m_tree->Branch("event_id", &m_data.event_id);
  m_tree->Branch("event_zvertex", &m_data.event_zvertex);
  m_tree->Branch("event_centrality", &m_data.event_centrality);
  m_tree->Branch("sepd_totalcharge", &m_data.sepd_totalcharge);
  m_tree->Branch("sepd_channel", &m_data.sepd_channel);
  m_tree->Branch("sepd_charge", &m_data.sepd_charge);
  m_tree->Branch("sepd_phi", &m_data.sepd_phi);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int sEPD_TreeGen::process_event_check(PHCompositeNode *topNode)
{
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");

  if (!vertexmap)
  {
    std::cout << PHWHERE << "GlobalVertexMap Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!vertexmap->empty())
  {
    GlobalVertex *vtx = vertexmap->begin()->second;
    m_data.event_zvertex = vtx->get_z();
  }

  MinimumBiasInfo *m_mb_info = findNode::getClass<MinimumBiasInfo>(topNode, "MinimumBiasInfo");
  if (!m_mb_info)
  {
    std::cout << PHWHERE << "MinimumBiasInfo Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // skip event if not minimum bias
  if (!m_mb_info->isAuAuMinimumBias())
  {
    if (Verbosity() > 2)
    {
      std::cout << "Event: " << m_data.event_id << ", Not Min Bias, Skipping" << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // skip event if zvtx is too large
  if (std::abs(m_data.event_zvertex) >= m_cuts.m_zvtx_max)
  {
    if (Verbosity() > 2)
    {
      std::cout << "Event: " << m_data.event_id << ", Z: " << m_data.event_zvertex << " cm, Skipping" << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int sEPD_TreeGen::process_centrality(PHCompositeNode *topNode)
{
  CentralityInfo *centInfo = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  if (!centInfo)
  {
    std::cout << PHWHERE << "CentralityInfo Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_data.event_centrality = centInfo->get_centile(CentralityInfo::PROP::mbd_NS) * 100;

  // skip event if centrality is too peripheral
  if (!std::isfinite(m_data.event_centrality) || m_data.event_centrality < 0 || m_data.event_centrality >= m_cuts.m_cent_max)
  {
    if(Verbosity() > 2)
    {
        std::cout << "Event: " << m_data.event_id << ", Centrality: " << m_data.event_centrality << ", Skipping" << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int sEPD_TreeGen::process_sEPD(PHCompositeNode *topNode)
{
  TowerInfoContainer *towerinfosEPD = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_SEPD");
  if (!towerinfosEPD)
  {
    std::cout << PHWHERE << "TOWERINFO_CALIB_SEPD Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  EpdGeom *epdgeom = findNode::getClass<EpdGeom>(topNode, "TOWERGEOM_EPD");
  if (!epdgeom)
  {
    std::cout << PHWHERE << "TOWERGEOM_EPD Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // sepd
  unsigned int sepd_channels = towerinfosEPD->size();

  if(sepd_channels != QVecShared::sepd_channels)
  {
    if (Verbosity() > 2)
    {
      std::cout << "Event: " << m_data.event_id << ", SEPD Channels = " << sepd_channels << " != " << QVecShared::sepd_channels << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_data.sepd_totalcharge = 0;

  for (unsigned int channel = 0; channel < sepd_channels; ++channel)
  {
    unsigned int key = TowerInfoDefs::encode_epd(channel);

    TowerInfo *tower = towerinfosEPD->get_tower_at_channel(channel);

    double charge = tower->get_energy();
    bool isZS = tower->get_isZS();
    double phi = epdgeom->get_phi(key);

    // exclude ZS
    // exclude Nmips
    if (isZS || charge < m_cuts.m_sepd_charge_min)
    {
      continue;
    }

    m_data.sepd_channel.push_back(channel);
    m_data.sepd_charge.push_back(charge);
    m_data.sepd_phi.push_back(phi);

    m_data.sepd_totalcharge += charge;

    hSEPD_Charge->Fill(channel, charge);
  }

  h2SEPD_totalcharge_centrality->Fill(m_data.sepd_totalcharge, m_data.event_centrality);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int sEPD_TreeGen::process_event(PHCompositeNode *topNode)
{
  EventHeader *eventInfo = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if (!eventInfo)
  {
    std::cout << PHWHERE << "EventHeader Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_data.event_id = eventInfo->get_EvtSequence();

  if (Verbosity() > 1 && m_event % PROGRESS_PRINT_INTERVAL == 0)
  {
    std::cout << "Progress: " << m_event << ", Global: " << m_data.event_id << std::endl;
  }
  ++m_event;

  int ret = process_event_check(topNode);
  if (ret)
  {
    return ret;
  }

  ret = process_centrality(topNode);
  if (ret)
  {
    return ret;
  }

  ret = process_sEPD(topNode);
  if (ret)
  {
    return ret;
  }

  // Fill the TTree
  m_tree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int sEPD_TreeGen::ResetEvent([[maybe_unused]] PHCompositeNode *topNode)
{
  // Event
  m_data.event_id = -1;
  m_data.event_zvertex = 9999;
  m_data.event_centrality = 9999;

  // sEPD
  m_data.sepd_totalcharge = 0;
  m_data.sepd_channel.clear();
  m_data.sepd_charge.clear();
  m_data.sepd_phi.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int sEPD_TreeGen::End([[maybe_unused]] PHCompositeNode *topNode)
{
  std::cout << "sEPD_TreeGen::End" << std::endl;

  TFile output(m_outfile_name.c_str(), "recreate");
  output.cd();

  hSEPD_Charge->Write();
  h2SEPD_totalcharge_centrality->Write();

  output.Close();

  // TTree
  m_output->cd();
  m_tree->Write();
  m_output->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}
