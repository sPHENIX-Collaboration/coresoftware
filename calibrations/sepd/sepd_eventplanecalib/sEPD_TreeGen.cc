#include "sEPD_TreeGen.h"
#include "QVecDefs.h"
#include "EventPlaneData.h"

// -- c++
#include <iomanip>
#include <iostream>

// -- Calo
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

// -- Vtx
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

// -- MB
#include <calotrigger/MinimumBiasClassifier.h>
#include <calotrigger/MinimumBiasInfo.h>
#include <centrality/CentralityInfo.h>

// -- sEPD
#include <epd/EpdGeom.h>

// -- event
#include <ffaobjects/EventHeader.h>

// -- Fun4All
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

// -- Nodes
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

//____________________________________________________________________________..
sEPD_TreeGen::sEPD_TreeGen(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int sEPD_TreeGen::Init(PHCompositeNode *topNode)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Print("NODETREE");

  unsigned int bins_sepd_totalcharge{100};
  double sepd_totalcharge_low{0};
  double sepd_totalcharge_high{2e4};

  unsigned int bins_centrality{80};
  double centrality_low{-0.5};
  double centrality_high{79.5};

  hSEPD_Charge = new TProfile("hSEPD_Charge", "|z| < 10 cm and MB; Channel; Avg Charge", QVecShared::SEPD_CHANNELS, 0, QVecShared::SEPD_CHANNELS);
  hSEPD_Charge->Sumw2();

  h2SEPD_totalcharge_centrality = new TH2F("h2SEPD_totalcharge_centrality",
                                           "|z| < 10 cm and MB; sEPD Total Charge; Centrality [%]",
                                           bins_sepd_totalcharge, sepd_totalcharge_low, sepd_totalcharge_high,
                                           bins_centrality, centrality_low, centrality_high);

  se->registerHisto(hSEPD_Charge);
  se->registerHisto(h2SEPD_totalcharge_centrality);

  PHNodeIterator node_itr(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(node_itr.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cout << PHWHERE << "DST node missing, cannot attach EventPlaneData." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  EventPlaneData *evtdata = findNode::getClass<EventPlaneData>(topNode, "EventPlaneData");
  if (!evtdata)
  {
    evtdata = new EventPlaneData();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(evtdata, "EventPlaneData", "PHObject");
    dstNode->addNode(newNode);
  }

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

  if (vertexmap->empty())
  {
    if (Verbosity() > 1)
    {
      std::cout << PHWHERE << "GlobalVertexMap Empty, Skipping Event: " << m_data.event_id << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  GlobalVertex *vtx = vertexmap->begin()->second;
  double zvtx = vtx->get_z();

  MinimumBiasInfo *m_mb_info = findNode::getClass<MinimumBiasInfo>(topNode, "MinimumBiasInfo");
  if (!m_mb_info)
  {
    std::cout << PHWHERE << "MinimumBiasInfo Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // skip event if not minimum bias
  if (!m_mb_info->isAuAuMinimumBias())
  {
    if (Verbosity() > 1)
    {
      std::cout << "Event: " << m_data.event_id << ", Not Min Bias, Skipping" << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // skip event if zvtx is too large
  if (std::abs(zvtx) >= m_cuts.m_zvtx_max)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Event: " << m_data.event_id << ", Z: " << zvtx << " cm, Skipping" << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_evtdata->set_event_zvertex(zvtx);

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

  double cent = centInfo->get_centile(CentralityInfo::PROP::mbd_NS) * 100;

  // skip event if centrality is too peripheral
  if (!std::isfinite(cent) || cent < 0 || cent >= m_cuts.m_cent_max)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Event: " << m_data.event_id << ", Centrality: " << cent << ", Skipping" << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_evtdata->set_event_centrality(cent);
  m_data.event_centrality = cent;

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

  if(sepd_channels != QVecShared::SEPD_CHANNELS)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Event: " << m_data.event_id << ", SEPD Channels = " << sepd_channels << " != " << QVecShared::SEPD_CHANNELS << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  double sepd_totalcharge = 0;

  for (unsigned int channel = 0; channel < sepd_channels; ++channel)
  {
    TowerInfo *tower = towerinfosEPD->get_tower_at_channel(channel);

    if (!tower)
    {
      if (Verbosity() > 1)
      {
        std::cout << PHWHERE << "Null SEPD tower at channel " << channel << std::endl;
      }
      continue;
    }

    double charge = tower->get_energy();
    bool isZS = tower->get_isZS();

    // exclude ZS
    // exclude Nmips
    if (isZS || charge < m_cuts.m_sepd_charge_min)
    {
      continue;
    }

    m_evtdata->set_sepd_charge(channel, charge);

    sepd_totalcharge += charge;

    hSEPD_Charge->Fill(channel, charge);
  }

  m_evtdata->set_sepd_totalcharge(sepd_totalcharge);
  h2SEPD_totalcharge_centrality->Fill(sepd_totalcharge, m_data.event_centrality);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void sEPD_TreeGen::Print([[maybe_unused]] const std::string &what) const
{
  // Only execute if Verbosity is high enough
  if (Verbosity() <= 2) return;

  std::cout << "\n============================================================" << std::endl;
  std::cout << "sEPD_TreeGen::Print -> Event Data State" << std::endl;

  if (!m_evtdata)
  {
    std::cout << " [WARNING] m_evtdata is null." << std::endl;
    return;
  }

  // Verbosity > 2: Print basic scalars
  std::cout << "  Event ID:          " << m_evtdata->get_event_id() << std::endl;
  std::cout << "  Z-Vertex:          " << m_evtdata->get_event_zvertex() << " cm" << std::endl;
  std::cout << "  Centrality:        " << m_evtdata->get_event_centrality() << " %" << std::endl;
  std::cout << "  sEPD Total Charge: " << m_evtdata->get_sepd_totalcharge() << std::endl;

  // Verbosity > 3: Print channel arrays
  if (Verbosity() > 3)
  {
    std::cout << "  Active Towers (Charge > 0):" << std::endl;
    for (int i = 0; i < QVecShared::SEPD_CHANNELS; ++i)
    {
      double charge = m_evtdata->get_sepd_charge(i);
      if (charge > 0)
      {
        std::cout << "    Channel: " << std::setw(3) << i
                  << " | Charge: " << std::fixed << std::setprecision(4) << charge << std::endl;
      }
    }
  }
  std::cout << "============================================================\n" << std::endl;
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

  if (Verbosity() && m_event % PROGRESS_PRINT_INTERVAL == 0)
  {
    std::cout << "Progress: " << m_event << ", Global: " << m_data.event_id << std::endl;
  }
  ++m_event;

  m_evtdata = findNode::getClass<EventPlaneData>(topNode, "EventPlaneData");
  if (!m_evtdata)
  {
    std::cout << PHWHERE << "EventPlaneData Node missing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_evtdata->set_event_id(m_data.event_id);

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
  
  Print();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int sEPD_TreeGen::ResetEvent([[maybe_unused]] PHCompositeNode *topNode)
{
  // Event
  m_data.event_id = -1;
  m_data.event_centrality = 9999;

  // DST
  if (m_evtdata)
  {
    m_evtdata->Reset();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int sEPD_TreeGen::End([[maybe_unused]] PHCompositeNode *topNode)
{
  std::cout << "sEPD_TreeGen::End" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
