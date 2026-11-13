#include "sEPD_TreeGen.h"
#include "QVecDefs.h"
#include "EventPlaneData.h"

// -- Calo
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

// -- Vtx
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

// -- MB
#include <calotrigger/MinimumBiasInfo.h>
#include <calotrigger/TriggerAnalyzer.h>
#include <calotrigger/TriggerRunInfo.h>

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
#include <phool/recoConsts.h>
#include <phool/RunnumberRange.h>

// -- ROOT
#include <TProfile.h>
#include <TH2.h>

// -- c++
#include <iomanip>
#include <iostream>

//____________________________________________________________________________..
sEPD_TreeGen::sEPD_TreeGen(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
// Defined here so unique_ptr<TriggerAnalyzer> sees the complete type
sEPD_TreeGen::~sEPD_TreeGen() = default;

//____________________________________________________________________________..
int sEPD_TreeGen::Init(PHCompositeNode *topNode)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  if (Verbosity() > 0)
  {
    se->Print("NODETREE");
  }

  m_triggerAnalyzer = std::make_unique<TriggerAnalyzer>();

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
int sEPD_TreeGen::InitRun(PHCompositeNode *topNode)
{
  TriggerRunInfo *triggerruninfo = findNode::getClass<TriggerRunInfo>(topNode, "TriggerRunInfo");
  if (!triggerruninfo)
  {
    std::cout << PHWHERE << "TriggerRunInfo Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // MBD trigger names for each run period (gl1_triggernames)
  const int runnumber = recoConsts::instance()->get_IntFlag("RUNNUMBER");
  if (runnumber >= RunnumberRange::RUN2AUAU_FIRST && runnumber <= RunnumberRange::RUN2AUAU_LAST)
  {
    // No 150 cm trigger in early Run 2 Au+Au, so use the 60 cm trigger present for the full run period
    m_trig_mbd_narrow = "MBD N&S >= 2, vtx < 10 cm";
    m_trig_mbd_wide = "MBD N&S >= 2, vtx < 60 cm";
  }
  else if (runnumber >= RunnumberRange::RUN3AUAU_FIRST && runnumber <= RunnumberRange::RUN3AUAU_LAST)
  {
    m_trig_mbd_narrow = "MBD N&S >= 2, vtx < 10 cm";
    m_trig_mbd_wide = "MBD N&S >= 2, vtx < 150 cm";
  }
  else if (runnumber >= RunnumberRange::RUN3PP_FIRST)
  {
    // Names changed from N&S >= 2 to N&S >= 1 starting with Run 3 p+p
    m_trig_mbd_narrow = "MBD N&S >= 1, vtx < 10 cm";
    m_trig_mbd_wide = "MBD N&S >= 1, vtx < 150 cm";
  }
  else
  {
    std::cout << PHWHERE << "Run " << runnumber << " is outside the supported run periods, no MBD trigger names defined." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // getTriggerBitByName returns 0 for an unknown name, so confirm the round trip
  auto resolve_bit = [triggerruninfo](const std::string &name) -> int
  {
    int bit = static_cast<int>(triggerruninfo->getTriggerBitByName(name));
    return (triggerruninfo->getTriggerName(bit) == name) ? bit : -1;
  };

  m_bit_mbd_narrow = resolve_bit(m_trig_mbd_narrow);
  m_bit_mbd_wide = resolve_bit(m_trig_mbd_wide);

  if (m_bit_mbd_narrow < 0 || m_bit_mbd_wide < 0)
  {
    std::cout << PHWHERE << "Run " << runnumber << ": MBD trigger not found in TriggerRunInfo: "
              << "\"" << m_trig_mbd_narrow << "\" -> " << m_bit_mbd_narrow << ", "
              << "\"" << m_trig_mbd_wide << "\" -> " << m_bit_mbd_wide << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (Verbosity() > 0)
  {
    std::cout << "sEPD_TreeGen::InitRun - Run " << runnumber << ": \"" << m_trig_mbd_narrow << "\" -> bit " << m_bit_mbd_narrow
              << ", \"" << m_trig_mbd_wide << "\" -> bit " << m_bit_mbd_wide << std::endl;
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

  // MBD Trigger
  // decodeTriggers leaves the previous event's trigger bits in place on failure
  if (m_triggerAnalyzer->decodeTriggers(topNode))
  {
    std::cout << PHWHERE << "Trigger decoding failed (TriggerRunInfo or GL1 packet missing), doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  bool mbd_trigger_fire = m_triggerAnalyzer->didTriggerFire(m_bit_mbd_narrow) || m_triggerAnalyzer->didTriggerFire(m_bit_mbd_wide);
  bool pass_MB = mbd_trigger_fire && m_mb_info->isAuAuMinimumBias();

  // skip event if not minimum bias
  if (!pass_MB)
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

  // skip event if centrality is bad or too peripheral
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
  TowerInfoContainer *towerinfosEPD = findNode::getClass<TowerInfoContainer>(topNode, m_inputNode);
  if (!towerinfosEPD)
  {
    std::cout << PHWHERE << m_inputNode << " Node missing, doing nothing." << std::endl;
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
  if (Verbosity() <= 2)
  {
    return;
  }

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

  return Fun4AllReturnCodes::EVENT_OK;
}
