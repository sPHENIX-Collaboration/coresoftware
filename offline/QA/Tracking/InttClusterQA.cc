#include "InttClusterQA.h"
#include <intt/CylinderGeomIntt.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/InttDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrDefs.h>

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TH1.h>
#include <TH2.h>

#include <boost/format.hpp>

//____________________________________________________________________________..
InttClusterQA::InttClusterQA(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int InttClusterQA::InitRun(PHCompositeNode * /*unused*/)
{
  for (auto &layer : {0, 1, 2, 3})
  {
    if (layer < 2)
    {
      m_layerLadderMap.insert(std::make_pair(layer, 12));
    }
    else
    {
      m_layerLadderMap.insert(std::make_pair(layer, 16));
    }
  }

  createHistos();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InttClusterQA::process_event(PHCompositeNode *topNode)
{
  auto clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!clusterContainer)
  {
    std::cout << PHWHERE << "No cluster container, bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto trkrHitSetContainer = findNode::getClass<TrkrHitSetContainerv1>(topNode, "TRKR_HITSET");
  if (!trkrHitSetContainer)
  {
    std::cout << PHWHERE << "No trkrhitset container, bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tGeometry)
  {
    std::cout << PHWHERE << "No acts geometry on node tree, bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  for (auto &hsk : clusterContainer->getHitSetKeys(TrkrDefs::TrkrId::inttId))
  {
    int numclusters = 0;
    auto range = clusterContainer->getClusters(hsk);
    auto layer = TrkrDefs::getLayer(hsk);
    auto ladderphiid = InttDefs::getLadderPhiId(hsk);
    auto sensor = InttDefs::getLadderZId(hsk);
    auto h_clusSize = dynamic_cast<TH1F *>(hm->getHisto((boost::format("%sclusterSize") % getHistoPrefix()).str()));

    if (m_sensorInfo)
    {  
      auto h = dynamic_cast<TH2 *>(hm->getHisto( boost::str(boost::format("%sncluspersensor%i_%i_%i") %getHistoPrefix() %(((int) layer) - 3) %((int) ladderphiid) %((int) sensor)).c_str()));

      for (auto iter = range.first; iter != range.second; ++iter)
      {
        // const auto cluskey = iter->first;
        const auto cluster = iter->second;
        h->Fill(cluster->getLocalY(), cluster->getLocalX());
        m_totalClusters++;
        numclusters++;
      }
      m_nclustersPerSensor[((int) layer) - 3][(int) ladderphiid][(int) sensor] += numclusters;
    }
    else
    {
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        const auto cluster = iter->second;
        h_clusSize->Fill(cluster->getSize()); 
      }
    }
  }
 
  TrkrHitSetContainer::ConstRange hitsetrange = trkrHitSetContainer->getHitSets(TrkrDefs::TrkrId::inttId);

  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first; hitsetitr != hitsetrange.second; ++hitsetitr)
  {
    int sensor_hits = hitsetitr->second->size();
    float sensor_occupancy = (float) sensor_hits / (128.*26.);
    auto h_occupancy = dynamic_cast<TH1F *>(hm->getHisto((boost::format("%ssensorOccupancy") % getHistoPrefix()).str()));
    h_occupancy->Fill(100.*sensor_occupancy);
  }

  m_event++;

  return Fun4AllReturnCodes::EVENT_OK;
}
int InttClusterQA::EndRun(const int /*runnumber*/)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

std::string InttClusterQA::getHistoPrefix() const
{
  return std::string("h_") + Name() + std::string("_");
}

void InttClusterQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  auto h_occupancy = new TH1F((boost::format("%ssensorOccupancy") % getHistoPrefix()).str().c_str(),"INTT Sensor Occupancy",100,0,5);
  h_occupancy->GetXaxis()->SetTitle("Sensor Occupancy [%]");
  h_occupancy->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_occupancy);
  auto h_clusSize = new TH1F((boost::format("%sclusterSize") % getHistoPrefix()).str().c_str(),"INTT Cluster Size",20,-0.5,19.5);
  h_clusSize->GetXaxis()->SetTitle("Cluster Size");
  h_clusSize->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_clusSize);

  if (m_sensorInfo)
  {
    for (const auto &[layer, ladders] : m_layerLadderMap)
    {
      for (int ladder = 0; ladder < ladders; ladder++)
      {
        //! 4 sensor on each ladder
        for (int sensor = 0; sensor < 4; sensor++)
        {
          auto h = new TH2F( boost::str(boost::format("%sncluspersensor%i_%i_%i") %getHistoPrefix()
	    			        %layer %ladder %sensor).c_str(),
                            "INTT clusters per sensor", 100, -5, 5, 1000, -1, 1);
          h->GetXaxis()->SetTitle("Local z [cm]");
          h->GetYaxis()->SetTitle("Local rphi [cm]");
          hm->registerHisto(h);
        }
      }
    }
  }

  return;
}
