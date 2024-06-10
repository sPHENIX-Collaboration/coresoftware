
#include "MvtxClusterQA.h"

#include <mvtx/CylinderGeom_Mvtx.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrClusterHitAssoc.h>
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
MvtxClusterQA::MvtxClusterQA(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int MvtxClusterQA::InitRun(PHCompositeNode *topNode)
{
  auto geomContainer = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  if (!geomContainer)
  {
    std::cout << PHWHERE
              << " CYLINDERGEOM_MVTX  node not found on node tree"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  for (auto &layer : {0, 1, 2})
  {
    auto layergeom = dynamic_cast<CylinderGeom_Mvtx *>(geomContainer->GetLayerGeom(layer));
    if (!layergeom)
    {
      std::cout << PHWHERE << "Did not get layergeom for layer "
                << layer << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    int nstaves = layergeom->get_N_staves();
    m_layerStaveMap.insert(std::make_pair(layer, nstaves));
  }

  createHistos();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MvtxClusterQA::process_event(PHCompositeNode *topNode)
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

  for (auto &hsk : clusterContainer->getHitSetKeys(TrkrDefs::TrkrId::mvtxId))
  {
    int numclusters = 0;
    auto range = clusterContainer->getClusters(hsk);
    auto layer = TrkrDefs::getLayer(hsk);
    auto stave = MvtxDefs::getStaveId(hsk);
    auto chip = MvtxDefs::getChipId(hsk);
    auto h_clusSize = dynamic_cast<TH1F *>(hm->getHisto((boost::format("%sclusterSize") % getHistoPrefix()).str()));

    if (m_chipInfo)
    {
      auto h = dynamic_cast<TH2 *>(hm->getHisto((boost::format("%snclusperchip%i_%i_%i") % getHistoPrefix() % (int) layer % (int) stave % (int) chip).str()));
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        // const auto cluskey = iter->first;
        const auto cluster = iter->second;
        h->Fill(cluster->getLocalY(), cluster->getLocalX());
        h_clusSize->Fill(cluster->getSize());
        m_totalClusters++;
        numclusters++;
      }
      m_nclustersPerChip[(int) layer][(int) stave][(int) chip] += numclusters; 
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

  TrkrHitSetContainer::ConstRange hitsetrange = trkrHitSetContainer->getHitSets(TrkrDefs::TrkrId::mvtxId);

  auto h_occupancy = dynamic_cast<TH1F *>(hm->getHisto((boost::format("%schipOccupancy") % getHistoPrefix()).str()));
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first; hitsetitr != hitsetrange.second; ++hitsetitr)
  {
    int chip_hits = hitsetitr->second->size();
    float chip_occupancy = (float) chip_hits / (512*1024);
    chip_occupancy = 100*chip_occupancy;
    h_occupancy->Fill(chip_occupancy);
  }

  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}
int MvtxClusterQA::EndRun(const int /*runnumber*/)
{
  
  return Fun4AllReturnCodes::EVENT_OK;
}
//____________________________________________________________________________..

std::string MvtxClusterQA::getHistoPrefix() const
{
  return std::string("h_") + Name() + std::string("_");
}

void MvtxClusterQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);
 
  auto h_occupancy = new TH1F((boost::format("%schipOccupancy") % getHistoPrefix()).str().c_str(),"MVTX Chip Occupancy",60,0,0.6); 
  h_occupancy->GetXaxis()->SetTitle("Chip Occupancy [%]");
  h_occupancy->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_occupancy);
  auto h_clusSize = new TH1F((boost::format("%sclusterSize") % getHistoPrefix()).str().c_str(),"MVTX Cluster Size",50,-0.5,49.5); 
  h_clusSize->GetXaxis()->SetTitle("Cluster Size");
  h_clusSize->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_clusSize);

  if (m_chipInfo)
  {
    for (const auto &[layer, nstave] : m_layerStaveMap)
    {
      for (int stave = 0; stave < nstave; stave++)
      {
        //! 9 chips on each stave
        for (int chip = 0; chip < 9; chip++)
        {
          auto h = new TH2F((boost::format("%snclusperchip%i_%i_%i") % getHistoPrefix() % layer % stave % chip).str().c_str(),
                            "MVTX clusters per chip", 2000, -2, 2, 2000, -1, 1);
          h->GetXaxis()->SetTitle("Local z [cm]");
          h->GetYaxis()->SetTitle("Local rphi [cm]");
          hm->registerHisto(h);
        }
      }
    }
  }

  return;
}
