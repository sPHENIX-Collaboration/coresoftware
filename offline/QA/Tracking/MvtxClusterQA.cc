
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

  auto h_clusSize = dynamic_cast<TH1F *>(hm->getHisto((boost::format("%sclusterSize") % getHistoPrefix()).str()));
  auto h_clusPhi_incl = dynamic_cast<TH1F *>(hm->getHisto((boost::format("%sclusterPhi_incl") % getHistoPrefix()).str()));
  auto h_clusPhi_l0 = dynamic_cast<TH1F *>(hm->getHisto((boost::format("%sclusterPhi_l0") % getHistoPrefix()).str()));
  auto h_clusPhi_l1 = dynamic_cast<TH1F *>(hm->getHisto((boost::format("%sclusterPhi_l1") % getHistoPrefix()).str()));
  auto h_clusPhi_l2 = dynamic_cast<TH1F *>(hm->getHisto((boost::format("%sclusterPhi_l2") % getHistoPrefix()).str()));
  auto h_clusZ_clusPhi_l0 = dynamic_cast<TH2F *>(hm->getHisto((boost::format("%sclusterZ_clusPhi_l0") % getHistoPrefix()).str()));
  auto h_clusZ_clusPhi_l1 = dynamic_cast<TH2F *>(hm->getHisto((boost::format("%sclusterZ_clusPhi_l1") % getHistoPrefix()).str()));
  auto h_clusZ_clusPhi_l2 = dynamic_cast<TH2F *>(hm->getHisto((boost::format("%sclusterZ_clusPhi_l2") % getHistoPrefix()).str()));

  for (auto &hsk : clusterContainer->getHitSetKeys(TrkrDefs::TrkrId::mvtxId))
  {
    int numclusters = 0;
    auto range = clusterContainer->getClusters(hsk);
    auto layer = TrkrDefs::getLayer(hsk);
    auto stave = MvtxDefs::getStaveId(hsk);
    auto chip = MvtxDefs::getChipId(hsk);

    if (m_chipInfo)
    {
      auto h = dynamic_cast<TH2 *>(hm->getHisto((boost::format("%snclusperchip%i_%i_%i") % getHistoPrefix() % (int) layer % (int) stave % (int) chip).str()));
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        const auto cluskey = iter->first;
        const auto cluster = iter->second;
        auto globalpos = tGeometry->getGlobalPosition(cluskey, cluster);
        auto phi = atan2(globalpos(1), globalpos(0));
        auto clayer = TrkrDefs::getLayer(cluskey);
        h->Fill(cluster->getLocalY(), cluster->getLocalX());
        h_clusSize->Fill(cluster->getSize());
        h_clusPhi_incl->Fill(phi);
        if (clayer == 0)
        {
          h_clusPhi_l0->Fill(phi);
          h_clusZ_clusPhi_l0->Fill(globalpos(2), phi);
        }
        else if (clayer == 1)
        {
          h_clusPhi_l1->Fill(phi);
          h_clusZ_clusPhi_l1->Fill(globalpos(2), phi);
        }
        else if (clayer == 2)
        {
          h_clusPhi_l2->Fill(phi);
          h_clusZ_clusPhi_l2->Fill(globalpos(2), phi);
        }
        m_totalClusters++;
        numclusters++;
      }
      m_nclustersPerChip[(int) layer][(int) stave][(int) chip] += numclusters; 
    }
    else
    {
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        const auto cluskey = iter->first;
        const auto cluster = iter->second;
        auto globalpos = tGeometry->getGlobalPosition(cluskey, cluster);
        auto phi = atan2(globalpos(1), globalpos(0));
        auto clayer = TrkrDefs::getLayer(cluskey);
        h_clusSize->Fill(cluster->getSize());
        h_clusPhi_incl->Fill(phi);
        if (clayer == 0)
        {
          h_clusPhi_l0->Fill(phi);
          h_clusZ_clusPhi_l0->Fill(globalpos(2), phi);
        }
        else if (clayer == 1)
        {
          h_clusPhi_l1->Fill(phi);
          h_clusZ_clusPhi_l1->Fill(globalpos(2), phi);
        }
        else if (clayer == 2)
        {
          h_clusPhi_l2->Fill(phi);
          h_clusZ_clusPhi_l2->Fill(globalpos(2), phi);
        }
      }
    }
  }

  TrkrHitSetContainer::ConstRange hitsetrange = trkrHitSetContainer->getHitSets(TrkrDefs::TrkrId::mvtxId);

  auto h_occupancy = dynamic_cast<TH1F *>(hm->getHisto((boost::format("%schipOccupancy") % getHistoPrefix()).str()));
  auto h_strobe = dynamic_cast<TH1I *>(hm->getHisto((boost::format("%sstrobeTiming") % getHistoPrefix()).str()));
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first; hitsetitr != hitsetrange.second; ++hitsetitr)
  {
    int chip_hits = hitsetitr->second->size();
    float chip_occupancy = (float) chip_hits / (512*1024);
    chip_occupancy = 100*chip_occupancy;
    h_occupancy->Fill(chip_occupancy);
    int strobe = MvtxDefs::getStrobeId(hitsetitr->first);
    for (int i = 0; i < chip_hits; i++) 
    {
      h_strobe->Fill(strobe);
    }
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
  auto h_clusPhi_incl = new TH1F((boost::format("%sclusterPhi_incl") % getHistoPrefix()).str().c_str(),"MVTX Cluster Phi",320,-3.2,3.2);
  h_clusPhi_incl->GetXaxis()->SetTitle("Cluster (layer 0+1+2) #phi [rad]");
  h_clusPhi_incl->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_clusPhi_incl);
  auto h_clusPhi_l0 = new TH1F((boost::format("%sclusterPhi_l0") % getHistoPrefix()).str().c_str(),"MVTX Cluster Phi",320,-3.2,3.2);
  h_clusPhi_l0->GetXaxis()->SetTitle("Cluster (layer 0) #phi [rad]");
  h_clusPhi_l0->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_clusPhi_l0);
  auto h_clusPhi_l1 = new TH1F((boost::format("%sclusterPhi_l1") % getHistoPrefix()).str().c_str(),"MVTX Cluster Phi",320,-3.2,3.2);
  h_clusPhi_l1->GetXaxis()->SetTitle("Cluster (layer 1) #phi [rad]");
  h_clusPhi_l1->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_clusPhi_l1);
  auto h_clusPhi_l2 = new TH1F((boost::format("%sclusterPhi_l2") % getHistoPrefix()).str().c_str(),"MVTX Cluster Phi",320,-3.2,3.2);
  h_clusPhi_l2->GetXaxis()->SetTitle("Cluster (layer 2) #phi [rad]");
  h_clusPhi_l2->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_clusPhi_l2);
  auto h_clusZ_clusPhi_l0 = new TH2F((boost::format("%sclusterZ_clusPhi_l0") % getHistoPrefix()).str().c_str(),"MVTX Cluster Z vs Phi",300,-15,15,350,-3.5,3.5);
  h_clusZ_clusPhi_l0->GetXaxis()->SetTitle("Cluster (layer 2) Z [cm]");
  h_clusZ_clusPhi_l0->GetYaxis()->SetTitle("Cluster (layer 0) #phi [rad]");
  hm->registerHisto(h_clusZ_clusPhi_l0);
  auto h_clusZ_clusPhi_l1 = new TH2F((boost::format("%sclusterZ_clusPhi_l1") % getHistoPrefix()).str().c_str(),"MVTX Cluster Z vs Phi",300,-15,15,350,-3.5,3.5);
  h_clusZ_clusPhi_l1->GetXaxis()->SetTitle("Cluster (layer 2) Z [cm]");
  h_clusZ_clusPhi_l1->GetYaxis()->SetTitle("Cluster (layer 1) #phi [rad]");
  hm->registerHisto(h_clusZ_clusPhi_l1);
  auto h_clusZ_clusPhi_l2 = new TH2F((boost::format("%sclusterZ_clusPhi_l2") % getHistoPrefix()).str().c_str(),"MVTX Cluster Z vs Phi",300,-15,15,350,-3.5,3.5);
  h_clusZ_clusPhi_l2->GetXaxis()->SetTitle("Cluster (layer 2) Z [cm]");
  h_clusZ_clusPhi_l2->GetYaxis()->SetTitle("Cluster (layer 2) #phi [rad]");
  hm->registerHisto(h_clusZ_clusPhi_l2);
  auto h_strobe = new TH1I((boost::format("%sstrobeTiming") % getHistoPrefix()).str().c_str(),"MVTX Strobe Timing per Hit",32,-15,16); 
  h_strobe->GetXaxis()->SetTitle("Strobe BCO - GL1 BCO");
  h_strobe->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_strobe);

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
