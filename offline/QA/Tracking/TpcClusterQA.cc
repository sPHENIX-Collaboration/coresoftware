#include "TpcClusterQA.h"

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <tpc/TpcDistortionCorrectionContainer.h>
#include <tpc/TpcGlobalPositionWrapper.h>

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

namespace
{
  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track)
  {
    std::vector<TrkrDefs::cluskey> out;
    for (const auto& seed : {track->get_silicon_seed(), track->get_tpc_seed()})
    {
      if (seed)
      {
        std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
      }
    }

    return out;
  }
}

//____________________________________________________________________________..
TpcClusterQA::TpcClusterQA(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int TpcClusterQA::InitRun(PHCompositeNode *topNode)
{
  // find tpc geometry
  auto geomContainer =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geomContainer)
  {
    std::cout << PHWHERE << " unable to find DST node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  // TPC has 3 regions, inner, mid and outer
  std::vector<int> region_layer_low = {7, 23, 39};
  std::vector<int> region_layer_high = {22, 38, 54};

  // make a layer to region multimap
  const auto range = geomContainer->get_begin_end();
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    m_layers.insert(iter->first);

    for (int region = 0; region < 3; ++region)
    {
      if (iter->first >= region_layer_low[region] && iter->first <= region_layer_high[region])
      {
        m_layerRegionMap.insert(std::make_pair(iter->first, region));
      }
    }
  }
  
  if (m_residQA)
  {
    m_dccStatic = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerStatic");
    if (m_dccStatic)
    {
      std::cout << PHWHERE << "  found static TPC distortion correction container" << std::endl;
    }
    m_dccAverage = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerAverage");
    if (m_dccAverage)
    {
      std::cout << PHWHERE << "  found average TPC distortion correction container" << std::endl;
    }
    m_dccFluctuation = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerFluctuation");
    if (m_dccFluctuation)
    {
      std::cout << PHWHERE << "  found fluctuation TPC distortion correction container" << std::endl;
    }
  }

  createHistos();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcClusterQA::process_event(PHCompositeNode *topNode)
{ 
  auto clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!clusterContainer)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  auto geomContainer =
              findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  auto hitmap = findNode::getClass<TrkrHitSetContainer>(topNode,"TRKR_HITSET");
  if(!hitmap)
    {
      std::cout << PHWHERE << "No hitmap found, bailing" << std::endl;
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


  TH2 *h_totalclusters = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "stotal_clusters")));
  TH2 *h_clusterssector = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "ncluspersector")));
  TH2F *h_hitpositions = dynamic_cast<TH2F *>(hm->getHisto(std::string(getHistoPrefix() +"hit_positions")));
  TH1F *h_hitzpositions_side0 = dynamic_cast<TH1F *>(hm->getHisto(std::string(getHistoPrefix() +"hitz_positions_side0")));
  TH1F *h_hitzpositions_side1 = dynamic_cast<TH1F *>(hm->getHisto(std::string(getHistoPrefix() +"hitz_positions_side1")));

  struct HistoList
  {
    TH1 *crphisize_side0 = nullptr;
    TH1 *crphisize_side1 = nullptr;
    TH1 *czsize = nullptr;
    TH1 *crphierr = nullptr;
    TH1 *czerr = nullptr;
    TH1 *cedge = nullptr;
    TH1 *coverlap = nullptr;
    TH1 *cxposition_side0 = nullptr;
    TH1 *cxposition_side1 = nullptr;
    TH1 *cyposition_side0 = nullptr;
    TH1 *cyposition_side1 = nullptr;
    TH1 *czposition_side0 = nullptr;
    TH1 *czposition_side1 = nullptr;
  };

  int hitsetkeynum = 0;
  using HistoMap = std::map<int, HistoList>;
  HistoMap histos;

  for (auto &region : {0, 1, 2})
  {
    HistoList hist;
    hist.crphisize_side0 = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sphisize_side0_%i") % getHistoPrefix() % region).str()));
    hist.crphisize_side1 = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sphisize_side1_%i") % getHistoPrefix() % region).str()));
    hist.czsize = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%szsize_%i") % getHistoPrefix() % region).str()));
    hist.crphierr = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%srphi_error_%i") % getHistoPrefix() % region).str()));
    hist.czerr = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sz_error_%i") % getHistoPrefix() % region).str()));
    hist.cedge = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sclusedge_%i") % getHistoPrefix() % region).str()));
    hist.coverlap = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sclusoverlap_%i") % getHistoPrefix() % region).str()));
    
      
    hist.cxposition_side0 = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sclusxposition_side0_%i") % getHistoPrefix() % region).str()));
    hist.cxposition_side1 = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sclusxposition_side1_%i") % getHistoPrefix() % region).str()));
        
    hist.cyposition_side0 = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sclusyposition_side0_%i") % getHistoPrefix() % region).str()));
    hist.cyposition_side1 = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sclusyposition_side1_%i") % getHistoPrefix() % region).str()));
      
    hist.czposition_side0 = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%scluszposition_side0_%i") % getHistoPrefix() % region).str()));
    hist.czposition_side1 = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%scluszposition_side1_%i") % getHistoPrefix() % region).str()));
      
    histos.insert(std::make_pair(region, hist));
  }
  auto fill = [](TH1 *h, float val)
  { if (h) { h->Fill(val); 
} };

TrkrHitSetContainer::ConstRange all_hitsets = hitmap->getHitSets(TrkrDefs::TrkrId::tpcId);
 for (TrkrHitSetContainer::ConstIterator hitsetiter = all_hitsets.first;
      hitsetiter != all_hitsets.second;
      ++hitsetiter)
   {
    auto hitsetkey = hitsetiter->first;
    TrkrHitSet* hitset = hitsetiter->second;
    if(TrkrDefs::getTrkrId(hitsetkey) != TrkrDefs::TrkrId::tpcId)
      {
	continue;
      }
    int hitlayer = TrkrDefs::getLayer(hitsetkey);
    //auto sector = TpcDefs::getSectorId(hitsetkey);
    auto m_side = TpcDefs::getSide(hitsetkey);
    auto hitrangei = hitset->getHits();
    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
         hitr != hitrangei.second;
         ++hitr)
    { 
      auto hitkey = hitr->first;
      //auto hit = hitr->second;
      //auto adc = hit->getAdc();
      auto hitpad = TpcDefs::getPad(hitkey);
      auto m_hittbin = TpcDefs::getTBin(hitkey);
      //Check TrackResiduals.cc
      auto geoLayer = geomContainer->GetLayerCellGeom(hitlayer);
      auto phi = geoLayer->get_phicenter(hitpad);
      auto radius = geoLayer->get_radius();
      auto m_hitgx = radius * std::cos(phi);
      auto m_hitgy = radius * std::sin(phi);
      float AdcClockPeriod = geoLayer->get_zstep();
      float m_zdriftlength = m_hittbin * tGeometry->get_drift_velocity() * AdcClockPeriod;
      double NZBinsSide = 249;  // physical z bins per TPC side
      double tdriftmax = AdcClockPeriod * NZBinsSide;
      auto m_hitgz = (tdriftmax * tGeometry->get_drift_velocity()) - m_zdriftlength;
      if (m_side == 0)
      {
          m_hitgz *= -1;
      }
      //geoLayer->identify(std::cout);
      h_hitpositions->Fill(m_hitgx,m_hitgy);
      if(m_side==0){h_hitzpositions_side0->Fill(m_hitgz);}
      if(m_side==1){h_hitzpositions_side1->Fill(m_hitgz);}
    }

   }
  float nclusperevent[24] = {0};
  for (auto &hsk : clusterContainer->getHitSetKeys(TrkrDefs::TrkrId::tpcId))
  {
    int numclusters = 0;
    auto range = clusterContainer->getClusters(hsk);
    int sector = TpcDefs::getSectorId(hsk);
    int side = TpcDefs::getSide(hsk);
    if (side > 0)
    {
      sector += 12;
    }
    for (auto iter = range.first; iter != range.second; ++iter)
    { 
      const auto cluskey = iter->first;
      const auto cluster = iter->second; //auto cluster = clusters->findCluster(key);
      auto glob = tGeometry->getGlobalPosition(cluskey, cluster);
      auto sclusgx = glob.x();
      auto sclusgy = glob.y();
      auto sclusgz = glob.z();
        
        
      const auto it = m_layerRegionMap.find(TrkrDefs::getLayer(cluskey));
      int region = it->second;
      const auto hiter = histos.find(region);
      if (hiter == histos.end())
      {
        continue;
      }
      fill(hiter->second.czsize, cluster->getZSize());
      fill(hiter->second.crphierr, cluster->getRPhiError());
      fill(hiter->second.czerr, cluster->getZError());
      fill(hiter->second.cedge, cluster->getEdge());
      fill(hiter->second.coverlap, cluster->getOverlap());
        
      if(side==0){
        fill(hiter->second.crphisize_side0, cluster->getPhiSize());
        fill(hiter->second.cxposition_side0, sclusgx);
        fill(hiter->second.cyposition_side0, sclusgy);
        fill(hiter->second.czposition_side0, sclusgz);
      }
      if(side==1){
        fill(hiter->second.crphisize_side1, cluster->getPhiSize());
        fill(hiter->second.cxposition_side1, sclusgx);
        fill(hiter->second.cyposition_side1, sclusgy);
        fill(hiter->second.czposition_side1, sclusgz);
      }
        
      numclusters++;
    }

    nclusperevent[sector] += numclusters;
    h_totalclusters->Fill(hitsetkeynum, numclusters);
    m_totalClusters += numclusters;
    hitsetkeynum++;
  }
  for (int i = 0; i < 24; i++)
  {
    h_clusterssector->Fill(i, nclusperevent[i]);
    m_clustersPerSector[i] += nclusperevent[i];
  }
  
  if (m_residQA)
  {
    auto tpcseedmap = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
    auto trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
    auto clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    auto geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    auto svtxseedmap = findNode::getClass<TrackSeedContainer>(topNode, "SvtxTrackSeedContainer");

    if (!tpcseedmap or !trackmap or !clustermap or !svtxseedmap or !geometry)
    {
      std::cout << "Missing node, can't continue" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;;
    }

    TH1I *h_ntpc = dynamic_cast<TH1I *>(hm->getHisto(std::string(getHistoPrefix() +"ntpc")));

    struct PhiHistoList
    {
      TH1 *cphisize1pT_side0 = nullptr;
      TH1 *cphisize1pT_side1 = nullptr;
      TH1 *cphisizegeq1pT_side0 = nullptr;
      TH1 *cphisizegeq1pT_side1 = nullptr;
    };

    using PhiHistoMap = std::map<int, PhiHistoList>;
    PhiHistoMap phihistos;
    for (auto &region : {0, 1, 2})
    {
      PhiHistoList phihist;
    
      phihist.cphisize1pT_side0 = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sclusphisize1pT_side0_%i") % getHistoPrefix() % region).str()));
      phihist.cphisize1pT_side1 = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sclusphisize1pT_side1_%i") % getHistoPrefix() % region).str()));
    
      phihist.cphisizegeq1pT_side0 = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sclusphisizegeq1pT_side0_%i") % getHistoPrefix() % region).str()));
      phihist.cphisizegeq1pT_side1 = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sclusphisizegeq1pT_side1_%i") % getHistoPrefix() % region).str()));
    
      phihistos.insert(std::make_pair(region, phihist));
    }

    std::set<unsigned int> tpc_seed_ids;
    for (const auto& [key, track] : *trackmap)
    {
      if (!track)
      {
        continue;
      }
      m_px = track->get_px();
      m_py = track->get_py();
      m_pt = std::sqrt(m_px*m_px + m_py*m_py); 
  
      m_ntpc = 0;
      m_region.clear();
      m_clusgz.clear();
      m_cluslayer.clear();
      m_clusphisize.clear();
      m_cluszsize.clear();
      for (const auto& ckey : get_cluster_keys(track))
      {
        TrkrCluster* cluster = clustermap->findCluster(ckey);
        Acts::Vector3 clusglob;
        if (TrkrDefs::getTrkrId(ckey) == TrkrDefs::tpcId)
        {
          clusglob = TpcGlobalPositionWrapper::getGlobalPositionDistortionCorrected(ckey, cluster, geometry, track->get_crossing(), 
                                                                                    m_dccStatic, m_dccAverage, m_dccFluctuation); //NEED TO DEFINE THESE
        }
        else
        {
          clusglob = geometry->getGlobalPosition(ckey, cluster);
        }
        switch (TrkrDefs::getTrkrId(ckey))
        {      
          case TrkrDefs::tpcId:
            m_ntpc++;
            break; 
        }
        const auto it = m_layerRegionMap.find(TrkrDefs::getLayer(ckey));
        int region = it->second;
        m_region.push_back(region);
        m_clusgz.push_back(clusglob.z());
        m_cluslayer.push_back(TrkrDefs::getLayer(ckey));
        m_clusphisize.push_back(cluster->getPhiSize());
        m_cluszsize.push_back(cluster->getZSize()); 
      }
      
      if (m_pt > 0.8)
      {
	h_ntpc->Fill(m_ntpc);
      }
     
      int nClus = m_cluslayer.size(); 
      for (int cl = 0; cl < nClus; cl++) 
      {
        if (m_pt > 0.8 && m_ntpc > 25)
        {
          if (m_clusphisize[cl] == 1 && m_cluszsize[cl] > 1)
          {
            if (m_clusgz[cl] < 0.)
            {
              const auto hiter = phihistos.find(m_region[cl]);
              if (hiter == phihistos.end())
              {
                continue;
              }
              fill(hiter->second.cphisize1pT_side0, m_pt); 
            }
            else if (m_clusgz[cl] > 0.)
            {
              const auto hiter = phihistos.find(m_region[cl]);
              if (hiter == phihistos.end())
              {
                continue;
              }
              fill(hiter->second.cphisize1pT_side1, m_pt); 
            } 
          }
          if (m_clusphisize[cl] >= 1 && m_cluszsize[cl] > 1)
          {
            if (m_clusgz[cl] < 0.)
            {
              const auto hiter = phihistos.find(m_region[cl]);
              if (hiter == phihistos.end())
              {
                continue;
              }
              fill(hiter->second.cphisizegeq1pT_side0, m_pt); 
            }
            else if (m_clusgz[cl] > 0.)
            {
              const auto hiter = phihistos.find(m_region[cl]);
              if (hiter == phihistos.end())
              {
                continue;
              }
              fill(hiter->second.cphisizegeq1pT_side1, m_pt); 
            } 
          }
        } 
      }
    }  
  }   
  
  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcClusterQA::EndRun(const int /*runnumber*/)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

std::string TpcClusterQA::getHistoPrefix() const
{
  return std::string("h_") + Name() + std::string("_");
}
void TpcClusterQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);
 
  {
    auto h = new TH2F(std::string(getHistoPrefix() + "ncluspersector").c_str(),
                      "TPC Clusters per event per sector", 24, 0, 24, 1000, 0, 1000);
    h->GetXaxis()->SetTitle("Sector number");
    h->GetYaxis()->SetTitle("Clusters per event");
    hm->registerHisto(h);
  }
  for (auto &region : {0, 1, 2})
  {
    {
      auto h = new TH1F((boost::format("%sphisize_side0_%i") % getHistoPrefix() % region).str().c_str(),
                        (boost::format("TPC (side 0) cluster #phi size region_%i") % region).str().c_str(), 10, 0, 10);
      h->GetXaxis()->SetTitle("Cluster #phi_{size}");
      hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%sphisize_side1_%i") % getHistoPrefix() % region).str().c_str(),
                        (boost::format("TPC (side 1) cluster #phi size region_%i") % region).str().c_str(), 10, 0, 10);
      h->GetXaxis()->SetTitle("Cluster #phi_{size}");
      hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%szsize_%i") % getHistoPrefix() % region).str().c_str(),
                        (boost::format("TPC cluster z size region_%i") % region).str().c_str(), 10, 0, 10);
      h->GetXaxis()->SetTitle("Cluster z_{size}");
      hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%srphi_error_%i") % getHistoPrefix() % region).str().c_str(),
                        (boost::format("TPC r#Delta#phi error region_%i") % region).str().c_str(), 100, 0, 0.075);
      h->GetXaxis()->SetTitle("r#Delta#phi error [cm]");
      hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%sz_error_%i") % getHistoPrefix() % region).str().c_str(),
                        (boost::format("TPC z error region_%i") % region).str().c_str(), 100, 0, 0.18);
      h->GetXaxis()->SetTitle("z error [cm]");
      hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%sclusedge_%i") % getHistoPrefix() % region).str().c_str(),
                        (boost::format("TPC hits on edge region_%i") % region).str().c_str(), 30, 0, 30);
      h->GetXaxis()->SetTitle("Cluster edge");
      hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%sclusoverlap_%i") % getHistoPrefix() % region).str().c_str(),
                        (boost::format("TPC clus overlap region_%i") % region).str().c_str(), 30, 0, 30);
      h->GetXaxis()->SetTitle("Cluster overlap");
      hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%sclusxposition_side0_%i") % getHistoPrefix() % region).str().c_str(),
                          (boost::format("TPC cluster x position side 0 region_%i") % region).str().c_str(), 210*2, -105, 105);
        h->GetXaxis()->SetTitle("x (cm)");
        hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%sclusxposition_side1_%i") % getHistoPrefix() % region).str().c_str(),
                            (boost::format("TPC cluster x position side 1 region_%i") % region).str().c_str(), 210*2, -105, 105);
      h->GetXaxis()->SetTitle("x (cm)");
      hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%sclusyposition_side0_%i") % getHistoPrefix() % region).str().c_str(),
                            (boost::format("TPC cluster y position side 0 region_%i") % region).str().c_str(), 210*2, -105, 105);
      h->GetXaxis()->SetTitle("y (cm)");
      hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%sclusyposition_side1_%i") % getHistoPrefix() % region).str().c_str(),
                              (boost::format("TPC cluster y position side 1 region_%i") % region).str().c_str(), 210*2, -105, 105);
      h->GetXaxis()->SetTitle("y (cm)");
      hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%scluszposition_side0_%i") % getHistoPrefix() % region).str().c_str(),
                            (boost::format("TPC cluster z position side 0 region_%i") % region).str().c_str(), 210*2, -105, 105);
      h->GetXaxis()->SetTitle("z (cm)");
      hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%scluszposition_side1_%i") % getHistoPrefix() % region).str().c_str(),
                              (boost::format("TPC cluster z position side 1 region_%i") % region).str().c_str(), 210*2, -105, 105);
      h->GetXaxis()->SetTitle("z (cm)");
      hm->registerHisto(h);
    }
    if (m_residQA)
    {
      auto h = new TH1F((boost::format("%sclusphisize1pT_side0_%i") % getHistoPrefix() % region).str().c_str(),
                              (boost::format("TPC Cluster Phi Size == 1, side 0, region_%i") % region).str().c_str(), 4, 0.8, 3.2);
      h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      hm->registerHisto(h);
    }
    if (m_residQA)
    {
      auto h = new TH1F((boost::format("%sclusphisize1pT_side1_%i") % getHistoPrefix() % region).str().c_str(),
                              (boost::format("TPC Cluster Phi Size == 1, side 1, region_%i") % region).str().c_str(), 4, 0.8, 3.2);
      h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      hm->registerHisto(h);
    }
    if (m_residQA)
    {
      auto h = new TH1F((boost::format("%sclusphisizegeq1pT_side0_%i") % getHistoPrefix() % region).str().c_str(),
                              (boost::format("TPC Cluster Phi Size >= 1, side 0, region_%i") % region).str().c_str(), 4, 0.8, 3.2);
      h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      hm->registerHisto(h);
    }
    if (m_residQA)
    {
      auto h = new TH1F((boost::format("%sclusphisizegeq1pT_side1_%i") % getHistoPrefix() % region).str().c_str(),
                              (boost::format("TPC Cluster Phi Size >= 1, side 1, region_%i") % region).str().c_str(), 4, 0.8, 3.2);
      h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      hm->registerHisto(h);
    }
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "stotal_clusters").c_str(),
                      "TPC clusters per hitsetkey", 1152, 0, 1152, 10000, 0, 10000);
    h->GetXaxis()->SetTitle("Hitsetkey number");
    h->GetYaxis()->SetTitle("Number of clusters");
    hm->registerHisto(h);
  }
  
    
  {
    auto h = new TH2F(std::string(getHistoPrefix()+"hit_positions").c_str(),
                                           "Histogram of hit x y positions", 160, 0, 80, 160, 0, 80);
    h->GetXaxis()->SetTitle("x (cm)");
    h->GetYaxis()->SetTitle("y (cm)");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1F(std::string(getHistoPrefix()+"hitz_positions_side0").c_str(),
                                             "Histogram of hit z positions side=0", 105*4, -105, 105);
    h->GetXaxis()->SetTitle("z (cm)");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1F(std::string(getHistoPrefix()+"hitz_positions_side1").c_str(),
                                             "Histogram of hit z positions side=1", 105*4, -105, 105);
    h->GetXaxis()->SetTitle("z (cm)");
    hm->registerHisto(h);
  }
  if (m_residQA)
  {
    auto h = new TH1I(std::string(getHistoPrefix()+"ntpc").c_str(),
                                             "Clusters per Track Full Detector", 50, 0, 50);
    h->GetXaxis()->SetTitle("nClusters/Track");
    hm->registerHisto(h);
  }
  return;
}
