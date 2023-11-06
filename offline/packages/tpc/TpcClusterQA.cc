
#include "TpcClusterQA.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>

#include <qautils/QAHistManagerDef.h>

#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>

#include <qautils/QAUtil.h>
#include <qautils/QAHistManagerDef.h>

#include <TH1F.h>
#include <TH2F.h>

//____________________________________________________________________________..
TpcClusterQA::TpcClusterQA(const std::string &name):
 SubsysReco(name)
{
}

//____________________________________________________________________________..
TpcClusterQA::~TpcClusterQA()
{
}

//____________________________________________________________________________..
int TpcClusterQA::Init(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
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
        m_layerRegionMap.insert(std::make_pair(iter->first, region));
    }
  }

  createHistos();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcClusterQA::process_event(PHCompositeNode* topNode)
{


  auto clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if(!clusterContainer)
    {
      std::cout << PHWHERE << "No cluster container, bailing"<<std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  auto tGeometry = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
  if(!tGeometry)
    {
      std::cout << PHWHERE << "No acts geometry on node tree, bailing"<<std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  
  TH2* h_totalclusters = dynamic_cast<TH2*>(hm->getHisto(Form("%stotal_clusters", getHistoPrefix().c_str())));
  
  struct HistoList
  {
    TH1 *crphisize = nullptr;
    TH1 *czsize = nullptr;
    TH1 *crphierr = nullptr;
    TH1 *czerr = nullptr;
    TH1 *cedge = nullptr;
    TH1 *coverlap = nullptr;
  };

  int hitsetkeynum = 0;
  using HistoMap = std::map<int, HistoList>;
  HistoMap histos;
  
  for(auto& region : {0,1,2})
    {
      HistoList hist;
      hist.crphisize = dynamic_cast<TH1*>(hm->getHisto(Form("%sphisize_%i", getHistoPrefix().c_str(), region)));
      hist.czsize = dynamic_cast<TH1*>(hm->getHisto(Form("%szsize_%i", getHistoPrefix().c_str(), region)));
      hist.crphierr = dynamic_cast<TH1*>(hm->getHisto(Form("%srphi_error_%i", getHistoPrefix().c_str(), region)));
      hist.czerr = dynamic_cast<TH1*>(hm->getHisto(Form("%sz_error_%i", getHistoPrefix().c_str(), region)));
      hist.cedge = dynamic_cast<TH1*>(hm->getHisto(Form("%sclusedge_%i", getHistoPrefix().c_str(), region)));
      hist.coverlap = dynamic_cast<TH1*>(hm->getHisto(Form("%sclusoverlap_%i", getHistoPrefix().c_str(), region)));
   
      histos.insert(std::make_pair(region, hist));
    }
  auto fill = [](TH1* h, float val) { if (h) h->Fill(val); };

  for(auto& hsk : clusterContainer->getHitSetKeys(TrkrDefs::TrkrId::tpcId))
    {
      int numclusters = 0;
      
      auto range = clusterContainer->getClusters(hsk);
      for(auto iter = range.first; iter != range.second; ++iter)
	{
	  const auto cluskey = iter->first;
	  const auto cluster = iter->second;
	  const auto it = m_layerRegionMap.find(TrkrDefs::getLayer(cluskey));
	  int region = it->second;
	  const auto hiter = histos.find(region);
	  if(hiter == histos.end()) 
	    {
	      continue;
	    }

	  fill(hiter->second.crphisize, cluster->getPhiSize());
	  fill(hiter->second.czsize, cluster->getZSize());
	  fill(hiter->second.crphierr, cluster->getRPhiError());
	  fill(hiter->second.czerr, cluster->getZError());
	  fill(hiter->second.cedge, cluster->getEdge());
	  fill(hiter->second.coverlap, cluster->getOverlap());
	  
	  numclusters++;
	}

      h_totalclusters->Fill(hitsetkeynum, numclusters);
      
      hitsetkeynum++;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________..
int TpcClusterQA::End(PHCompositeNode*)
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

  for(auto& region : {0,1,2})
    {
      {
	auto h = new TH1F(Form("%sphisize_%i", getHistoPrefix().c_str(), region),
			  Form("TPC cluster #phi size region_%i", region), 10,0,10);
	h->GetXaxis()->SetTitle("Cluster #phi_{size}");
	hm->registerHisto(h);
      }
      {
	auto h = new TH1F(Form("%szsize_%i", getHistoPrefix().c_str(), region),
			  Form("TPC cluster z size region_%i", region), 10,0,10);
	h->GetXaxis()->SetTitle("Cluster z_{size}");
	hm->registerHisto(h);
      }
      {
	auto h = new TH1F(Form("%srphi_error_%i", getHistoPrefix().c_str(), region), 
			  Form("TPC r#Delta#phi error region_%i", region), 100, 0, 0.075);
      h->GetXaxis()->SetTitle("r#Delta#phi error [cm]");
      hm->registerHisto(h);
      }
      {
	auto h = new TH1F(Form("%sz_error_%i", getHistoPrefix().c_str(), region), 
			  Form("TPC z error region_%i", region), 100, 0, 0.18);
	h->GetXaxis()->SetTitle("z error [cm]");
	hm->registerHisto(h);
      }
      {
	auto h = new TH1F(Form("%sphisize_%i", getHistoPrefix().c_str(), region),
			  Form("TPC #phi size region_%i", region), 30,0,30);
	h->GetXaxis()->SetTitle("Cluster #phi_{size}");
	hm->registerHisto(h);
      }     
      {
	auto h = new TH1F(Form("%sclusedge_%i", getHistoPrefix().c_str(), region),
			  Form("TPC hits on edge_%i", region), 30,0,30);
	h->GetXaxis()->SetTitle("Cluster edge");
	hm->registerHisto(h);
      }   
	  {
	auto h = new TH1F(Form("%sclusoverlap_%i", getHistoPrefix().c_str(), region),
			  Form("TPC clus overlap_%i", region), 30,0,30);
	h->GetXaxis()->SetTitle("Cluster overlap");
	hm->registerHisto(h);
      }   
    }
  
  {
    auto h = new TH2F(Form("%stotal_clusters", getHistoPrefix().c_str()),
		      "TPC clusters per hitsetkey", 1152,0,1152, 10000,0,10000);
    h->GetXaxis()->SetTitle("Hitsetkey number");
    h->GetYaxis()->SetTitle("Number of clusters");
    hm->registerHisto(h);
  }


  return;
}
