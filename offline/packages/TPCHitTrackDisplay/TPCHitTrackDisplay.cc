#include "TPCHitTrackDisplay.h"

#include <centrality/CentralityInfo.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>

#include <phgeom/PHGeomUtility.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TVector3.h>

#include <format>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

/*************************************************************/
/*              TPC Event Display Generator                  */
/*         Thomas Marshall,Aditya Dash,Ejiro Umaka           */
/*rosstom@ucla.edu,aditya55@physics.ucla.edu,eumaka1@bnl.gov */
/*************************************************************/

//----------------------------------------------------------------------------//
//-- Constructor:
//--  simple initialization
//----------------------------------------------------------------------------//

TPCHitTrackDisplay::TPCHitTrackDisplay(const std::string& name)
  : SubsysReco(name)
  , m_cut_ADC(75.0)
  , m_trackless_clusters(true)
  , _fileName(name)
{
}

//----------------------------------------------------------------------------//
//-- process_event():
//--   Call user instructions for every event.
//--   This function contains the analysis structure.
//----------------------------------------------------------------------------//

int TPCHitTrackDisplay::process_event(PHCompositeNode* topNode)
{
  _event++;
  SimulationOut(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

// Produce json file from raw TPC csv data
//  Will adjust this soon, once TrkrHitSetContainer is updated
/*
void TPCHitTrackDisplay::TPCRawOut(PHCompositeNode *topNode)
{

    std::cout << _event << std::endl;

    outfile1 << "\"TRACKS\": {" << std::endl;
    outfile1 <<"\""<<"INNERTRACKER"<<"\": [";

    bool first = true;

    std::stringstream spts;
    float c = std::numeric_limits<float>::quiet_NaN();

    float x = 0.; float y = 0.; float z = 0.;
    float px = 0.; float py = 0.; float pz = 0.;
    TVector3 pos; TVector3 mom;

    for

    for (SvtxTrackMap::Iter iter = trackmap->begin(); iter != trackmap->end(); ++iter)
    {
      SvtxTrack* track = iter->second;
      px = track->get_px();
      py = track->get_py();
      pz = track->get_pz();
      mom.SetX(px); mom.SetY(py); mom.SetZ(pz);
      c = track->get_charge();


      std::vector<TrkrDefs::cluskey> clusters;
      auto siseed = track->get_silicon_seed();
      if(siseed)
      {
        for (auto iter = siseed->begin_cluster_keys(); iter != siseed->end_cluster_keys(); ++iter)
            {
                  TrkrDefs::cluskey cluster_key = *iter;
                  clusters.push_back(cluster_key);
            }
        }

      auto tpcseed = track->get_tpc_seed();
      if(tpcseed)
      {
        for (auto iter = tpcseed->begin_cluster_keys(); iter != tpcseed->end_cluster_keys(); ++iter)
            {
                  TrkrDefs::cluskey cluster_key = *iter;
                  clusters.push_back(cluster_key);
            }
      }


     for(unsigned int iclus = 0; iclus < clusters.size(); ++iclus)
      {
        TrkrDefs::cluskey cluster_key = clusters[iclus];
        TrkrCluster* cluster = clusterContainer->findCluster(cluster_key);
        if(!cluster) continue;

        Acts::Vector3 globalClusterPosition = geometry->getGlobalPosition(cluster_key, cluster);
        x = globalClusterPosition(0);
        y = globalClusterPosition(1);
        z = globalClusterPosition(2);
        pos.SetX(x); pos.SetY(y); pos.SetZ(z);

        if (first)
        {
            first = false;
        }

        else
        spts << ",";
        spts << "[";
        spts << pos.x();
        spts << ",";
        spts << pos.y();
        spts << ",";
        spts << pos.z();
        spts << "]";

        first = true;
    }

     outfile1
        << std::format("{ \"pt\": {}, \"e\": {}, \"p\": {}, \"c\": {}, \"pdgcode \": {}, \"pts\":[ {} ]},",
          mom.Pt(), mom.PseudoRapidity(), mom.Phi(), c, _pdgcode, spts.str() ) << std::endl;
           spts.str("");
}


    outfile1 << "]" << std::endl;
    outfile1 << "}" << std::endl;



return;

}
*/

// write json file from simulation data for silicon and tpc clusters/tracks
void TPCHitTrackDisplay::SimulationOut(PHCompositeNode* topNode)
{
  std::string fname = _fileName + "event" + std::to_string(_event) + ".json";
  std::fstream outfile1(fname, std::ios_base::out);

  std::cout << _event << std::endl;

  TrkrClusterContainer* clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  CentralityInfo* cent = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  if (!cent)
  {
    std::cout << " ERROR -- can't find CentralityInfo node" << std::endl;
    return;
  }

  int cent_index = cent->get_centile(CentralityInfo::PROP::bimp);

  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  ActsGeometry* geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!geometry)
  {
    std::cout << PHWHERE << "No Acts geometry on node tree. Can't  continue."
              << std::endl;
  }

  std::cout << "This event has centrality: \t" << cent_index << std::endl;

  // Header information for the json file

  outfile1 << "{\n    \"EVENT\": {\n        \"runid\": 1, \n        \"evtid\": 1, \n        \"time\": 0, \n        \"type\": \"Collision\", \n        \"s_nn\": 0, \n        \"B\": 3.0,\n        \"pv\": [0,0,0]  \n    },\n"
           << std::endl;

  outfile1 << "    \"META\": {\n       \"HITS\": {\n          \"INNERTRACKER\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 5,\n              \"color\": 16777215\n              } \n          },\n"
           << std::endl;
  outfile1 << "          \"TRACKHITS\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 5,\n              \"transparent\": 0.5,\n              \"color\": 16777215\n              } \n          },\n"
           << std::endl;
  outfile1 << "    \"JETS\": {\n        \"type\": \"JET\",\n        \"options\": {\n            \"rmin\": 0,\n            \"rmax\": 78,\n            \"emin\": 0,\n            \"emax\": 30,\n            \"color\": 16777215,\n            \"transparent\": 0.5 \n        }\n    }\n        }\n    }\n," << std::endl;
  outfile1 << "    \"HITS\": {\n        \"CEMC\":[{\"eta\": 0, \"phi\": 0, \"e\": 0}\n            ],\n        \"HCALIN\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n            ],\n        \"HCALOUT\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n \n            ],\n\n"
           << std::endl;
  outfile1 << "    \"TRACKHITS\": [\n\n ";

  std::vector<TrkrDefs::cluskey> usedClusters;

  // Fill usedClusters with cluster keys that are associated with a reconstructed track
  for (auto& iter : *trackmap)
  {
    if (!m_trackless_clusters)
    {
      break;
    }
    SvtxTrack* track = iter.second;
    std::vector<TrkrDefs::cluskey> clusters;
    auto *siseed = track->get_silicon_seed();
    if (siseed)
    {
      for (auto iter1 = siseed->begin_cluster_keys(); iter1 != siseed->end_cluster_keys(); ++iter1)
      {
        TrkrDefs::cluskey cluster_key = *iter1;
        usedClusters.push_back(cluster_key);
      }
    }

    auto *tpcseed = track->get_tpc_seed();
    if (tpcseed)
    {
      for (auto iter1 = tpcseed->begin_cluster_keys(); iter1 != tpcseed->end_cluster_keys(); ++iter1)
      {
        TrkrDefs::cluskey cluster_key = *iter1;
        usedClusters.push_back(cluster_key);
      }
    }
  }

  bool firstClus = true;
  int counter = 0;
  std::stringstream spts;
  float c = std::numeric_limits<float>::quiet_NaN();

  float x = 0.;
  float y = 0.;
  float z = 0.;
  float px = 0.;
  float py = 0.;
  float pz = 0.;
  TVector3 pos;
  TVector3 mom;

  // iterate over all clusters and write out only those without an associated track
  for (const auto& hitsetkey : clusterContainer->getHitSetKeys())
  {
    if (!m_trackless_clusters)
    {
      break;  // cut on whether or not to display trackless clusters
    }
    auto range = clusterContainer->getClusters(hitsetkey);
    for (auto iter = range.first; iter != range.second; ++iter)
    {
      counter = counter + 1;
      TrkrDefs::cluskey cluster_key = iter->first;
      bool clusFromTrack = false;
      for (unsigned long trackClusKey : usedClusters)
      {
        if (trackClusKey == cluster_key)
        {
          clusFromTrack = true;
          break;
        }
      }
      if (clusFromTrack)
      {
        continue;
      }
      TrkrCluster* cluster = iter->second;
      if (!cluster)
      {
        continue;
      }

      unsigned int clusADC = cluster->getAdc();
      if (clusADC < m_cut_ADC)
      {
        continue;  // ADC check on cluster to not over-saturate display
      }

      Acts::Vector3 globalClusterPosition = geometry->getGlobalPosition(cluster_key, cluster);
      x = globalClusterPosition(0);
      y = globalClusterPosition(1);
      z = globalClusterPosition(2);
      pos.SetX(x);
      pos.SetY(y);
      pos.SetZ(z);

      if (firstClus)
      {
        firstClus = false;
      }
      else
      {
        spts << ",";
      }

      spts << "{ \"x\": ";
      spts << pos.x();
      spts << ", \"y\": ";
      spts << pos.y();
      spts << ", \"z\": ";
      spts << pos.z();
      spts << ", \"e\": ";
      spts << clusADC;
      spts << "}";

      outfile1 << std::format("{}", spts.str());
      spts.str("");
    }
  }

  outfile1 << "],\n    \"JETS\": [\n         ]\n    }," << std::endl;

  outfile1 << "\"TRACKS\": {" << std::endl;
  outfile1 << "\""
           << "INNERTRACKER"
           << "\": [";

  firstClus = true;
  bool firstTrack = true;

  // tracking

  c = std::numeric_limits<float>::quiet_NaN();

  x = 0.;
  y = 0.;
  z = 0.;
  px = 0.;
  py = 0.;
  pz = 0.;

  // iterate over tracks and write to file all tracks and their associated clusters
  for (auto& iter : *trackmap)
  {
    SvtxTrack* track = iter.second;
    px = track->get_px();
    py = track->get_py();
    pz = track->get_pz();
    mom.SetX(px);
    mom.SetY(py);
    mom.SetZ(pz);
    c = track->get_charge();

    std::vector<TrkrDefs::cluskey> clusters;
    auto *siseed = track->get_silicon_seed();
    if (siseed)
    {
      for (auto iter1 = siseed->begin_cluster_keys(); iter1 != siseed->end_cluster_keys(); ++iter1)
      {
        TrkrDefs::cluskey cluster_key = *iter1;
        clusters.push_back(cluster_key);
      }
    }

    auto *tpcseed = track->get_tpc_seed();
    if (tpcseed)
    {
      for (auto iter1 = tpcseed->begin_cluster_keys(); iter1 != tpcseed->end_cluster_keys(); ++iter1)
      {
        TrkrDefs::cluskey cluster_key = *iter1;
        clusters.push_back(cluster_key);
      }
    }

    for (unsigned long cluster_key : clusters)
    {
      TrkrCluster* cluster = clusterContainer->findCluster(cluster_key);
      if (!cluster)
      {
        continue;
      }

      Acts::Vector3 globalClusterPosition = geometry->getGlobalPosition(cluster_key, cluster);
      x = globalClusterPosition(0);
      y = globalClusterPosition(1);
      z = globalClusterPosition(2);
      pos.SetX(x);
      pos.SetY(y);
      pos.SetZ(z);

      if (firstClus)
      {
        firstClus = false;
      }

      else
      {
        spts << ",";
      }
      spts << "[";
      spts << pos.x();
      spts << ",";
      spts << pos.y();
      spts << ",";
      spts << pos.z();
      spts << "]";
    }

    firstClus = true;
    if (firstTrack)
    {
      outfile1 << std::format("{{ \"pt\": {}, \"e\": {}, \"p\": {}, \"c\": {}, \"pdgcode \": {}, \"centrality \": {}, \"pts\": [ {} ] }}",
		  mom.Pt(), mom.PseudoRapidity(), mom.Phi(), c, _pdgcode, cent_index, spts.str())
          << std::endl;
      spts.str("");
      firstTrack = false;
    }
    else
    {
      outfile1 << std::format(",{{ \"pt\": {}, \"e\": {}, \"p\": {}, \"c\": {}, \"pdgcode\": {}, \"centrality\": {}, \"pts\": [ {} ] }}",
			      mom.Pt(), mom.PseudoRapidity(), mom.Phi(), c, _pdgcode, cent_index, spts.str())
          << std::endl;
      spts.str("");
    }
  }

  outfile1 << "]" << std::endl;
  outfile1 << "}" << std::endl;
  outfile1 << "}" << std::endl;

  usedClusters.clear();

  return;
}
