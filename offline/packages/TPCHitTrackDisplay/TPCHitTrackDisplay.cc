#include "TPCHitTrackDisplay.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllServer.h>

#include <centrality/CentralityInfo.h>
#include <centrality/CentralityInfov1.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Particle.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase/TrkrCluster.h>

#include <phgeom/PHGeomUtility.h>

#include <TTree.h>
#include <TH2D.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TMath.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>
#include <jetbase/JetContainerv1.h>
#include <g4main/PHG4Utils.h>

#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>
#include <boost/format.hpp>
#include <boost/math/special_functions/sign.hpp>

/*************************************************************/
/*              TPC Event Display Generator                  */
/*         Thomas Marshall,Aditya Dash,Ejiro Umaka           */
/*rosstom@ucla.edu,aditya55@physics.ucla.edu,eumaka1@bnl.gov */
/*************************************************************/


using namespace std;

//----------------------------------------------------------------------------//
//-- Constructor:
//--  simple initialization
//----------------------------------------------------------------------------//

TPCHitTrackDisplay::TPCHitTrackDisplay(const string &name) :
  SubsysReco(name)
  , m_cut_ADC(75.0)
  , m_trackless_clusters(true)
  , _pdgcode(0)
  , _fileName(name)
{
	//initialize
	_event = 0;
}

//----------------------------------------------------------------------------//
//-- Init():
//--   Intialize all histograms, trees, and ntuples
//----------------------------------------------------------------------------//
int TPCHitTrackDisplay::Init(PHCompositeNode *topNode) {
    
	return Fun4AllReturnCodes::EVENT_OK;
}

int TPCHitTrackDisplay::InitRun(PHCompositeNode *topNode)
{
	return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//
//-- process_event():
//--   Call user instructions for every event.
//--   This function contains the analysis structure.
//----------------------------------------------------------------------------//

int TPCHitTrackDisplay::process_event(PHCompositeNode *topNode)
{
	_event++;
	SimulationOut(topNode);
	
       return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//
//-- End():
//--   End method, wrap everything up
//----------------------------------------------------------------------------//

int TPCHitTrackDisplay::EndRun(const int runnumber)
{

    return Fun4AllReturnCodes::EVENT_OK;

}

//Produce json file from raw TPC csv data
// Will adjust this soon, once TrkrHitSetContainer is updated
/*
void TPCHitTrackDisplay::TPCRawOut(PHCompositeNode *topNode)
{

    cout << _event << endl;

    outfile1 << "\"TRACKS\": {" << endl;
    outfile1 <<"\""<<"INNERTRACKER"<<"\": [";
    
    bool first = true;
 
    stringstream spts;
    float c = NAN;
    
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
        << (boost::format(
          "{ \"pt\": %1%, \"e\": %2%, \"p\": %3%, \"c\": %4%, \"pdgcode \": %5%, \"pts\":[ %6% ]},")
           % mom.Pt() % mom.PseudoRapidity() % mom.Phi() % c % _pdgcode % spts.str() ) << endl;
           spts.clear();
           spts.str("");
}
  

    outfile1 << "]" << endl;
    outfile1 << "}" << endl;
 
    
    
return;

}
*/

// write json file from simulation data for silicon and tpc clusters/tracks
void TPCHitTrackDisplay::SimulationOut(PHCompositeNode *topNode)
{
    stringstream fname;
    fname << _fileName << "event" << _event << ".json";
    fstream outfile1(fname.str(), ios_base::out); 
    
    cout << _event << endl;

    TrkrClusterContainer *clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

    CentralityInfov1 *cent = findNode::getClass<CentralityInfov1>(topNode, "CentralityInfo");
    if (!cent)
    {
      std::cout << " ERROR -- can't find CentralityInfo node" << std::endl;
      return;
    }

    int cent_index = cent->get_centile(CentralityInfo::PROP::bimp);


    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

    ActsGeometry *geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    if(!geometry)
    {
      std::cout << PHWHERE << "No Acts geometry on node tree. Can't  continue."
              << std::endl;
    }

    std::cout<< "This event has centrality: \t" << cent_index << std::endl;
   

    // Header information for the json file
 
    outfile1 << "{\n    \"EVENT\": {\n        \"runid\": 1, \n        \"evtid\": 1, \n        \"time\": 0, \n        \"type\": \"Collision\", \n        \"s_nn\": 0, \n        \"B\": 3.0,\n        \"pv\": [0,0,0]  \n    },\n" << endl;

    outfile1 << "    \"META\": {\n       \"HITS\": {\n          \"INNERTRACKER\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 5,\n              \"color\": 16777215\n              } \n          },\n" << endl;
    outfile1 << "          \"TRACKHITS\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 5,\n              \"transparent\": 0.5,\n              \"color\": 16777215\n              } \n          },\n" << endl;
    outfile1 << "    \"JETS\": {\n        \"type\": \"JET\",\n        \"options\": {\n            \"rmin\": 0,\n            \"rmax\": 78,\n            \"emin\": 0,\n            \"emax\": 30,\n            \"color\": 16777215,\n            \"transparent\": 0.5 \n        }\n    }\n        }\n    }\n," << endl;
    outfile1 << "    \"HITS\": {\n        \"CEMC\":[{\"eta\": 0, \"phi\": 0, \"e\": 0}\n            ],\n        \"HCALIN\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n            ],\n        \"HCALOUT\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n \n            ],\n\n" << endl;
    outfile1 << "    \"TRACKHITS\": [\n\n "; 

    std::vector<TrkrDefs::cluskey> usedClusters; 

    // Fill usedClusters with cluster keys that are associated with a reconstructed track
    for (SvtxTrackMap::Iter iter = trackmap->begin(); iter != trackmap->end(); ++iter)
    {
      if (!m_trackless_clusters) break;
      SvtxTrack* track = iter->second;
      std::vector<TrkrDefs::cluskey> clusters;
      auto siseed = track->get_silicon_seed();
      if(siseed)
      {
        for (auto iter = siseed->begin_cluster_keys(); iter != siseed->end_cluster_keys(); ++iter)
            {
                  TrkrDefs::cluskey cluster_key = *iter;
 	          usedClusters.push_back(cluster_key);
            }
      }
        
      auto tpcseed = track->get_tpc_seed();
      if(tpcseed)
      {
        for (auto iter = tpcseed->begin_cluster_keys(); iter != tpcseed->end_cluster_keys(); ++iter)
            {
                  TrkrDefs::cluskey cluster_key = *iter;
	          usedClusters.push_back(cluster_key);
            }
      }
    }   

    bool firstClus = true;
    int counter = 0;
    stringstream spts;
    float c = NAN;
    
    float x = 0.; float y = 0.; float z = 0.;
    float px = 0.; float py = 0.; float pz = 0.;
    TVector3 pos; TVector3 mom;
   
    // iterate over all clusters and write out only those without an associated track 
    for(const auto& hitsetkey:clusterContainer->getHitSetKeys())
    {
      if (!m_trackless_clusters) break; // cut on whether or not to display trackless clusters
      auto range = clusterContainer->getClusters(hitsetkey);
      for( auto iter = range.first; iter != range.second; ++iter )
      {
        counter = counter + 1;
        TrkrDefs::cluskey cluster_key = iter->first;
        bool clusFromTrack = false;
        for(unsigned int iclus = 0; iclus < usedClusters.size(); ++iclus)
        {
          TrkrDefs::cluskey trackClusKey = usedClusters[iclus];
          if (trackClusKey == cluster_key) 
          {
            clusFromTrack = true;
            break;
          }
        } 
        if (clusFromTrack) continue;
        TrkrCluster* cluster = iter->second;
        if(!cluster) continue;

        unsigned int clusADC = cluster->getAdc();
        if (clusADC < m_cut_ADC) continue;      // ADC check on cluster to not over-saturate display 
 
        Acts::Vector3 globalClusterPosition = geometry->getGlobalPosition(cluster_key, cluster);
        x = globalClusterPosition(0);
        y = globalClusterPosition(1);
        z = globalClusterPosition(2);
        pos.SetX(x); pos.SetY(y); pos.SetZ(z);
               
        if (firstClus) firstClus = false;
        else spts << ","; 

        spts << "{ \"x\": ";
        spts << pos.x();
        spts << ", \"y\": ";
        spts << pos.y();
        spts << ", \"z\": ";
        spts << pos.z();
        spts << ", \"e\": ";
        spts << clusADC;
        spts << "}";
    
        outfile1 << (boost::format("%1%") % spts.str());
        spts.clear();
        spts.str("");
      }
    }

    outfile1 << "],\n    \"JETS\": [\n         ]\n    }," << endl;

    outfile1 << "\"TRACKS\": {" << endl;
    outfile1 <<"\""<<"INNERTRACKER"<<"\": [";   
 
    firstClus = true;
    bool firstTrack = true;
   
    //tracking

    c = NAN;
    
    x = 0.; y = 0.; z = 0.;
    px = 0.; py = 0.; pz = 0.;

    // iterate over tracks and write to file all tracks and their associated clusters
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
                
        if (firstClus)
        {
            firstClus = false;
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

     }
   
     firstClus = true;
     if (firstTrack)
     {
       outfile1
        << (boost::format(
          "{ \"pt\": %1%, \"e\": %2%, \"p\": %3%, \"c\": %4%, \"pdgcode \": %5%, \"centrality \": %6%, \"pts\":[ %7% ]}")
           % mom.Pt() % mom.PseudoRapidity() % mom.Phi() % c % _pdgcode % cent_index % spts.str() ) << endl;
           spts.clear();
           spts.str("");
           firstTrack = false;
     }
     else
     {
       outfile1
        << (boost::format(
          ",{ \"pt\": %1%, \"e\": %2%, \"p\": %3%, \"c\": %4%, \"pdgcode \": %5%, \"centrality \": %6%, \"pts\":[ %7% ]}")
           % mom.Pt() % mom.PseudoRapidity() % mom.Phi() % c % _pdgcode % cent_index % spts.str() ) << endl;
           spts.clear();
           spts.str("");
     }
    }
    
    
    outfile1 << "]" << endl;
    outfile1 << "}" << endl; 
    outfile1 << "}" << endl; 
    
    usedClusters.clear();
  
    return;

}

