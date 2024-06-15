#include "TrackerEventDisplay.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrackFitUtils.h>

#include <trackbase_historic/TrackSeedContainer_v1.h>

#include <trackbase/TpcDefs.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/TrkrHitSet.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TVector3.h>
#include <boost/format.hpp>
#include <boost/math/special_functions/sign.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>  // for shared_ptr
#include <set>     // for _Rb_tree_cons...
#include <utility>
#include <vector>

using namespace std;

TrackerEventDisplay::TrackerEventDisplay(const string& /*name*/, const string& filename, const string& runnumber, const string& date)
  : SubsysReco("TrackerEventDisplay")
  , _hit(true)
  , _cluster(false)
  , _ievent(0)
  , _filename(filename)
  , _runnumber(runnumber)
  , _date(date)
{
}

TrackerEventDisplay::~TrackerEventDisplay()
{
}

int TrackerEventDisplay::Init(PHCompositeNode* /*topNode*/)
{
  _ievent = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

int TrackerEventDisplay::InitRun(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrackerEventDisplay::process_event(PHCompositeNode* topNode)
{
  makeJsonFile(topNode);
  ++_ievent;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrackerEventDisplay::End(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1)
  {
    cout << "========================= TrackerEventDisplay::End() ============================" << endl;
    cout << " " << _ievent << " events of output written to: " << _filename << endl;
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void TrackerEventDisplay::makeJsonFile(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    cout << "TrackerEventDisplay::makeJsonFile() entered" << endl;
  }

  ActsGeometry* tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tgeometry)
  {
    std::cout << "No Acts geometry on node tree. Can't  continue."
              << std::endl;
    return;
  }

  PHG4TpcCylinderGeomContainer* geom_container =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return;
  }

  //--------------------
  // fill the Hit json
  //--------------------

  if (_hit)
  {
    bool firstHit = true;   

    outdata.open((_filename+"_event"+std::to_string(_ievent)+"_hits.json").c_str(),std::ofstream::out | std::ofstream::trunc);
  
    if( !outdata ) 
    { // file couldn't be opened
      cerr << "ERROR: file could not be opened" << endl;
      exit(1);
    }
  
    outdata << "{\n    \"EVENT\": {\n        \"runid\":" << _runnumber << ", \n        \"evtid\": 1, \n        \"time\": 0, \n \"timeStr\": \"2023-08-23, 15:23:30 EST\", \n       \"type\": \"Cosmics\", \n        \"s_nn\": 0, \n        \"B\": 0.0,\n        \"pv\": [0,0,0],\n  \"runstats\": [ \n  \"sPHENIX Time Projection Chamber\", \"" << _date << ", Run " << _runnumber << " - Event " << _ievent << "\", \"All Hits in Event\"] \n   },\n" << endl;

    outdata << "    \"META\": {\n       \"HITS\": {\n          \"INNERTRACKER\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 2,\n              \"color\": 16777215\n              } \n          },\n" << endl;
    outdata << "          \"TRACKHITS\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 2,\n              \"transparent\": 0.5,\n              \"color\": 16777215\n              } \n          },\n" << endl;
    outdata << "    \"JETS\": {\n        \"type\": \"JET\",\n        \"options\": {\n            \"rmin\": 0,\n            \"rmax\": 78,\n            \"emin\": 0,\n            \"emax\": 30,\n            \"color\": 16777215,\n            \"transparent\": 0.5 \n        }\n    }\n        }\n    }\n," << endl;
    outdata << "    \"HITS\": {\n        \"CEMC\":[{\"eta\": 0, \"phi\": 0, \"e\": 0}\n            ],\n        \"HCALIN\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n            ],\n        \"HCALOUT\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n \n            ],\n\n" << endl;
    outdata << "    \"TRACKHITS\": [\n\n ";
    
    auto m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

    if (Verbosity() >= 1)
    {
      cout << "Filling hit json" << endl;
    }
    // need things off of the DST...
    TrkrHitSetContainer* hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
 
    if (hitmap)
    {
      TrkrHitSetContainer::ConstRange all_hitsets = hitmap->getHitSets();
      for (TrkrHitSetContainer::ConstIterator iter = all_hitsets.first;
           iter != all_hitsets.second;
           ++iter)
      {
        const TrkrDefs::hitsetkey hitset_key = iter->first;
        TrkrHitSet* hitset = iter->second;
	      // get all hits for this hitset
        TrkrHitSet::ConstRange hitrangei = hitset->getHits();
        for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
             hitr != hitrangei.second;
             ++hitr)
        {
          TrkrDefs::hitkey hit_key = hitr->first;
          TrkrHit* hit = hitr->second;
          //float event = _ievent;
          //float hitID = hit_key;
          //float e = hit->getEnergy();
          float adc = hit->getAdc();
          float layer_local = TrkrDefs::getLayer(hitset_key);
          //float sector = TpcDefs::getSectorId(hitset_key);
          float side = TpcDefs::getSide(hitset_key);
          //float cellID = 0;
          //float ecell = hit->getAdc();

          float phibin = NAN;
          float tbin = NAN;
          //float phi = NAN;
	  float phi_center = NAN;
          float x = NAN;
          float y = NAN;
          float z = NAN;
        
          if (TrkrDefs::getTrkrId(hitset_key) == TrkrDefs::TrkrId::tpcId)
          {
            PHG4TpcCylinderGeom* GeoLayer_local = geom_container->GetLayerCellGeom(layer_local);
 	    double radius = GeoLayer_local->get_radius();
            phibin = (float) TpcDefs::getPad(hit_key);
            tbin = (float) TpcDefs::getTBin(hit_key);
            //phi = GeoLayer_local->get_phicenter(phibin);

            double zdriftlength = tbin * m_tGeometry->get_drift_velocity() * AdcClockPeriod;
            // convert z drift length to z position in the TPC
            //		cout << " tbin: " << tbin << " vdrift " <<m_tGeometry->get_drift_velocity() << " l drift: " << zdriftlength  <<endl;
            unsigned short NTBins = (unsigned short) GeoLayer_local->get_zbins();
            double m_tdriftmax = AdcClockPeriod * NTBins / 2.0;
            double clusz = (m_tdriftmax * m_tGeometry->get_drift_velocity()) - zdriftlength;
            if (side == 0)
            {
              clusz = -clusz;
            }
            z = clusz;
            phi_center = GeoLayer_local->get_phicenter(phibin);
	    x = radius * cos(phi_center);
	    y = radius * sin(phi_center);
   
            stringstream spts;

            if (firstHit) firstHit = false;
            else spts << ",";

            spts << "{ \"x\": ";
            spts << x;
            spts << ", \"y\": ";
            spts << y;
            spts << ", \"z\": ";
            spts << z;
            spts << ", \"e\": ";
            spts << adc;
            spts << "}";
     
            outdata << (boost::format("%1%") % spts.str());
            spts.clear();
            spts.str("");
          }
        }
      }
    }
    outdata << "],\n    \"JETS\": [\n         ]\n    }," << endl;
    outdata << "\"TRACKS\": {" << endl;
    outdata <<"\""<<"INNERTRACKER"<<"\": [";
    outdata << "]" << endl;
    outdata << "}" << endl;
    outdata << "}" << endl;
    outdata.close();
  }

  //------------------------
  // fill the Cluster JSON
  //------------------------

  if (_cluster)
  {
    bool firstHit = true;   
    
    outdata.open((_filename+"_event"+std::to_string(_ievent)+"_clusters.json").c_str(),std::ofstream::out | std::ofstream::trunc);
  
    if( !outdata ) 
    { // file couldn't be opened
      cerr << "ERROR: file could not be opened" << endl;
      exit(1);
    }
  
    outdata << "{\n    \"EVENT\": {\n        \"runid\":" << _runnumber << ", \n        \"evtid\": 1, \n        \"time\": 0, \n \"timeStr\": \"2023-08-23, 15:23:30 EST\", \n       \"type\": \"Cosmics\", \n        \"s_nn\": 0, \n        \"B\": 0.0,\n        \"pv\": [0,0,0],\n  \"runstats\": [ \n  \"sPHENIX Time Projection Chamber\", \"" << _date << ", Run " << _runnumber << " - Event " << _ievent << "\", \"All Clusters in Event\"] \n   },\n" << endl;

    outdata << "    \"META\": {\n       \"HITS\": {\n          \"INNERTRACKER\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 2,\n              \"color\": 16777215\n              } \n          },\n" << endl;
    outdata << "          \"TRACKHITS\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 2,\n              \"transparent\": 0.5,\n              \"color\": 16777215\n              } \n          },\n" << endl;
    outdata << "    \"JETS\": {\n        \"type\": \"JET\",\n        \"options\": {\n            \"rmin\": 0,\n            \"rmax\": 78,\n            \"emin\": 0,\n            \"emax\": 30,\n            \"color\": 16777215,\n            \"transparent\": 0.5 \n        }\n    }\n        }\n    }\n," << endl;
    outdata << "    \"HITS\": {\n        \"CEMC\":[{\"eta\": 0, \"phi\": 0, \"e\": 0}\n            ],\n        \"HCALIN\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n            ],\n        \"HCALOUT\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n \n            ],\n\n" << endl;
    outdata << "    \"TRACKHITS\": [\n\n ";
    
    if (Verbosity() > 1)
    {
      cout << "Filling cluster json " << endl;
    }
    
    // need things off of the DST...
    TrkrClusterContainer* clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
    if (!clustermap)
    {
      clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    }

    if (Verbosity() > 1)
    {
      if (clustermap != nullptr)
      {
        cout << "got clustermap" << endl;
      }
      else
      {
        cout << "no clustermap" << endl;
      }
    }

    if (clustermap)
    {
      for (const auto& hitsetkey : clustermap->getHitSetKeys())
      {
        auto range = clustermap->getClusters(hitsetkey);
        for (auto iter = range.first; iter != range.second; ++iter)
        {
          TrkrDefs::cluskey cluster_key = iter->first;
          TrkrCluster* cluster = clustermap->findCluster(cluster_key);
         
          Acts::Vector3 cglob;
          cglob = tgeometry->getGlobalPosition(cluster_key, cluster);
          float x = cglob(0);
          float y = cglob(1);
          float z = cglob(2);
          float adc = cluster->getAdc();

          stringstream spts;  
  
          if (firstHit) firstHit = false;
          else spts << ",";

          spts << "{ \"x\": ";
          spts << x;
          spts << ", \"y\": ";
          spts << y;
          spts << ", \"z\": ";
          spts << z;
          spts << ", \"e\": ";
          spts << adc;
          spts << "}";

          outdata << (boost::format("%1%") % spts.str());
          spts.clear();
          spts.str(""); 
        }
      }
    }
    outdata << "],\n    \"JETS\": [\n         ]\n    }," << endl;
    outdata << "\"TRACKS\": {" << endl;
    outdata <<"\""<<"INNERTRACKER"<<"\": [";
    outdata << "]" << endl;
    outdata << "}" << endl;
    outdata << "}" << endl;
    outdata.close();
  }
  return;
}
