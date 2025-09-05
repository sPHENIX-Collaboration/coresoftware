#include "TrackerEventDisplay.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSetContainer.h>

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
#include <limits>
#include <map>
#include <memory>  // for shared_ptr
#include <set>     // for _Rb_tree_cons...
#include <utility>
#include <vector>

TrackerEventDisplay::TrackerEventDisplay(const std::string& /*name*/, const std::string& filename, const std::string& runnumber, const std::string& date)
  : SubsysReco("TrackerEventDisplay")
  , _filename(filename)
  , _runnumber(runnumber)
  , _date(date)
{
}

int TrackerEventDisplay::Init(PHCompositeNode* /*topNode*/)
{
  _ievent = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}
int TrackerEventDisplay::InitRun(PHCompositeNode *topNode)
{
  auto geom =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  AdcClockPeriod = geom->GetFirstLayerCellGeom()->get_zstep();

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
    std::cout << "========================= TrackerEventDisplay::End() ============================" << std::endl;
    std::cout << " " << _ievent << " events of output written to: " << _filename << std::endl;
    std::cout << "===========================================================================" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void TrackerEventDisplay::makeJsonFile(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "TrackerEventDisplay::makeJsonFile() entered" << std::endl;
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

    outdata.open((_filename + "_event" + std::to_string(_ievent) + "_hits.json").c_str(), std::ofstream::out | std::ofstream::trunc);

    if (!outdata)
    {  // file couldn't be opened
      std::cout << "ERROR: file could not be opened" << std::endl;
      exit(1);
    }

    outdata << "{\n    \"EVENT\": {\n        \"runid\":" << _runnumber << ", \n        \"evtid\": 1, \n        \"time\": 0, \n \"timeStr\": \"2023-08-23, 15:23:30 EST\", \n       \"type\": \"Cosmics\", \n        \"s_nn\": 0, \n        \"B\": 0.0,\n        \"pv\": [0,0,0],\n  \"runstats\": [ \n  \"sPHENIX Time Projection Chamber\", \"" << _date << ", Run " << _runnumber << " - Event " << _ievent << "\", \"All Hits in Event\"] \n   },\n"
            << std::endl;

    outdata << "    \"META\": {\n       \"HITS\": {\n          \"INNERTRACKER\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 2,\n              \"color\": 16777215\n              } \n          },\n"
            << std::endl;
    outdata << "          \"TRACKHITS\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 2,\n              \"transparent\": 0.5,\n              \"color\": 16777215\n              } \n          },\n"
            << std::endl;
    outdata << "    \"JETS\": {\n        \"type\": \"JET\",\n        \"options\": {\n            \"rmin\": 0,\n            \"rmax\": 78,\n            \"emin\": 0,\n            \"emax\": 30,\n            \"color\": 16777215,\n            \"transparent\": 0.5 \n        }\n    }\n        }\n    }\n," << std::endl;
    outdata << "    \"HITS\": {\n        \"CEMC\":[{\"eta\": 0, \"phi\": 0, \"e\": 0}\n            ],\n        \"HCALIN\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n            ],\n        \"HCALOUT\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n \n            ],\n\n"
            << std::endl;
    outdata << "    \"TRACKHITS\": [\n\n ";

    auto m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

    if (Verbosity() >= 1)
    {
      std::cout << "Filling hit json" << std::endl;
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
          // float event = _ievent;
          // float hitID = hit_key;
          // float e = hit->getEnergy();
          float adc = hit->getAdc();
          float layer_local = TrkrDefs::getLayer(hitset_key);
          // float sector = TpcDefs::getSectorId(hitset_key);
          float side = TpcDefs::getSide(hitset_key);
          // float cellID = 0;
          // float ecell = hit->getAdc();

          float phibin = std::numeric_limits<float>::quiet_NaN();
          float tbin = std::numeric_limits<float>::quiet_NaN();
          // float phi = std::numeric_limits<float>::quiet_NaN();
          float phi_center = std::numeric_limits<float>::quiet_NaN();
          float x = std::numeric_limits<float>::quiet_NaN();
          float y = std::numeric_limits<float>::quiet_NaN();
          float z = std::numeric_limits<float>::quiet_NaN();

          if (TrkrDefs::getTrkrId(hitset_key) == TrkrDefs::TrkrId::tpcId)
          {
            PHG4TpcCylinderGeom* GeoLayer_local = geom_container->GetLayerCellGeom(layer_local);
            double radius = GeoLayer_local->get_radius();
            phibin = (float) TpcDefs::getPad(hit_key);
            tbin = (float) TpcDefs::getTBin(hit_key);
            // phi = GeoLayer_local->get_phicenter(phibin);

            double zdriftlength = tbin * m_tGeometry->get_drift_velocity() * AdcClockPeriod;
            // convert z drift length to z position in the TPC
            //		std::cout << " tbin: " << tbin << " vdrift " <<m_tGeometry->get_drift_velocity() << " l drift: " << zdriftlength  <<std::endl;
            unsigned short NTBins = (unsigned short) GeoLayer_local->get_zbins();
            double m_tdriftmax = AdcClockPeriod * NTBins / 2.0;
            double clusz = (m_tdriftmax * m_tGeometry->get_drift_velocity()) - zdriftlength;
            if (side == 0)
            {
              clusz = -clusz;
            }
            z = clusz;
            phi_center = GeoLayer_local->get_phicenter(phibin);
            x = radius * std::cos(phi_center);
            y = radius * std::sin(phi_center);

            std::stringstream spts;

            if (firstHit)
            {
              firstHit = false;
            }
            else
            {
              spts << ",";
            }

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
    outdata << "],\n    \"JETS\": [\n         ]\n    }," << std::endl;
    outdata << "\"TRACKS\": {" << std::endl;
    outdata << "\""
            << "INNERTRACKER"
            << "\": [";
    outdata << "]" << std::endl;
    outdata << "}" << std::endl;
    outdata << "}" << std::endl;
    outdata.close();
  }

  //------------------------
  // fill the Cluster JSON
  //------------------------

  if (_cluster)
  {
    bool firstHit = true;

    outdata.open((_filename + "_event" + std::to_string(_ievent) + "_clusters.json").c_str(), std::ofstream::out | std::ofstream::trunc);

    if (!outdata)
    {  // file couldn't be opened
      std::cout << "ERROR: file could not be opened" << std::endl;
      exit(1);
    }

    outdata << "{\n    \"EVENT\": {\n        \"runid\":" << _runnumber << ", \n        \"evtid\": 1, \n        \"time\": 0, \n \"timeStr\": \"2023-08-23, 15:23:30 EST\", \n       \"type\": \"Cosmics\", \n        \"s_nn\": 0, \n        \"B\": 0.0,\n        \"pv\": [0,0,0],\n  \"runstats\": [ \n  \"sPHENIX Time Projection Chamber\", \"" << _date << ", Run " << _runnumber << " - Event " << _ievent << "\", \"All Clusters in Event\"] \n   },\n"
            << std::endl;

    outdata << "    \"META\": {\n       \"HITS\": {\n          \"INNERTRACKER\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 2,\n              \"color\": 16777215\n              } \n          },\n"
            << std::endl;
    outdata << "          \"TRACKHITS\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 2,\n              \"transparent\": 0.5,\n              \"color\": 16777215\n              } \n          },\n"
            << std::endl;
    outdata << "    \"JETS\": {\n        \"type\": \"JET\",\n        \"options\": {\n            \"rmin\": 0,\n            \"rmax\": 78,\n            \"emin\": 0,\n            \"emax\": 30,\n            \"color\": 16777215,\n            \"transparent\": 0.5 \n        }\n    }\n        }\n    }\n," << std::endl;
    outdata << "    \"HITS\": {\n        \"CEMC\":[{\"eta\": 0, \"phi\": 0, \"e\": 0}\n            ],\n        \"HCALIN\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n            ],\n        \"HCALOUT\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n \n            ],\n\n"
            << std::endl;
    outdata << "    \"TRACKHITS\": [\n\n ";

    if (Verbosity() > 1)
    {
      std::cout << "Filling cluster json " << std::endl;
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
        std::cout << "got clustermap" << std::endl;
      }
      else
      {
        std::cout << "no clustermap" << std::endl;
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

          std::stringstream spts;

          if (firstHit)
          {
            firstHit = false;
          }
          else
          {
            spts << ",";
          }

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
    outdata << "],\n    \"JETS\": [\n         ]\n    }," << std::endl;
    outdata << "\"TRACKS\": {" << std::endl;
    outdata << "\""
            << "INNERTRACKER"
            << "\": [";
    outdata << "]" << std::endl;
    outdata << "}" << std::endl;
    outdata << "}" << std::endl;
    outdata.close();
  }
  return;
}
