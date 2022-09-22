/*!
 *  \file PHHybridSeeding.cc
 *  \brief Track Seeding using STAR "CA" algorithm and ALICE simplified Kalman filter
 *  \detail 
 *  \author Michael Peters & Christof Roland
 */

#include "PHHybridSeeding.h"
#include "GPUTPCTrackLinearisation.h"
#include "GPUTPCTrackParam.h"

#include "../PHTpcTracker/externals/kdfinder.hpp"
#include "../PHTpcTracker/PHTpcTrackerUtil.h"

// trackbase_historic includes
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_v3.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>  // for getLayer, clu...

// sPHENIX Geant4 includes
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

// sPHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHTimer.h>  // for PHTimer
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

//ROOT includes
#include <TVector3.h>  // for TVector3
#include <TFile.h>
#include <TNtuple.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <utility>  // for pair, make_pair
#include <vector>
#include <algorithm> // for find
#include <unordered_set>

// forward declarations
class PHCompositeNode;

#define _DEBUG_

#if defined(_DEBUG_)
#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#else
#define LogDebug(exp) (void)0
#endif

#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

//end

typedef std::vector<TrkrDefs::cluskey> keylist;
using namespace std;

PHHybridSeeding::PHHybridSeeding(
    const string &name,
    double maxSinPhi,
    double fieldDir,
    double search_radius1,
    double search_angle1,
    size_t min_track_size1,
    double search_radius2,
    double search_angle2,
    size_t min_track_size2,
    size_t nthreads
    )
  : PHTrackSeeding(name)
  , _max_sin_phi(maxSinPhi)
  , _fieldDir(fieldDir)
  , _search_radius1(search_radius1)
  , _search_angle1(search_angle1)
  , _min_track_size1(min_track_size1)
  , _search_radius2(search_radius2)
  , _search_angle2(search_angle2)
  , _min_track_size2(min_track_size2)
  , _nthreads(nthreads)
{
}

int PHHybridSeeding::InitializeGeometry(PHCompositeNode *topNode)
{
  PHG4TpcCylinderGeomContainer *cellgeos = findNode::getClass<
      PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  PHG4CylinderGeomContainer *laddergeos = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  PHG4CylinderGeomContainer *mapsladdergeos = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");

  //_nlayers_seeding = _seeding_layer.size();
  //_radii.assign(_nlayers_seeding, 0.0);
  map<double, int> radius_layer_map;

  _radii_all.assign(60, 0.0);
  _layer_ilayer_map.clear();
  _layer_ilayer_map_all.clear();
  if (cellgeos)
  {
    PHG4TpcCylinderGeomContainer::ConstRange layerrange =
        cellgeos->get_begin_end();
    for (PHG4TpcCylinderGeomContainer::ConstIterator layeriter =
             layerrange.first;
         layeriter != layerrange.second; ++layeriter)
    {
      radius_layer_map.insert(
          make_pair(layeriter->second->get_radius(),
                    layeriter->second->get_layer()));
    }
  }

  if (laddergeos)
  {
    PHG4CylinderGeomContainer::ConstRange layerrange =
        laddergeos->get_begin_end();
    for (PHG4CylinderGeomContainer::ConstIterator layeriter =
             layerrange.first;
         layeriter != layerrange.second; ++layeriter)
    {
      radius_layer_map.insert(
          make_pair(layeriter->second->get_radius(),
                    layeriter->second->get_layer()));
    }
  }

  if (mapsladdergeos)
  {
    PHG4CylinderGeomContainer::ConstRange layerrange =
        mapsladdergeos->get_begin_end();
    for (PHG4CylinderGeomContainer::ConstIterator layeriter =
             layerrange.first;
         layeriter != layerrange.second; ++layeriter)
    {
      radius_layer_map.insert(
          make_pair(layeriter->second->get_radius(),
                    layeriter->second->get_layer()));
    }
  }
  for (map<double, int>::iterator iter = radius_layer_map.begin();
       iter != radius_layer_map.end(); ++iter)
  {
    _layer_ilayer_map_all.insert(make_pair(iter->second, _layer_ilayer_map_all.size()));

    /*if (std::find(_seeding_layer.begin(), _seeding_layer.end(),
                  iter->second) != _seeding_layer.end())
    {
      _layer_ilayer_map.insert(make_pair(iter->second, ilayer));
      ++ilayer;
      }*/
  }
  if (cellgeos)
  {
    PHG4TpcCylinderGeomContainer::ConstRange begin_end =
        cellgeos->get_begin_end();
    PHG4TpcCylinderGeomContainer::ConstIterator miter = begin_end.first;
    for (; miter != begin_end.second; ++miter)
    {
      PHG4TpcCylinderGeom *geo = miter->second;
      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
          geo->get_radius() + 0.5 * geo->get_thickness();

      /*if (_layer_ilayer_map.find(geo->get_layer()) != _layer_ilayer_map.end())
      {
        _radii[_layer_ilayer_map[geo->get_layer()]] =
            geo->get_radius();
	    }*/
    }
  }

  if (laddergeos)
  {
    PHG4CylinderGeomContainer::ConstRange begin_end =
        laddergeos->get_begin_end();
    PHG4CylinderGeomContainer::ConstIterator miter = begin_end.first;
    for (; miter != begin_end.second; ++miter)
    {
      PHG4CylinderGeom *geo = miter->second;
      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
          geo->get_radius() + 0.5 * geo->get_thickness();

      /*if (_layer_ilayer_map.find(geo->get_layer()) != _layer_ilayer_map.end())
      {
        _radii[_layer_ilayer_map[geo->get_layer()]] = geo->get_radius();
	}*/
    }
  }

  if (mapsladdergeos)
  {
    PHG4CylinderGeomContainer::ConstRange begin_end =
        mapsladdergeos->get_begin_end();
    PHG4CylinderGeomContainer::ConstIterator miter = begin_end.first;
    for (; miter != begin_end.second; ++miter)
    {
      PHG4CylinderGeom *geo = miter->second;

      //if(geo->get_layer() > (int) _radii.size() ) continue;

      //			if (Verbosity() >= 2)
      //				geo->identify();

      //TODO
      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
          geo->get_radius();

      /*if (_layer_ilayer_map.find(geo->get_layer()) != _layer_ilayer_map.end())
      {
        _radii[_layer_ilayer_map[geo->get_layer()]] = geo->get_radius();
	}*/
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHHybridSeeding::Process(PHCompositeNode */*topNode*/)
{
  // wipe previous vertex coordinates
  _vertex_x.clear();
  _vertex_y.clear();
  _vertex_z.clear();
  _vertex_xerr.clear();
  _vertex_yerr.clear();
  _vertex_zerr.clear();
  _vertex_ids.clear();
  // fill new vertex coordinates
  for(map<unsigned int, SvtxVertex*>::iterator iter = _vertex_map->begin(); iter != _vertex_map->end(); ++iter)
  {
    SvtxVertex* v = dynamic_cast<SvtxVertex*>(iter->second->CloneMe());
    _vertex_x.push_back(v->get_x());
    _vertex_y.push_back(v->get_y());
    _vertex_z.push_back(v->get_z());
    _vertex_xerr.push_back(sqrt(v->get_error(0,0)));
    _vertex_yerr.push_back(sqrt(v->get_error(1,1)));
    _vertex_zerr.push_back(sqrt(v->get_error(2,2)));
    _vertex_ids.push_back(v->get_id());
  }
  if(Verbosity()>1) cout << "vertices:\n";
  for(size_t i=0;i<_vertex_x.size();i++)
  {
    if(Verbosity()>1) cout << "(" << _vertex_x[i] << "," << _vertex_y[i] << "," << _vertex_z[i] << ")\n";
  }
  vector<vector<double>> kdhits(PHTpcTrackerUtil::convert_clusters_to_hits(_cluster_map,_hitsets));
  vector<vector<double> > unused_hits;
  vector<vector<vector<double> > > kdtracks;

  bool print_stats = (Verbosity()>0);

  kdtracks = kdfinder::find_tracks_iterative<double>(kdhits, unused_hits,
                                                     _search_radius1, /* max distance in cm*/
                                                     _search_angle1, /* triplet angle */
                                                     _min_track_size1, /* min hits to keep track */
                                                     // first iteration
                                                     _search_radius2, 
                                                     _search_angle2, 
                                                     _min_track_size2,                                                                           // second iteration params
                                                     _nthreads,
                                                     print_stats);
  if(Verbosity()>0) cout << "n_kdtracks: " << kdtracks.size() << "\n";
  vector<keylist> clusterLists;
  for(auto track : kdtracks)
  {
    keylist k;
    for(auto kdhit : track)
    {
      // see PHTpcTracker/PHTpcTrackerUtil.cc; this recovers the cluster key, apparently
       k.push_back(*((int64_t*)&kdhit[3])); 
    }
    clusterLists.push_back(k);
  }
  for(auto& clist : clusterLists)
  {
    if(Verbosity()>1) cout << "front layer: " << (int)TrkrDefs::getLayer(clist.front()) << " back layer: " << (int)TrkrDefs::getLayer(clist.back());
    if(clist.size()>1 && ((int)TrkrDefs::getLayer(clist.front()))<((int)TrkrDefs::getLayer(clist.back())))
    {
      if(Verbosity()>1) cout << "reversing\n";
      std::reverse(clist.begin(),clist.end());
    }
    if(Verbosity()>1) cout << "final layer order:\n";
    for(auto c : clist) if(Verbosity()>1) cout << (int)TrkrDefs::getLayer(c) << endl;
  }
  for(auto clist : clusterLists)
  {
    if(Verbosity()>1) cout << "layers:\n";
    for(auto c : clist)
    {
      if(Verbosity()>1) cout << (int)TrkrDefs::getLayer(c) << endl;
    }
  }
  vector<SvtxTrack_v3> seeds = fitter->ALICEKalmanFilter(clusterLists,false);
  if(Verbosity()>0) cout << "nseeds: " << seeds.size() << "\n";
  publishSeeds(seeds);
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHHybridSeeding::publishSeeds(vector<SvtxTrack_v3> seeds)
{
  for(size_t i=0;i<seeds.size();i++)
  {
    _track_map->insert(&(seeds[i]));
  }
}

int PHHybridSeeding::Setup(PHCompositeNode *topNode)
{
  if(Verbosity()>0) cout << "Called Setup" << endl;
  if(Verbosity()>0) cout << "topNode:" << topNode << endl;
  PHTrackSeeding::Setup(topNode);

  _vertex_map = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_vertex_map)
    {
      cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  auto surfmaps = findNode::getClass<ActsSurfaceMaps>(topNode, "ActsSurfaceMaps");
  if(!surfmaps)
    {
      std::cout << "No Acts surface maps, bailing." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  auto tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode,"ActsTrackingGeometry");
  if(!tGeometry)
    {
      std::cout << "No Acts tracking geometry, exiting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  InitializeGeometry(topNode);
  fitter = std::make_shared<ALICEKF>(topNode,_cluster_map,surfmaps, tGeometry,
				     _fieldDir,_min_fit_track_size,_max_sin_phi,Verbosity()); 
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHHybridSeeding::End()
{
  if(Verbosity()>0) cout << "Called End " << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
