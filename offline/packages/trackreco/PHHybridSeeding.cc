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
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>  // for getLayer, clu...
#include <trackbase/TrkrClusterHitAssoc.h>

// sPHENIX Geant4 includes
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

// sPHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHTimer.h>  // for PHTimer
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

//ROOT includes
#include <TVector3.h>  // for TVector3

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

#if defined(_DEBUG_)
#define LogDebug(exp) if(Verbosity()>0) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#else
#define LogDebug(exp) (void)0
#endif

#define LogError(exp) if(Verbosity()>0) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) if(Verbosity()>0) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

//end

typedef std::vector<TrkrDefs::cluskey> keylist;
using namespace std;

PHHybridSeeding::PHHybridSeeding(
    const string &name,
    unsigned int min_clusters_per_track,
    float cluster_z_error,
    float cluster_alice_y_error,
    float maxSinPhi,
    float Bz,
    float search_radius1,
    float search_angle1,
    size_t min_track_size1,
    float search_radius2,
    float search_angle2,
    size_t min_track_size2,
    size_t nthreads
    )
  : PHTrackSeeding(name)
  , _min_clusters_per_track(min_clusters_per_track)
  , _cluster_z_error(cluster_z_error)
  , _cluster_alice_y_error(cluster_alice_y_error)
  , _max_sin_phi(maxSinPhi)
  , _Bz(Bz)
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
  PHG4CylinderCellGeomContainer *cellgeos = findNode::getClass<
      PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  PHG4CylinderGeomContainer *laddergeos = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  PHG4CylinderGeomContainer *mapsladdergeos = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");

  //_nlayers_seeding = _seeding_layer.size();
  //_radii.assign(_nlayers_seeding, 0.0);
  map<float, int> radius_layer_map;

  _radii_all.assign(60, 0.0);
  _layer_ilayer_map.clear();
  _layer_ilayer_map_all.clear();
  if (cellgeos)
  {
    PHG4CylinderCellGeomContainer::ConstRange layerrange =
        cellgeos->get_begin_end();
    for (PHG4CylinderCellGeomContainer::ConstIterator layeriter =
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
  for (map<float, int>::iterator iter = radius_layer_map.begin();
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
    PHG4CylinderCellGeomContainer::ConstRange begin_end =
        cellgeos->get_begin_end();
    PHG4CylinderCellGeomContainer::ConstIterator miter = begin_end.first;
    for (; miter != begin_end.second; ++miter)
    {
      PHG4CylinderCellGeom *geo = miter->second;
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

int PHHybridSeeding::Process(PHCompositeNode *topNode)
{
  vector<vector<double>> kdhits(PHTpcTrackerUtil::convert_clusters_to_hits(_cluster_map));
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
  for(auto clist : clusterLists)
  {
    if(clist.size()>1 && TrkrDefs::getLayer(clist[0])<TrkrDefs::getLayer(clist[1])) std::reverse(clist.begin(),clist.end());
  }
  vector<SvtxTrack_v1> seeds = ALICEKalmanFilter(clusterLists,true);
  if(Verbosity()>0) cout << "nseeds: " << seeds.size() << "\n";
  publishSeeds(seeds);
  return Fun4AllReturnCodes::EVENT_OK;
}

bool PHHybridSeeding::checknan(float val, std::string name, int num)
{
  if(std::isnan(val))
  {
    std::cout << "WARNING: " << name << " is NaN for seed " << num << ". Aborting this seed.\n";
  }
  return std::isnan(val);
}

vector<SvtxTrack_v1> PHHybridSeeding::ALICEKalmanFilter(vector<keylist> trackSeedKeyLists,bool use_nhits_limit)
{
  vector<SvtxTrack_v1> seeds_vector;
  int nseeds = 0;
  if(Verbosity()>0) std::cout << "min clusters per track: " << _min_clusters_per_track << "\n";
  for(vector<keylist>::iterator trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
  {
    if(trackKeyChain->size()<2) continue;
    if(use_nhits_limit && trackKeyChain->size() < _min_clusters_per_track) continue;
    // get starting cluster from key
    TrkrCluster* startCluster = _cluster_map->findCluster(trackKeyChain->at(0));
    // Transform sPHENIX coordinates into ALICE-compatible coordinates
    double x0 = startCluster->getPosition(0);
    double y0 = startCluster->getPosition(1);
    double z0 = startCluster->getPosition(2);
    LogDebug("Initial (x,y,z): (" << x0 << "," << y0 << "," << z0 << ")" << endl);
    // ALICE x coordinate = distance from beampipe
    float alice_x0 = sqrt(x0*x0+y0*y0);
    float alice_y0 = 0;
    float alice_z0 = z0;
    // Initialize track and linearisation
    GPUTPCTrackParam trackSeed;
    trackSeed.InitParam();
    trackSeed.SetX(alice_x0);
    trackSeed.SetY(alice_y0);
    trackSeed.SetZ(alice_z0);
    float x = x0;
    float y = y0;
    #if defined(_DEBUG_)
    float z = z0;
    float alice_x = sqrt(x0*x0+y0*y0);
    #endif
    float trackCartesian_x = 0.;
    float trackCartesian_y = 0.;
    // float trackCartesian_z = 0.;
    // Pre-set momentum-based parameters to improve numerical stability
    TrkrCluster* SecondCluster = _cluster_map->findCluster(trackKeyChain->at(1));
    float second_x = SecondCluster->getPosition(0);
    float second_y = SecondCluster->getPosition(1);
    float second_z = SecondCluster->getPosition(2);
    float second_alice_x = sqrt(second_x*second_x+second_y*second_y);
    float delta_alice_x = second_alice_x - alice_x0;
    float first_phi = atan2(y0,x0);
    float second_alice_y = (second_x/cos(first_phi)-second_y/sin(first_phi))/(sin(first_phi)/cos(first_phi)+cos(first_phi)/sin(first_phi));
    float init_SinPhi = second_alice_y / sqrt(delta_alice_x*delta_alice_x + second_alice_y*second_alice_y);
    float delta_z = second_z - z0;
    float init_DzDs = delta_z / sqrt(delta_alice_x*delta_alice_x + second_alice_y*second_alice_y + delta_z*delta_z);
    trackSeed.SetSinPhi(init_SinPhi);
    LogDebug("Set initial SinPhi to " << init_SinPhi << endl);
    trackSeed.SetDzDs(init_DzDs);
    LogDebug("Set initial DzDs to " << init_DzDs << endl);
    GPUTPCTrackLinearisation trackLine(trackSeed);

    LogDebug(endl << endl << "------------------------" << endl << "seed size: " << trackKeyChain->size() << endl << endl << endl);
    int cluster_ctr = 1;
    bool aborted = false;
    // starting at second cluster, perform track propagation
    for(keylist::iterator clusterkey = next(trackKeyChain->begin()); clusterkey != trackKeyChain->end(); ++clusterkey)
    {
      LogDebug("cluster " << cluster_ctr << " -> " << cluster_ctr + 1 << endl);
      LogDebug("this cluster (x,y,z) = (" << x << "," << y << "," << z << ")" << endl);
      // get cluster from key
      TrkrCluster* nextCluster = _cluster_map->findCluster(*clusterkey);
      // find ALICE x-coordinate
      float nextCluster_x = nextCluster->getPosition(0);
      float nextCluster_y = nextCluster->getPosition(1);
      float nextCluster_z = nextCluster->getPosition(2);
      float nextAlice_x = sqrt(nextCluster_x*nextCluster_x+nextCluster_y*nextCluster_y);
      // rotate track coordinates to match orientation of next cluster
      float newPhi = atan2(nextCluster_y,nextCluster_x);
      LogDebug("new phi = " << newPhi << endl);
      float oldPhi = atan2(y,x);
      LogDebug("old phi = " << oldPhi << endl);
      float alpha = newPhi - oldPhi;
      LogDebug("alpha = " << alpha << endl);
      if(!trackSeed.Rotate(alpha,trackLine,_max_sin_phi))
      {
        LogWarning("Rotate failed! Aborting for this seed...\n");
        aborted = true;
        break;
      }
      LogDebug("track coordinates (ALICE) after rotation: (" << trackSeed.GetX() << "," << trackSeed.GetY() << "," << trackSeed.GetZ() << ")" << endl);
      LogDebug("Transporting from " << alice_x << " to " << nextAlice_x << "...");
      if(!trackSeed.TransportToX(nextAlice_x,trackLine,_Bz,_max_sin_phi))
      {
        LogWarning("Transport failed! Aborting for this seed...\n");
        aborted = true;
        break;
      }
      // convert ALICE coordinates to sPHENIX cartesian coordinates, for debugging
      float predicted_alice_x = trackSeed.GetX();
      LogDebug("new track ALICE x = " << trackSeed.GetX() << endl);
      float predicted_alice_y = trackSeed.GetY();
      LogDebug("new track ALICE y = " << trackSeed.GetY() << endl);
      // float predicted_z = trackSeed.GetZ();
      LogDebug("new track z = " << trackSeed.GetZ() << endl);
      float cos_phi = x/sqrt(x*x+y*y);
      LogDebug("cos_phi = " << cos_phi << endl);
      float sin_phi = y/sqrt(x*x+y*y);
      LogDebug("sin phi = " << sin_phi << endl);
      trackCartesian_x = predicted_alice_x*cos_phi-predicted_alice_y*sin_phi;
      trackCartesian_y = predicted_alice_x*sin_phi+predicted_alice_y*cos_phi;
      // trackCartesian_z = predicted_z;
      LogDebug("Track transported to (x,y,z) = (" << trackCartesian_x << "," << trackCartesian_y << "," << trackCartesian_z << ")" << endl);
      LogDebug("Next cluster is at (x,y,z) = (" << nextCluster_x << "," << nextCluster_y << "," << nextCluster_z << ")" << endl);
      LogDebug("track coordinates (ALICE) after rotation: (" << trackSeed.GetX() << "," << trackSeed.GetY() << "," << trackSeed.GetZ() << ")" << endl);
      //float nextCluster_alice_y = (nextCluster_x/cos(newPhi) - nextCluster_y/sin(newPhi))/(tan(newPhi)+1./tan(newPhi));
      float nextCluster_alice_y = 0.;
      LogDebug("next cluster ALICE y = " << nextCluster_alice_y << endl);
      float y2_error = _cluster_alice_y_error*_cluster_alice_y_error;
      float z2_error = _cluster_z_error*_cluster_z_error;
      LogDebug("track ALICE SinPhi = " << trackSeed.GetSinPhi() << endl);
      // Apply Kalman filter
      if(!trackSeed.Filter(nextCluster_alice_y,nextCluster_z,y2_error,z2_error,_max_sin_phi))
      {
	      if (Verbosity() >= 1)
	      LogError("Kalman filter failed for seed " << nseeds << "! Aborting for this seed..." << endl);
        aborted = true;
        break;
      }
      
      x = nextCluster_x;
      y = nextCluster_y;
      #if defined(_DEBUG_)
      z = nextCluster_z;
      alice_x = nextAlice_x;
      #endif
      ++cluster_ctr;
    }
    if(!trackSeed.TransportToX(0.,trackLine,_Bz,_max_sin_phi))
    {
      LogWarning("Transport failed! Aborting for this seed...\n");
      aborted = true;
      break;
    }

/*    
    if(cluster_ctr!=1 && !trackSeed.CheckNumericalQuality())
    {
      cout << "ERROR: Track seed failed numerical quality check before conversion to sPHENIX coordinates! Skipping this one.\n";
      aborted = true;
      continue;
    } 
*/    
    //    pt:z:dz:phi:dphi:c:dc
    // Fill NT with track parameters
    // float StartEta = -log(tan(atan(z0/sqrt(x0*x0+y0*y0))));
    if(aborted) continue;
    float track_pt = fabs( 1./(trackSeed.GetQPt()));
    if(checknan(track_pt,"pT",nseeds)) continue;
    float track_pterr = sqrt(trackSeed.GetErr2QPt())/(trackSeed.GetQPt()*trackSeed.GetQPt());
    if(checknan(track_pterr,"pT err",nseeds)) continue;
    LogDebug("Track pterr = " << track_pterr << endl);
    float track_z = trackSeed.GetZ();
    if(checknan(track_z,"z",nseeds)) continue;
    float track_zerr = sqrt(trackSeed.GetErr2Z());
    if(checknan(track_zerr,"zerr",nseeds)) continue;
    float track_phi = atan2(trackCartesian_y,trackCartesian_x);
    if(checknan(track_phi,"phi",nseeds)) continue;
    float last_cluster_phierr = _cluster_map->findCluster(trackKeyChain->back())->getPhiError();
    // phi error assuming error in track radial coordinate is zero
    float track_phierr = sqrt(pow(last_cluster_phierr,2)+(pow(trackSeed.GetX(),2)*trackSeed.GetErr2Y()) / 
      pow(pow(trackSeed.GetX(),2)+pow(trackSeed.GetY(),2),2));
    if(checknan(track_phierr,"phierr",nseeds)) continue;
    LogDebug("Track phi = " << track_phi << endl);
    LogDebug("Track phierr = " << track_phierr << endl);
    float track_curvature = trackSeed.GetKappa(_Bz);
    if(checknan(track_curvature,"curvature",nseeds)) continue;
    float track_curverr = sqrt(trackSeed.GetErr2QPt())*_Bz;
    if(checknan(track_curverr,"curvature error",nseeds)) continue;
    SvtxTrack_v1 track;
    track.set_id(nseeds);
    for (unsigned int j = 0; j < trackKeyChain->size(); ++j)
    {
      track.insert_cluster_key(trackKeyChain->at(j));
    }
    track.set_ndf(trackSeed.GetNDF());
    int track_charge = 0;
    if(trackSeed.GetQPt()<0) track_charge = -1 * _fieldDir;
    else track_charge = 1 * _fieldDir;
    track.set_charge(track_charge);
    TrkrCluster *cl = _cluster_map->findCluster(trackKeyChain->at(0));
    track.set_x(cl->getX());  //track.set_x(cl->getX());
    track.set_y(cl->getY());  //track.set_y(cl->getY());
    track.set_z(cl->getZ());  //track.set_z(cl->getZ());
    float s = sin(track_phi);
    float c = cos(track_phi);
    float p = trackSeed.GetSinPhi();
    if(checknan(p,"ALICE sinPhi",nseeds)) continue;
    float d = trackSeed.GetDzDs();
    if(checknan(d,"ALICE dz/ds",nseeds)) continue;
    float pY = track_pt*p;
    float pX = sqrt(track_pt*track_pt-pY*pY);
    track.set_px(pX*c-pY*s);
    track.set_py(pX*s+pY*c);
    track.set_pz(track_pt * trackSeed.GetDzDs()); 
    const float* cov = trackSeed.GetCov();
    bool cov_nan = false;
    for(int i=0;i<15;i++)
    {
      if(checknan(cov[i],"covariance element "+std::to_string(i),nseeds)) cov_nan = true;
    }
    if(cov_nan) continue;
    // make this into an actual Eigen matrix
    Eigen::Matrix<float,5,5> ecov;
    ecov(0,0)=cov[0];
    ecov(0,1)=cov[1];
    ecov(0,2)=cov[2];
    ecov(0,3)=cov[3];
    ecov(0,4)=cov[4];
    ecov(1,1)=cov[5];
    ecov(1,2)=cov[6];
    ecov(1,3)=cov[7];
    ecov(1,4)=cov[8];
    ecov(2,2)=cov[9];
    ecov(2,3)=cov[10];
    ecov(2,4)=cov[11];
    ecov(3,3)=cov[12];
    ecov(3,4)=cov[13];
    ecov(4,4)=cov[14];
    // symmetrize
    ecov(1,0)=ecov(0,1);
    ecov(2,0)=ecov(0,2);
    ecov(3,0)=ecov(0,3);
    ecov(4,0)=ecov(0,4);
    ecov(2,1)=ecov(1,2);
    ecov(3,1)=ecov(1,3);
    ecov(4,1)=ecov(1,4);
    ecov(3,2)=ecov(2,3);
    ecov(4,2)=ecov(2,4);
    ecov(4,3)=ecov(3,4);
    // make rotation matrix based on the following:
    // x = X*cos(track_phi) - Y*sin(track_phi)
    // y = X*sin(track_phi) + Y*cos(track_phi)
    // z = Z
    // pY = pt*sinphi
    // pX = sqrt(pt**2 - pY**2)
    // px = pX*cos(track_phi) - pY*sin(track_phi)
    // py = pX*sin(track_phi) + pY*cos(track_phi)
    // pz = pt*(dz/ds)
    Eigen::Matrix<float,6,5> J;
    J(0,0) = -s; // dx/dY
    J(0,1) = 0.; // dx/dZ
    J(0,2) = 0.; // dx/d(sinphi)
    J(0,3) = 0.; // dx/d(dz/ds)
    J(0,4) = 0.; // dx/d(Q/pt)

    J(1,0) = c;  // dy/dY
    J(1,1) = 0.; // dy/dZ
    J(1,2) = 0.; // dy/d(sinphi)
    J(1,3) = 0.; // dy/d(dz/ds)
    J(1,4) = 0.; // dy/d(Q/pt)

    J(2,0) = 0.; // dz/dY
    J(2,1) = 1.; // dz/dZ
    J(2,2) = 0.; // dz/d(sinphi)
    J(2,3) = 0.; // dz/d(dz/ds)
    J(2,4) = 0.; // dz/d(Q/pt)

    J(3,0) = 0.; // dpx/dY
    J(3,1) = 0.; // dpx/dZ
    J(3,2) = -track_pt*(p*c/sqrt(1-p*p)+s); // dpx/d(sinphi)
    J(3,3) = 0.; // dpx/d(dz/ds)
    J(3,4) = track_pt*track_pt*track_charge*(p*s-c*sqrt(1-p*p)); // dpx/d(Q/pt)

    J(4,0) = 0.; // dpy/dY
    J(4,1) = 0.; // dpy/dZ
    J(4,2) = track_pt*(c-p*s/sqrt(1-p*p)); // dpy/d(sinphi)
    J(4,3) = 0.; // dpy/d(dz/ds)
    J(4,4) = -track_pt*track_pt*track_charge*(p*c+s*sqrt(1-p*p)); // dpy/d(Q/pt)

    J(5,0) = 0.; // dpz/dY
    J(5,1) = 0.; // dpz/dZ
    J(5,2) = 0.; // dpz/d(sinphi)
    J(5,3) = track_pt; // dpz/d(dz/ds)
    J(5,4) = -track_pt*track_pt*track_charge*d; // dpz/d(Q/pt)
    bool cov_rot_nan = false;
    for(int i=0;i<6;i++)
    {
      for(int j=0;j<5;j++)
      {
        if(checknan(J(i,j),"covariance rotator element ("+std::to_string(i)+","+std::to_string(j)+")",nseeds))
        {
          cov_rot_nan = true;
          continue;
        }
      }
    }
    if(cov_rot_nan) continue;

    // the heavy lifting happens here
    Eigen::Matrix<float,6,6> scov = J*ecov*J.transpose();
    
    // fill SvtxTrack covariance matrix with results
    for(int i=0;i<6;i++)
    {
      for(int j=0;j<6;j++)
      {
        track.set_error(i, j, scov(i,j));
      }
    }
/*
    // Proceed with the absolutely hellish coordinate transformation of the covariance matrix.
    // Derived from:
    // 1) Taking the Jacobian of the conversion from (Y,Z,SinPhi,DzDs,Q/Pt) to (x,y,z,px,py,pz)
    // 2) Computing (Jacobian)*(ALICE covariance matrix)*(transpose of Jacobian)
    track.set_error(0, 0, cov[0]*s*s);
    track.set_error(0, 1, -cov[0]*c*s);
    track.set_error(0, 2, -cov[1]*s);
    track.set_error(0, 3, cov[2]*s*s/q-cov[4]*s*(-c/(q*q)+p*s/(q*q)));
    track.set_error(0, 4, -cov[2]*c*s/q-cov[4]*s*(-c*p/(q*q)-s/(q*q)));
    track.set_error(0, 5, cov[4]*d*s/(q*q)-cov[3]*s/q);
    track.set_error(1, 1, cov[0]*c*c);
    track.set_error(1, 2, cov[1]*c);
    track.set_error(1, 3, -cov[2]*c*s/q+cov[4]*c*(-c/(q*q)+p*s/(q*q)));
    track.set_error(1, 4, cov[2]*c*c/q+cov[4]*c*(-c*p/(q*q)-s/(q*q)));
    track.set_error(1, 5, cov[4]*d*c/(q*q)+cov[3]*c/q);
    track.set_error(2, 2, cov[5]);
    track.set_error(2, 3, -cov[6]*s/q+cov[8]*(-c/(q*q)+p*s/(q*q)));
    track.set_error(2, 4, cov[6]*c/q+cov[8]*(-c*p/(q*q)-s/(q*q)));
    track.set_error(2, 5, -cov[8]*d/(q*q)+cov[7]/q);
    track.set_error(3, 3, cov[9]*s*s/(q*q)-cov[11]*(-c/(q*q*q)+p*s/(q*q*q)) + (-c/(q*q)+p*s/(q*q))*(-cov[11]*s/q+cov[14]*(-c/(q*q)+p*s/(q*q))));
    track.set_error(3, 4, -cov[9]*c*s/(q*q)+cov[11]*(-c/(q*q*q)+p*s/(q*q*q)) + (-c*p/(q*q)-s/(q*q))*(-cov[11]*s/q+cov[14]*(-c/(q*q)+p*s/(q*q))));
    track.set_error(3, 5, -cov[10]*s/(q*q)+cov[13]/q*(-c/(q*q)+p*s/(q*q))-d/(q*q)*(-cov[11]*s/q+cov[14]*(-c/(q*q)+p*s/(q*q))));
    track.set_error(4, 4, c/q*(c/q*cov[9]+cov[11]*(-c*p/(q*q)-s/(q*q)))+(-c*p/(q*q)-s/(q*q))*(c/q*cov[11]+cov[14]*(-c*p/(q*q)-s/(q*q))));
    track.set_error(4, 5, cov[10]*c/(q*q)+cov[13]/q*(-c*p/(q*q)-s/(q*q))-d/(q*q)*(c/q*cov[11]+cov[14]*(-c*p/(q*q)-s/(q*q))));
    track.set_error(5, 5, -d/(q*q)*(-d*cov[14]/(q*q)+cov[13]/q)-d*cov[13]/(q*q*q)+cov[12]/(q*q));
    // symmetrize covariance
    track.set_error(1, 0, track.get_error(0, 1));
    track.set_error(2, 0, track.get_error(0, 2));
    track.set_error(3, 0, track.get_error(0, 3));
    track.set_error(4, 0, track.get_error(0, 4));
    track.set_error(5, 0, track.get_error(0, 5));
    track.set_error(2, 1, track.get_error(1, 2));
    track.set_error(3, 1, track.get_error(1, 3));
    track.set_error(4, 1, track.get_error(1, 4));
    track.set_error(5, 1, track.get_error(1, 5));
    track.set_error(3, 2, track.get_error(2, 3));
    track.set_error(4, 2, track.get_error(2, 4));
    track.set_error(5, 2, track.get_error(2, 5));
    track.set_error(4, 3, track.get_error(3, 4));
    track.set_error(5, 3, track.get_error(3, 5));
    track.set_error(5, 4, track.get_error(4, 5));
*/

    if(!covIsPosDef(track))
    {
      repairCovariance(track);
    }


    seeds_vector.push_back(track);
    ++nseeds;
  }
  if(Verbosity()>0) cout << "number of seeds: " << nseeds << "\n";
  return seeds_vector;
}

void PHHybridSeeding::publishSeeds(vector<SvtxTrack_v1> seeds)
{
  for(size_t i=0;i<seeds.size();i++)
  {
    _track_map->insert(&(seeds[i]));
  }
}

Eigen::Matrix<float,6,6> PHHybridSeeding::getEigenCov(SvtxTrack_v1 &track)
{
  Eigen::Matrix<float,6,6> cov;
  for(int i=0;i<6;i++)
  {
    for(int j=0;j<6;j++)
    {
      cov(i,j) = track.get_error(i,j);
    }
  }
  return cov;
}

bool PHHybridSeeding::covIsPosDef(SvtxTrack_v1 &track)
{
  // put covariance matrix into Eigen container
  Eigen::Matrix<float,6,6> cov = getEigenCov(track);
  // attempt Cholesky decomposition
  Eigen::LLT<Eigen::Matrix<float,6,6>> chDec(cov);
  // if Cholesky decomposition does not exist, matrix is not positive definite
  return (chDec.info() != Eigen::NumericalIssue);
}

void PHHybridSeeding::repairCovariance(SvtxTrack_v1 &track)
{
  // find closest positive definite matrix
  Eigen::Matrix<float,6,6> cov = getEigenCov(track);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<float,6,6>> solver(cov);
  Eigen::Matrix<float,6,1> D = solver.eigenvalues();
  Eigen::Matrix<float,6,6> Q = solver.eigenvectors();
  Eigen::Matrix<float,6,1> Dp = D.cwiseMax(1e-15);
  Eigen::Matrix<float,6,6> Z = Q*Dp.asDiagonal()*Q.transpose();
  // updates covariance matrix
  for(int i=0;i<6;i++)
  {
    for(int j=0;j<6;j++)
    {
      track.set_error(i,j,Z(i,j));
    }
  }
}

int PHHybridSeeding::Setup(PHCompositeNode *topNode)
{
  if(Verbosity()>0) cout << "Called Setup" << endl;
  if(Verbosity()>0) cout << "topNode:" << topNode << endl;
  PHTrackSeeding::Setup(topNode);
  InitializeGeometry(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHHybridSeeding::End()
{
  if(Verbosity()>0) cout << "Called End " << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
