#include "PHSimpleKFProp.h"

#include "AssocInfoContainer.h"
#include "nanoflann.hpp"
#include "GPUTPCTrackParam.h"
#include "GPUTPCTrackLinearisation.h"

#include <Eigen/Core>
#include <Eigen/Dense>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxTrack_v2.h>

#include <phfield/PHField.h>
#include <phfield/PHFieldUtility.h>
#include <Geant4/G4SystemOfUnits.hh>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrDefs.h>

#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>                       // for PHWHERE

#include <iostream>                            // for operator<<, basic_ostream
#include <vector>

//#define _DEBUG_

#if defined(_DEBUG_)
#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#else
#define LogDebug(exp) (void)0
#endif

#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

using namespace std;

PHSimpleKFProp::PHSimpleKFProp(const std::string& name)
  : PHTrackPropagating(name)
  , _field_map(nullptr)
{
  //cout << "created PHSimpleKFProp\n";
}
/*
int PHSimpleKFProp::InitRun(PHCompositeNode* topNode)
{
  cout << "initrun started\n";
  return Setup(topNode);
}

int PHSimpleKFProp::process_event(PHCompositeNode* topNode)
{
  cout << "process_event started\n";
  return Process();
}

int PHSimpleKFProp::End(PHCompositeNode* topNode)
{
  End();
  return Fun4AllReturnCodes::EVENT_OK;
}
*/
int PHSimpleKFProp::End()
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSimpleKFProp::Setup(PHCompositeNode* topNode)
{
  int ret = get_nodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  fitter = std::make_shared<ALICEKF>(topNode,_cluster_map,_fieldDir,_min_clusters_per_track,_max_sin_phi,Verbosity());
  fitter->useConstBField(_use_const_field);
  fitter->useFixedClusterError(_use_fixed_clus_err);
  fitter->setFixedClusterError(0,_fixed_clus_err.at(0));
  fitter->setFixedClusterError(1,_fixed_clus_err.at(1));
  fitter->setFixedClusterError(2,_fixed_clus_err.at(2));
  _field_map = PHFieldUtility::GetFieldMapNode(nullptr,topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

double PHSimpleKFProp::get_Bz(double x, double y, double z)
{
  if(_use_const_field) return 1.4;
  double p[4] = {x*cm,y*cm,z*cm,0.*cm};
  double bfield[3];
  _field_map->GetFieldValue(p,bfield);
  return bfield[2]/tesla;
}

int PHSimpleKFProp::get_nodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  if(_use_truth_clusters)
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_TRUTH");
  else
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  if (!_cluster_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _vertex_map = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_vertex_map)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name);
  if (!_track_map)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << _track_map_name << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _assoc_container = findNode::getClass<AssocInfoContainer>(topNode, "AssocInfoContainer");
  if (!_assoc_container)
  {
    cerr << PHWHERE << " ERROR: Can't find AssocInfoContainer." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  PHG4CylinderCellGeomContainer *geom_container =
      findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cerr << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  for(int i=7;i<=54;i++)
  {
    radii.push_back(geom_container->GetLayerCellGeom(i)->get_radius());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSimpleKFProp::Process()
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
  MoveToFirstTPCCluster();
  PrepareKDTrees();
  std::vector<std::vector<TrkrDefs::cluskey>> new_chains;
  std::vector<SvtxTrack> unused_tracks;
  for(SvtxTrackMap::Iter track_it = _track_map->begin(); track_it != _track_map->end(); )
  {
    // if not a TPC track, ignore
    bool is_tpc = false;
    std::vector<TrkrDefs::cluskey> ckeys;
    SvtxTrack* track = track_it->second;
    std::copy(track->begin_cluster_keys(),track->end_cluster_keys(),std::back_inserter(ckeys));
    for(size_t i=0;i<ckeys.size();i++)
    {
      if(TrkrDefs::getLayer(ckeys[i])>=7) is_tpc = true;
    }
    if(is_tpc)
    {
      if(Verbosity()>0) cout << "is tpc track" << endl;
      new_chains.push_back(PropagateTrack(track));
      ++track_it;
    }
    else
    {
      if(Verbosity()>0) cout << "is NOT tpc track" << endl;
      unused_tracks.push_back(*track);
      ++track_it;
    }
  }
  _track_map->Reset();
  std::vector<SvtxTrack_v2> ptracks = fitter->ALICEKalmanFilter(new_chains,true);
  publishSeeds(ptracks);
  publishSeeds(unused_tracks);
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHSimpleKFProp::PrepareKDTrees()
{
  //***** convert clusters to kdhits, and divide by layer
  std::vector<std::vector<std::vector<double> > > kdhits;
  kdhits.resize(56);
  if (!_cluster_map)
  {
    std::cout << "WARNING: (tracking.PHTpcTrackerUtil.convert_clusters_to_hits) cluster map is not provided" << endl;
    return;
  }
  auto hitsetrange = _hitsets->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (auto hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {
    auto range = _cluster_map->getClusters(hitsetitr->first);
    for( auto it = range.first; it != range.second; ++it )
    {
      TrkrDefs::cluskey cluskey = it->first;
      TrkrCluster* cluster = it->second;
      int layer = TrkrDefs::getLayer(cluskey);
      std::vector<double> kdhit(4);
      kdhit[0] = cluster->getPosition(0);
      kdhit[1] = cluster->getPosition(1);
      kdhit[2] = cluster->getPosition(2);
      uint64_t key = cluster->getClusKey();
      std::memcpy(&kdhit[3], &key, sizeof(key));
    
      //      HINT: way to get original uint64_t value from double:
      //
      //      LOG_DEBUG("tracking.PHTpcTrackerUtil.convert_clusters_to_hits")
      //        << "orig: " << cluster->getClusKey() << ", readback: " << (*((int64_t*)&kdhit[3]));

      kdhits[layer].push_back(kdhit);
    }
  }
  _ptclouds.resize(kdhits.size());
  _kdtrees.resize(kdhits.size());
  for(size_t l=0;l<kdhits.size();++l)
  {
    if(Verbosity()>0) cout << "l: " << l << endl;
    _ptclouds[l] = std::make_shared<KDPointCloud<double>>();
    _ptclouds[l]->pts.resize(kdhits[l].size());
    if(Verbosity()>0) cout << "resized to " << kdhits[l].size() << endl;
    for(size_t i=0;i<kdhits[l].size();++i)
    {
      _ptclouds[l]->pts[i] = kdhits[l][i];
    }
    _kdtrees[l] = std::make_shared<nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, KDPointCloud<double>>,
                                        KDPointCloud<double>, 3>>(3,*_ptclouds[l],nanoflann::KDTreeSingleIndexAdaptorParams(10));
    _kdtrees[l]->buildIndex();
  }
}

std::vector<TrkrDefs::cluskey> PHSimpleKFProp::PropagateTrack(SvtxTrack* track)
{
  // extract cluster list
  std::vector<TrkrDefs::cluskey> ckeys;
  //ckeys.resize(track->size_cluster_keys());
  std::copy(track->begin_cluster_keys(),track->end_cluster_keys(),std::back_inserter(ckeys));
  if(ckeys.size()>1 && ((int)TrkrDefs::getLayer(ckeys.front()))>((int)TrkrDefs::getLayer(ckeys.back())))
  {
    std::reverse(ckeys.begin(),ckeys.end());
  } 
  // get track parameters
  GPUTPCTrackParam kftrack;
  kftrack.InitParam();
  float track_phi = atan2(track->get_py(),track->get_px());
  kftrack.SetQPt(track->get_charge()/track->get_pt());
  float track_pY = -track->get_px()*sin(track_phi)+track->get_py()*cos(track_phi);
  kftrack.SetSinPhi(track_pY/track->get_pt());
  kftrack.SetDzDs(track->get_pz()/track->get_pt());
  Eigen::Matrix<double,6,6> xyzCov;
  for(int i=0;i<6;i++)
  {
     for(int j=0;j<6;j++)
     {
       xyzCov(i,j) = track->get_error(i,j);
     }
  }
  // Y = y
  // Z = z
  // SinPhi = py/sqrt(px^2+py^2)
  // DzDs = pz/sqrt(px^2+py^2)
  // QPt = 1/sqrt(px^2+py^2)

  double track_px = track->get_px();
  double track_py = track->get_py();
  double track_pz = track->get_pz();

  Eigen::Matrix<double,6,5> Jrot;
  Jrot(0,0) = 0; // dY/dx
  Jrot(1,0) = 1; // dY/dy
  Jrot(2,0) = 0; // dY/dz
  Jrot(3,0) = 0; // dY/dpx
  Jrot(4,0) = 0; // dY/dpy
  Jrot(5,0) = 0; // dY/dpz
  
  Jrot(0,1) = 0; // dZ/dx
  Jrot(1,1) = 0; // dZ/dy
  Jrot(2,1) = 1; // dZ/dz
  Jrot(3,1) = 0; // dZ/dpx
  Jrot(4,1) = 0; // dZ/dpy
  Jrot(5,1) = 0; // dZ/dpz

  Jrot(0,2) = 0; // dSinPhi/dx
  Jrot(1,2) = 0; // dSinPhi/dy
  Jrot(2,2) = 0; // dSinPhi/dz
  Jrot(3,2) = -track_py*track_px/pow(track_py*track_py+track_px*track_px,3./2.); // dSinPhi/dpx
  Jrot(4,2) = track_px*track_px/pow(track_py*track_py+track_px*track_px,3./2.); // dSinPhi/dpy
  Jrot(5,2) = 0; // dSinPhi/dpz

  Jrot(0,3) = 0; // dDzDs/dx
  Jrot(1,3) = 0; // dDzDs/dy
  Jrot(2,3) = 0; // dDzDs/dz
  Jrot(3,3) = -track_px*track_pz/pow(track_px*track_px+track_py*track_py,3./2.); // dDzDs/dpx
  Jrot(4,3) = -track_py*track_pz/pow(track_px*track_px+track_py*track_py,3./2.); // dDzDs/dpy
  Jrot(5,3) = 1./sqrt(track_px*track_px+track_py*track_py); // dDzDs/dpz

  Jrot(0,4) = 0; // dQPt/dx
  Jrot(1,4) = 0; // dQPt/dy
  Jrot(2,4) = 0; // dQPt/dz
  Jrot(3,4) = -track_px/pow(track_px*track_px+track_py*track_py,3./2.); // dQPt/dpx
  Jrot(4,4) = -track_py/pow(track_px*track_px+track_py*track_py,3./2.); // dQPt/dpy
  Jrot(5,4) = 0; // dQPt/dpz

  Eigen::Matrix<double,5,5> kfCov = Jrot.transpose()*xyzCov*Jrot;

  int ctr = 0;
  for(int i=0;i<5;i++)
  {
    for(int j=0;j<5;j++)
    {
      if(i>=j)
      {
        kftrack.SetCov(ctr,kfCov(i,j));
        ctr++;
      }
    }
  }

  std::vector<TrkrDefs::cluskey> propagated_track;

  // setup ALICE track model on first cluster
  TrkrCluster* firstclus = _cluster_map->findCluster(ckeys[0]);
//  TrkrCluster* lastclus = _cluster_map->findCluster(ckeys.back());
  double fx = firstclus->getX();
  double fy = firstclus->getY();
//  double fz = firstclus->getZ();
  double fphi = atan2(fy,fx);
//  double lx = lastclus->getX();
//  double ly = lastclus->getY();
//  double lz = lastclus->getZ();
//  double lphi = atan2(ly,lx);
  GPUTPCTrackLinearisation kfline(kftrack);
//  kftrack.Rotate(fphi-lphi,kfline,10.);
  kftrack.SetX(track->get_x()*cos(fphi)+track->get_y()*sin(fphi));
  kftrack.SetY(-track->get_x()*sin(fphi)+track->get_y()*cos(fphi));
  kftrack.SetZ(track->get_z());
  if(Verbosity()>0)
  {
    cout << "initial track params:" << endl;
    cout << "X: " << kftrack.GetX() << endl;
    cout << "Y: " << kftrack.GetY() << endl;
    cout << "Z: " << kftrack.GetZ() << endl;
    cout << "SinPhi: " << kftrack.GetSinPhi() << endl;
    cout << "DzDs: " << kftrack.GetDzDs() << endl;
    cout << "QPt: " << kftrack.GetQPt() << endl;  
  }
//  kftrack.Rotate(fphi,kfline,10.);
//  double fX = fx*cos(fphi)+fy*sin(fphi);
//  double oldphi = track_phi;
//  double track_x = kftrack.GetX()*cos(oldphi)-kftrack.GetY()*sin(oldphi);
//  double track_y = kftrack.GetX()*sin(oldphi)+kftrack.GetX()*cos(oldphi);
//  kftrack.TransportToX(fX,kfline,_Bzconst*get_Bz(track_x,track_y,kftrack.GetZ()),10.);
//  double new_phi = atan2(track_y,track_x);
//  kftrack.Rotate(new_phi-oldphi,kfline,10.);
//  kftrack.SetX(fx*cos(fphi)+fy*sin(fphi));
//  kftrack.SetY(-fx*sin(fphi)+fy*cos(fphi));
//  kftrack.SetZ(fz);

  // copy ALICE track model for later
  GPUTPCTrackParam kftrack_up(kftrack);
  GPUTPCTrackLinearisation kfline_up(kftrack_up);


  // get layer for each cluster
  std::vector<unsigned int> layers;
  if(Verbosity()>0) cout << "cluster layers:" << endl;
  for(size_t i=0;i<ckeys.size();i++)
  {
    layers.push_back(TrkrDefs::getLayer(ckeys[i]));
    if(Verbosity()>0) cout << layers[i] << endl;
  }

  double old_phi = fphi;
  unsigned int old_layer = TrkrDefs::getLayer(ckeys[0]);
  if(Verbosity()>0) cout << "first layer: " << old_layer << endl;

  propagated_track.push_back(ckeys[0]);
  // first, propagate downward
  for(unsigned int l=old_layer+1;l<=54;l++)
  {
    if(Verbosity()>0) cout << "layer " << l << ":" << endl;
    // check to see whether layer is already occupied by at least one cluster
    // choosing the last one first (clusters organized from inside out)
    bool layer_filled = false;
    TrkrDefs::cluskey next_ckey = 0;
    for(int k=ckeys.size()-1; k>=0; k--)
    {
      if(layer_filled) continue;
      if(layers[k]==l)
      {
        layer_filled = true;
        next_ckey = ckeys[k];
      }
    }
    // if layer is already occupied, reset track parameters to last cluster in layer
    if(layer_filled)
    {
      if(Verbosity()>0) cout << "layer is filled" << endl;
      TrkrCluster* nc = _cluster_map->findCluster(next_ckey);
      double cx = nc->getX();
      double cy = nc->getY();
      double cz = nc->getZ();
      double cphi = atan2(cy,cx);
      double alpha = cphi-old_phi;
      double tX = kftrack.GetX();
      double tY = kftrack.GetY();
      double tx = tX*cos(old_phi)-tY*sin(old_phi);
      double ty = tX*sin(old_phi)+tY*cos(old_phi);
      double tz = kftrack.GetZ();
      GPUTPCTrackParam::GPUTPCTrackFitParam fp;
      kftrack.CalculateFitParameters(fp);
      kftrack.Rotate(alpha,kfline,10.);
      kftrack.TransportToX(cx*cos(cphi)+cy*sin(cphi),kfline,_Bzconst*get_Bz(tx,ty,tz),10.);
      if(Verbosity()>0) cout << "track position: (" << tx << ", " << ty << ", " << tz << ")" << endl;
      if(Verbosity()>0) cout << "cluster position: (" << cx << ", " << cy << ", " << cz << ")" << endl;
      if(Verbosity()>0) cout << "new track position: (" << kftrack.GetX()*cos(cphi)-kftrack.GetY()*sin(cphi) << ", " << kftrack.GetX()*sin(cphi)+kftrack.GetY()*cos(cphi) << ", " << kftrack.GetZ() << ")" << endl;
      kftrack.SetX(cx*cos(cphi)+cy*sin(cphi));
      kftrack.SetY(-cx*sin(cphi)+cy*cos(cphi));
      kftrack.SetZ(cz);
      propagated_track.push_back(next_ckey);
      old_phi = cphi;
    }
    // if layer is not occupied, search for the nearest available cluster to projected track position
    else
    {
      if(Verbosity()>0) cout << "layer not filled" << endl;
      double tX = kftrack.GetX();
      double tY = kftrack.GetY();
      double tx = tX*cos(old_phi)-tY*sin(old_phi);
      double ty = tX*sin(old_phi)+tY*cos(old_phi);
      double tz = kftrack.GetZ();
      GPUTPCTrackParam::GPUTPCTrackFitParam fp;
      kftrack.CalculateFitParameters(fp);
      kftrack.TransportToX(radii[l-7],kfline,_Bzconst*get_Bz(tx,ty,tz),10.);
      if(Verbosity()>0) cout << "transported to " << radii[l-7] << "\n";
      if(Verbosity()>0) cout << "track position: (" << tx << ", " << ty << ", " << tz << ")" << endl;
      double query_pt[3] = {tx, ty, tz};
      std::vector<long unsigned int> index_out(1);
      std::vector<double> distance_out(1);
      _kdtrees[l]->knnSearch(&query_pt[0],1,&index_out[0],&distance_out[0]);
      if(Verbosity()>0) cout << "index_out: " << index_out[0] << endl;
      if(Verbosity()>0) cout << "squared_distance_out: " << distance_out[0] << endl;
      if(Verbosity()>0) cout << "solid_angle_dist: " << atan2(sqrt(distance_out[0]),radii[l-7]) << endl;
      std::vector<double> point = _ptclouds[l]->pts[index_out[0]];
      TrkrDefs::cluskey closest_ckey = (*((int64_t*)&point[3]));
      TrkrCluster* cc = _cluster_map->findCluster(closest_ckey);                         
      double ccX = cc->getX();
      double ccY = cc->getY();
      double ccphi = atan2(ccY,ccX);
      if(Verbosity()>0) cout << "cluster position: (" << ccX << ", " << ccY << ", " << cc->getZ() << ")" << endl;
      if(Verbosity()>0) cout << "cluster X: " << ccX*cos(ccphi)+ccY*sin(ccphi) << endl;
      if(atan2(sqrt(distance_out[0]),radii[l-7])<_max_dist)
      {
        propagated_track.push_back(closest_ckey);
/*        TrkrCluster* cc = _cluster_map->findCluster(closest_ckey);
        double ccX = cc->getX();
        cout << "cluster X: " << ccX << endl;
        double ccY = cc->getY();
        double ccx = ccX*cos(old_phi)-ccY*sin(old_phi);
        double ccy = ccX*sin(old_phi)+ccY*cos(old_phi);
        cout << "cluster position: (" << ccx << ", " << ccy << ", " << cc->getZ() << ")" << endl;
*/ 
        
        double alpha = ccphi-old_phi;
        kftrack.Rotate(alpha,kfline,10.);
        kftrack.SetX(ccX*cos(ccphi)+ccY*sin(ccphi));
        kftrack.SetY(-ccX*sin(ccphi)+ccY*cos(ccphi));
        kftrack.SetZ(cc->getZ());
        double ccaY = -ccX*sin(ccphi)+ccY*cos(ccphi);
        double ccerrY = cc->getError(0,0)*sin(ccphi)*sin(ccphi)+cc->getError(0,1)*sin(ccphi)*cos(ccphi)+cc->getError(1,1)*cos(ccphi)*cos(ccphi);
        double ccerrZ = cc->getError(2,2);
        kftrack.Filter(ccaY,cc->getZ(),ccerrY,ccerrZ,_max_sin_phi);
        if(Verbosity()>0) cout << "added cluster" << endl;
        old_phi = ccphi;
      }
    }
    old_layer = l;
  }
  // reverse propagated track so we push_back onto the correct end
  std::reverse(propagated_track.begin(),propagated_track.end());
  old_layer = TrkrDefs::getLayer(ckeys[0]);
  // then, propagate upward
  for(unsigned int l=old_layer-1;l>=7;l--)
  {
    if(Verbosity()>0) cout << "layer " << l << ":" << endl;
    // check to see whether layer is already occupied by at least one cluster
    // choosing the first one first (clusters organized from outside in)
    bool layer_filled = false;
    TrkrDefs::cluskey next_ckey = 0;
    for(size_t k=0; k<ckeys.size(); k++)
    {
      if(layer_filled) continue;
      if(layers[k]==l)
      {
        layer_filled = true;
        next_ckey = ckeys[k];
      }
    }
    // if layer is already occupied, reset track parameters to last cluster in layer
    if(layer_filled)
    {
      if(Verbosity()>0) cout << "layer is filled" << endl;
      TrkrCluster* nc = _cluster_map->findCluster(next_ckey);
      double cx = nc->getX();
      double cy = nc->getY();
      double cz = nc->getZ();
      double cphi = atan2(cy,cx);
      double alpha = cphi-old_phi;
      double tX = kftrack_up.GetX();
      double tY = kftrack_up.GetY();
      double tx = tX*cos(old_phi)-tY*sin(old_phi);
      double ty = tX*sin(old_phi)+tY*cos(old_phi);
      double tz = kftrack_up.GetZ();
      GPUTPCTrackParam::GPUTPCTrackFitParam fp;
      kftrack_up.CalculateFitParameters(fp);
      kftrack_up.Rotate(alpha,kfline_up,10.);
      kftrack_up.TransportToX(cx*cos(cphi)+cy*sin(cphi),kfline_up,_Bzconst*get_Bz(tx,ty,tz),10.);
      if(Verbosity()>0) cout << "transported to " << radii[l-7] << "\n";
      if(Verbosity()>0) cout << "track position: (" << kftrack_up.GetX()*cos(cphi)-kftrack_up.GetY()*sin(cphi) << ", " << kftrack_up.GetX()*sin(cphi)+kftrack_up.GetY()*cos(cphi) << ", " << kftrack_up.GetZ() << ")" << endl;
      if(Verbosity()>0) cout << "cluster position: (" << cx << ", " << cy << ", " << cz << ")" << endl;
      kftrack_up.SetX(cx*cos(cphi)+cy*sin(cphi));
      kftrack_up.SetY(-cx*sin(cphi)+cy*cos(cphi));
      kftrack_up.SetZ(cz);
      propagated_track.push_back(next_ckey);
      old_phi = cphi;
    }
    // if layer is not occupied, search for the nearest available cluster to projected track position
    else
    {
      if(Verbosity()>0) cout << "layer not filled" << endl;
      double tX = kftrack_up.GetX();
      double tY = kftrack_up.GetY();
      double tx = tX*cos(old_phi)-tY*sin(old_phi);
      double ty = tX*sin(old_phi)+tY*cos(old_phi);
      double tz = kftrack_up.GetZ();
//      GPUTPCTrackParam::GPUTPCTrackFitParam fp;
//      kftrack_up.CalculateFitParameters(fp);
      kftrack_up.TransportToX(radii[l-7],kfline_up,_Bzconst*get_Bz(tx,ty,tz),10.);
      if(Verbosity()>0) cout << "transported to " << radii[l-7] << "\n";
      if(Verbosity()>0) cout << "track position: (" << kftrack_up.GetX()*cos(old_phi)-kftrack_up.GetY()*sin(old_phi) << ", " << kftrack_up.GetX()*sin(old_phi)+kftrack_up.GetY()*cos(old_phi) << ", " << kftrack_up.GetZ() << ")" << endl;
      double query_pt[3] = {tx, ty, tz};
      std::vector<long unsigned int> index_out(1);
      std::vector<double> distance_out(1);
      _kdtrees[l]->knnSearch(&query_pt[0],1,&index_out[0],&distance_out[0]);
      if(Verbosity()>0) cout << "index_out: " << index_out[0] << endl;
      if(Verbosity()>0) cout << "squared_distance_out: " << distance_out[0] << endl;
      if(Verbosity()>0) cout << "solid_angle_dist: " << atan2(sqrt(distance_out[0]),radii[l-7]) << endl;
      std::vector<double> point = _ptclouds[l]->pts[index_out[0]];
      TrkrDefs::cluskey closest_ckey = (*((int64_t*)&point[3]));
      TrkrCluster* cc = _cluster_map->findCluster(closest_ckey);
      double ccX = cc->getX();
      double ccY = cc->getY();
      if(Verbosity()>0) cout << "cluster position: (" << ccX << ", " << ccY << ", " << cc->getZ() << ")" << endl;
      double ccphi = atan2(ccY,ccX);
      if(Verbosity()>0) cout << "cluster X: " << ccX*cos(ccphi)+ccY*sin(ccphi) << endl;
      double alpha = ccphi-old_phi;
      if(atan2(sqrt(distance_out[0]),radii[l-7])<_max_dist)
      {
        propagated_track.push_back(closest_ckey);
/*        TrkrCluster* cc = _cluster_map->findCluster(closest_ckey);
        double ccX = cc->getX();
        cout << "cluster X: " << ccX << endl;
        double ccY = cc->getY();
        double ccx = ccX*cos(old_phi)-ccY*sin(old_phi);
        double ccy = ccX*sin(old_phi)+ccY*cos(old_phi);
        cout << "cluster position: (" << ccx << ", " << ccy << ", " << cc->getZ() << ")" << endl;
        double ccphi = atan2(ccy,ccx);
        double alpha = ccphi-old_phi;
*/
        kftrack_up.Rotate(alpha,kfline_up,10.);
        double ccaY = -ccX*sin(ccphi)+ccY*cos(ccphi);
        double ccerrY = cc->getError(0,0)*sin(ccphi)*sin(ccphi)+cc->getError(0,1)*sin(ccphi)*cos(ccphi)+cc->getError(1,1)*cos(ccphi)*cos(ccphi);
        double ccerrZ = cc->getError(2,2);
        kftrack_up.Filter(ccaY,cc->getZ(),ccerrY,ccerrZ,_max_sin_phi);
        kftrack_up.SetX(ccX*cos(ccphi)+ccY*sin(ccphi));
        kftrack_up.SetY(-ccX*sin(ccphi)+ccY*cos(ccphi));
        kftrack_up.SetZ(cc->getZ());
        if(Verbosity()>0) cout << "added cluster" << endl;
        old_phi = ccphi;
      }
    }
    old_layer = l;
  }
  return propagated_track;
}

void PHSimpleKFProp::publishSeeds(vector<SvtxTrack_v2> seeds)
{
  for(size_t i=0;i<seeds.size();i++)
  {
    _track_map->insert(&(seeds[i]));
  }
}

void PHSimpleKFProp::publishSeeds(vector<SvtxTrack> seeds)
{
  for(size_t i=0;i<seeds.size();i++)
  {
    _track_map->insert(&(seeds[i]));
  }
}

void PHSimpleKFProp::MoveToVertex()
{
  // _track_map contains the TPC seed track stubs
  // We want to associate these TPC track seeds with a collision vertex
  // Then we add the collision vertex position as the track seed position

  // All we need is to project the TPC clusters in Z to the beam line.


  if(Verbosity() > 0)
    cout << PHWHERE << " TPC track map size " << _track_map->size()  << endl;
  /*
 // We remember the original size of the TPC track map here
  const unsigned int original_track_map_lastkey = _track_map->end()->first;
  */

  // loop over the TPC track seeds
  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      /*
      // we may add tracks to the map, so we stop at the last original track
      if(phtrk_iter->first >= original_track_map_lastkey)  break;
      */      
      SvtxTrack* _tracklet_tpc = phtrk_iter->second;
      
      if (Verbosity() >= 1)
	{
	  std::cout
	    << __LINE__
	    << ": Processing seed itrack: " << phtrk_iter->first
	    << ": nhits: " << _tracklet_tpc-> size_cluster_keys()
	    << ": phi: " << _tracklet_tpc->get_phi()
	    << ": eta: " << _tracklet_tpc->get_eta()
	    << endl;
	}

      // get the tpc track seed cluster positions in z and r

      // Get the outermost TPC clusters for this tracklet
      std::map<unsigned int, TrkrCluster*> tpc_clusters_map;
      std::vector<TrkrCluster*> clusters;
      
      for (SvtxTrack::ConstClusterKeyIter key_iter = _tracklet_tpc->begin_cluster_keys();
	   key_iter != _tracklet_tpc->end_cluster_keys();
	   ++key_iter)
	{
	  TrkrDefs::cluskey cluster_key = *key_iter;
	  unsigned int layer = TrkrDefs::getLayer(cluster_key);

	  //if(layer < _min_tpc_layer) continue;
	  //if(layer >= _max_tpc_layer) continue;

	  // get the cluster
	  TrkrCluster *tpc_clus =  _cluster_map->findCluster(cluster_key);

	  tpc_clusters_map.insert(std::make_pair(layer, tpc_clus));
	  clusters.push_back(tpc_clus);

	  if(Verbosity() > 5) 
	    std::cout << "  TPC cluster in layer " << layer << " with position " << tpc_clus->getX() 
		      << "  " << tpc_clus->getY() << "  " << tpc_clus->getZ() << " clusters.size() " << tpc_clusters_map.size() << std::endl;
	}


      // need at least 3 clusters to fit a circle
      if(tpc_clusters_map.size() < 3)
	{
	  if(Verbosity() > 3) std::cout << PHWHERE << "  -- skip this tpc tracklet, not enough clusters " << std::endl; 
	  continue;  // skip to the next TPC tracklet
	}

      /*
      // fit a circle to the clusters
      double R, X0, Y0;
      CircleFitByTaubin(clusters, R, X0, Y0);
      if(Verbosity() > 10) std::cout << " Fitted circle has R " << R << " X0 " << X0 << " Y0 " << Y0 << std::endl;
      // toss tracks for which the fitted circle could not have come from the vertex
      if(R < 40.0) continue;
      */

      // get the straight line representing the z trajectory in the form of z vs radius
      double A = 0; double B = 0;
      line_fit_clusters(clusters, A, B);
      if(Verbosity() > 5) std::cout << " Fitted line has A " << A << " B " << B << std::endl;

      // Project this TPC tracklet  to the beam line and store the projections
      //bool skip_tracklet = false;
      // z projection is unique
      double _z_proj = B;
      
      // Find the nearest collision vertex
      // if find multiple close matches, duplicate the track?

      int trackVertexId = 9999;
      double dz = 9999.;	  
      for(SvtxVertexMap::Iter viter = _vertex_map->begin();
	  viter != _vertex_map->end();
	  ++viter)
	{
	  auto vertexKey = viter->first;
	  auto vertex = viter->second;
	  if(Verbosity() > 100)
	    vertex->identify();

	  const double vertexZ = vertex->get_z();
	  
	  if(fabs(_z_proj - vertexZ) < dz )
	    {
	      dz = fabs(_z_proj - vertexZ);
	      trackVertexId = vertexKey;
	    }	  
	}  // end loop over collision vertices

      _tracklet_tpc->set_vertex_id(trackVertexId);
      auto vertex = _vertex_map->find(trackVertexId)->second;
      vertex->insert_track(phtrk_iter->first);

      // set the track position to the vertex position
      _tracklet_tpc->set_x(vertex->get_x());
      _tracklet_tpc->set_y(vertex->get_y());
      _tracklet_tpc->set_z(vertex->get_z());

      if(Verbosity() > 1)
	{
	  std::cout << "TPC seed track " << phtrk_iter->first << " matched to vertex " << trackVertexId << endl; 
	} 

      // Finished association of track with vertex
      // Now we modify the track parameters

      // Repeat the z line fit including the vertex position, get theta, update pz
      std::vector<std::pair<double, double>> points;
      double r_vertex = sqrt(vertex->get_x()*vertex->get_x() + vertex->get_y()*vertex->get_y());
      double z_vertex = vertex->get_z();
      points.push_back(make_pair(r_vertex, z_vertex));
      for (unsigned int i=0; i<clusters.size(); ++i)
	{
	  double z = clusters[i]->getZ();
	  double r = sqrt(pow(clusters[i]->getX(),2) + pow(clusters[i]->getY(), 2));
	  
	  points.push_back(make_pair(r,z));
	}
      
      line_fit(points, A, B);
      if(Verbosity() > 5) 
	std::cout << " Fitted line including vertex has A " << A << " B " << B << std::endl;      

      // extract the track theta
      double track_angle = atan(A);  // referenced to 90 degrees

      //  update pz of track
      double pt_track = _tracklet_tpc->get_pt();
      double ptrack = sqrt(pt_track*pt_track + _tracklet_tpc->get_pz()*_tracklet_tpc->get_pz());
      double pz_new = ptrack * sin(track_angle);
      if(Verbosity() > 5)
	std::cout << " Original pz = " << _tracklet_tpc->get_pz() << " new pz " << pz_new << " track angle " << track_angle << std::endl;
      _tracklet_tpc->set_pz(pz_new);
      if(Verbosity() > 5)
	std::cout << "       new eta " <<  _tracklet_tpc->get_eta() << std::endl;
      // total momentum is now a bit different because pt was not changed - OK - we measure pt from bend, pz from dz/dr

      // make circle fit including vertex as point
      std::vector<std::pair<double, double>> cpoints;
      double x_vertex = vertex->get_x();
      double y_vertex = vertex->get_y();
      cpoints.push_back(std::make_pair(x_vertex, y_vertex));
      for (unsigned int i=0; i<clusters.size(); ++i)
	{
	  double x = clusters[i]->getX();
	  double y = clusters[i]->getY();	  
	  cpoints.push_back(make_pair(x, y));
	}
      double R, X0, Y0;
      CircleFitByTaubin(cpoints, R, X0, Y0);
      if(Verbosity() > 5) 
	std::cout << " Fitted circle has R " << R << " X0 " << X0 << " Y0 " << Y0 << std::endl;

      //  could take new pT from radius of circle - we choose to keep the seed pT

      // We want the angle of the tangent relative to the positive x axis
      // start with the angle of the radial line from vertex to circle center
      double dx = X0 - x_vertex;
      double dy = Y0 - y_vertex;
      double phi= atan2(dy,dx);
      //std::cout << "x_vertex " << x_vertex << " y_vertex " << y_vertex << " X0 " << X0 << " Y0 " << Y0 << " angle " << phi * 180 / 3.14159 << std::endl; 
      // convert to the angle of the tangent to the circle
      // we need to know if the track proceeds clockwise or CCW around the circle
      double dx0 = cpoints[0].first - X0;
      double dy0 = cpoints[0].second - Y0;
      double phi0 = atan2(dy0, dx0);
      double dx1 = cpoints[1].first - X0;
      double dy1 = cpoints[1].second - Y0;
      double phi1 = atan2(dy1, dx1);
      double dphi = phi1 - phi0;

      if(Verbosity() > 5) 
	{
	  int charge = _tracklet_tpc->get_charge();   // needed for diagnostic output only
	  std::cout << " charge " << charge << " phi0 " << phi0*180.0 / M_PI << " phi1 " << phi1*180.0 / M_PI << " dphi " << dphi*180.0 / M_PI << std::endl;
	}

      // whether we add or subtract 90 degrees depends on the track propagation direction determined above
      if(dphi < 0)
	phi += M_PI / 2.0;  
      else
	phi -= M_PI / 2.0;  
      if(Verbosity() > 5) 
	std::cout << " input track phi " << _tracklet_tpc->get_phi() * 180.0 / M_PI << " new phi " << phi * 180 / M_PI << std::endl;  

      // update px, py of track
      double px_new = pt_track * cos(phi);
      double py_new = pt_track * sin(phi);
      if(Verbosity() > 5)
	std::cout << " input track px " << _tracklet_tpc->get_px()  << " new px " << px_new << " input py " << _tracklet_tpc->get_py() << " new py " << py_new << std::endl;

      // update track on node tree
      _tracklet_tpc->set_px(px_new);
      _tracklet_tpc->set_py(py_new);
      
    }  // end loop over TPC track seeds
}

void PHSimpleKFProp::MoveToFirstTPCCluster()
{
  for(auto trkmap_iter = _track_map->begin();
      trkmap_iter != _track_map->end();
      ++trkmap_iter)
  {
    SvtxTrack* track = trkmap_iter->second;
    // check to see if track is, in fact, not in the TPC
    double track_x = track->get_x();
    double track_y = track->get_y();
    if(sqrt(track_x*track_x+track_y*track_y)>10.)
    {
      std::cout << "WARNING: attempting to move track to TPC which is already in TPC! Aborting for this track." << std::endl;
      continue;
    }
    // get cluster keys
    std::vector<TrkrDefs::cluskey> ckeys;
    std::copy(track->begin_cluster_keys(),track->end_cluster_keys(),std::back_inserter(ckeys));
    std::vector<TrkrCluster*> tpc_clusters;
    // extract TPC clusters
    for(auto ckey : ckeys)
    {
      if(TrkrDefs::getLayer(ckey)>=7) tpc_clusters.push_back(_cluster_map->findCluster(ckey));
    }
    // get circle fit for TPC clusters plus vertex
    std::vector<std::pair<double,double>> pts;
    for(auto cluster : tpc_clusters)
    {
      double x = cluster->getX();
      double y = cluster->getY();
      pts.push_back(std::make_pair(x,y));
    }
    SvtxVertex* vertex = _vertex_map->get(track->get_vertex_id());
    double vertex_x = vertex->get_x();
    double vertex_y = vertex->get_y();
    pts.push_back(std::make_pair(vertex_x,vertex_y));
    double R = 0;
    double xc = 0;
    double yc = 0;
    CircleFitByTaubin(pts,R,xc,yc);
    // want angle of tangent to circle at innermost (i.e. last) cluster
    size_t inner_index;
    if(TrkrDefs::getLayer(ckeys[0])>TrkrDefs::getLayer(ckeys.back()))
    {
      inner_index = ckeys.size()-1;
    }
    else
    {
      inner_index = 0;
    }
    double cluster_x = pts[inner_index].first;
    double cluster_y = pts[inner_index].second;
    double dy = cluster_y-yc;
    double dx = cluster_x-xc;
    double phi = atan2(dy,dx);
    double dx0 = pts[0].first - xc;     
    double dy0 = pts[0].second - yc;
    double phi0 = atan2(dy0, dx0);
    double dx1 = pts[1].first - xc;
    double dy1 = pts[1].second - yc;
    double phi1 = atan2(dy1, dx1);
    double dphi = phi1 - phi0;
    if(dphi < 0)
      phi += M_PI / 2.0;
    else
      phi -= M_PI / 2.0;
    // rotate track momentum vector (pz stays the same)
    double pt = track->get_pt();
    track->set_px(pt*cos(phi));
    track->set_py(pt*sin(phi));
    // set track position
    track->set_x(tpc_clusters[inner_index]->getX());
    track->set_y(tpc_clusters[inner_index]->getY());
    track->set_z(tpc_clusters[inner_index]->getZ());
  }
}

void  PHSimpleKFProp::line_fit(std::vector<std::pair<double,double>> points, double &a, double &b)
{
  // copied from: https://www.bragitoff.com
  // we want to fit z vs radius
  
   double xsum=0,x2sum=0,ysum=0,xysum=0;                //variables for sums/sigma of xi,yi,xi^2,xiyi etc
   for (unsigned int i=0; i<points.size(); ++i)
    {
      //double z = clusters[i]->getZ();
      //double r = sqrt(pow(clusters[i]->getX(),2) + pow(clusters[i]->getY(), 2));
      double r = points[i].first;
      double z = points[i].second;

      xsum=xsum+r;                        //calculate sigma(xi)
      ysum=ysum+z;                        //calculate sigma(yi)
      x2sum=x2sum+pow(r,2);                //calculate sigma(x^2i)
      xysum=xysum+r*z;                    //calculate sigma(xi*yi)
    }
   a=(points.size()*xysum-xsum*ysum)/(points.size()*x2sum-xsum*xsum);            //calculate slope
   b=(x2sum*ysum-xsum*xysum)/(x2sum*points.size()-xsum*xsum);            //calculate intercept

   if(Verbosity() > 10)
     {
       for (unsigned int i=0;i<points.size(); ++i)
	 {
	   double r = points[i].first;
	   double z_fit = a * r + b;                    //to calculate z(fitted) at given r points
	   std::cout << " r " << r << " z " << points[i].second << " z_fit " << z_fit << std::endl; 
	 } 
     }

    return;
}   

void  PHSimpleKFProp::line_fit_clusters(std::vector<TrkrCluster*> clusters, double &a, double &b)
{
  std::vector<std::pair<double,double>> points;
  
   for (unsigned int i=0; i<clusters.size(); ++i)
     {
       double z = clusters[i]->getZ();
       double r = sqrt(pow(clusters[i]->getX(),2) + pow(clusters[i]->getY(), 2));

       points.push_back(make_pair(r,z));
     }

   line_fit(points, a, b);

    return;
}

void PHSimpleKFProp::CircleFitByTaubin (std::vector<std::pair<double,double>> points, double &R, double &X0, double &Y0)
/*  
      Circle fit to a given set of data points (in 2D)
      This is an algebraic fit, due to Taubin, based on the journal article
      G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
                  Space Curves Defined By Implicit Equations, With 
                  Applications To Edge And Range Image Segmentation",
                  IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
*/
{
  int iter,IterMAX=99;
  
  double Mz,Mxy,Mxx,Myy,Mxz,Myz,Mzz,Cov_xy,Var_z;
  double A0,A1,A2,A22,A3,A33;
  double x,y;
  double DET,Xcenter,Ycenter;
  
  // Compute x- and y- sample means   
  double meanX = 0;
  double meanY = 0;
  double weight = 0;
  for(unsigned int i = 0; i < points.size(); ++i)
    {
      meanX += points[i].first;
      meanY += points[i].second;
      weight++;
    }
  meanX /= weight;
  meanY /= weight;

  //     computing moments 
  
  Mxx=Myy=Mxy=Mxz=Myz=Mzz=0.;
  
  for (unsigned int i=0; i<points.size(); i++)
    {
      double Xi = points[i].first - meanX;   //  centered x-coordinates
      double Yi = points[i].second - meanY;   //  centered y-coordinates
      double Zi = Xi*Xi + Yi*Yi;
      
      Mxy += Xi*Yi;
      Mxx += Xi*Xi;
      Myy += Yi*Yi;
      Mxz += Xi*Zi;
      Myz += Yi*Zi;
      Mzz += Zi*Zi;
    }
  Mxx /= weight;
  Myy /= weight;
  Mxy /= weight;
  Mxz /= weight;
  Myz /= weight;
  Mzz /= weight;
  
  //  computing coefficients of the characteristic polynomial
  
  Mz = Mxx + Myy;
  Cov_xy = Mxx*Myy - Mxy*Mxy;
  Var_z = Mzz - Mz*Mz;
  A3 = 4*Mz;
  A2 = -3*Mz*Mz - Mzz;
  A1 = Var_z*Mz + 4*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
  A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
  A22 = A2 + A2;
  A33 = A3 + A3 + A3;
  
  //    finding the root of the characteristic polynomial
  //    using Newton's method starting at x=0  
  //    (it is guaranteed to converge to the right root)
  
  for (x=0.,y=A0,iter=0; iter<IterMAX; iter++)  // usually, 4-6 iterations are enough
    {
      double Dy = A1 + x*(A22 + A33*x);
      double xnew = x - y/Dy;
      if ((xnew == x)||(!isfinite(xnew))) break;
      double ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3));
      if (fabs(ynew)>=fabs(y))  break;
      x = xnew;  y = ynew;
    }
  
  //  computing parameters of the fitting circle
  
  DET = x*x - x*Mz + Cov_xy;
  Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2;
  Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2;
  
  //  assembling the output
  
  X0 = Xcenter + meanX;
  Y0 = Ycenter + meanY;
  R = sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz);
}
