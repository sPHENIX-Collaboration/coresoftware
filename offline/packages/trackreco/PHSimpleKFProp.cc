/*!
 *  \file PHSimpleKFProp.cc
 *  \brief		kalman filter based propagator
 *  \author Michael Peters & Christof Roland
 */

#include "PHSimpleKFProp.h"

#include "ALICEKF.h"
#include "nanoflann.hpp"
#include "GPUTPCTrackParam.h"
#include "GPUTPCTrackLinearisation.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <phfield/PHField.h>
#include <phfield/PHFieldUtility.h>

#include <phool/getClass.h>
#include <phool/phool.h>                       // for PHWHERE

// tpc distortion correction
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterIterationMapv1.h>

#include <Geant4/G4SystemOfUnits.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

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

// anonymous namespace for local functions
namespace
{
  // square
  template<class T> inline constexpr T square( const T& x ) { return x*x; }


  void CircleFitByTaubin ( const std::vector<Acts::Vector3>& points, double &R, double &X0, double &Y0)
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
    for( const auto& point:points )
    {
      meanX += point(0);
      meanY += point(1);
      weight++;
    }
    meanX /= weight;
    meanY /= weight;
    
    //     computing moments 
    
    Mxx=Myy=Mxy=Mxz=Myz=Mzz=0.;
    for( const auto& point:points )
    {
      double Xi = point(0) - meanX;   //  centered x-coordinates
      double Yi = point(1) - meanY;   //  centered y-coordinates
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
      if ((xnew == x)||(!std::isfinite(xnew))) break;
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
    R = std::sqrt(square(Xcenter) + square(Ycenter));
  }

//   void findRoot(const double R, const double X0, const double Y0, double& x, double& y)
//   {
//     /**
//     * We need to determine the closest point on the circle to the origin
//     * since we can't assume that the track originates from the origin
//     * The eqn for the circle is (x-X0)^2+(y-Y0)^2=R^2 and we want to
//     * minimize d = sqrt((0-x)^2+(0-y)^2), the distance between the
//     * origin and some (currently, unknown) point on the circle x,y.
//     *
//     * Solving the circle eqn for x and substituting into d gives an eqn for
//     * y. Taking the derivative and setting equal to 0 gives the following
//     * two solutions. We take the smaller solution as the correct one, as
//     * usually one solution is wildly incorrect (e.g. 1000 cm)
//     */
// 
//     double miny = (sqrt(pow(X0, 2) * pow(R, 2) * pow(Y0, 2) + pow(R, 2)
//       * pow(Y0,4)) + pow(X0,2) * Y0 + pow(Y0, 3))
//       / (pow(X0, 2) + pow(Y0, 2));
// 
//     double miny2 = (-sqrt(pow(X0, 2) * pow(R, 2) * pow(Y0, 2) + pow(R, 2)
//       * pow(Y0,4)) + pow(X0,2) * Y0 + pow(Y0, 3))
//       / (pow(X0, 2) + pow(Y0, 2));
// 
//     double minx = std::sqrt(square(R) - square(miny - Y0)) + X0;
//     double minx2 = -std::sqrt(square(R) - square(miny2 - Y0)) + X0;
// 
//     /// Figure out which of the two roots is actually closer to the origin
//     if(fabs(minx) < fabs(minx2))
//       x = minx;
//     else
//       x = minx2;
// 
//     if(fabs(miny) < fabs(miny2))
//       y = miny;
//     else
//       y = miny2;
// 
//   }

}

using keylist = std::vector<TrkrDefs::cluskey>;

PHSimpleKFProp::PHSimpleKFProp(const std::string& name)
  : SubsysReco(name)
{}

int PHSimpleKFProp::End(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSimpleKFProp::InitRun(PHCompositeNode* topNode)
{
  
  int ret = get_nodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  fitter = std::make_unique<ALICEKF>(topNode,_cluster_map,_fieldDir,
				     _min_clusters_per_track,_max_sin_phi,Verbosity());
  fitter->useConstBField(_use_const_field);
  fitter->useFixedClusterError(_use_fixed_clus_err);
  fitter->setFixedClusterError(0,_fixed_clus_err.at(0));
  fitter->setFixedClusterError(1,_fixed_clus_err.at(1));
  fitter->setFixedClusterError(2,_fixed_clus_err.at(2));
  _field_map = PHFieldUtility::GetFieldMapNode(nullptr,topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

double PHSimpleKFProp::get_Bz(double x, double y, double z) const
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
    std::cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!_track_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find TrackSeedContainer" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if(!_hitsets)
    {
      std::cerr << PHWHERE << "No hitset container on node tree. Bailing."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  PHG4CylinderCellGeomContainer *geom_container =
      findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cerr << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  _surfmaps = findNode::getClass<ActsSurfaceMaps>(topNode, "ActsSurfaceMaps");
  if(!_surfmaps)
    {
      std::cout << "No Acts surface maps, exiting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  _tgeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  if(!_tgeometry)
    {
      std::cout << "No Acts tracking geometry, exiting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
   
  // tpc distortion correction
  m_dcc = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainer");
  if( m_dcc )
  { std::cout << "PHSimpleKFProp::InitRun - found TPC distortion correction container" << std::endl; }


  for(int i=7;i<=54;i++)
  {
    radii.push_back(geom_container->GetLayerCellGeom(i)->get_radius());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSimpleKFProp::process_event(PHCompositeNode* topNode)
{
  if(_n_iteration!=0){
    _iteration_map = findNode::getClass<TrkrClusterIterationMapv1>(topNode, "CLUSTER_ITERATION_MAP");
    if (!_iteration_map){
      std::cerr << PHWHERE << "Cluster Iteration Map missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  
  if(Verbosity()>0) std::cout << "starting Process" << std::endl;
  const auto globalPositions = PrepareKDTrees();

  std::vector<std::vector<TrkrDefs::cluskey>> keylist = getKeyList();
  std::vector<GPUTPCTrackParam> aliceSeeds = fitter->ALICEKalmanFilter(keylist, true, globalPositions);

  if(Verbosity()>0) std::cout << "prepared KD trees" << std::endl;
  
  if(Verbosity()>0) std::cout << "moved tracks into TPC" << std::endl;
  std::vector<std::vector<TrkrDefs::cluskey>> new_chains;
  std::vector<TrackSeed> unused_tracks;
  unsigned int trackid = 0;
  for(TrackSeedContainer::Iter track_it = _track_map->begin(); 
      track_it != _track_map->end(); ++track_it )
  {
    // if not a TPC track, ignore
    TrackSeed* track = *track_it;
    const bool is_tpc = std::any_of(
      track->begin_cluster_keys(),
      track->end_cluster_keys(),
      []( const TrkrDefs::cluskey& key ) { return TrkrDefs::getTrkrId(key) == TrkrDefs::tpcId; } );

    if(is_tpc)
    {
      MoveToFirstTPCCluster(track, globalPositions);
      GPUTPCTrackParam kftrack = aliceSeeds.at(trackid);
      std::cout << "propagating track id " << trackid << std::endl;
      std::cout << "with track params : " <<kftrack.GetX() << ", " 
	      << kftrack.GetY() << ", " << kftrack.GetZ() << ", " 
	      << kftrack.GetQPt() << ", " << kftrack.GetSinPhi() 
	      << ", " << kftrack.GetDzDs() << std::endl;
      if(Verbosity()>0) std::cout << "is tpc track" << std::endl;
      new_chains.push_back(PropagateTrack(track, kftrack, globalPositions));
    }
    else
    {
      // this is bad: it copies the track to its base class, which is essentially nothing
      if(Verbosity()>0) std::cout << "is NOT tpc track" << std::endl;
      unused_tracks.push_back(*track);
    }

    trackid++;
  }
  std::cout << "Finished tpc track loop"<<std::endl;
  _track_map->Reset();
  std::vector<std::vector<TrkrDefs::cluskey>> clean_chains = RemoveBadClusters(new_chains, globalPositions); 
  std::vector<GPUTPCTrackParam> ptracks = fitter->ALICEKalmanFilter(clean_chains,true, globalPositions);
  publishSeeds(clean_chains, ptracks, globalPositions);
  publishSeeds(unused_tracks);
  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::Vector3 PHSimpleKFProp::getGlobalPosition( TrkrCluster* cluster ) const
{
  // get global position from Acts transform
  auto globalpos = m_transform.getGlobalPosition(cluster,
    _surfmaps,
    _tgeometry);

  // check if TPC distortion correction are in place and apply
  if( m_dcc ) { globalpos = m_distortionCorrection.get_corrected_position( globalpos, m_dcc ); }

  return globalpos;
}

PositionMap PHSimpleKFProp::PrepareKDTrees()
{
  PositionMap globalPositions;
  //***** convert clusters to kdhits, and divide by layer
  std::vector<std::vector<std::vector<double> > > kdhits;
  kdhits.resize(58);
  if (!_cluster_map)
  {
    std::cout << "WARNING: (tracking.PHTpcTrackerUtil.convert_clusters_to_hits) cluster map is not provided" << std::endl;
    return globalPositions;
  }

  auto hitsetrange = _hitsets->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (auto hitsetitr = hitsetrange.first; hitsetitr != hitsetrange.second; ++hitsetitr)
  {
    auto range = _cluster_map->getClusters(hitsetitr->first);
    for (TrkrClusterContainer::ConstIterator it = range.first; it != range.second; ++it)
    {
      TrkrDefs::cluskey cluskey = it->first;
      TrkrCluster* cluster = it->second;
      if(!cluster) continue;
      if(_n_iteration!=0){
        if(_iteration_map != NULL ){
          //	  std::cout << "map exists entries: " << _iteration_map->size() << std::endl;
          if(_iteration_map->getIteration(cluskey)>0){ 
            //std::cout << "hit used, continue" << std::endl;
            continue; // skip hits used in a previous iteration
          }
        }
      }

      const Acts::Vector3 globalpos_d = getGlobalPosition(cluster);
      const Acts::Vector3 globalpos = { (float) globalpos_d.x(), (float) globalpos_d.y(), (float) globalpos_d.z()};
      globalPositions.insert(std::make_pair(cluskey, globalpos));

      int layer = TrkrDefs::getLayer(cluskey);
      std::vector<double> kdhit(4);
      kdhit[0] = globalpos_d.x();
      kdhit[1] = globalpos_d.y();
      kdhit[2] = globalpos_d.z();
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
    if(Verbosity()>0) std::cout << "l: " << l << std::endl;
    _ptclouds[l] = std::make_shared<KDPointCloud<double>>();
    _ptclouds[l]->pts.resize(kdhits[l].size());
    if(Verbosity()>0) std::cout << "resized to " << kdhits[l].size() << std::endl;
    for(size_t i=0;i<kdhits[l].size();++i)
    {
      _ptclouds[l]->pts[i] = kdhits[l][i];
    }
    _kdtrees[l] = std::make_shared<nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, KDPointCloud<double>>, KDPointCloud<double>, 3>>(3,*(_ptclouds[l]),nanoflann::KDTreeSingleIndexAdaptorParams(10));
    _kdtrees[l]->buildIndex();
  }

  return globalPositions;
}

void PHSimpleKFProp::MoveToFirstTPCCluster( const TrackSeed* track, const PositionMap& globalPositions )
{
  double track_x = track->get_x();
  double track_y = track->get_y();
  if(sqrt(track_x*track_x+track_y*track_y)>10.)
    {
      if(Verbosity()>0) std::cout << "WARNING: attempting to move track to TPC which is already in TPC! Aborting for this track." << std::endl;
      std::cout << "Track x,y " << track_x << ", " << track_y << std::endl;
      continue;
    }
  
  // get cluster keys
  std::vector<TrkrDefs::cluskey> ckeys;
  std::copy(track->begin_cluster_keys(),track->end_cluster_keys(),std::back_inserter(ckeys));
  
  std::vector<Acts::Vector3> trkGlobPos;
  for(const auto& ckey : ckeys)
    {
      if(TrkrDefs::getTrkrId(ckey) == TrkrDefs::tpcId )
        { trkGlobPos.push_back(globalPositions.at(ckey)); }
    }
  
  // get circle fit for TPC clusters plus vertex
  //double R = 1. / fabs(track->get_qOverR());
  float xc = track->get_X0();
  float yc = track->get_Y0();
  
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
  double cluster_x = trkGlobPos.at(inner_index)(0);
  double cluster_y = trkGlobPos.at(inner_index)(1);
  double dy = cluster_y-yc;
  double dx = cluster_x-xc;
  double phi = atan2(dy,dx);
  double dx0 = trkGlobPos.at(0)(0) - xc;     
  double dy0 = trkGlobPos.at(0)(1) - yc;
  double phi0 = atan2(dy0, dx0);
  double dx1 = trkGlobPos.at(1)(0) - xc;
  double dy1 = trkGlobPos.at(1)(1) - yc;
  double phi1 = atan2(dy1, dx1);
  double dphi = phi1 - phi0;
  if(dphi < 0)
    phi += M_PI / 2.0;
  else
    phi -= M_PI / 2.0;
  // rotate track momentum vector (pz stays the same)
  //double pt = track->get_pt();
  //track->set_px(pt*cos(phi));
  //track->set_py(pt*sin(phi));
  // set track position
  //track->set_x(trkGlobPos.at(0)(0));
  //track->set_y(trkGlobPos.at(0)(1));
  //track->set_z(trkGlobPos.at(0)(2));
  

}

std::vector<TrkrDefs::cluskey> PHSimpleKFProp::PropagateTrack(TrackSeed* track, GPUTPCTrackParam& kftrack, const PositionMap& globalPositions) const
{
  // extract cluster list
  std::vector<TrkrDefs::cluskey> ckeys;
  std::copy(track->begin_cluster_keys(),track->end_cluster_keys(),std::back_inserter(ckeys));
  if(ckeys.size()>1 && ((int)TrkrDefs::getLayer(ckeys.front()))>((int)TrkrDefs::getLayer(ckeys.back())))
  {
    std::reverse(ckeys.begin(),ckeys.end());
  } 
  // get track parameters
  
  float track_phi = atan2(track->get_py(),track->get_px());
 
  // Y = y
  // Z = z
  // SinPhi = py/sqrt(px^2+py^2)
  // DzDs = pz/sqrt(px^2+py^2)
  // QPt = 1/sqrt(px^2+py^2)
  
  std::vector<TrkrDefs::cluskey> propagated_track;

  if(Verbosity()>0)
  {
    std::cout << "initial track params:" << std::endl;
    std::cout << "X: " << kftrack.GetX() << std::endl;
    std::cout << "Y: " << kftrack.GetY() << std::endl;
    std::cout << "Z: " << kftrack.GetZ() << std::endl;
    std::cout << "SinPhi: " << kftrack.GetSinPhi() << std::endl;
    std::cout << "DzDs: " << kftrack.GetDzDs() << std::endl;
    std::cout << "QPt: " << kftrack.GetQPt() << std::endl;  
  }

  GPUTPCTrackLinearisation kfline(kftrack);

  // get layer for each cluster
  std::vector<unsigned int> layers;
  if(Verbosity()>0) std::cout << "cluster layers:" << std::endl;
  std::transform( ckeys.begin(), ckeys.end(), std::back_inserter( layers ), []( const TrkrDefs::cluskey& key ) { return TrkrDefs::getLayer(key); } );

  double old_phi = track_phi;
  std::cout << "old phi " << old_phi << std::endl;
  unsigned int old_layer = TrkrDefs::getLayer(ckeys[0]);
  if(Verbosity()>0) std::cout << "first layer: " << old_layer << std::endl;

  propagated_track.push_back(ckeys[0]);
  // first, propagate downward
  for(unsigned int l=old_layer+1;l<=54;l++)
  {
    if(std::isnan(kftrack.GetX()) ||
       std::isnan(kftrack.GetY()) ||
       std::isnan(kftrack.GetZ())) continue;
    if(fabs(kftrack.GetZ())>105.) continue;
    if(Verbosity()>0) std::cout << "\nlayer " << l << ":" << std::endl;
    // check to see whether layer is already occupied by at least one cluster
    // choosing the last one first (clusters organized from inside out)
    bool layer_filled = false;
    TrkrDefs::cluskey next_ckey = 0;
    for(int k=layers.size()-1; k>=0; k--)
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
      if(Verbosity()>0) std::cout << "layer is filled" << std::endl;
      TrkrCluster* nc = _cluster_map->findCluster(next_ckey);
      auto globalpos = globalPositions.at(next_ckey);
      double cx = globalpos(0);
      double cy = globalpos(1);
      double cz = globalpos(2);
      double cphi = atan2(cy,cx);
      double cxerr = sqrt(fitter->getClusterError(nc,globalpos,0,0));
      double cyerr = sqrt(fitter->getClusterError(nc,globalpos,1,1));
      double czerr = sqrt(fitter->getClusterError(nc,globalpos,2,2));
      double alpha = cphi-old_phi;
      double tX = kftrack.GetX();
      double tY = kftrack.GetY();
      double tx = tX*cos(old_phi)-tY*sin(old_phi);
      double ty = tX*sin(old_phi)+tY*cos(old_phi);
      double tz = kftrack.GetZ();
//      GPUTPCTrackParam::GPUTPCTrackFitParam fp;
//      kftrack.CalculateFitParameters(fp);
      if(Verbosity()>0) std::cout << "track position: (" << tx << ", " << ty << ", " << tz << ")" << std::endl;
      kftrack.Rotate(alpha,kfline,10.);
      kftrack.TransportToX(cx*cos(cphi)+cy*sin(cphi),kfline,_Bzconst*get_Bz(tx,ty,tz),10.);
      if(std::isnan(kftrack.GetX()) ||
       std::isnan(kftrack.GetY()) ||
       std::isnan(kftrack.GetZ())) continue;
      tX = kftrack.GetX();
      tY = kftrack.GetY();
      tx = tX*cos(cphi)-tY*sin(cphi);
      ty = tX*sin(cphi)+tY*cos(cphi);
      tz = kftrack.GetZ();
      double tYerr = sqrt(kftrack.GetCov(0));
      double tzerr = sqrt(kftrack.GetCov(5));
      double txerr = fabs(tYerr*sin(cphi));
      double tyerr = fabs(tYerr*cos(cphi));
      if(Verbosity()>0) std::cout << "cluster position: (" << cx << ", " << cy << ", " << cz << ")" << std::endl;
      if(Verbosity()>0) std::cout << "cluster position errors: (" << cxerr << ", " << cyerr << ", " << czerr << ")" << std::endl;
      if(Verbosity()>0) std::cout << "new track position: (" << kftrack.GetX()*cos(cphi)-kftrack.GetY()*sin(cphi) << ", " << kftrack.GetX()*sin(cphi)+kftrack.GetY()*cos(cphi) << ", " << kftrack.GetZ() << ")" << std::endl;
      if(Verbosity()>0) std::cout << "track position errors: (" << txerr << ", " << tyerr << ", " << tzerr << ")" << std::endl;
      if(Verbosity()>0) std::cout << "distance: " << sqrt(square(kftrack.GetX()*cos(cphi)-kftrack.GetY()*sin(cphi)-cx)+square(kftrack.GetX()*sin(cphi)+kftrack.GetY()*cos(cphi)-cy)+square(kftrack.GetZ()-cz)) << std::endl;
      if(fabs(tx-cx)<_max_dist*sqrt(txerr*txerr+cxerr*cxerr) &&
         fabs(ty-cy)<_max_dist*sqrt(tyerr*tyerr+cyerr*cyerr) &&
         fabs(tz-cz)<_max_dist*sqrt(tzerr*tzerr+czerr*czerr))
      {
        if(Verbosity()>0) std::cout << "Kept cluster" << std::endl;
        propagated_track.push_back(next_ckey);
//      kftrack.SetX(cx*cos(cphi)+cy*sin(cphi));
//      kftrack.SetY(-cx*sin(cphi)+cy*cos(cphi));
//      kftrack.SetZ(cz);
      }
      else
      {
        if(Verbosity()>0)
        {
          std::cout << "Rejected cluster" << std::endl;
          std::cout << "x: " << fabs(tx-cx) << " vs. " << _max_dist*sqrt(txerr*txerr+cxerr*cxerr) << std::endl;
          std::cout << "y: " << fabs(ty-cy) << " vs. " << _max_dist*sqrt(tyerr*tyerr+cyerr*cyerr) << std::endl;
          std::cout << "z: " << fabs(tz-cz) << " vs. " << _max_dist*sqrt(tzerr*tzerr+czerr*czerr) << std::endl;
        }
        kftrack.SetNDF(kftrack.GetNDF()-2);
        //ckeys.erase(std::remove(ckeys.begin(),ckeys.end(),next_ckey),ckeys.end());
      }
      old_phi = cphi;
    }
    // if layer is not occupied, search for the nearest available cluster to projected track position
    else
    {
      if(Verbosity()>0) std::cout << "layer not filled" << std::endl;
      // get current track coordinates to extract B field from map
      double tX = kftrack.GetX();
      double tY = kftrack.GetY();
      double tx = tX*cos(old_phi)-tY*sin(old_phi);
      double ty = tX*sin(old_phi)+tY*cos(old_phi);
      double tz = kftrack.GetZ();
//      GPUTPCTrackParam::GPUTPCTrackFitParam fp;
//      kftrack.CalculateFitParameters(fp);
      kftrack.TransportToX(radii[l-7],kfline,_Bzconst*get_Bz(tx,ty,tz),10.);
      if(std::isnan(kftrack.GetX()) ||
       std::isnan(kftrack.GetY()) ||
       std::isnan(kftrack.GetZ())) continue;
      // update track coordinates after transport
      tX = kftrack.GetX();
      tY = kftrack.GetY();
      double tYerr = sqrt(kftrack.GetCov(0));
      double tzerr = sqrt(kftrack.GetCov(5));
      tx = tX*cos(old_phi)-tY*sin(old_phi);
      ty = tX*sin(old_phi)+tY*cos(old_phi);
      tz = kftrack.GetZ();
      double query_pt[3] = {tx, ty, tz};

      if(m_dcc)
	{
	  // The distortion corrected cluster positions in globalPos are not at the layer radius
	  // We want to project to the radius appropriate for the globalPos values
	  // Get the distortion correction for the projection point, and calculate the radial increment

	  double proj_radius = sqrt(tx*tx+ty*ty);
	  if(proj_radius > 78.0 || abs(tz) > 105.5) continue;   // projection is bad, no cluster will be found

	  Acts::Vector3 proj_pt(tx,ty,tz);
	  if(Verbosity() > 2) 
	    std::cout << " call distortion correction for layer " << l  << " tx " << tx << " ty " << ty << " tz " << tz << " radius " << proj_radius << std::endl;
	  proj_pt = m_distortionCorrection.get_corrected_position( proj_pt, m_dcc ); 
	  // this point is meaningless, except that it gives us an estimate of the corrected radius of a point measured in this layer
	  double radius = sqrt(proj_pt[0]*proj_pt[0] + proj_pt[1]*proj_pt[1]);
	  // now project the track to that radius
	  if(Verbosity() > 2) 
	    std::cout << " call transport again for layer " << l  << " x " << proj_pt[0] << " y " << proj_pt[1] << " z " << proj_pt[2] 
		      << " radius " << radius << std::endl;
	  kftrack.TransportToX(radius,kfline,_Bzconst*get_Bz(tx,ty,tz),10.);
	  if(std::isnan(kftrack.GetX()) ||
	     std::isnan(kftrack.GetY()) ||
	     std::isnan(kftrack.GetZ())) continue;
	  tX = kftrack.GetX();
	  tY = kftrack.GetY();
	  tYerr = sqrt(kftrack.GetCov(0));
	  tzerr = sqrt(kftrack.GetCov(5));
	  tx = tX*cos(old_phi)-tY*sin(old_phi);
	  ty = tX*sin(old_phi)+tY*cos(old_phi);
	  tz = kftrack.GetZ();
	  query_pt[0] = tx;
	  query_pt[1] = ty;
	  query_pt[2] = tz; 
	}

      double txerr = fabs(tYerr*sin(old_phi));
      double tyerr = fabs(tYerr*cos(old_phi));
      if(Verbosity()>0) std::cout << "transported to " << radii[l-7] << "\n";
      if(Verbosity()>0) std::cout << "track position: (" << tx << ", " << ty << ", " << tz << ")" << std::endl;
      if(Verbosity()>0) std::cout << "track position error: (" << txerr << ", " << tyerr << ", " << tzerr << ")" << std::endl;

//      size_t ret_index;
//      double out_dist_sqr;
//      nanoflann::KNNResultSet<double> resultSet(1);
//      resultSet.init(&ret_index,&out_dist_sqr);
//      _kdtrees[l]->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));
      std::vector<long unsigned int> index_out(1);
      std::vector<double> distance_out(1);
      int n_results = _kdtrees[l]->knnSearch(&query_pt[0],1,&index_out[0],&distance_out[0]);
      if(Verbosity()>0) std::cout << "index_out: " << index_out[0] << std::endl;
      if(Verbosity()>0) std::cout << "squared_distance_out: " << distance_out[0] << std::endl;
      if(Verbosity()>0) std::cout << "solid_angle_dist: " << atan2(sqrt(distance_out[0]),radii[l-7]) << std::endl;
      if(n_results==0) continue;
      std::vector<double> point = _ptclouds[l]->pts[index_out[0]];
      TrkrDefs::cluskey closest_ckey = (*((int64_t*)&point[3]));
      TrkrCluster* cc = _cluster_map->findCluster(closest_ckey);
      auto ccglob = globalPositions.at(closest_ckey);
      double ccX = ccglob(0);
      double ccY = ccglob(1);
      double ccZ = ccglob(2);

      /*
      // alternatively:
      if(m_dcc)
	{
	  // The distortion corrected cluster positions in globalPos are not at the layer radius
	  // We want to project to the radius appropriate for the globalPos values
	  // Get the radius from the nearest associated cluster found above
	  Acts::Vector3 proj_pt(ccX, ccY, ccZ);
	  double radius = sqrt(proj_pt[0]*proj_pt[0] + proj_pt[1]*proj_pt[1]);

	  // now project the track to that radius
	  if(Verbosity() > 2) 
	    std::cout << " call transport again for layer " << l  << " x " << proj_pt[0] << " y " << proj_pt[1] << " z " << proj_pt[2] 
		      << " radius " << radius << std::endl;
	  kftrack.TransportToX(radius,kfline,_Bzconst*get_Bz(ccX,ccY,ccZ),10.);
	  if(std::isnan(kftrack.GetX()) ||
	     std::isnan(kftrack.GetY()) ||
	     std::isnan(kftrack.GetZ())) continue;
	  tX = kftrack.GetX();
	  tY = kftrack.GetY();
	  tYerr = sqrt(kftrack.GetCov(0));
	  tzerr = sqrt(kftrack.GetCov(5));
	  tx = tX*cos(old_phi)-tY*sin(old_phi);
	  ty = tX*sin(old_phi)+tY*cos(old_phi);
	  tz = kftrack.GetZ();
	  query_pt[0] = tx;
	  query_pt[1] = ty;
	  query_pt[2] = tz; 

	  n_results = _kdtrees[l]->knnSearch(&query_pt[0],1,&index_out[0],&distance_out[0]);
	  if(Verbosity()>0) std::cout << "index_out: " << index_out[0] << std::endl;
	  if(Verbosity()>0) std::cout << "squared_distance_out: " << distance_out[0] << std::endl;
	  if(Verbosity()>0) std::cout << "solid_angle_dist: " << atan2(sqrt(distance_out[0]),radii[l-7]) << std::endl;
	  if(n_results==0) continue;
	  point = _ptclouds[l]->pts[index_out[0]];
	  closest_ckey = (*((int64_t*)&point[3]));
	  cc = _cluster_map->findCluster(closest_ckey);
	  ccglob = globalPositions.at(closest_ckey);
	  ccX = ccglob(0);
	  ccY = ccglob(1);
	  ccZ = ccglob(2);
	}
      */      

      double cxerr = sqrt(fitter->getClusterError(cc,ccglob,0,0));
      double cyerr = sqrt(fitter->getClusterError(cc,ccglob,1,1));
      double czerr = sqrt(fitter->getClusterError(cc,ccglob,2,2));
      double ccphi = atan2(ccY,ccX);
      if(Verbosity()>0) std::cout << "cluster position: (" << ccX << ", " << ccY << ", " << ccZ << ")" << std::endl;
      if(Verbosity()>0) std::cout << "cluster position error: (" << cxerr << ", " << cyerr << ", " << czerr << ")" << std::endl;
      if(Verbosity()>0) std::cout << "cluster X: " << ccX*cos(ccphi)+ccY*sin(ccphi) << std::endl;
      if(fabs(tx-ccX)<_max_dist*sqrt(txerr*txerr+cxerr*cxerr) &&
         fabs(ty-ccY)<_max_dist*sqrt(tyerr*tyerr+cyerr*cyerr) &&
         fabs(tz-ccZ)<_max_dist*sqrt(tzerr*tzerr+czerr*czerr))
      {
	std::cout << "pushing back"<<std::endl;
        propagated_track.push_back(closest_ckey);
        layers.push_back(TrkrDefs::getLayer(closest_ckey));
	std::cout << "pushed back"<<std::endl;
/*        TrkrCluster* cc = _cluster_map->findCluster(closest_ckey);
        double ccX = cc->getX();
        std::cout << "cluster X: " << ccX << std::endl;
        double ccY = cc->getY();
        double ccx = ccX*cos(old_phi)-ccY*sin(old_phi);
        double ccy = ccX*sin(old_phi)+ccY*cos(old_phi);
        std::cout << "cluster position: (" << ccx << ", " << ccy << ", " << cc->getZ() << ")" << std::endl;
*/ 
        
        double alpha = ccphi-old_phi;
	std::cout << "rotating kftrack"<<std::endl;
        kftrack.Rotate(alpha,kfline,10.);
//        kftrack.SetX(ccX*cos(ccphi)+ccY*sin(ccphi));
//        kftrack.SetY(-ccX*sin(ccphi)+ccY*cos(ccphi));
//        kftrack.SetZ(cc->getZ());
	std::cout << "get cluserr"<<std::endl;
        double ccaY = -ccX*sin(ccphi)+ccY*cos(ccphi);
        double ccerrY = fitter->getClusterError(cc,ccglob,0,0)*sin(ccphi)*sin(ccphi)+fitter->getClusterError(cc,ccglob,0,1)*sin(ccphi)*cos(ccphi)+fitter->getClusterError(cc,ccglob,1,1)*cos(ccphi)*cos(ccphi);
        double ccerrZ = fitter->getClusterError(cc,ccglob,2,2);
	std::cout << "kftrack filter"<<std::endl;
        kftrack.Filter(ccaY,ccZ,ccerrY,ccerrZ,_max_sin_phi);
        if(Verbosity()>0) std::cout << "added cluster" << std::endl;
        old_phi = ccphi;
      }
    }
    old_layer = l;
  }
//  old_layer = TrkrDefs::getLayer(ckeys[0]);
//  std::reverse(ckeys.begin(),ckeys.end());
  layers.clear();
  if(Verbosity()>0) std::cout << "\nlayers after outward propagation:" << std::endl;
  for(int i=0;i<propagated_track.size();i++)
  {
    layers.push_back(TrkrDefs::getLayer(propagated_track[i]));
    if(Verbosity()>0) std::cout << layers[i] << std::endl;
  }
  // then, propagate upward
  for(unsigned int l=old_layer-1;l>=7;l--)
  {
    if(std::isnan(kftrack.GetX()) ||
       std::isnan(kftrack.GetY()) ||
       std::isnan(kftrack.GetZ())) continue;
    if(Verbosity()>0) std::cout << "\nlayer " << l << ":" << std::endl;
    // check to see whether layer is already occupied by at least one cluster
    // choosing the first one first (clusters organized from outside in)
    bool layer_filled = false;
    TrkrDefs::cluskey next_ckey = 0;
    for(size_t k=0; k<layers.size(); k++)
    {
      if(layer_filled) continue;
      if(layers[k]==l)
      {
        layer_filled = true;
        next_ckey = propagated_track[k];
      }
    }
    // if layer is already occupied, reset track parameters to last cluster in layer
    if(layer_filled)
    {
      if(Verbosity()>0) std::cout << "layer is filled" << std::endl;
      auto ncglob = globalPositions.at(next_ckey);
      double cx = ncglob(0);
      double cy = ncglob(1);
      double cz = ncglob(2);
      double cphi = atan2(cy,cx);
      double alpha = cphi-old_phi;
      double tX = kftrack.GetX();
      double tY = kftrack.GetY();
      double tx = tX*cos(old_phi)-tY*sin(old_phi);
      double ty = tX*sin(old_phi)+tY*cos(old_phi);
      double tz = kftrack.GetZ();
//      GPUTPCTrackParam::GPUTPCTrackFitParam fp;
//      kftrack_up.CalculateFitParameters(fp);
      kftrack.Rotate(alpha,kfline,10.);
      kftrack.TransportToX(cx*cos(cphi)+cy*sin(cphi),kfline,_Bzconst*get_Bz(tx,ty,tz),10.);
      if(std::isnan(kftrack.GetX()) ||
       std::isnan(kftrack.GetY()) ||
       std::isnan(kftrack.GetZ())) continue;
      if(Verbosity()>0) std::cout << "transported to " << radii[l-7] << "\n";
      if(Verbosity()>0) std::cout << "track position: (" << kftrack.GetX()*cos(cphi)-kftrack.GetY()*sin(cphi) << ", " << kftrack.GetX()*sin(cphi)+kftrack.GetY()*cos(cphi) << ", " << kftrack.GetZ() << ")" << std::endl;
      if(Verbosity()>0) std::cout << "cluster position: (" << cx << ", " << cy << ", " << cz << ")" << std::endl;
//      kftrack.SetX(cx*cos(cphi)+cy*sin(cphi));
//      kftrack.SetY(-cx*sin(cphi)+cy*cos(cphi));
//      kftrack.SetZ(cz);
//      propagated_track.push_back(next_ckey);
      old_phi = cphi;
    }
    // if layer is not occupied, search for the nearest available cluster to projected track position
    else
    {
      if(Verbosity()>0) std::cout << "layer not filled" << std::endl;
      double tX = kftrack.GetX();
      double tY = kftrack.GetY();
      double tx = tX*cos(old_phi)-tY*sin(old_phi);
      double ty = tX*sin(old_phi)+tY*cos(old_phi);
      double tz = kftrack.GetZ();
//      GPUTPCTrackParam::GPUTPCTrackFitParam fp;
//      kftrack_up.CalculateFitParameters(fp);
      kftrack.TransportToX(radii[l-7],kfline,_Bzconst*get_Bz(tx,ty,tz),10.);
      if(std::isnan(kftrack.GetX()) ||
       std::isnan(kftrack.GetY()) ||
       std::isnan(kftrack.GetZ())) continue;
      tX = kftrack.GetX();
      tY = kftrack.GetY();
      double tYerr = sqrt(kftrack.GetCov(0));
      double tzerr = sqrt(kftrack.GetCov(5));
      tx = tX*cos(old_phi)-tY*sin(old_phi);
      ty = tX*sin(old_phi)+tY*cos(old_phi);
      tz = kftrack.GetZ();
      double query_pt[3] = {tx, ty, tz};

      // Now look for the nearest cluster to this projection point (tx,ty,tz), which is at the nominal layer radius
      if(m_dcc)
	{
	  // The distortion corrected cluster positions in globalPos are not at the layer radius
	  // We want to project to the radius appropriate for the globalPos values
	  // Get the distortion correction for the projection point, and calculate the radial increment
	  double proj_radius = sqrt(tx*tx+ty*ty);
	  if(proj_radius > 78.0 || abs(tz) > 105.5) continue;   // projection is bad, no cluster will be found

	  Acts::Vector3 proj_pt(tx,ty,tz);
	  if(Verbosity() > 2)
	    std::cout << " call distortion correction for layer " << l  << " tx " << tx << " ty " << ty << " tz " << tz << " radius " << proj_radius << std::endl;
	  proj_pt = m_distortionCorrection.get_corrected_position( proj_pt, m_dcc ); 
	  // this point is meaningless, except that it givs us an estimate of the corrected radius of a point measured in this layer
	  double radius = sqrt(proj_pt[0]*proj_pt[0] + proj_pt[1]*proj_pt[1]);

	  // now project the track to that radius
	  if(Verbosity() > 2)
	    std::cout << " call transport again for layer " << l  << " x " << proj_pt[0] << " y " << proj_pt[1] << " z " << proj_pt[2] 
		      << " radius " << radius << std::endl;
	  kftrack.TransportToX(radius,kfline,_Bzconst*get_Bz(tx,ty,tz),10.);
	  if(std::isnan(kftrack.GetX()) ||
	     std::isnan(kftrack.GetY()) ||
	     std::isnan(kftrack.GetZ())) continue;
	  tX = kftrack.GetX();
	  tY = kftrack.GetY();
	  tYerr = sqrt(kftrack.GetCov(0));
	  tzerr = sqrt(kftrack.GetCov(5));
	  tx = tX*cos(old_phi)-tY*sin(old_phi);
	  ty = tX*sin(old_phi)+tY*cos(old_phi);
	  tz = kftrack.GetZ();
	  query_pt[0] = tx;
	  query_pt[1] = ty;
	  query_pt[2] = tz; 
	}

      double txerr = fabs(tYerr*sin(old_phi));
      double tyerr = fabs(tYerr*cos(old_phi));
      if(Verbosity()>0) std::cout << "transported to " << radii[l-7] << "\n";
      if(Verbosity()>0) std::cout << "track position: (" << kftrack.GetX()*cos(old_phi)-kftrack.GetY()*sin(old_phi) << ", " << kftrack.GetX()*sin(old_phi)+kftrack.GetY()*cos(old_phi) << ", " << kftrack.GetZ() << ")" << std::endl;
      if(Verbosity()>0) std::cout << "track position errors: (" << txerr << ", " << tyerr << ", " << tzerr << ")" << std::endl;

      std::vector<long unsigned int> index_out(1);
      std::vector<double> distance_out(1);
      int n_results = _kdtrees[l]->knnSearch(&query_pt[0],1,&index_out[0],&distance_out[0]);
      if(Verbosity()>0) std::cout << "index_out: " << index_out[0] << std::endl;
      if(Verbosity()>0) std::cout << "squared_distance_out: " << distance_out[0] << std::endl;
      if(Verbosity()>0) std::cout << "solid_angle_dist: " << atan2(sqrt(distance_out[0]),radii[l-7]) << std::endl;
      if(n_results==0) continue;
      std::vector<double> point = _ptclouds[l]->pts[index_out[0]];
      TrkrDefs::cluskey closest_ckey = (*((int64_t*)&point[3]));
      TrkrCluster* cc = _cluster_map->findCluster(closest_ckey);
      auto ccglob2 = globalPositions.at(closest_ckey);
      double ccX = ccglob2(0);
      double ccY = ccglob2(1);
      double ccZ = ccglob2(2);

      /*
      // alternatively:
      if(m_dcc)
	{
	  // The distortion corrected cluster positions in globalPos are not at the layer radius
	  // We want to project to the radius appropriate for the globalPos values
	  // Get the radius from the global position of the nearest cluster found above 
	  Acts::Vector3 proj_pt(ccX, ccY, ccZ);
	  double radius = sqrt(proj_pt[0]*proj_pt[0] + proj_pt[1]*proj_pt[1]);

	  // now project the track to that radius
	  if(Verbosity() > 2) 
	    std::cout << " call transport again for layer " << l  << " x " << proj_pt[0] << " y " << proj_pt[1] << " z " << proj_pt[2] 
		      << " radius " << radius << std::endl;
	  kftrack.TransportToX(radius,kfline,_Bzconst*get_Bz(ccX,ccY,ccZ),10.);
	  if(std::isnan(kftrack.GetX()) ||
	     std::isnan(kftrack.GetY()) ||
	     std::isnan(kftrack.GetZ())) continue;
	  tX = kftrack.GetX();
	  tY = kftrack.GetY();
	  tYerr = sqrt(kftrack.GetCov(0));
	  tzerr = sqrt(kftrack.GetCov(5));
	  tx = tX*cos(old_phi)-tY*sin(old_phi);
	  ty = tX*sin(old_phi)+tY*cos(old_phi);
	  tz = kftrack.GetZ();
	  query_pt[0] = tx;
	  query_pt[1] = ty;
	  query_pt[2] = tz; 

	  n_results = _kdtrees[l]->knnSearch(&query_pt[0],1,&index_out[0],&distance_out[0]);
	  if(Verbosity()>0) std::cout << "index_out: " << index_out[0] << std::endl;
	  if(Verbosity()>0) std::cout << "squared_distance_out: " << distance_out[0] << std::endl;
	  if(Verbosity()>0) std::cout << "solid_angle_dist: " << atan2(sqrt(distance_out[0]),radii[l-7]) << std::endl;
	  if(n_results==0) continue;
	  point = _ptclouds[l]->pts[index_out[0]];
	  closest_ckey = (*((int64_t*)&point[3]));
	  cc = _cluster_map->findCluster(closest_ckey);
	  ccglob2 = globalPositions.at(closest_ckey);
	  ccX = ccglob2(0);
	  ccY = ccglob2(1);
	  ccZ = ccglob2(2);
	}
      */

      double cxerr = sqrt(fitter->getClusterError(cc,ccglob2,0,0));
      double cyerr = sqrt(fitter->getClusterError(cc,ccglob2,1,1));
      double czerr = sqrt(fitter->getClusterError(cc,ccglob2,2,2));
      if(Verbosity()>0) std::cout << "cluster position: (" << ccX << ", " << ccY << ", " << ccZ << ")" << std::endl;
      double ccphi = atan2(ccY,ccX);
      if(Verbosity()>0) std::cout << "cluster position errors: (" << cxerr << ", " << cyerr << ", " << czerr << ")" << std::endl;
      if(Verbosity()>0) std::cout << "cluster X: " << ccX*cos(ccphi)+ccY*sin(ccphi) << std::endl;
      double alpha = ccphi-old_phi;
      if(fabs(tx-ccX)<_max_dist*sqrt(txerr*txerr+cxerr*cxerr) &&
         fabs(ty-ccY)<_max_dist*sqrt(tyerr*tyerr+cyerr*cyerr) &&
         fabs(tz-ccZ)<_max_dist*sqrt(tzerr*tzerr+czerr*czerr))
      {
        propagated_track.push_back(closest_ckey);
        layers.push_back(TrkrDefs::getLayer(closest_ckey));
/*        TrkrCluster* cc = _cluster_map->findCluster(closest_ckey);
        double ccX = cc->getX();
        std::cout << "cluster X: " << ccX << std::endl;
        double ccY = cc->getY();
        double ccx = ccX*cos(old_phi)-ccY*sin(old_phi);
        double ccy = ccX*sin(old_phi)+ccY*cos(old_phi);
        std::cout << "cluster position: (" << ccx << ", " << ccy << ", " << cc->getZ() << ")" << std::endl;
        double ccphi = atan2(ccy,ccx);
        double alpha = ccphi-old_phi;
*/
        kftrack.Rotate(alpha,kfline,10.);
        double ccaY = -ccX*sin(ccphi)+ccY*cos(ccphi);
        double ccerrY = fitter->getClusterError(cc,ccglob2,0,0)*sin(ccphi)*sin(ccphi)+fitter->getClusterError(cc,ccglob2,1,0)*sin(ccphi)*cos(ccphi)+fitter->getClusterError(cc,ccglob2,1,1)*cos(ccphi)*cos(ccphi);
        double ccerrZ = fitter->getClusterError(cc,ccglob2,2,2);
        kftrack.Filter(ccaY,ccZ,ccerrY,ccerrZ,_max_sin_phi);
//        kftrack.SetX(ccX*cos(ccphi)+ccY*sin(ccphi));
//        kftrack.SetY(-ccX*sin(ccphi)+ccY*cos(ccphi));
//        kftrack.SetZ(cc->getZ());
        if(Verbosity()>0) std::cout << "added cluster" << std::endl;
        old_phi = ccphi;
      }
    }
    old_layer = l;
  }
  std::sort(propagated_track.begin(),propagated_track.end(),
            [](TrkrDefs::cluskey a, TrkrDefs::cluskey b)
            {return (TrkrDefs::getLayer(a)<TrkrDefs::getLayer(b));});
  return propagated_track;
}

std::vector<keylist> PHSimpleKFProp::RemoveBadClusters(const std::vector<keylist>& chains, const PositionMap& globalPositions) const
{
  if(Verbosity()>0) std::cout << "removing bad clusters" << std::endl;
  std::vector<keylist> clean_chains;
  for(const keylist& chain : chains)
  {
    if(chain.size()<3) continue;
    keylist clean_chain;

    std::vector<std::pair<double,double>> xy_pts;
    std::vector<std::pair<double,double>> rz_pts;
    for(const TrkrDefs::cluskey& ckey : chain)
    {
      auto global = globalPositions.at(ckey);
      xy_pts.push_back(std::make_pair(global(0),global(1)));
      float r = sqrt(global(0)*global(0) + global(1)*global(1));
      rz_pts.push_back(std::make_pair(r,global(2)));
    }
    if(Verbosity()>0) std::cout << "chain size: " << chain.size() << std::endl;
    double A;
    double B;
    double R;
    double X0;
    double Y0;
    fitter->CircleFitByTaubin(xy_pts,R,X0,Y0);
    fitter->line_fit(rz_pts,A,B);
    std::vector<double> xy_resid = fitter->GetCircleClusterResiduals(xy_pts,R,X0,Y0);
    std::vector<double> rz_resid = fitter->GetLineClusterResiduals(rz_pts,A,B);
    for(size_t i=0;i<chain.size();i++)
    {
      if(xy_resid[i]>_xy_outlier_threshold) continue;
      clean_chain.push_back(chain[i]);
    }
    clean_chains.push_back(clean_chain);
    if(Verbosity()>0) std::cout << "pushed clean chain with " << clean_chain.size() << " clusters" << std::endl;
  }
  return clean_chains;
}

void PHSimpleKFProp::publishSeeds(const std::vector<std::vector<TrkrDefs::cluskey>>& chains, const std::vector<GPUTPCTrackParam>& seeds, const PositionMap& globalPositions)
{
  unsigned int trackid = 0;
  for(const auto& aliceSeed : seeds)
    {
      TrackSeed* sphenixSeed = _track_map->get(trackid);
      double track_pt = fabs(1./aliceSeed.GetQPt());
      double track_pterr = sqrt(aliceSeed.GetErr2QPt())/(aliceSeed.GetQPt()*aliceSeed.GetQPt());
      for(const auto ckey : chains.at(trackid))
	{
	  if(sphenixSeed->find_cluster_key(ckey) == sphenixSeed->end_cluster_keys())
	    { sphenixSeed->insert_cluster_key(ckey); }
	}
      // If Kalman filter doesn't do its job (happens often with short seeds), use the circle-fit estimate as the central value
      if(sphenixSeed->size_cluster_keys()<10) 
	{ track_pt = sphenixSeed->get_pt(); }

      if(fitter->checknan(track_pt, "pT", trackid)) continue;
      if(fitter->checknan(track_pterr,"pT err", trackid)) continue;

      auto lckey = *(sphenixSeed->end_cluster_keys());
      auto lcluster = _cluster_map->findCluster(lckey);
      const auto& lclusterglob = globalPositions.at(lckey);

      double track_phi = atan2(lclusterglob(1), lclusterglob(0));
      double track_x = aliceSeed.GetX()*cos(track_phi)-aliceSeed.GetY()*sin(track_phi);
      double track_y = aliceSeed.GetX()*sin(track_phi)+aliceSeed.GetY()*cos(track_phi);
      double track_z = aliceSeed.GetZ();
      if(fitter->checknan(track_z,"z",trackid)) continue;
      double track_zerr = sqrt(aliceSeed.GetErr2Z());
      if(fitter->checknan(track_zerr,"zerr",trackid)) continue;

    
      const float lclusterrad = sqrt(lclusterglob(0)*lclusterglob(0) + lclusterglob(1)*lclusterglob(1));
      double last_cluster_phierr = lcluster->getRPhiError() / lclusterrad;
      // phi error assuming error in track radial coordinate is zero
      double track_phierr = sqrt(pow(last_cluster_phierr,2)+(pow(aliceSeed.GetX(),2)*aliceSeed.GetErr2Y()) / 
				 pow(pow(aliceSeed.GetX(),2)+pow(aliceSeed.GetY(),2),2));
      if(fitter->checknan(track_phierr,"phierr",trackid)) continue;
      double track_curvature = aliceSeed.GetKappa(fitter->get_Bzconst()*fitter->get_Bz(track_x,track_y,track_z));
    if(fitter->checknan(track_curvature,"curvature",trackid)) continue;
    double track_curverr = sqrt(aliceSeed.GetErr2QPt())*fitter->get_Bzconst()*fitter->get_Bz(track_x,track_y,track_z);
    if(fitter->checknan(track_curverr,"curvature error",trackid)) continue;

    float track_pz = track_pt*aliceSeed.GetDzDs();
    float track_p = sqrt(track_pt*track_pt+track_pz*track_pz);
    float qoverPt = aliceSeed.GetQPt();
    float eta = atanh(track_pz / track_p);
    float theta = 2*atan(exp(-1*eta));

    sphenixSeed->set_Z0(track_z);
    // Update circle center first
    sphenixSeed->circleFitByTaubin(_cluster_map, _surfmaps, _tgeometry,0,58);
    // Then update with better pT estimate
    sphenixSeed->set_slope(1./tan(theta));
    sphenixSeed->set_qOverR(qoverPt * 100 / (0.3*1.4));

      trackid++;
    }

}

void PHSimpleKFProp::publishSeeds(const std::vector<TrackSeed>& seeds)
{
  for( const auto& seed:seeds )
  { _track_map->insert(&seed); }
}


std::vector<std::vector<TrkrDefs::cluskey>> PHSimpleKFProp::getKeyList()
{
  std::vector<std::vector<TrkrDefs::cluskey>> keylist;
  for(auto& seed : *_track_map)
    {
      std::vector<TrkrDefs::cluskey> dumvec;
      for(TrackSeed::ConstClusterKeyIter iter = seed->begin_cluster_keys();
	  iter != seed->end_cluster_keys();
	  ++iter)
	{
	  dumvec.push_back(*iter);
	}

      keylist.push_back(dumvec);
    }

  return keylist;
}
