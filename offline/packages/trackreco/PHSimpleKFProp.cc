/*!
 *  \file PHSimpleKFProp.cc
 *  \brief		kalman filter based propagator
 *  \author Michael Peters & Christof Roland
 */

#include "PHSimpleKFProp.h"
#include "PHGhostRejection.h"
#include "ALICEKF.h"
#include "nanoflann.hpp"
#include "GPUTPCTrackParam.h"
#include "GPUTPCTrackLinearisation.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <phfield/PHField.h>
#include <phfield/PHFieldUtility.h>
#include <phfield/PHFieldConfig.h>
#include <phfield/PHFieldConfigv1.h>

#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>                       // for PHWHERE

// tpc distortion correction
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/ActsGeometry.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed_v1.h>

#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
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
  PHFieldConfigv1 fcfg;
  fcfg.set_field_config(PHFieldConfig::FieldConfigTypes::Field3DCartesian);
  char *calibrationsroot = getenv("CALIBRATIONROOT");
  assert(calibrationsroot);
  auto magField = std::string(calibrationsroot) +
    std::string("/Field/Map/sphenix3dtrackingmapxyz.root"); 
  fcfg.set_filename(magField);
  //  fcfg.set_rescale(1);
  _field_map = std::unique_ptr<PHField>(PHFieldUtility::BuildFieldMap(&fcfg));

  fitter = std::make_unique<ALICEKF>(topNode,_cluster_map,_field_map.get(), _fieldDir,
				     _min_clusters_per_track,_max_sin_phi,Verbosity());
  fitter->useConstBField(_use_const_field);
  fitter->useFixedClusterError(_use_fixed_clus_err);
  fitter->set_cluster_version( m_cluster_version );
  std::cout << "simple kf cluster version " << m_cluster_version << std::endl;
  fitter->setFixedClusterError(0,_fixed_clus_err.at(0));
  fitter->setFixedClusterError(1,_fixed_clus_err.at(1));
  fitter->setFixedClusterError(2,_fixed_clus_err.at(2));
  //  _field_map = PHFieldUtility::GetFieldMapNode(nullptr,topNode);
  // m_Cache = magField->makeCache(m_tGeometry->magFieldContext);
  return Fun4AllReturnCodes::EVENT_OK;
}

double PHSimpleKFProp::get_Bz(double x, double y, double z) const
{
  if(_use_const_field) return 1.4;
  double p[4] = {x*cm,y*cm,z*cm,0.*cm};
  double bfield[3];
  _field_map->GetFieldValue(p,bfield);
  /*  Acts::Vector3 loc(0,0,0);
  int mfex = (magField != nullptr);
  int tgex = (m_tGeometry != nullptr);
  std::cout << " getting acts field " << mfex << " " << tgex << std::endl;
  auto Cache =  m_tGeometry->magField->makeCache(m_tGeometry->magFieldContext);
  auto bf = m_tGeometry->magField->getField(loc,Cache);
  if(bf.ok())
    {
      Acts::Vector3 val = bf.value();
      std::cout << "bz big: " <<  bfield[2]/tesla << " bz acts: " << val(2) << std::endl;
    }
  */
  return bfield[2]/tesla;
}

int PHSimpleKFProp::get_nodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  _tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if(!_tgeometry)
    {
      std::cout << "No Acts tracking geometry, exiting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
   
  // tpc distortion correction
  m_dcc = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainer");
  if( m_dcc )
  { std::cout << "PHSimpleKFProp::InitRun - found TPC distortion correction container" << std::endl; }

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
    std::cerr << PHWHERE << " ERROR: Can't find TrackSeedContainer " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  PHG4TpcCylinderGeomContainer *geom_container =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
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

int PHSimpleKFProp::process_event(PHCompositeNode* topNode)
{
  if(_n_iteration!=0){
    _iteration_map = findNode::getClass<TrkrClusterIterationMapv1>(topNode, "CLUSTER_ITERATION_MAP");
    if (!_iteration_map){
      std::cerr << PHWHERE << "Cluster Iteration Map missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  
  PHTimer timer("KFPropTimer");
  
  timer.stop();
  timer.restart();

  if(Verbosity()>0) std::cout << "starting Process" << std::endl;
  PositionMap globalPositions = PrepareKDTrees();
  if(Verbosity()>0) std::cout << "prepared KD trees" << std::endl;

  std::vector<std::vector<TrkrDefs::cluskey>> new_chains;
  std::vector<TrackSeed_v1> unused_tracks;
  for(int track_it = 0; track_it != _track_map->size(); ++track_it )
  {
    if(Verbosity()>0) std::cout << "TPC seed " << track_it << std::endl;
    // if not a TPC track, ignore
    TrackSeed* track = _track_map->get(track_it);
    const bool is_tpc = std::any_of(
      track->begin_cluster_keys(),
      track->end_cluster_keys(),
      []( const TrkrDefs::cluskey& key ) { return TrkrDefs::getTrkrId(key) == TrkrDefs::tpcId; } );

    if(is_tpc)
    {
      std::vector<std::vector<TrkrDefs::cluskey>> keylist;
      std::vector<TrkrDefs::cluskey> dumvec;
      std::map<TrkrDefs::cluskey, Acts::Vector3> trackClusPositions;
      for(TrackSeed::ConstClusterKeyIter iter = track->begin_cluster_keys();
	  iter != track->end_cluster_keys();
	  ++iter)
	{
	  dumvec.push_back(*iter);
	  auto pos = globalPositions.at(*iter);
	  trackClusPositions.insert(std::make_pair(*iter,pos));
	}

      /// Can't circle fit a seed with less than 3 clusters, skip it
      if(dumvec.size() < 3)
	{ continue; }

      keylist.push_back(dumvec);
    
      /// This will by definition return a single pair with each vector 
      /// in the pair length 1 corresponding to the seed info
      std::vector<float> trackChi2;
      timer.stop();
      timer.restart();
      
      auto seedpair = fitter->ALICEKalmanFilter(keylist, false, 
						trackClusPositions, trackChi2);

      timer.stop();
      if(Verbosity() > 3)
	{ std::cout << "single track ALICEKF time " << timer.elapsed()
		<< std::endl;
	}
      timer.restart();
      /// circle fit back to update track parameters
      track->circleFitByTaubin(trackClusPositions, 7, 55);
      track->lineFit(trackClusPositions, 7, 55);
      timer.stop();
      if(Verbosity() > 3)
	{ std::cout << "single track circle fit time " << timer.elapsed() << std::endl; }
      if(seedpair.first.size() == 0 || seedpair.second.size() == 0)
	{ continue; }
      if(Verbosity()>0) std::cout << "is tpc track" << std::endl;

      timer.stop();
      timer.restart();

      if(Verbosity()>0) std::cout << "propagate first round" << std::endl;
      auto preseed = PropagateTrack(track,seedpair.second.at(0),globalPositions);

      if(Verbosity()>0) std::cout << "preseed size " << preseed.size() << std::endl;

      if(preseed.size()>40){
        new_chains.push_back(preseed);
        continue;
      }

      std::vector<std::vector<TrkrDefs::cluskey>> kl;
      kl.push_back(preseed);

      if(Verbosity()>0) std::cout << "kl size " << kl.size() << std::endl;
      std::vector<float> pretrackChi2;
      auto prepair = fitter->ALICEKalmanFilter(kl,false,globalPositions,pretrackChi2);
      if(prepair.first.size()==0 || prepair.second.size()==0) continue;

      auto pretrack = prepair.first.at(0);
      std::vector<TrkrDefs::cluskey> dumvec2;
      std::map<TrkrDefs::cluskey, Acts::Vector3> pretrackClusPositions;
      for(TrackSeed::ConstClusterKeyIter iter = pretrack.begin_cluster_keys();
          iter != pretrack.end_cluster_keys();
          ++iter)
        {
          dumvec2.push_back(*iter);
          auto pos = globalPositions.at(*iter);
          pretrackClusPositions.insert(std::make_pair(*iter,pos));
        }

      pretrack.circleFitByTaubin(pretrackClusPositions, 7, 55);
      pretrack.lineFit(pretrackClusPositions, 7, 55);

      new_chains.push_back(PropagateTrack(&pretrack, prepair.second.at(0),
					  globalPositions));
      timer.stop();
      auto propagatetime = timer.elapsed();
      if(Verbosity() > 3)
	{ std::cout << "propagate track time " << propagatetime << std::endl; }
    }
    else
    {
      // this is bad: it copies the track to its base class, which is essentially nothing
      if(Verbosity()>0) std::cout << "is NOT tpc track" << std::endl;
      unused_tracks.push_back(*track);
    }
  }
  
  _track_map->Reset();
  timer.stop();
  timer.restart();

  std::vector<std::vector<TrkrDefs::cluskey>> clean_chains = RemoveBadClusters(new_chains, globalPositions); 
  timer.stop();

  timer.stop();
  timer.restart();
  std::vector<float> trackChi2;
  auto seeds = fitter->ALICEKalmanFilter(clean_chains, true, globalPositions,
					 trackChi2);
  timer.stop();
  auto alicekftime = timer.elapsed();
  if(Verbosity() > 0 )
    { std::cout << "full alice kf time all tracks " << alicekftime << std::endl; }
  timer.stop();
  timer.restart();
  publishSeeds(seeds.first, globalPositions);

  timer.stop();
  auto circlefittime = timer.elapsed();
  if(Verbosity() > 0)
    { std::cout << "circle fit all tracks time " << circlefittime << std::endl; }
  publishSeeds(unused_tracks);

  /// Remove tracks that are duplicates from the KFProp
  PHGhostRejection rejector(Verbosity());
  rejector.positionMap(globalPositions);
  rejector.trackSeedContainer(_track_map);
  timer.stop();
  timer.restart();
  rejector.rejectGhostTracks(trackChi2);
  timer.stop();
  if(Verbosity() > 2)
    { std::cout << "ghost rejection time " << timer.elapsed() << std::endl; }

  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::Vector3 PHSimpleKFProp::getGlobalPosition( TrkrDefs::cluskey key, TrkrCluster* cluster ) const
{
  // get global position from Acts transform
  auto globalpos = _tgeometry->getGlobalPosition(key, cluster);

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

  for(const auto& hitsetkey:_cluster_map->getHitSetKeys(TrkrDefs::TrkrId::tpcId))
  {
    auto range = _cluster_map->getClusters(hitsetkey);
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

      const Acts::Vector3 globalpos_d = getGlobalPosition( cluskey, cluster);
      const Acts::Vector3 globalpos = { (float) globalpos_d.x(), (float) globalpos_d.y(), (float) globalpos_d.z()};
      globalPositions.insert(std::make_pair(cluskey, globalpos));

      int layer = TrkrDefs::getLayer(cluskey);
      std::vector<double> kdhit(4);
      kdhit[0] = globalpos_d.x();
      kdhit[1] = globalpos_d.y();
      kdhit[2] = globalpos_d.z();
      uint64_t key = cluskey;
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


std::vector<TrkrDefs::cluskey> PHSimpleKFProp::PropagateTrack(TrackSeed* track, Eigen::Matrix<double,6,6>& xyzCov, const PositionMap& globalPositions) const
{
  // extract cluster list
 
  std::vector<TrkrDefs::cluskey> ckeys;
  std::copy(track->begin_cluster_keys(),track->end_cluster_keys(),std::back_inserter(ckeys));
 
  if(ckeys.size()>1 && ((int)TrkrDefs::getLayer(ckeys.front()))>((int)TrkrDefs::getLayer(ckeys.back())))
  {
    std::reverse(ckeys.begin(),ckeys.end());
  } 

  double track_x = track->get_x();
  double track_y = track->get_y();
  double track_z = track->get_z();

  if(Verbosity()>0) std::cout << "track (x,y,z) = (" << track_x << ", " << track_y << ", " << track_z << ")" << std::endl;

  double track_px = track->get_px(_cluster_map,_tgeometry);
  double track_py = track->get_py(_cluster_map,_tgeometry);
  double track_pz = track->get_pz();

  if(Verbosity()>0) std::cout << "track (px,py,pz) = (" << track_px << ", " << track_py << ", " << track_pz << ")" << std::endl;

  // Move to first tpc cluster layer if necessary
//  if(sqrt(track_x*track_x+track_y*track_y)<10.)
//    {
      if(Verbosity()>0) std::cout << "WARNING: moving track into TPC" << std::endl;
      std::vector<Acts::Vector3> trkGlobPos;
      for(const auto& ckey : ckeys)
	{
	  if(TrkrDefs::getTrkrId(ckey) == TrkrDefs::tpcId )
	    { trkGlobPos.push_back(globalPositions.at(ckey)); }
	}
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

      double xc = track->get_X0();
      double yc = track->get_Y0();
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

      double pt = sqrt(track_px*track_px+track_py*track_py);
      // rotate track momentum vector (pz stays the same)
      track_px = pt * cos(phi);
      track_py = pt * sin(phi);
      track_x = trkGlobPos.at(0)(0);
      track_y = trkGlobPos.at(0)(1);
      track_z = trkGlobPos.at(0)(2);

      if(Verbosity()>0) std::cout << "new track (x,y,z) = " << track_x << ", " << track_y << ", " << track_z << ")" << std::endl;
      if(Verbosity()>0) std::cout << "new track (px,py,pz) = " << track_px << ", " << track_py << ", " << track_pz << ")" << std::endl;

//    }

  double track_pt = sqrt(track_px*track_px + track_py*track_py);
  if(Verbosity()>0) std::cout << "track pt: " << track_pt << std::endl;

  // get track parameters
  GPUTPCTrackParam kftrack;
  kftrack.InitParam();
  float track_phi = atan2(track_y,track_x);
  if(Verbosity()>0) std::cout << "track phi: " << track_phi << std::endl;
  kftrack.SetQPt(track->get_charge()/track_pt);

  float track_pX = track_px*cos(track_phi)+track_py*sin(track_phi);
  float track_pY = -track_px * sin(track_phi) + track_py * cos(track_phi);

  kftrack.SetSignCosPhi(track_pX/track_pt);
  kftrack.SetSinPhi(track_pY/track_pt);
  kftrack.SetDzDs(-track_pz/track_pt);

  // Y = y
  // Z = z
  // SinPhi = py/sqrt(px^2+py^2)
  // DzDs = pz/sqrt(px^2+py^2)
  // QPt = 1/sqrt(px^2+py^2)

  const double track_pt3 = std::pow( track_pt, 3. );
  
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
  Jrot(3,2) = -track_py*track_px/track_pt3; // dSinPhi/dpx
  Jrot(4,2) = track_px*track_px/track_pt3; // dSinPhi/dpy
  Jrot(5,2) = 0; // dSinPhi/dpz

  Jrot(0,3) = 0; // dDzDs/dx
  Jrot(1,3) = 0; // dDzDs/dy
  Jrot(2,3) = 0; // dDzDs/dz
  Jrot(3,3) = -track_px*track_pz/track_pt3; // dDzDs/dpx
  Jrot(4,3) = -track_py*track_pz/track_pt3; // dDzDs/dpy
  Jrot(5,3) = 1./track_pt; // dDzDs/dpz

  Jrot(0,4) = 0; // dQPt/dx
  Jrot(1,4) = 0; // dQPt/dy
  Jrot(2,4) = 0; // dQPt/dz
  Jrot(3,4) = -track_px/track_pt3; // dQPt/dpx
  Jrot(4,4) = -track_py/track_pt3; // dQPt/dpy
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

  kftrack.SetX(track_x*cos(track_phi)+track_y*sin(track_phi));
  kftrack.SetY(-track_x*sin(track_phi)+track_y*cos(track_phi));
  kftrack.SetZ(track_z);

  if(Verbosity()>0)
  {
    std::cout << "initial track params:" << std::endl;
    std::cout << "X: " << kftrack.GetX() << std::endl;
    std::cout << "Y: " << kftrack.GetY() << std::endl;
    std::cout << "Z: " << kftrack.GetZ() << std::endl;
    std::cout << "SinPhi: " << kftrack.GetSinPhi() << std::endl;
    std::cout << "DzDs: " << kftrack.GetDzDs() << std::endl;
    std::cout << "QPt: " << kftrack.GetQPt() << std::endl;  
    std::cout << "cov: " << std::endl;
    for(int i=0; i<15; i++) std::cout << kftrack.GetCov(i) << ", ";
    std::cout << std::endl;

  }

  // get layer for each cluster
  std::vector<unsigned int> layers;
  std::transform( ckeys.begin(), ckeys.end(), std::back_inserter( layers ), []( const TrkrDefs::cluskey& key ) { return TrkrDefs::getLayer(key); } );

  double old_phi = track_phi;
  unsigned int old_layer = TrkrDefs::getLayer(ckeys[0]);
  if(Verbosity()>0) std::cout << "first layer: " << old_layer << std::endl;
  if(Verbosity()>0) std::cout << "cluster (x,y,z) = (" << globalPositions.at(ckeys[0])(0) << ", " << globalPositions.at(ckeys[0])(1) << ", " << globalPositions.at(ckeys[0])(2) << ")" << std::endl;
  propagated_track.push_back(ckeys[0]);

  GPUTPCTrackLinearisation kfline(kftrack);
  GPUTPCTrackParam::GPUTPCTrackFitParam fp;
  kftrack.CalculateFitParameters(fp);

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
      double cxerr = sqrt(fitter->getClusterError(nc,next_ckey,globalpos,0,0));
      double cyerr = sqrt(fitter->getClusterError(nc,next_ckey,globalpos,1,1));
      double czerr = sqrt(fitter->getClusterError(nc,next_ckey,globalpos,2,2));
      double alpha = cphi-old_phi;
      double tX = kftrack.GetX();
      double tY = kftrack.GetY();
      double tx = tX*cos(old_phi)-tY*sin(old_phi);
      double ty = tX*sin(old_phi)+tY*cos(old_phi);
      double tz = kftrack.GetZ();
//      GPUTPCTrackParam::GPUTPCTrackFitParam fp;
//      kftrack.CalculateFitParameters(fp);
      if(Verbosity()>0) std::cout << "track position: (" << tx << ", " << ty << ", " << tz << ")" << std::endl;
      kftrack.Rotate(alpha/2.,kfline,10.);
      double newX = cx*cos(cphi)+cy*sin(cphi);
      kftrack.TransportToXWithMaterial((tX+newX)/2.,kfline,fp,_Bzconst*get_Bz(tx,ty,tz),10.);
      kftrack.Rotate(alpha/2.,kfline,10.);
      kftrack.TransportToXWithMaterial(cx*cos(cphi)+cy*sin(cphi),kfline,fp,_Bzconst*get_Bz(tx,ty,tz),10.);
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
      if(1)//fabs(tx-cx)<_max_dist*sqrt(txerr*txerr+cxerr*cxerr) &&
         //fabs(ty-cy)<_max_dist*sqrt(tyerr*tyerr+cyerr*cyerr) &&
         //fabs(tz-cz)<_max_dist*sqrt(tzerr*tzerr+czerr*czerr))
      {
        if(Verbosity()>0) std::cout << "Kept cluster" << std::endl;
        propagated_track.push_back(next_ckey);
//        kftrack.SetX(cx*cos(cphi)+cy*sin(cphi));
//        kftrack.SetY(-cx*sin(cphi)+cy*cos(cphi));
//        kftrack.SetZ(cz);
      }
      else
      {
        if(Verbosity()>0)
        {
          std::cout << "Rejected cluster" << std::endl;
          std::cout << "dx: " << fabs(tx-cx) << " when max = " << _max_dist*sqrt(txerr*txerr+cxerr*cxerr) << std::endl;
          std::cout << "dy: " << fabs(ty-cy) << " when max = " << _max_dist*sqrt(tyerr*tyerr+cyerr*cyerr) << std::endl;
          std::cout << "dz: " << fabs(tz-cz) << " when max = " << _max_dist*sqrt(tzerr*tzerr+czerr*czerr) << std::endl;
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
      double newX = radii[l-7];
      double alpha = atan2(ty,tx)-old_phi;
      kftrack.Rotate(alpha/2.,kfline,10.);
      kftrack.TransportToXWithMaterial((tX+newX)/2.,kfline,fp,_Bzconst*get_Bz(tx,ty,tz),10.);
      kftrack.Rotate(alpha/2.,kfline,10.);
      kftrack.TransportToXWithMaterial(radii[l-7],kfline,fp,_Bzconst*get_Bz(tx,ty,tz),10.);
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
          double newX = radius;
          double alpha = atan2(ty,tx)-old_phi;
          kftrack.Rotate(alpha/2.,kfline,10.);
          kftrack.TransportToXWithMaterial((tX+newX)/2.,kfline,fp,_Bzconst*get_Bz(tx,ty,tz),10.);
          kftrack.Rotate(alpha/2.,kfline,10.);
	  kftrack.TransportToXWithMaterial(radius,kfline,fp,_Bzconst*get_Bz(tx,ty,tz),10.);
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

      double cxerr = sqrt(fitter->getClusterError(cc,closest_ckey,ccglob,0,0));
      double cyerr = sqrt(fitter->getClusterError(cc,closest_ckey,ccglob,1,1));
      double czerr = sqrt(fitter->getClusterError(cc,closest_ckey,ccglob,2,2));
      double ccphi = atan2(ccY,ccX);
      if(Verbosity()>0) std::cout << "cluster position: (" << ccX << ", " << ccY << ", " << ccZ << ")" << std::endl;
      if(Verbosity()>0) std::cout << "cluster position error: (" << cxerr << ", " << cyerr << ", " << czerr << ")" << std::endl;
      if(Verbosity()>0) std::cout << "cluster X: " << ccX*cos(ccphi)+ccY*sin(ccphi) << std::endl;
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
*/ 
        
        double alpha = ccphi-old_phi;
        kftrack.Rotate(alpha,kfline,10.);
//        kftrack.SetX(ccX*cos(ccphi)+ccY*sin(ccphi));
//        kftrack.SetY(-ccX*sin(ccphi)+ccY*cos(ccphi));
//        kftrack.SetZ(cc->getZ());
        double ccaY = -ccX*sin(ccphi)+ccY*cos(ccphi);
        double ccerrY = fitter->getClusterError(cc,closest_ckey,ccglob,0,0)*sin(ccphi)*sin(ccphi)+fitter->getClusterError(cc,closest_ckey,ccglob,0,1)*sin(ccphi)*cos(ccphi)+fitter->getClusterError(cc,closest_ckey,ccglob,1,1)*cos(ccphi)*cos(ccphi);
        double ccerrZ = fitter->getClusterError(cc,closest_ckey,ccglob,2,2);
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
//      kftrack.CalculateFitParameters(fp);
      kftrack.Rotate(alpha/2.,kfline,10.);
      double newX = cx*cos(cphi)+cy*sin(cphi);
      kftrack.TransportToXWithMaterial((tX+newX)/2.,kfline,fp,_Bzconst*get_Bz(tx,ty,tz),10.);
      kftrack.Rotate(alpha/2.,kfline,10.);
      kftrack.TransportToXWithMaterial(cx*cos(cphi)+cy*sin(cphi),kfline,fp,_Bzconst*get_Bz(tx,ty,tz),10.);
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
//      kftrack.CalculateFitParameters(fp);
      double newX = radii[l-7];
      double alpha = atan2(ty,tx)-old_phi;
      kftrack.Rotate(alpha/2.,kfline,10.);
      kftrack.TransportToXWithMaterial((tX+newX)/2.,kfline,fp,_Bzconst*get_Bz(tx,ty,tz),10.);
      kftrack.Rotate(alpha/2.,kfline,10.);
      kftrack.TransportToXWithMaterial(radii[l-7],kfline,fp,_Bzconst*get_Bz(tx,ty,tz),10.);
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
          double newX = radius;
          double alpha = atan2(ty,tx)-old_phi;
          kftrack.Rotate(alpha/2.,kfline,10.);
          kftrack.TransportToXWithMaterial((tX+newX)/2.,kfline,fp,_Bzconst*get_Bz(tx,ty,tz),10.);
          kftrack.Rotate(alpha/2.,kfline,10.);
	  kftrack.TransportToXWithMaterial(radius,kfline,fp,_Bzconst*get_Bz(tx,ty,tz),10.);
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

      double cxerr = sqrt(fitter->getClusterError(cc,closest_ckey,ccglob2,0,0));
      double cyerr = sqrt(fitter->getClusterError(cc,closest_ckey,ccglob2,1,1));
      double czerr = sqrt(fitter->getClusterError(cc,closest_ckey,ccglob2,2,2));
      if(Verbosity()>0) std::cout << "cluster position: (" << ccX << ", " << ccY << ", " << ccZ << ")" << std::endl;
      double ccphi = atan2(ccY,ccX);
      if(Verbosity()>0) std::cout << "cluster position errors: (" << cxerr << ", " << cyerr << ", " << czerr << ")" << std::endl;
      if(Verbosity()>0) std::cout << "cluster X: " << ccX*cos(ccphi)+ccY*sin(ccphi) << std::endl;
      double alpha2 = ccphi-old_phi;
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
        kftrack.Rotate(alpha2,kfline,10.);
        double ccaY = -ccX*sin(ccphi)+ccY*cos(ccphi);
        double ccerrY = fitter->getClusterError(cc,closest_ckey,ccglob2,0,0)*sin(ccphi)*sin(ccphi)+fitter->getClusterError(cc,closest_ckey,ccglob2,1,0)*sin(ccphi)*cos(ccphi)+fitter->getClusterError(cc,closest_ckey,ccglob2,1,1)*cos(ccphi)*cos(ccphi);
        double ccerrZ = fitter->getClusterError(cc,closest_ckey,ccglob2,2,2);
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

    //fit a circle through x,y coordinates and calculate residuals
    const auto [R, X0, Y0] = TrackFitUtils::circle_fit_by_taubin( xy_pts );
    const std::vector<double> xy_resid = fitter->GetCircleClusterResiduals(xy_pts,R,X0,Y0);

    // fit a line through r,z coordinates and calculate residuals
    const auto [A, B] = TrackFitUtils::line_fit( rz_pts );
    const std::vector<double> rz_resid = fitter->GetLineClusterResiduals(rz_pts,A,B);
    
    for(size_t i=0;i<chain.size();i++)
    {
      //if(xy_resid[i]>_xy_outlier_threshold) continue;
      clean_chain.push_back(chain[i]);
    }
    clean_chains.push_back(clean_chain);
    if(Verbosity()>0) std::cout << "pushed clean chain with " << clean_chain.size() << " clusters" << std::endl;
  }
  return clean_chains;
}



void PHSimpleKFProp::publishSeeds(std::vector<TrackSeed_v1>& seeds, PositionMap& positions)
{
  for(auto& seed: seeds )
  { 
    /// The ALICEKF gives a better charge determination at high pT
    int q = seed.get_charge();
    seed.circleFitByTaubin(positions,7,55);
    seed.lineFit(positions,7,55);
    seed.set_qOverR(fabs(seed.get_qOverR()) * q);
    _track_map->insert(&seed); 

  }
}

void PHSimpleKFProp::publishSeeds(const std::vector<TrackSeed_v1>& seeds)
{
  for( const auto& seed:seeds )
  { _track_map->insert(&seed); }
}
