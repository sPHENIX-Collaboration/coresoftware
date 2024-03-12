/*!
 *  \file PrelimDistortionCorrection.cc
 *  \brief	Makes preliminary TPC distortion corrections in pp running and replaces TPC seed parameters with new fits
 *  \author Tony Frawley
 */

#include "PrelimDistortionCorrection.h"
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

PrelimDistortionCorrection::PrelimDistortionCorrection(const std::string& name)
  : SubsysReco(name)
{}

int PrelimDistortionCorrection::End(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PrelimDistortionCorrection::InitRun(PHCompositeNode* topNode)
{
  
  int ret = get_nodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) { return ret;
}
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
  fitter->setConstBField(_const_field);
  fitter->useFixedClusterError(_use_fixed_clus_err);
  fitter->setFixedClusterError(0,_fixed_clus_err.at(0));
  fitter->setFixedClusterError(1,_fixed_clus_err.at(1));
  fitter->setFixedClusterError(2,_fixed_clus_err.at(2));
  //  _field_map = PHFieldUtility::GetFieldMapNode(nullptr,topNode);
  // m_Cache = magField->makeCache(m_tGeometry->magFieldContext);

  return Fun4AllReturnCodes::EVENT_OK;
}

double PrelimDistortionCorrection::get_Bz(double x, double y, double z) const
{
  if(_use_const_field) { return _const_field;
}
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

int PrelimDistortionCorrection::get_nodes(PHCompositeNode* topNode)
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
  m_dcc = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerStatic");
  if( m_dcc )
  { 
    std::cout << "PrelimDistortionCorrection::InitRun - found TPC distortion correction container: TpcDistortionContainerStatic" << std::endl; 
  }
  else
    {
      std::cout << "PrelimDistortionCorrection::InitRun - failed to find TPC distortion correction container: TpcDistortionContainerStatic" << std::endl; 
    }

  if(_use_truth_clusters) {
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_TRUTH");
  } else {
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
}

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

int PrelimDistortionCorrection::process_event(PHCompositeNode* /*topNode*/)
{
  if(!_pp_mode) 
    {
      std::cout << "PrelimDistortionCorrection called but not in pp_mode, do nothing!" << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  
  PHTimer timer("PrelimDistortionCorrectionTimer");
  timer.stop();
  timer.restart();

  if(Verbosity()>0) { std::cout << "starting PrelimDistortionCorrection process_event" << std::endl;
}

  // These accumulate trackseeds as we loop over tracks, keylist is a vector of vectors of seed cluskeys
  // The seeds are all given to the fitter at once
    std::vector<std::vector<TrkrDefs::cluskey>> keylist;
    PositionMap correctedOffsetTrackClusPositions;   //  this is an std::map<TrkrDefs::cluskey, Acts::Vector3>

  for(int track_it = 0; track_it != _track_map->size(); ++track_it )
  {
    if(Verbosity()>0) { std::cout << "TPC seed " << track_it << std::endl; }
    // if not a TPC track, ignore
    TrackSeed* track = _track_map->get(track_it);
    if(!track) 
      {
	continue;
      }

    if(Verbosity() > 0)
      {
	std::cout << "Input seed pars for " << track_it
		  << " q " << track->get_charge()
		  << " qOverR " << fabs(track->get_qOverR()) * track->get_charge()
		  << " X0 " << track->get_x()
		  << " Y0 " << track->get_y()
		  << " Z0 " << track->get_z()
		  << " eta " << track->get_eta()
		  << std::endl;
      }

    const bool is_tpc = std::any_of(
      track->begin_cluster_keys(),
      track->end_cluster_keys(),
      []( const TrkrDefs::cluskey& key ) { return TrkrDefs::getTrkrId(key) == TrkrDefs::tpcId; } );

    if(is_tpc)
    {
      // We want to move all clusters in this seed to point to Z0 = 0
      float Z0 = track->get_Z0();
      float offset_Z = 0.0 - Z0;

      if(Verbosity() > 0)
	{
	  std::cout << "  processing seed with offset_Z " << offset_Z
		    << " eta " << track->get_eta()
		    << " x " << track->get_x() 
		    << " y " << track->get_y() 
		    << " z " << track->get_z() 
		    << std::endl;
	}

      // We want to make  distortion corrections to all clusters in this seed after offsetting the z values
      std::vector<TrkrDefs::cluskey> dumvec;
      for(TrackSeed::ConstClusterKeyIter iter = track->begin_cluster_keys();
	  iter != track->end_cluster_keys();
	  ++iter)
	{
	  TrkrDefs::cluskey cluskey = *iter;
	  TrkrCluster *cluster = _cluster_map->findCluster(cluskey);

	  Acts::Vector3 pos = getGlobalPosition(cluskey, cluster);
	  Acts::Vector3 offsetpos(pos(0), pos(1), pos(2) + offset_Z);
	  // Distortion correct the offset positions
	  if( m_dcc ) { offsetpos = m_distortionCorrection.get_corrected_position( offsetpos, m_dcc ); }

	  // now move the distortion corrected cluster back by offset_Z, to preserve the z measurement info
	  Acts::Vector3 corrpos(offsetpos(0), offsetpos(1), offsetpos(2) - offset_Z);
	  correctedOffsetTrackClusPositions.insert(std::make_pair(cluskey,corrpos));
	  dumvec.push_back(cluskey);

	  if(Verbosity() > 0)
	    {
	      std::cout << " cluskey " << cluskey << " input pos " << pos(0) << "  " << pos(1) << "  " << pos(2) 
			<< "   corr. pos " << corrpos(0) << "  " << corrpos(1) << "  " << corrpos(2) << std::endl
			<< "distortion correction " <<  corrpos(0) - pos(0) << "  " << corrpos(1) - pos(1) << "  " << corrpos(2) - pos(2)
			<< std::endl;

	    }
	}

      /// Can't circle fit a seed with less than 3 clusters, skip it
      if(dumvec.size() < 3)
	{ continue; }
      keylist.push_back(dumvec);

      if(Verbosity() > 0)
	{
	  std::cout << "Added  input seed " << track_it << "  becomes output seed " << keylist.size() - 1 << std::endl;
	}

    } // end if TPC seed

  }  // end loop over tracks

  // reset the seed map for the TPC
  _track_map->Reset();
      
  // refit the corrected clusters 
  std::vector<float> trackChi2;
  auto seeds = fitter->ALICEKalmanFilter(keylist, true, correctedOffsetTrackClusPositions,
					 trackChi2);
  
  // update the seed parameters on the node tree
  // This calls circle fit and line fit (setting x, y, z, phi, eta) and sets qOverR explicitly using q from KF
  publishSeeds(seeds.first, correctedOffsetTrackClusPositions);

  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::Vector3 PrelimDistortionCorrection::getGlobalPosition( TrkrDefs::cluskey key, TrkrCluster* cluster ) const
{
  // get global position from Acts transform

  auto globalpos = _tgeometry->getGlobalPosition(key, cluster);

  // check if TPC distortion correction are in place and apply if this is a triggered event (ie. crossing is known)
  if(!_pp_mode)
    {
      if( m_dcc ) { globalpos = m_distortionCorrection.get_corrected_position( globalpos, m_dcc ); }
    }
  
  return globalpos;
}

PositionMap PrelimDistortionCorrection::PrepareKDTrees()
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
      if(!cluster) { continue;
}
      if(_n_iteration!=0){
        if(_iteration_map != nullptr ){
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
    if(Verbosity()>0) { std::cout << "l: " << l << std::endl;
}
    _ptclouds[l] = std::make_shared<KDPointCloud<double>>();
    _ptclouds[l]->pts.resize(kdhits[l].size());
    if(Verbosity()>0) { std::cout << "resized to " << kdhits[l].size() << std::endl;
}
    for(size_t i=0;i<kdhits[l].size();++i)
    {
      _ptclouds[l]->pts[i] = kdhits[l][i];
    }
    _kdtrees[l] = std::make_shared<nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, KDPointCloud<double>>, KDPointCloud<double>, 3>>(3,*(_ptclouds[l]),nanoflann::KDTreeSingleIndexAdaptorParams(10));
    _kdtrees[l]->buildIndex();
  }

  return globalPositions;
}


void PrelimDistortionCorrection::publishSeeds(std::vector<TrackSeed_v1>& seeds, PositionMap& positions)
{
  int seed_index = 0;
  for(auto& seed: seeds )
  { 
    /// The ALICEKF gives a better charge determination at high pT
    int q = seed.get_charge();
    seed.circleFitByTaubin(positions, 7, 55);
    seed.lineFit(positions, 7, 55);
    
    seed.set_qOverR(fabs(seed.get_qOverR()) * q);
    _track_map->insert(&seed); 

    if(Verbosity() > 0)
      {
	std::cout << "Publishing seed " << seed_index
		  << " q " << q
		  << " qOverR " << fabs(seed.get_qOverR()) * q 
		  << " X0 " << seed.get_x()
		  << " Y0 " << seed.get_y()
		  << " Z0 " << seed.get_z()
		  << " eta " << seed.get_eta()
		  << std::endl;
      }

    seed_index++;
  }
}

void PrelimDistortionCorrection::publishSeeds(const std::vector<TrackSeed_v1>& seeds)
{
  for( const auto& seed:seeds )
  { _track_map->insert(&seed); }
}
