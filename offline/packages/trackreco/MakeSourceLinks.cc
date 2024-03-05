#include "MakeSourceLinks.h"

#include <trackbase/TrkrCluster.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/ActsSourceLink.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/alignmentTransformationContainer.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/SvtxTrackState_v2.h>
#include <trackbase_historic/TrackSeed.h>

#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <Acts/EventData/ParticleHypothesis.hpp>
#include <Acts/EventData/SourceLink.hpp>

#include <phool/PHTimer.h>
#include <phool/phool.h>

namespace
{
 template <class T>
  inline T square(const T& x)
  {
    return x * x;
  }

}

void MakeSourceLinks::initialize(PHG4TpcCylinderGeomContainer* cellgeo)
{
  // get the TPC layer radii from the geometry object  
  if (cellgeo)
  {
    _clusterMover.initialize_geometry(cellgeo);
  }

}

  //___________________________________________________________________________________
SourceLinkVec MakeSourceLinks::getSourceLinks(TrackSeed* track,
						  ActsTrackFittingAlgorithm::MeasurementContainer& measurements,
						  TrkrClusterContainer*  clusterContainer,
						  ActsGeometry* tGeometry,
					          const TpcDistortionCorrectionContainer* dcc_static,
					          const TpcDistortionCorrectionContainer* dcc_average,
					          const TpcDistortionCorrectionContainer* dcc_fluctuation,
						  alignmentTransformationContainer* transformMapTransient,
						   std::set< Acts::GeometryIdentifier>& transient_id_set,
						  short int crossing)
{
  if(m_verbosity > 1) { std::cout << "Entering MakeSourceLinks::getSourceLinks " << std::endl; }

  SourceLinkVec sourcelinks;

  if (m_pp_mode && crossing == SHRT_MAX)
  {
    // Need to skip this in the pp case, for AuAu it should not happen
    if(m_verbosity > 1) 
      { std::cout << "Seed has no crossing, and in pp mode: skip this seed" << std::endl; }

    return sourcelinks;
  }

  PHTimer SLTrackTimer("SLTrackTimer");
  SLTrackTimer.stop();
  SLTrackTimer.restart();

  // loop over all clusters
  std::vector<TrkrDefs::cluskey> cluster_vec;

  for (auto clusIter = track->begin_cluster_keys();
       clusIter != track->end_cluster_keys();
       ++clusIter)
  {
    auto key = *clusIter;
    auto cluster = clusterContainer->findCluster(key);
    if (!cluster)
    {
      if (m_verbosity > 0)
        {std::cout << "MakeSourceLinks: Failed to get cluster with key " << key << " for track seed" << std::endl;}
      else
        {std::cout << "MakeSourceLinks: Key " << key << " for track seed " << std::endl;}
      continue;
    }

     /// Make a safety check for clusters that couldn't be attached  to a surface
    auto surf = tGeometry->maps().getSurface(key, cluster);
    if (!surf)
    {
      continue;
    }

    const unsigned int trkrid = TrkrDefs::getTrkrId(key);
    const unsigned int side = TpcDefs::getSide(key);

    if(m_verbosity > 1) { std::cout << "    Cluster key " << key << " trkrid " << trkrid << " crossing " << crossing << std::endl; }

    // For the TPC, cluster z has to be corrected for the crossing z offset, distortion, and TOF z offset
    // we do this by modifying the fake surface transform, to move the cluster to the corrected position
    if (trkrid == TrkrDefs::tpcId)
    {
      Acts::Vector3 global = tGeometry->getGlobalPosition(key, cluster);
      Acts::Vector3 global_in = global;

      // make all corrections to global position of TPC cluster
      float z = _clusterCrossingCorrection.correctZ(global[2], side, crossing);
      global[2] = z;

      // apply distortion corrections
      if (dcc_static)
      {
        global = _distortionCorrection.get_corrected_position(global, dcc_static);
      }
      if (dcc_average)
      {
        global = _distortionCorrection.get_corrected_position(global, dcc_average);
      }
      if (dcc_fluctuation)
      {
        global = _distortionCorrection.get_corrected_position(global, dcc_fluctuation);
      }

      if(m_verbosity > 2)
	{
	  std::cout << " global_in " << global_in(0) << "  " << global_in(1) << "  " << global_in(2) 
		    << " corr glob " << global(0) << "  " << global(1) << "  " << global(2) << std::endl
		    << " crossing z correction " << z - global_in(2) 
		    << " distortion correction " << global(0)-global_in(0) << "  " << global(1) - global_in(1) << "  " << global(2) - z 
		    << std::endl;
	}

      // Make an afine transform that implements the correction as a translation 
      auto correction_translation = (global - global_in) * 10.0;  // need mm
      Acts::Vector3 correction_rotation(0,0,0);   // null rotation
      Acts::Transform3 tcorr = tGeometry->makeAffineTransform(correction_rotation, correction_translation);

      auto this_surf = tGeometry->maps().getSurface(key, cluster);
      Acts::GeometryIdentifier id = this_surf->geometryId();

      auto check_cluster = clusterContainer->findCluster(key);
      Acts::Vector2 check_local2d = tGeometry->getLocalCoords(key, check_cluster) * Acts::UnitConstants::cm; // need mm
      Acts::Vector3 check_local3d (check_local2d(0), check_local2d(1), 0);
      Acts::GeometryContext temp_transient_geocontext;
      temp_transient_geocontext =  transformMapTransient;
      Acts::Vector3 check_before_pos_surf = this_surf->localToGlobal( temp_transient_geocontext,
				  check_local2d,
				  Acts::Vector3(1,1,1));
      if(m_verbosity > 2)
	{ 
      std::cout << "Check global from transient transform BEFORE via surface method " << check_before_pos_surf(0)/10.0 << "  " 
		<< "  " << check_before_pos_surf(1)/10.0 << "  " << check_before_pos_surf(2)/10.0 << std::endl;
	}
      // replace the the default alignment transform with the corrected one
      auto ctxt = tGeometry->geometry().getGeoContext();
      alignmentTransformationContainer* transformMap = ctxt.get<alignmentTransformationContainer*>();
      auto corrected_transform = tcorr * transformMap->getTransform(id);
      transformMapTransient->replaceTransform(id, corrected_transform);
      transient_id_set.insert(id);

      Acts::Vector3 check_after_pos_surf = this_surf->localToGlobal( temp_transient_geocontext,
				  check_local2d,
				  Acts::Vector3(1,1,1));
      if(m_verbosity > 2)
	{
      std::cout << "Check global from transient transform AFTER via surface method " << check_after_pos_surf(0)/10.0 << "  " 
		<< "  " << check_after_pos_surf(1)/10.0 << "  " << check_after_pos_surf(2)/10.0 << std::endl;
	}
      //      std::cout << " Print ideal transform matrix from surface object: " << std::endl
      //	<< surf->transform(tGeometry->geometry().getGeoContext()).matrix() << std::endl;
      //      std::cout << " Print transient transform matrix from surface object: " << std::endl
      //	<< surf->transform(temp_transient_geocontext).matrix() << std::endl;

    }  // end TPC specific treatment
    
    // corrected TPC transforms are installed, capture the cluster key    
    cluster_vec.push_back(key);

  }  // end loop over clusters here

  Acts::GeometryContext transient_geocontext;
  transient_geocontext =  transformMapTransient;
  
  // loop over cluster_vec and make source links
  for(auto& cluskey : cluster_vec)
    {
      if (m_ignoreLayer.find(TrkrDefs::getLayer(cluskey)) != m_ignoreLayer.end())
	{
	  if (m_verbosity > 3)
	    {
	      std::cout << PHWHERE << "skipping cluster in layer "
			<< (unsigned int) TrkrDefs::getLayer(cluskey) << std::endl;
	    }
	  continue;
	}      

      // get local coordinates (TPC time needs conversion to cm)
      auto cluster = clusterContainer->findCluster(cluskey);
      Acts::Vector2 localPos = tGeometry->getLocalCoords(cluskey, cluster);   //  cm

      Surface surf = tGeometry->maps().getSurface(cluskey, cluster);
      
      Acts::ActsVector<2> loc;
      loc[Acts::eBoundLoc0] = localPos(0) * Acts::UnitConstants::cm;   // mm
      loc[Acts::eBoundLoc1] = localPos(1) * Acts::UnitConstants::cm;
      
      std::array<Acts::BoundIndices, 2> indices = 
        {Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundLoc1};
      Acts::ActsSquareMatrix<2> cov = Acts::ActsSquareMatrix<2>::Zero();

      // get errors
      Acts::Vector3 global = tGeometry->getGlobalPosition(cluskey, cluster);      
      double clusRadius = sqrt(global[0] * global[0] + global[1] * global[1]);
      auto para_errors = _ClusErrPara.get_clusterv5_modified_error(cluster, clusRadius, cluskey);
      cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = para_errors.first * Acts::UnitConstants::cm2;
      cov(Acts::eBoundLoc0, Acts::eBoundLoc1) = 0;
      cov(Acts::eBoundLoc1, Acts::eBoundLoc0) = 0;
      cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = para_errors.second * Acts::UnitConstants::cm2;
      
      ActsSourceLink::Index index = measurements.size();
      
      SourceLink sl(surf->geometryId(), index, cluskey);
      Acts::SourceLink actsSL{sl};
      Acts::Measurement<Acts::BoundIndices, 2> meas(actsSL, indices, loc, cov);
      if (m_verbosity > 3)
	{
	  std::cout << "source link " << sl.index() << ", loc : "
		    << loc.transpose() << std::endl
		    << ", cov : " << cov.transpose() << std::endl
		    << " geo id " << sl.geometryId() << std::endl;
	  std::cout << "Surface : " << std::endl;
	  //surf.get()->toStream(tGeometry->geometry().getGeoContext(), std::cout);
	  surf.get()->toStream(transient_geocontext, std::cout);
	  std::cout << std::endl;
	  std::cout << "Corrected surface transform:" << std::endl;
	  std::cout << transformMapTransient->getTransform(surf->geometryId()).matrix() << std::endl;
	  std::cout << "Cluster error " << cluster->getRPhiError() << " , " << cluster->getZError() << std::endl;
	  std::cout << "For key " << cluskey << " with local pos " << std::endl
		    << localPos(0) << ", " << localPos(1)
		    << std::endl << std::endl;
	}
      
      sourcelinks.push_back(actsSL);
      measurements.push_back(meas);
    }
 
  SLTrackTimer.stop();
  auto SLTime = SLTrackTimer.get_accumulated_time();
 
  if (m_verbosity > 1)
    {
      std::cout << "PHActsTrkFitter Source Links generation time:  "
              << SLTime << std::endl;
    }
  return sourcelinks;
}

void MakeSourceLinks::resetTransientTransformMap(
						  alignmentTransformationContainer* transformMapTransient,
						  std::set< Acts::GeometryIdentifier>& transient_id_set,
						  ActsGeometry* tGeometry )
{
  if(m_verbosity > 2) { std::cout << "Resetting TransientTransformMap with transient_id_set size " << transient_id_set.size() << std::endl; }

  // loop over modifiedTransformSet and replace transient elements modified for the last track with the default transforms
  for(auto& id : transient_id_set)
    {
      auto ctxt = tGeometry->geometry().getGeoContext();
      alignmentTransformationContainer* transformMap = ctxt.get<alignmentTransformationContainer*>();
      auto transform = transformMap->getTransform(id);
      transformMapTransient->replaceTransform(id, transform);
    }
  transient_id_set.clear();
}


//___________________________________________________________________________________
SourceLinkVec MakeSourceLinks::getSourceLinksClusterMover(
							  TrackSeed* track,
							  ActsTrackFittingAlgorithm::MeasurementContainer& measurements,
							  TrkrClusterContainer*  clusterContainer,
							  ActsGeometry* tGeometry,
							  const TpcDistortionCorrectionContainer* dcc_static,
							  const TpcDistortionCorrectionContainer* dcc_average,
							  const TpcDistortionCorrectionContainer* dcc_fluctuation,
							  short int crossing
							  )
{
  if(m_verbosity > 1) { std::cout << "Entering MakeSourceLinks::getSourceLinksClusterMover  for seed  " 
				  << " with crossing " << crossing
				  << std::endl; }

  SourceLinkVec sourcelinks;

  if (m_pp_mode && crossing == SHRT_MAX)
  {
    // Need to skip this in the pp case, for AuAu it should not happen
    if(m_verbosity > 1) 
      { std::cout << "Seed has no crossing, and in pp mode: skip this seed" << std::endl; }

    return sourcelinks;
  }

  PHTimer SLTrackTimer("SLTrackTimer");
  SLTrackTimer.stop();
  SLTrackTimer.restart();

  // loop over all clusters
  std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> global_raw;

  for (auto clusIter = track->begin_cluster_keys();
       clusIter != track->end_cluster_keys();
       ++clusIter)
  {
    auto key = *clusIter;
    auto cluster = clusterContainer->findCluster(key);
    if (!cluster)
    {
      if (m_verbosity > 0)
        {std::cout << "Failed to get cluster with key " << key << " for track " << track << std::endl;}
      else
        {std::cout << "PHActsTrkFitter :: Key: " << key << " for track " << track << std::endl;}
      continue;
    }

    //    auto subsurfkey = cluster->getSubSurfKey();

    /// Make a safety check for clusters that couldn't be attached
    /// to a surface
    auto surf = tGeometry->maps().getSurface(key, cluster);
    if (!surf)
    {
      continue;
    }

    const unsigned int trkrid = TrkrDefs::getTrkrId(key);
    const unsigned int side = TpcDefs::getSide(key);

    if(m_verbosity > 1) { std::cout << "    Cluster key " << key << " trkrid " << trkrid << " crossing " << crossing << std::endl; }

    // For the TPC, cluster z has to be corrected for the crossing z offset, distortion, and TOF z offset
    // we do this locally here and do not modify the cluster, since the cluster may be associated with multiple silicon tracks
    Acts::Vector3 global = tGeometry->getGlobalPosition(key, cluster);
    Acts::Vector3 global_in = global;
    
    if (trkrid == TrkrDefs::tpcId)
      {
	// make all corrections to global position of TPC cluster
	float z = _clusterCrossingCorrection.correctZ(global[2], side, crossing);
	global[2] = z;

	// apply distortion corrections
	if (dcc_static)
	  {
	    global = _distortionCorrection.get_corrected_position(global, dcc_static);
	  }
	if (dcc_average)
	  {
	    global = _distortionCorrection.get_corrected_position(global, dcc_average);
	  }
	if (dcc_fluctuation)
	  {
	    global = _distortionCorrection.get_corrected_position(global, dcc_fluctuation);
	  }

        if (m_verbosity > 1)
	  {
	    std::cout << " global_in " << global_in(0) << "  " << global_in(1) << "  " << global_in(2) 
		      << " corr glob (unmoved) " << global(0) << "  " << global(1) << "  " << global(2) << std::endl
		      << " crossing z correction " << z - global_in(2) 
		      << " distortion correction " << global(0)-global_in(0) << "  " << global(1) - global_in(1) << "  " << global(2) - z 
		      << std::endl;
	  }
      }

    if (m_verbosity > 1)
      {
	std::cout << " AGAIN: global_in " << global_in(0) << "  " << global_in(1) << "  " << global_in(2) 
		  << " corr glob " << global(0) << "  " << global(1) << "  " << global(2) << std::endl
		  << std::endl;
      }
    
    // add the global positions to a vector to give to the cluster mover
    global_raw.emplace_back(std::make_pair(key, global));
    
  }  // end loop over clusters here

  // move the cluster positions back to the original readout surface
  auto global_moved = _clusterMover.processTrack(global_raw);

if (m_verbosity > 1)
    {
      std::cout << "Cluster global positions after mover puts them on readout surface:" << std::endl;
    }

  // loop over global positions returned by cluster mover
  for(auto& [cluskey, global] : global_moved)
  {
    if (m_ignoreLayer.find(TrkrDefs::getLayer(cluskey)) != m_ignoreLayer.end())
    {
      if (m_verbosity > 3)
      {
        std::cout << PHWHERE << "skipping cluster in layer "
                  << (unsigned int) TrkrDefs::getLayer(cluskey) << std::endl;
      }
      continue;
    }

    auto cluster = clusterContainer->findCluster(cluskey);
    Surface surf = tGeometry->maps().getSurface(cluskey, cluster);

    //    std::cout << " Initial surface " <<  std::endl;
    //surf.get()->toStream(tGeometry->geometry().getGeoContext(), std::cout);
    //std::cout << std::endl;

    // if this is a TPC cluster, the crossing correction may have moved it across the central membrane, check the surface
    auto trkrid = TrkrDefs::getTrkrId(cluskey);
    if (trkrid == TrkrDefs::tpcId)
    {
      TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);
      TrkrDefs::subsurfkey new_subsurfkey = 0;
      surf = tGeometry->get_tpc_surface_from_coords(hitsetkey, global, new_subsurfkey);
    }

    if (!surf)
    {
      if(m_verbosity > 2) { std::cout << " Failed to find surface for cluskey " << cluskey << std::endl; }
      continue;
    }

    //std::cout << " Moved surface " <<  std::endl;
    //surf.get()->toStream(tGeometry->geometry().getGeoContext(), std::cout);
    //std::cout << std::endl;

    // get local coordinates
    Acts::Vector2 localPos;
    global *= Acts::UnitConstants::cm;  // we want mm for transformations

    Acts::Vector3 normal = surf->normal(tGeometry->geometry().getGeoContext());
    auto local = surf->globalToLocal(tGeometry->geometry().getGeoContext(),
                                     global, normal);

    if (local.ok())
    {
      localPos = local.value() / Acts::UnitConstants::cm;
    }
    else
    {
      if(m_verbosity > 2) { std::cout << "Taking manual calculation for global to local " << std::endl; }

      /// otherwise take the manual calculation for the TPC
      // doing it this way just avoids the bounds check that occurs in the surface class method
      Acts::Vector3 loct = surf->transform(tGeometry->geometry().getGeoContext()).inverse() * global ;  // global is in mm
      loct /= Acts::UnitConstants::cm;

      localPos(0) = loct(0);
      localPos(1) = loct(1);
    }

    if (m_verbosity > 2)
    {
      std::cout << "Cluster " << cluskey << " cluster global after mover: " << global(0)/10.0 << "  " << global(1)/10.0 << "  " << global(2)/10.0 << std::endl;
      std::cout << " stored: cluster local X " << cluster->getLocalX() << " cluster local Y " << cluster->getLocalY() << std::endl;
      Acts::Vector2 localTest = tGeometry->getLocalCoords(cluskey, cluster);   //  cm  
      std::cout << " localTest from getLocalCoords: " << localTest(0) << "  " << localTest(1) << std::endl;
      std::cout << " new from inverse transform of cluster global after mover:     local X " << localPos(0) << " new       local Y " << localPos(1) << std::endl;
      Acts::Vector3 globalTest = surf->localToGlobal(tGeometry->geometry().getGeoContext(), localTest*Acts::UnitConstants::cm, normal);
      std::cout << " global from localTest: " << globalTest(0)/10.0 << "  " << globalTest(1)/10.0 << "  " << globalTest(2)/10.0 << std::endl; 
      Acts::Vector3 globalNew = surf->localToGlobal(tGeometry->geometry().getGeoContext(), localPos*Acts::UnitConstants::cm, normal);
      std::cout << " global from new local: " << globalNew(0)/10.0 << "  " << globalNew(1)/10.0 << "  " << globalNew(2)/10.0 << std::endl; 
    }

    Acts::ActsVector<2> loc;
    loc[Acts::eBoundLoc0] = localPos(0) * Acts::UnitConstants::cm;
    loc[Acts::eBoundLoc1] = localPos(1) * Acts::UnitConstants::cm;
    std::array<Acts::BoundIndices, 2> indices = 
      {Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundLoc1};
    
    Acts::ActsSquareMatrix<2> cov = Acts::ActsSquareMatrix<2>::Zero();

    double clusRadius = sqrt(global[0] * global[0] + global[1] * global[1]);
    auto para_errors = _ClusErrPara.get_clusterv5_modified_error(cluster, clusRadius, cluskey);
    cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = para_errors.first * Acts::UnitConstants::cm2;
    cov(Acts::eBoundLoc0, Acts::eBoundLoc1) = 0;
    cov(Acts::eBoundLoc1, Acts::eBoundLoc0) = 0;
    cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = para_errors.second * Acts::UnitConstants::cm2;

    ActsSourceLink::Index index = measurements.size();

    SourceLink sl(surf->geometryId(), index, cluskey);
    Acts::SourceLink actsSL{sl};
    Acts::Measurement<Acts::BoundIndices, 2> meas(actsSL, indices, loc, cov);
    if (m_verbosity > 3)
    {
      std::cout << "source link " << sl.index() << ", loc : "
                << loc.transpose() << std::endl
                << ", cov : " << cov.transpose() << std::endl
                << " geo id " << sl.geometryId() << std::endl;
      std::cout << "Surface : " << std::endl;
      surf.get()->toStream(tGeometry->geometry().getGeoContext(), std::cout);
      std::cout << std::endl;
      std::cout << "Cluster error " << cluster->getRPhiError() << " , " << cluster->getZError() << std::endl;
      std::cout << "For key " << cluskey << " with local pos " << std::endl
                << localPos(0) << ", " << localPos(1)
                << std::endl;
    }

    sourcelinks.push_back(actsSL);
    measurements.push_back(meas);
  }

  SLTrackTimer.stop();
  auto SLTime = SLTrackTimer.get_accumulated_time();

  if (m_verbosity > 1)
    {
      std::cout << "PHActsTrkFitter Source Links generation time:  "
              << SLTime << std::endl;
    }
  return sourcelinks;
}
