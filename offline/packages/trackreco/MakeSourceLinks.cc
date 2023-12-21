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

#include <Acts/EventData/ParticleHypothesis.hpp>
#include <Acts/EventData/SourceLink.hpp>

#include <phool/PHTimer.h>
#include <phool/phool.h>

namespace
{
  /// get radius from coordinates
  template<class T> T radius(const T& x, const T& y)
  { return std::sqrt(square(x) + square(y));}

}

  //___________________________________________________________________________________
SourceLinkVec MakeSourceLinks::getSourceLinks(TrackSeed* track,
						  ActsTrackFittingAlgorithm::MeasurementContainer& measurements,
						  TrkrClusterContainer*  clusterContainer,
						  ActsGeometry* tGeometry,
						  alignmentTransformationContainer* transformMapTransient,
						   std::set< Acts::GeometryIdentifier> transient_id_set,
						  short int crossing)
{
  SourceLinkVec sourcelinks;

  if (m_pp_mode && crossing == SHRT_MAX)
  {
    // Need to skip this in the pp case, for AuAu it should not happen
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
        std::cout << "Failed to get cluster with key " << key << " for track seed" << std::endl;
      else
        std::cout << "PHActsTrkFitter :: Key: " << key << " for track seed " << std::endl;
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
      if (_dcc_static)
      {
        global = _distortionCorrection.get_corrected_position(global, _dcc_static);
      }
      if (_dcc_average)
      {
        global = _distortionCorrection.get_corrected_position(global, _dcc_average);
      }
      if (_dcc_fluctuation)
      {
        global = _distortionCorrection.get_corrected_position(global, _dcc_fluctuation);
      }

      // Make an afine transform that implements the correction as a translation 
      auto correction_translation = (global - global_in) * 10.0;  // need mm
      Acts::Vector3 correction_rotation(0,0,0);   // null rotation
      Acts::Transform3 tcorr = tGeometry->makeAffineTransform(correction_rotation, correction_translation);

      auto this_surf = tGeometry->maps().getSurface(key, cluster);
      Acts::GeometryIdentifier id = this_surf->geometryId();

      // replace the the default alignment transform with the corrected one
      auto ctxt = tGeometry->geometry().getGeoContext();
      alignmentTransformationContainer* transformMap = ctxt.get<alignmentTransformationContainer*>();
      auto corrected_transform = tcorr * transformMap->getTransform(id);
      transformMapTransient->replaceTransform(id, corrected_transform);
      transient_id_set.insert(id);

      /*
      bool checks = true;
      if(checks)
	{
	  auto new_alignmentTransformTransient = m_alignmentTransformationMapTransient->getTransform(id);
	  std::cout << std::endl << " TransformTransient after is: " << std::endl <<   new_alignmentTransformTransient.matrix() << std::endl;

	  std::cout << "Get transform for surf id " << id << " key " << key << " crossing " << crossing << " trkrid " << trkrid  << std::endl
		    << " global_in " << global_in(0) << "  " << global_in(1) << "  " << global_in(2) << std::endl
		    << " global corrected " << global(0) << "  " << global(1) << "  " << global(2) << std::endl
		    << " correction (mm) " << correction_translation(0) << "  " <<  correction_translation(1)  << "  " <<  correction_translation(2)  
		    << std::endl;

	  Acts::Vector2 local = tGeometry->getLocalCoords(key, cluster) * Acts::UnitConstants::cm;	  
	  Eigen::Vector3d clusterLocalPosition(local(0), local(1), 0);  // follows the convention for the acts transform of local = (x,z,y)
	  Eigen::Vector3d finalCoords = corrected_transform * clusterLocalPosition / 10.0;
	  std::cout << " corrected coords  " << global(0) << "  " << global(1) << "  " << global(2) << std::endl;
	  std::cout << " corrected from transform " << finalCoords(0) << "  " << finalCoords(1) << "  " << finalCoords(2) << std::endl;
	}  // end checks
      */

    }  // end TPC specific treatment
    
    // corrected TPC transforms are installed, capture the cluster key    
    cluster_vec.push_back(key);

  }  // end loop over clusters here
  
  // loop over cluster_vec and make source links
  for (unsigned int i = 0; i < cluster_vec.size(); ++i)
    {
      TrkrDefs::cluskey cluskey = cluster_vec[i];

      if (m_ignoreLayer.find(TrkrDefs::getLayer(cluskey)) != m_ignoreLayer.end())
	{
	  if (m_verbosity > 3)
	    {
	      std::cout << PHWHERE << "skipping cluster in layer "
			<< (unsigned int) TrkrDefs::getLayer(cluskey) << std::endl;
	    }
	  continue;
	}      

      // get local coordinates
      auto cluster = clusterContainer->findCluster(cluskey);
      Acts::Vector2 localPos = tGeometry->getLocalCoords(cluskey, cluster);

      Surface surf = tGeometry->maps().getSurface(cluskey, cluster);
      
      Acts::ActsVector<2> loc;
      loc[Acts::eBoundLoc0] = localPos(0) * Acts::UnitConstants::cm;
      loc[Acts::eBoundLoc1] = localPos(1) * Acts::UnitConstants::cm;
      
      std::array<Acts::BoundIndices, 2> indices;
      indices[0] = Acts::BoundIndices::eBoundLoc0;
      indices[1] = Acts::BoundIndices::eBoundLoc1;
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
    std::cout << "PHActsTrkFitter Source Links generation time:  "
              << SLTime << std::endl;

  return sourcelinks;
}

void MakeSourceLinks::resetTransientTransformMap(
						  alignmentTransformationContainer* transformMapTransient,
						  std::set< Acts::GeometryIdentifier>& transient_id_set,
						  ActsGeometry* tGeometry )
{
  // loop over modifiedTransformSet and replace transient elements modified for the last track with the default transforms
  for(auto it = transient_id_set.begin(); it != transient_id_set.end(); ++it)
    {
      Acts::GeometryIdentifier id = *it;
      auto ctxt = tGeometry->geometry().getGeoContext();
      alignmentTransformationContainer* transformMap = ctxt.get<alignmentTransformationContainer*>();
      auto transform = transformMap->getTransform(id);
      // auto transient_transform = m_alignmentTransformationMapTransient->getTransform(id);
      //std::cout << "    transient transform for id " << id << " is " << std::endl << transient_transform.matrix() << std::endl;
      transformMapTransient->replaceTransform(id, transform);
    }
  transient_id_set.clear();
   }
