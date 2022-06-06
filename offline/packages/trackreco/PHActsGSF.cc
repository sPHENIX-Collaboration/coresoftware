#include "PHActsGSF.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <tpc/TpcDefs.h>

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>

#include <ActsExamples/EventData/Index.hpp>

//____________________________________________________________________________..
PHActsGSF::PHActsGSF(const std::string &name):
 SubsysReco(name)
{
}

//____________________________________________________________________________..
PHActsGSF::~PHActsGSF()
{
}

//____________________________________________________________________________..
int PHActsGSF::InitRun(PHCompositeNode *topNode)
{

  if(Verbosity() > 1)
    { std::cout << "PHActsGSF::InitRun begin" << std::endl; }

  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    { return Fun4AllReturnCodes::ABORTEVENT; }

  m_fitCfg.fit = ActsExamples::TrackFittingAlgorithm::makeGsfFitterFunction(
			       m_tGeometry->tGeometry,
			       m_tGeometry->magField);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHActsGSF::process_event(PHCompositeNode*)
{
  auto logLevel = Acts::Logging::FATAL;
  if(Verbosity() > 4) 
    { logLevel = Acts::Logging::VERBOSE; }

  auto logger = Acts::getDefaultLogger("PHActsGSF", logLevel);

  for(const auto& [key, track] : *m_trackMap)
    {
      auto pSurface = makePerigee(track);
      const auto seed = makeSeed(track, pSurface);

      ActsExamples::MeasurementContainer measurements;
      TrackSeed *tpcseed = track->get_tpc_seed();
      TrackSeed *silseed = track->get_silicon_seed();

      /// We only fit full sPHENIX tracks
      if(!silseed or !tpcseed)
	{ continue; }
      
      auto crossing = silseed->get_crossing();
      if(crossing == SHRT_MAX)
	{ continue; }

      auto sourceLinks = getSourceLinks(tpcseed, measurements);
      auto silSourceLinks = getSourceLinks(silseed, measurements);

      for(auto& siSL : silSourceLinks)
	{ sourceLinks.push_back(siSL); }

      /// Acts requires a wrapped vector, so we need to replace the
      /// std::vector contents with a wrapper vector to get the memory
      /// access correct
      std::vector<std::reference_wrapper<const SourceLink>>  wrappedSls;
      for(const auto& sl : sourceLinks)
	{ wrappedSls.push_back(std::cref(sl)); }

      ActsExamples::MeasurementCalibrator calibrator(measurements);
      ActsExamples::TrackFittingAlgorithm::GeneralFitterOptions options{m_tGeometry->geoContext,
	  m_tGeometry->magFieldContext,
	  m_tGeometry->calibContext,
	  calibrator, &(*pSurface), Acts::LoggerWrapper(*logger)};

      auto result = fitTrack(wrappedSls, seed, options);
      if(result.ok())
	{
	  const FitResult& output = result.value();
	  updateTrack(output, track);
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;
}


std::shared_ptr<Acts::PerigeeSurface> PHActsGSF::makePerigee(SvtxTrack* track) const
{
  SvtxVertex* vertex = m_vertexMap->get(track->get_vertex_id());
  
  Acts::Vector3 vertexpos(vertex->get_x() * Acts::UnitConstants::cm,
			  vertex->get_y() * Acts::UnitConstants::cm,
			  vertex->get_z() * Acts::UnitConstants::cm);
  
  return Acts::Surface::makeShared<Acts::PerigeeSurface>(
					  vertexpos);
}

ActsExamples::TrackParameters PHActsGSF::makeSeed(SvtxTrack* track,
						  std::shared_ptr<Acts::PerigeeSurface> psurf) const
{
  
  Acts::Vector4 fourpos(track->get_x() * Acts::UnitConstants::cm,
			track->get_y() * Acts::UnitConstants::cm,
			track->get_z() * Acts::UnitConstants::cm,
			10 * Acts::UnitConstants::ns);



  int charge = track->get_charge();
  Acts::Vector3 momentum(track->get_px(),
			 track->get_py(),
			 track->get_pz());

  ActsTransformations transformer;
  auto cov = transformer.rotateSvtxTrackCovToActs(track);

  return ActsExamples::TrackParameters::create(psurf, 
					       m_tGeometry->geoContext,
					       fourpos, 
					       momentum, 
					       charge / momentum.norm(),
					       cov).value();
}

SourceLinkVec PHActsGSF::getSourceLinks(TrackSeed* track,
					ActsExamples::MeasurementContainer& measurements)
{
  SourceLinkVec sls;
  // loop over all clusters
  std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> global_raw;

  ActsTransformations transformer;

  for (auto clusIter = track->begin_cluster_keys();
       clusIter != track->end_cluster_keys();
       ++clusIter)
    {
      auto key = *clusIter;
      auto cluster = m_clusterContainer->findCluster(key);
      if(!cluster)
	{
	  if(Verbosity() > 0) std::cout << "Failed to get cluster with key " << key << std::endl;
	  continue;
	}

      auto subsurfkey = cluster->getSubSurfKey();
      
      /// Make a safety check for clusters that couldn't be attached
      /// to a surface
      auto surf = transformer.getSurface(key, cluster, m_surfMaps);
      if(!surf)
	{ continue; }

      unsigned int trkrid = TrkrDefs::getTrkrId(key);
      unsigned int side = TpcDefs::getSide(key);

      // For the TPC, cluster z has to be corrected for the crossing z offset, distortion, and TOF z offset 
      // we do this locally here and do not modify the cluster, since the cluster may be associated with multiple silicon tracks  
      
      // transform to global coordinates for z correction 
      ActsTransformations transformer;
      auto global = transformer.getGlobalPosition(key, cluster,
						  m_surfMaps,
						  m_tGeometry);
      
      if(Verbosity() > 0)
	{
	  std::cout << " zinit " << global[2] << " xinit " << global[0] << " yinit " << global[1] << " side " << side << " crossing " << track->get_crossing() 
		    << " cluskey " << key << " subsurfkey " << subsurfkey << std::endl;
	}
      
      if(trkrid ==  TrkrDefs::tpcId)
	{	  
	  // make all corrections to global position of TPC cluster
	  float z = m_clusterCrossingCorrection.correctZ(global[2], side, track->get_crossing());
	  global[2] = z;
	  
	  // apply distortion corrections
	  if(m_dccStatic) { global = m_distortionCorrection.get_corrected_position( global, m_dccStatic ); }
	  if(m_dccAverage) { global = m_distortionCorrection.get_corrected_position( global, m_dccAverage ); }
	  if(m_dccFluctuation) { global = m_distortionCorrection.get_corrected_position( global, m_dccFluctuation ); }
	 
	  // add the global positions to a vector to give to the cluster mover
	  global_raw.push_back(std::make_pair(key, global));
	}
      else
	{
	  // silicon cluster, no corrections needed
	  global_raw.push_back(std::make_pair(key, global));
	}
    }	  // end loop over clusters here
  
  // move the cluster positions back to the original readout surface
  auto global_moved = m_clusterMover.processTrack(global_raw);
  
  // loop over global positions returned by cluster mover
  for(int i=0; i<global_moved.size(); ++i)
    {
      TrkrDefs::cluskey cluskey = global_moved[i].first;
      Acts::Vector3 global = global_moved[i].second;

      Surface surf;
      TrkrDefs::subsurfkey subsurfkey;

      unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);
      if(trkrid ==  TrkrDefs::tpcId)
	{      
	  // get the new surface corresponding to this global position
	  TrkrDefs::hitsetkey tpcHitSetKey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);
	  surf = transformer.get_tpc_surface_from_coords(tpcHitSetKey,
							 global,
							 m_surfMaps,
							 m_tGeometry,
							 subsurfkey);
	  
	  if(!surf)
	    {
	      /// If the surface can't be found, we can't track with it. So 
	      /// just continue and don't add the cluster to the source links
	      if(Verbosity() > 0) std::cout << PHWHERE << "Failed to find surface for cluster " << cluskey << std::endl;
	      continue;
	    }
	  
	  if(Verbosity() > 0)
	    {
	      unsigned int side = TpcDefs::getSide(cluskey);       
	      std::cout << "      global z corrected and moved " << global[2] << " xcorr " << global[0] << " ycorr " << global[1] << " side " << side << " crossing " << track->get_crossing()
			<< " cluskey " << cluskey << " subsurfkey " << subsurfkey << std::endl;
	    }
	}
      else
	{
	  // silicon cluster, no changes possible
	  auto cluster = m_clusterContainer->findCluster(cluskey);
	  subsurfkey = cluster->getSubSurfKey();
	  surf = transformer.getSurface(cluskey, cluster, m_surfMaps);	  
	}

      // get local coordinates
      Acts::Vector2 localPos;
      Acts::Vector3 normal = surf->normal(m_tGeometry->geoContext);
      auto local = surf->globalToLocal(m_tGeometry->geoContext,
				       global * Acts::UnitConstants::cm,
				       normal);
      
      if(local.ok())
	{
	  localPos = local.value() / Acts::UnitConstants::cm;
	}
      else
	{
	  /// otherwise take the manual calculation
	  Acts::Vector3 center = surf->center(m_tGeometry->geoContext)/Acts::UnitConstants::cm;
	  double clusRadius = sqrt(global[0]*global[0] + global[1]*global[1]);
	  double clusphi = atan2(global[1], global[0]);
	  double rClusPhi = clusRadius * clusphi;
	  double surfRadius = sqrt(center(0)*center(0) + center(1)*center(1));
	  double surfPhiCenter = atan2(center[1], center[0]);
	  double surfRphiCenter = surfPhiCenter * surfRadius;
	  double surfZCenter = center[2];
	  
	  localPos(0) = rClusPhi - surfRphiCenter;
	  localPos(1) = global[2] - surfZCenter; 
	}
      
      auto cluster = m_clusterContainer->findCluster(cluskey);
      if(Verbosity() > 0)
	{
	  std::cout << " cluster local X " << cluster->getLocalX() << " cluster local Y " << cluster->getLocalY() << std::endl;
	  std::cout << " new      local X " << localPos(0) << " new       local Y " << localPos(1) << std::endl;
	}
            
      Acts::ActsVector<2> loc;
      loc[Acts::eBoundLoc0] = localPos(0) * Acts::UnitConstants::cm;
      loc[Acts::eBoundLoc1] = localPos(1) * Acts::UnitConstants::cm;
      std::array<Acts::BoundIndices,2> indices;
      indices[0] = Acts::BoundIndices::eBoundLoc0;
      indices[1] = Acts::BoundIndices::eBoundLoc1;
      Acts::ActsSymMatrix<2> cov = Acts::ActsSymMatrix<2>::Zero();
      cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = 
	cluster->getActsLocalError(0,0) * Acts::UnitConstants::cm2;
      cov(Acts::eBoundLoc0, Acts::eBoundLoc1) =
	cluster->getActsLocalError(0,1) * Acts::UnitConstants::cm2;
      cov(Acts::eBoundLoc1, Acts::eBoundLoc0) = 
	cluster->getActsLocalError(1,0) * Acts::UnitConstants::cm2;
      cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = 
	cluster->getActsLocalError(1,1) * Acts::UnitConstants::cm2;
      ActsExamples::Index index = measurements.size();
      
      SourceLink sl(surf->geometryId(), index, cluskey);
      
      Acts::Measurement<Acts::BoundIndices,2> meas(sl, indices, loc, cov);
      if(Verbosity() > 3)
	{
	  std::cout << "source link " << sl.index() << ", loc : " 
		    << loc.transpose() << std::endl 
		    << ", cov : " << cov.transpose() << std::endl
		    << " geo id " << sl.geometryId() << std::endl;
	  std::cout << "Surface : " << std::endl;
	  surf.get()->toStream(m_tGeometry->geoContext, std::cout);
	  std::cout << std::endl;
	  std::cout << "Cluster error " << cluster->getRPhiError() << " , " << cluster->getZError() << std::endl;
	  std::cout << "For key " << cluskey << " with local pos " << std::endl
		    << localPos(0) << ", " << localPos(1)
		    << std::endl;
	}
      
      sls.push_back(sl);
      measurements.push_back(meas);
    }

  return sls;
}

ActsExamples::TrackFittingAlgorithm::TrackFitterResult PHActsGSF::fitTrack(
      const std::vector<std::reference_wrapper<const SourceLink>>& sourceLinks,
      const ActsExamples::TrackParameters& seed,
      const ActsExamples::TrackFittingAlgorithm::GeneralFitterOptions& options)
{
  return (*m_fitCfg.fit)(sourceLinks, seed, options);
}


void PHActsGSF::updateTrack(const FitResult& result, SvtxTrack* track)
{


}

//____________________________________________________________________________..
int PHActsGSF::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsGSF::getNodes(PHCompositeNode* topNode)
{
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  if(!m_trackMap)
    {
      std::cout << PHWHERE << " The input track map is not available. Exiting PHActsGSF" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if(!m_clusterContainer)
    {
      std::cout << PHWHERE << "The input cluster container is not available. Exiting PHActsGSF" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << PHWHERE << "The input Acts tracking geometry is not available. Exiting PHActsGSF" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_surfMaps = findNode::getClass<ActsSurfaceMaps>(topNode, "ActsSurfaceMaps");
  if(!m_surfMaps)
    {
      std::cout << PHWHERE << "The input Acts surface maps are not available. Exiting PHActsGSF" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

   // tpc distortion corrections
  m_dccStatic = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerStatic");
  if( m_dccStatic )
    { 
      std::cout << PHWHERE << "  found static TPC distortion correction container" << std::endl; 
    }

  m_dccAverage = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerAverage");
  if( m_dccAverage )
    { 
      std::cout << PHWHERE << "  found average TPC distortion correction container" << std::endl; 
    }

  m_dccFluctuation = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerFluctuation");
  if( m_dccFluctuation )
    { 
      std::cout << PHWHERE << "  found fluctuation TPC distortion correction container" << std::endl; 
    }
  
  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if(!m_vertexMap)
    {
      std::cout << PHWHERE << "Vertex map unavailable, exiting PHActsGSF" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}
