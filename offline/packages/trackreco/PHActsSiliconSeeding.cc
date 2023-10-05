#include "PHActsSiliconSeeding.h"

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrackFitUtils.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#include <intt/CylinderGeomIntt.h>

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeed_v1.h>
#include <trackbase/TrkrCluster.h>            
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
 
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#include <Acts/Seeding/BinnedSPGroup.hpp>
#pragma GCC diagnostic pop
#include <Acts/Seeding/InternalSeed.hpp>
#include <Acts/Seeding/InternalSpacePoint.hpp>
#include <Acts/Seeding/Seed.hpp>
#include <Acts/Seeding/SeedFilter.hpp>

namespace
{
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }
}

PHActsSiliconSeeding::PHActsSiliconSeeding(const std::string& name)
  : SubsysReco(name)
{}

int PHActsSiliconSeeding::Init(PHCompositeNode */*topNode*/)
{
  configureSeeder();
  configureSPGrid();
 
  Acts::SeedFilterConfig sfCfg = configureSeedFilter();
  sfCfg = sfCfg.toInternalUnits();
  
  if(Verbosity() > 5)
    {
      printSeedConfigs(sfCfg);
    }
  // vector containing the map of z bins in the top and bottom layers
  m_bottomBinFinder = std::make_shared<const Acts::BinFinder<SpacePoint>>(
      zBinNeighborsBottom, nphineighbors);
  m_topBinFinder = std::make_shared<const Acts::BinFinder<SpacePoint>>(
      zBinNeighborsTop, nphineighbors);

  m_seedFinderCfg.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
     Acts::SeedFilter<SpacePoint>(sfCfg));

  if(m_seedAnalysis)
    {    
      m_file = new TFile(std::string(Name() + ".root").c_str(),"recreate");
      createHistograms();
    }  
  
  return Fun4AllReturnCodes::EVENT_OK;
}
int PHActsSiliconSeeding::InitRun(PHCompositeNode *topNode)
{

  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    { return Fun4AllReturnCodes::ABORTEVENT; }
  
  if(createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    { return Fun4AllReturnCodes::ABORTEVENT; }
  
  auto beginend = m_geomContainerIntt->get_begin_end();
  int i = 0;
  for(auto iter = beginend.first; iter != beginend.second; ++iter)
    {
      m_nInttLayerRadii[i] = iter->second->get_radius();
      i++;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsSiliconSeeding::process_event(PHCompositeNode *topNode)
{
  if(_n_iteration>0){
    _iteration_map = findNode::getClass<TrkrClusterIterationMapv1>(topNode, "CLUSTER_ITERATION_MAP");
    if (!_iteration_map){
      std::cerr << PHWHERE << "Cluster Iteration Map missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  auto eventTimer = std::make_unique<PHTimer>("eventTimer");
  eventTimer->stop();
  eventTimer->restart();

  if(Verbosity() > 0)
    std::cout << "Processing PHActsSiliconSeeding event "
	      << m_event << std::endl;

  std::vector<const SpacePoint*> spVec;
  auto seedVector = runSeeder(spVec);

  eventTimer->stop();
  auto seederTime = eventTimer->get_accumulated_time();
  eventTimer->restart();
  
  makeSvtxTracks(seedVector);

  eventTimer->stop();
  auto circleFitTime = eventTimer->get_accumulated_time();

  for(auto sp : spVec)
    { delete sp; }
  spVec.clear();

  if(Verbosity() > 0)
    std::cout << "Finished PHActsSiliconSeeding process_event"
	      << std::endl;

  if(Verbosity() > 0)
    {
      std::cout << "PHActsSiliconSeeding Acts seed time "
		<< seederTime << std::endl;
      std::cout << "PHActsSiliconSeeding circle fit time "
		<< circleFitTime << std::endl;
      std::cout << "PHActsSiliconSeeding total event time " 
		<< circleFitTime + seederTime  << std::endl;
    }

  m_event++;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsSiliconSeeding::End(PHCompositeNode */*topNode*/)
{
  if(m_seedAnalysis)
    {
      writeHistograms();
    }

  if(Verbosity() > 1)
    {
      std::cout << "There were " << m_nBadInitialFits 
		<< " bad initial circle fits" << std::endl;
      std::cout << "There were " << m_nBadUpdates 
		<< " bad second circle fits" << std::endl;
    }
     
  return Fun4AllReturnCodes::EVENT_OK;
}

GridSeeds PHActsSiliconSeeding::runSeeder(std::vector<const SpacePoint*>& spVec)
{
  
  Acts::SeedFinder<SpacePoint> seedFinder(m_seedFinderCfg);
  
  /// Covariance converter functor needed by seed finder
  auto covConverter = 
    [=](const SpacePoint& sp, float zAlign, float rAlign, float sigmaError)
    -> std::pair<Acts::Vector3, Acts::Vector2> { 
       Acts::Vector3 position{sp.x(), sp.y(), sp.z()};
       Acts::Vector2 cov;
       cov[0] = (sp.m_varianceR + rAlign*rAlign) * sigmaError;
       cov[1] = (sp.m_varianceZ + zAlign*zAlign) * sigmaError;
       return std::make_pair(position, cov);
  };

  Acts::Extent rRangeSPExtent;

  spVec = getMvtxSpacePoints(rRangeSPExtent);

  if(m_seedAnalysis)
    { h_nInputMeas->Fill(spVec.size()); }
  
  auto grid = 
    Acts::SpacePointGridCreator::createGrid<SpacePoint>(m_gridCfg,
							m_gridOptions);
  
  auto spGroup = Acts::BinnedSPGroup<SpacePoint>(spVec.begin(),
						 spVec.end(),
						 covConverter,
						 m_bottomBinFinder,
						 m_topBinFinder,
						 std::move(grid), 
						 rRangeSPExtent,
						 m_seedFinderCfg,
						 m_seedFinderOptions);

  /// variable middle SP radial region of interest
  const Acts::Range1D<float> rMiddleSPRange(
	 std::floor(rRangeSPExtent.min(Acts::binR) / 2) * 2 + 1.5, 
	 std::floor(rRangeSPExtent.max(Acts::binR) / 2) * 2 - 1.5);

  GridSeeds seedVector;
  SeedContainer seeds;
  seeds.clear();
  decltype(seedFinder)::SeedingState state;
  state.spacePointData.resize(spVec.size(),
			      m_seedFinderCfg.useDetailedDoubleMeasurementInfo);
  for(const auto [bottom, middle, top] : spGroup)
    {
    
      seedFinder.createSeedsForGroup(m_seedFinderOptions,
				     state, spGroup.grid(),
				     std::back_inserter(seeds),
				     bottom,
				     middle,
				     top,
				     rMiddleSPRange);

    }
  
  seedVector.push_back(seeds);

  return seedVector;
}

void PHActsSiliconSeeding::makeSvtxTracks(GridSeeds& seedVector)
{
  int numSeeds = 0;
  int numGoodSeeds = 0;
  
  /// Loop over grid volumes
  for(auto& seeds : seedVector)
    {
      /// Loop over actual seeds in this grid volume
      for(auto& seed : seeds)
	{
	  if(Verbosity() > 1) {
	    std::cout << "Seed " << numSeeds << " has "
		      << seed.sp().size() << " measurements " 
		      << std::endl;
	  }

	  numSeeds++;

	  std::vector<TrkrDefs::cluskey> cluster_keys;

	  std::vector<Acts::Vector3> globalPositions;

	  std::map<TrkrDefs::cluskey, Acts::Vector3> positions;
	  auto trackSeed = std::make_unique<TrackSeed_v1>();

	  for(auto& spacePoint : seed.sp())
	    {
	      const auto& cluskey = spacePoint->Id();
	      cluster_keys.push_back(cluskey);
	  
	      trackSeed->insert_cluster_key(cluskey);
	      auto globalPosition = m_tGeometry->getGlobalPosition(
                     cluskey,
		     m_clusterMap->findCluster(cluskey));
	      globalPositions.push_back(globalPosition);      

	      positions.insert(std::make_pair(cluskey, globalPosition));
	      if(Verbosity() > 1) {
		std::cout << "Adding cluster with x,y "
			  << spacePoint->x() <<", " << spacePoint->y()
			  << " mm in detector " 
			  << TrkrDefs::getTrkrId(cluskey)
			  << std::endl;
	      }
	    }
	  
	  double z = seed.z() / Acts::UnitConstants::cm;
	  
	  auto fitTimer = std::make_unique<PHTimer>("trackfitTimer");
	  fitTimer->stop();
	  fitTimer->restart();

	  trackSeed->circleFitByTaubin(positions, 0, 8);
	  if(fabs(trackSeed->get_x()) > m_maxSeedPCA || 
	     fabs(trackSeed->get_y()) > m_maxSeedPCA)
	    { 
	      if(Verbosity() > 1)
		{ 
		  std::cout << "Large PCA seed " <<std::endl;
		  trackSeed->identify();
		}
	      m_nBadInitialFits++;
	      continue;
	    }
	 
	  trackSeed->lineFit(positions, 0, 8);
	  z = trackSeed->get_Z0();

	  fitTimer->stop();
	  auto circlefittime = fitTimer->get_accumulated_time();
	  fitTimer->restart();
	  
	  /// Project to INTT and find matches
	  int mvtxsize = globalPositions.size();
	  auto additionalClusters = findInttMatches(globalPositions, *trackSeed);

	  /// Add possible matches to cluster list to be parsed when
	  /// Svtx tracks are made
	  for(int newkey = 0; newkey < additionalClusters.size(); newkey++)
	    {
	      trackSeed->insert_cluster_key(additionalClusters[newkey]);
	      positions.insert(std::make_pair(additionalClusters[newkey],globalPositions[mvtxsize+newkey]));
	      
	      if(Verbosity() > 1)
		{ std::cout << "adding additional intt key " << additionalClusters[newkey] << std::endl; }
	    }

	  fitTimer->stop();
	  auto addClusters = fitTimer->get_accumulated_time();
	  fitTimer->restart();

	  //! Circle fit again to take advantage of INTT lever arm
	  trackSeed->circleFitByTaubin(positions, 0, 8);

	  if(Verbosity() > 0)
	    { std::cout << "find intt clusters time " << addClusters << std::endl; }

	  numGoodSeeds++;

	  /// The Acts z projection has better resolution than the circle fit
	  trackSeed->set_Z0(z);

	  if(Verbosity() > 1)
	    { 
	      std::cout << "Silicon seed id " << m_seedContainer->size() << std::endl;
	      std::cout << "seed phi, theta, eta : " 
			<< trackSeed->get_phi(m_clusterMap, m_tGeometry) << ", " << trackSeed->get_theta()
			<< ", " << trackSeed->get_eta() << std::endl;
	      trackSeed->identify(); 
	    }

	  m_seedContainer->insert(trackSeed.get());

	  fitTimer->stop();
	  auto svtxtracktime = fitTimer->get_accumulated_time();
	  if(Verbosity() > 0)
	    {
	      std::cout << "Intt fit time " << circlefittime << " and svtx time "
			<< svtxtracktime << std::endl;
	    }
	}
    }

  if(m_seedAnalysis)
    {
      h_nSeeds->Fill(numGoodSeeds);
      h_nActsSeeds->Fill(numSeeds);
    }
  if(Verbosity() > 1)
    {
      std::cout << "Total number of seeds found in " 
		<< seedVector.size() << " volume regions gives " 
		<< numSeeds << " Acts seeds " << std::endl;
    }

  return;
  
}
 
std::vector<TrkrDefs::cluskey> PHActsSiliconSeeding::findInttMatches(
		               std::vector<Acts::Vector3>& clusters,
			       TrackSeed& seed)
{
  const float R = fabs(1. / seed.get_qOverR());
  const float X0 = seed.get_X0();
  const float Y0 = seed.get_Y0();
  const float B = seed.get_Z0();
  const float m = seed.get_slope();

  double xProj[m_nInttLayers];
  double yProj[m_nInttLayers];
  double zProj[m_nInttLayers];
  
  /// Diagnostic 
  if(m_seedAnalysis)
    {
      for(const auto& glob : clusters)
	{
	  h_hits->Fill(glob(0), glob(1));
	  h_zhits->Fill(glob(2),
			std::sqrt(square(glob(0)) + square(glob(1))));
	  h_projHits->Fill(glob(0), glob(1));
	  h_zprojHits->Fill(glob(2),
			    sqrt(square(glob(0)) + square(glob(1))));
	}
    }

  /// Project the seed to the INTT to find matches
  for(int layer = 0; layer < m_nInttLayers; ++layer)
    {
      auto cci = TrackFitUtils::circle_circle_intersection(
                 m_nInttLayerRadii[layer],
		 R, X0, Y0);
      double xplus = std::get<0>(cci);
      double yplus = std::get<1>(cci);
      double xminus = std::get<2>(cci);
      double yminus = std::get<3>(cci);
  
      /// If there are no real solutions to the intersection, skip
      if(std::isnan(xplus))
	{
	  if(Verbosity() > 2)
	    {
	      std::cout << "Circle intersection calc failed, skipping" 
			<< std::endl;
	      std::cout << "layer radius " << m_nInttLayerRadii[layer] 
			<< " and circ rad " << R << " with center " << X0 
			<< ", " << Y0 << std::endl;
	    }
	  continue;
	}
      
      /// Figure out which solution is correct based on the position 
      /// of the last layer in the mvtx seed
      const auto lastclusglob = clusters.at(clusters.size() - 1);
      const double lastClusPhi = atan2(lastclusglob(1), lastclusglob(0));
      const double plusPhi = atan2(yplus, xplus);
      const double minusPhi = atan2(yminus, xminus);
  
      if(fabs(lastClusPhi - plusPhi) < fabs(lastClusPhi - minusPhi))
	{
	  xProj[layer] = xplus;
	  yProj[layer] = yplus;
	}
      else
	{
	  xProj[layer] = xminus;
	  yProj[layer] = yminus;
	}
      
      zProj[layer] = m * m_nInttLayerRadii[layer] + B;
      
      if(m_seedAnalysis) {
	h_projHits->Fill(xProj[layer], yProj[layer]);
	h_zprojHits->Fill(zProj[layer], std::sqrt(square(xProj[layer]) + 
					     square(yProj[layer])));
      }

      if(Verbosity() > 2)
	{
	  std::cout << "Projected point is : " << xProj[layer] << ", "
		    << yProj[layer] << ", " << zProj[layer] << std::endl;
	}
    }
  
  return matchInttClusters(clusters, xProj, yProj, zProj);
}

std::vector<TrkrDefs::cluskey> PHActsSiliconSeeding::matchInttClusters(
  std::vector<Acts::Vector3>& clusters,
  const double xProj[],
  const double yProj[],
  const double zProj[])
{
  std::vector<TrkrDefs::cluskey> matchedClusters;
  TrkrDefs::cluskey matchedClusterLay0;
  Acts::Vector3 matchedGlobPosLay0;
  TrkrDefs::cluskey matchedClusterLay1;
  Acts::Vector3 matchedGlobPosLay1;
  float minResidLay0 = std::numeric_limits<float>::max();
  float minResidLay1 = std::numeric_limits<float>::max();
  
  for(int inttlayer = 0; inttlayer < m_nInttLayers; inttlayer++)
    {
      const double projR = std::sqrt(square(xProj[inttlayer]) + square(yProj[inttlayer]));
      const double projPhi = std::atan2(yProj[inttlayer], xProj[inttlayer]);
      const double projRphi = projR * projPhi;

      for( const auto& hitsetkey : m_clusterMap->getHitSetKeys(TrkrDefs::TrkrId::inttId, inttlayer+3))
      {
	  double ladderLocation[3] = {0.,0.,0.};

	  // Add three to skip the mvtx layers for comparison
	  // to projections
	  auto layerGeom = dynamic_cast<CylinderGeomIntt*>
	    (m_geomContainerIntt->GetLayerGeom(inttlayer+3));
	  
	  auto surf = m_tGeometry->maps().getSiliconSurface(hitsetkey);
	  layerGeom->find_segment_center(surf, m_tGeometry, ladderLocation);
        
	  const double ladderphi = atan2(ladderLocation[1], ladderLocation[0]) + layerGeom->get_strip_phi_tilt();
	  const auto stripZSpacing = layerGeom->get_strip_z_spacing();
	  float dphi = ladderphi - projPhi;
	  if(dphi > M_PI)
	    { dphi -= 2. * M_PI; }
	  else if (dphi < -1 * M_PI)
	    { dphi += 2. * M_PI; }
	  
	  /// Check that the projection is within some reasonable amount of the segment
	  /// to reject e.g. looking at segments in the opposite hemisphere. This is about
	  /// the size of one intt segment (256 * 80 micron strips in a segment)
	  if(fabs(dphi) > 0.2)
	    { continue; }

	  TVector3 projectionLocal(0,0,0);
	  TVector3 projectionGlobal(xProj[inttlayer],yProj[inttlayer],zProj[inttlayer]);
	 
	  projectionLocal = layerGeom->get_local_from_world_coords(surf, 
								   m_tGeometry,
								   projectionGlobal);
        
	  auto range = m_clusterMap->getClusters(hitsetkey);	
	  for(auto clusIter = range.first; clusIter != range.second; ++clusIter )
	    {
	      const auto cluskey = clusIter->first;
	      const auto cluster = clusIter->second;
	      /// Diagnostic
	      if(m_seedAnalysis)
		{ 
		  const auto globalP = m_tGeometry->getGlobalPosition(
                        cluskey, cluster);
		  h_nInttProj->Fill(projectionLocal[1] - cluster->getLocalX(),
				    projectionLocal[2] - cluster->getLocalY()); 
		  h_hits->Fill(globalP(0), globalP(1));
		  h_zhits->Fill(globalP(2),
        std::sqrt(square(globalP(0))+square(globalP(1))));
		  
		  h_resids->Fill(zProj[inttlayer] - cluster->getLocalY(),
				 projRphi - cluster->getLocalX());
		}
	     
	      /// Z strip spacing is the entire strip, so because we use fabs
	      /// we divide by two
	      float rphiresid = fabs(projectionLocal[1] - cluster->getLocalX());
	      float zresid = fabs(projectionLocal[2] - cluster->getLocalY());
	      if(rphiresid < m_rPhiSearchWin and
		 zresid < stripZSpacing / 2.)
		{
		  
		  const auto globalPos = m_tGeometry->getGlobalPosition(
                    cluskey, cluster);
		  
		  if(inttlayer<2 && rphiresid < minResidLay0)
		    {
		      matchedClusterLay0 = cluskey;
		      /// Cache INTT global positions with seed
		      matchedGlobPosLay0 = globalPos;
		      minResidLay0 = rphiresid;
		    }
		  if(inttlayer>1 && rphiresid < minResidLay1)
		    {
		      matchedClusterLay1 = cluskey;
		      matchedGlobPosLay1 = globalPos;
		      minResidLay1 = rphiresid;
		    }
		    
		 		  	      
		  if(Verbosity() > 4)
		    {
		      std::cout << "Matched INTT cluster with cluskey " << cluskey 
				<< " and position " << globalPos.transpose() 
				<< std::endl << " with projections rphi "
				<< projRphi << " and inttclus rphi " << cluster->getLocalX()
				<< " and proj z " << zProj[inttlayer] << " and inttclus z "
				<< cluster->getLocalY() << " in layer " << inttlayer 
				<< " with search windows " << m_rPhiSearchWin 
				<< " in rphi and strip z spacing " << stripZSpacing 
				<< std::endl;
		    }
		}
	    }
	}  
    }

  if(minResidLay0 < std::numeric_limits<float>::max())
    {
      matchedClusters.push_back(matchedClusterLay0);
      clusters.push_back(matchedGlobPosLay0);
    }
  if(minResidLay1 < std::numeric_limits<float>::max())
    {
      matchedClusters.push_back(matchedClusterLay1);
      clusters.push_back(matchedGlobPosLay1);
    }

  if(m_seedAnalysis) {
    h_nMatchedClusters->Fill(matchedClusters.size());
  }

  return matchedClusters;
}

SpacePointPtr PHActsSiliconSeeding::makeSpacePoint(
  const Surface& surf,
  const TrkrDefs::cluskey key,
  TrkrCluster* clus)
{
  Acts::Vector2 localPos(clus->getLocalX() * Acts::UnitConstants::cm, 
			 clus->getLocalY() * Acts::UnitConstants::cm);
  Acts::Vector3 globalPos(0,0,0);
  Acts::Vector3 mom(1,1,1);
  globalPos = surf->localToGlobal(m_tGeometry->geometry().getGeoContext(),
				  localPos, mom);

  Acts::SymMatrix2 localCov = Acts::SymMatrix2::Zero();
  localCov(0,0) = pow(clus->getRPhiError(),2) * Acts::UnitConstants::cm2;
  localCov(1,1) = pow(clus->getZError(),2) * Acts::UnitConstants::cm2;
  
    
  float x = globalPos.x();
  float y = globalPos.y();
  float z = globalPos.z();
  float r = std::sqrt(x * x + y * y);

  /// The space point requires only the variance of the transverse and
  /// longitudinal position. Reduce computations by transforming the
  /// covariance directly from local to r/z.
  ///
  /// compute Jacobian from global coordinates to r/z
  ///
  ///         r = sqrt(x² + y²)
  /// dr/d{x,y} = (1 / sqrt(x² + y²)) * 2 * {x,y}
  ///             = 2 * {x,y} / r
  ///       dz/dz = 1 

  Acts::RotationMatrix3 rotLocalToGlobal =
    surf->referenceFrame(m_tGeometry->geometry().getGeoContext(), globalPos, mom);
  auto scale = 2 / std::hypot(x,y);
  Acts::ActsMatrix<2, 3> jacXyzToRhoZ = Acts::ActsMatrix<2, 3>::Zero();
  jacXyzToRhoZ(0, Acts::ePos0) = scale * x;
  jacXyzToRhoZ(0, Acts::ePos1) = scale * y;
  jacXyzToRhoZ(1, Acts::ePos2) = 1;
  // compute Jacobian from local coordinates to rho/z
  Acts::ActsMatrix<2, 2> jac =
    jacXyzToRhoZ *
    rotLocalToGlobal.block<3, 2>(Acts::ePos0, Acts::ePos0);
  // compute rho/z variance
  Acts::ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

  /*
   * From Acts v17 to v19 the scattering uncertainty value allowed was changed.
   * This led to a decrease in efficiency. To offset this, we scale the 
   * uncertainties by a tuned factor that gives the v17 performance
   * Track reconstruction is an art as much as it is a science...
   */
  SpacePointPtr spPtr(new SpacePoint{key, x, y, z, r,  surf->geometryId(), var[0]*m_uncfactor, var[1]*m_uncfactor});

  if(Verbosity() > 2)
    std::cout << "Space point has " 
	      << x << ", " << y << ", " << z << " with local coords "
	      << localPos.transpose() 
	      << " with rphi/z variances " << localCov(0,0) 
	      << ", " << localCov(1,1) << " and rotated variances "
	      << var[0] << ", " << var[1] 
	      << " and cluster key "
	      << key << " and geo id "
	      << surf->geometryId() << std::endl;
  
  return spPtr;

}

std::vector<const SpacePoint*> PHActsSiliconSeeding::getMvtxSpacePoints(Acts::Extent& rRangeSPExtent)
{
  std::vector<const SpacePoint*> spVec;
  unsigned int numSiliconHits = 0;
 
  for(const auto& hitsetkey:m_clusterMap->getHitSetKeys(TrkrDefs::TrkrId::mvtxId))
    {
      auto range = m_clusterMap->getClusters(hitsetkey);
      for( auto clusIter = range.first; clusIter != range.second; ++clusIter )
	{
	  const auto cluskey = clusIter->first;
	  const auto cluster = clusIter->second;
	  if(_iteration_map != NULL && _n_iteration >0 ){
	    //	    std::cout << "map exists entries: " << _iteration_map->size() << std::endl;
	    if(_iteration_map->getIteration(cluskey)>0){
	      continue; // skip hits used in a previous iteration
	    }
	  }

	  const auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);
	  const auto surface = m_tGeometry->maps().getSiliconSurface(hitsetkey);
	  if(!surface)
	    { continue; }

	  auto sp = makeSpacePoint(surface, cluskey, cluster).release();
	  spVec.push_back(sp);
	  rRangeSPExtent.extend({sp->x(), sp->y(), sp->z()});
	  numSiliconHits++;
	}
    }

  if(m_seedAnalysis)
    { h_nInputMvtxMeas->Fill(numSiliconHits); }

  if(Verbosity() > 1)
    std::cout << "Total number of silicon hits to seed find with is "
	      << numSiliconHits << std::endl;


  return spVec;
}

void PHActsSiliconSeeding::configureSPGrid()
{

  m_gridCfg.minPt = m_minSeedPt / m_gridFactor;
  m_gridCfg.rMax = m_rMax;
  m_gridCfg.zMax = m_zMax;
  m_gridCfg.zMin = m_zMin;
  m_gridCfg.deltaRMax = m_deltaRMax;
  m_gridCfg.cotThetaMax = m_cotThetaMax;
  m_gridCfg.impactMax = m_impactMax;
  m_gridCfg.phiBinDeflectionCoverage = m_numPhiNeighbors;

  m_gridOptions.bFieldInZ = m_bField;
  m_gridCfg = m_gridCfg.toInternalUnits();
  m_gridOptions = m_gridOptions.toInternalUnits();
}

Acts::SeedFilterConfig PHActsSiliconSeeding::configureSeedFilter()
{
  Acts::SeedFilterConfig config;
  config.maxSeedsPerSpM = m_maxSeedsPerSpM;
  return config;
}

void PHActsSiliconSeeding::configureSeeder()
{
  /// these are default values that used to be set in Acts
  m_seedFinderCfg.deltaRMinTopSP = 5 * Acts::UnitConstants::mm;
  m_seedFinderCfg.deltaRMaxTopSP = 270 * Acts::UnitConstants::mm;
  m_seedFinderCfg.deltaRMinBottomSP = 5 * Acts::UnitConstants::mm;
  m_seedFinderCfg.deltaRMaxBottomSP = 270 * Acts::UnitConstants::mm;

  /// Limiting location of measurements (e.g. detector constraints)
  m_seedFinderCfg.rMax = m_rMax;
  m_seedFinderCfg.rMin = m_rMin;
  m_seedFinderCfg.zMin = m_zMin;
  m_seedFinderCfg.zMax = m_zMax;

  /// Min/max distance between two measurements in one seed
  m_seedFinderCfg.deltaRMin = m_deltaRMin;
  m_seedFinderCfg.deltaRMax = m_deltaRMax;

  /// Limiting collision region in z
  m_seedFinderCfg.collisionRegionMin = -300. * Acts::UnitConstants::mm;
  m_seedFinderCfg.collisionRegionMax = 300. * Acts::UnitConstants::mm;
  m_seedFinderCfg.sigmaScattering = m_sigmaScattering;
  m_seedFinderCfg.maxSeedsPerSpM = m_maxSeedsPerSpM;
  m_seedFinderCfg.cotThetaMax = m_cotThetaMax;
  m_seedFinderCfg.minPt = m_minSeedPt;
  m_seedFinderOptions.bFieldInZ = m_bField;

  /// Average radiation length traversed per seed
  m_seedFinderCfg.radLengthPerSeed = 0.05;

  /// Maximum impact parameter must be smaller than rMin
  m_seedFinderCfg.impactMax = m_impactMax;

  m_seedFinderCfg.rMinMiddle = 3. * Acts::UnitConstants::cm;
  m_seedFinderCfg.rMaxMiddle = 3.7 * Acts::UnitConstants::cm;

  /// Configurations for dealing with misalignment
  m_seedFinderCfg.zAlign = m_zalign;
  m_seedFinderCfg.rAlign = m_ralign;
  m_seedFinderCfg.toleranceParam = m_tolerance;
  m_seedFinderCfg.maxPtScattering = m_maxPtScattering;
  m_seedFinderCfg.sigmaError = m_sigmaError;
  m_seedFinderCfg.helixcut = m_helixcut;

  m_seedFinderCfg = 
    m_seedFinderCfg.toInternalUnits().calculateDerivedQuantities();
  m_seedFinderOptions = m_seedFinderOptions.toInternalUnits().calculateDerivedQuantities(m_seedFinderCfg);
  
}

int PHActsSiliconSeeding::getNodes(PHCompositeNode *topNode)
{

  m_geomContainerIntt = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if(!m_geomContainerIntt)
    {
      std::cout << PHWHERE << "CYLINDERGEOM_INTT node not found on node tree"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }


  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if(!m_tGeometry)
    {
      std::cout << PHWHERE << "No ActsGeometry on node tree. Bailing."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }


  if(m_useTruthClusters)
    m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode,
							    "TRKR_CLUSTER_TRUTH");
  else
    m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode,
							     "TRKR_CLUSTER");

  if(!m_clusterMap)
    {
      std::cout << PHWHERE << "No cluster container on the node tree. Bailing."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  return Fun4AllReturnCodes::EVENT_OK;

}
int PHActsSiliconSeeding::createNodes(PHCompositeNode *topNode)
{
 
  PHNodeIterator iter(topNode);
  
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHActsSiliconSeeding::createNodes");
  }

  PHNodeIterator dstIter(dstNode);
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_seedContainer = findNode::getClass<TrackSeedContainer>(topNode, _track_map_name);
  if(!m_seedContainer)
    {
      m_seedContainer = new TrackSeedContainer_v1;
      PHIODataNode<PHObject> *trackNode = 
	new PHIODataNode<PHObject>(m_seedContainer, _track_map_name, "PHObject");
      svtxNode->addNode(trackNode);

    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsSiliconSeeding::writeHistograms()
{
  m_file->cd();
  h_nInttProj->Write();
  h_nMatchedClusters->Write();
  h_nMvtxHits->Write();
  h_nSeeds->Write();
  h_nActsSeeds->Write();
  h_nTotSeeds->Write();
  h_nInttHits->Write();
  h_nInputMeas->Write();
  h_nHits->Write();
  h_nInputMvtxMeas->Write();
  h_nInputInttMeas->Write();
  h_hits->Write();
  h_zhits->Write();
  h_projHits->Write();
  h_zprojHits->Write();
  h_resids->Write();
  m_file->Write();      
  m_file->Close(); 
}

void PHActsSiliconSeeding::createHistograms()
{
  h_nMatchedClusters = new TH1F("nMatchedClusters",";N_{matches}", 50,0,50);
  h_nInttProj = new TH2F("nInttProj",";l_{0}^{proj}-l_{0}^{clus} [cm]; l_{1}^{proj}-l_{1}^{clus} [cm]",
			 10000,-10,10,10000,-50,50);
  h_nMvtxHits = new TH1I("nMvtxHits",";N_{MVTX}",6,0,6);
  h_nInttHits = new TH1I("nInttHits",";N_{INTT}",80,0,80);
  h_nHits = new TH2I("nHits",";N_{MVTX};N_{INTT}",10,0,10,
		     80,0,80);
  h_nActsSeeds = new TH1I("nActsSeeds",";N_{Seeds}",400,0,400);
  h_nSeeds = new TH1I("nActsGoodSeeds",";N_{Seeds}",400,0,400);
  h_nTotSeeds = new TH1I("nTotSeeds",";N_{Seeds}",500,0,500);
  h_nInputMeas = new TH1I("nInputMeas",";N_{Meas}",2000,0,2000);
  h_nInputMvtxMeas = new TH1I("nInputMvtxMeas",";N_{meas}^{mvtx}",
			      150,0,150);
  h_nInputInttMeas = new TH1I("nInputInttMeas",";N_{meas}^{intt}",
			      150,0,150);
  h_hits = new TH2F("hits",";x [cm]; y [cm]",1000,-20,20,
		    1000,-20,20);
  h_zhits = new TH2F("zhits",";z [cm]; r [cm]",1000,-30,30,
		     1000,-30,30);
  h_projHits = new TH2F("projhits",";x [cm]; y [cm]",1000,-20,20,
			1000,-20,20);
  h_zprojHits = new TH2F("zprojhits",";z [cm]; r [cm]",1000,-30,30,
			 1000,-30,30);
  h_resids = new TH2F("resids",";z_{resid} [cm]; rphi_{resid} [cm]",
		      100,-1,1,100,-1,1);
}

double PHActsSiliconSeeding::normPhi2Pi(const double phi)
{
  double returnPhi = phi;
  if(returnPhi < 0)
    returnPhi += 2 * M_PI;
  return returnPhi;
}


void PHActsSiliconSeeding::largeGridSpacing(const bool spacing)
{
  if(!spacing)
    {
      m_gridFactor = 1.;
      m_rMax = 50.;
      m_cotThetaMax = 1.335647;
      m_maxSeedPCA = 0.1;
    }

}

void PHActsSiliconSeeding::printSeedConfigs(Acts::SeedFilterConfig& sfconfig)
{

  // helper function for when updating acts to ensure all seeding parameters
  // are actually the same as before

  std::cout << " --------------- SeedFilterConfig ------------------ " << std::endl;
  std::cout << "deltaInvHelixDiameter = " 
	    << sfconfig.deltaInvHelixDiameter << std::endl 
	    << "impactWeightFactor = " << sfconfig.impactWeightFactor 
	    << std::endl << "zOriginWeightFactor = " 
	    << sfconfig.zOriginWeightFactor << std::endl 
	    << "compatSeedWeight = " << sfconfig.compatSeedWeight 
	    << std::endl << "deltaRMin = " << sfconfig.deltaRMin 
	    << std::endl << "maxSeedsPerSpM = " 
	    << sfconfig.maxSeedsPerSpM << std::endl 
	    << "compatSeedLimit = " << sfconfig.compatSeedLimit 
	    << std::endl << "seedWeightIncrement = " 
	    << sfconfig.seedWeightIncrement << std::endl 
	    << "numSeedIncrement = " << sfconfig.numSeedIncrement 
	    << std::endl;

  std::cout << " --------------- SeedFinderConfig ------------------ " << std::endl;

  std::cout << "deltaRMinTopSp = " << m_seedFinderCfg.deltaRMinTopSP 
	    << std::endl << "deltaRMaxTopSP = " << m_seedFinderCfg.deltaRMaxTopSP 
	    << std::endl << "deltaRMinBottomSP = " << m_seedFinderCfg.deltaRMinBottomSP 
	    << std::endl << "deltaRMaxBottomSP = " << m_seedFinderCfg.deltaRMaxBottomSP 
	    << std::endl << "minpt = " << m_seedFinderCfg.minPt
	    << std::endl << "cotThetaMax = " << m_seedFinderCfg.cotThetaMax
	    << std::endl << "deltaRMin = " << m_seedFinderCfg.deltaRMin
	    << std::endl << "deltaRMax = " << m_seedFinderCfg.deltaRMax
	    << std::endl << "binSizeR = " << m_seedFinderCfg.binSizeR
	    << std::endl << "deltaRMiddleMinSPRange = " << m_seedFinderCfg.deltaRMiddleMinSPRange
	    << std::endl << "deltaRMiddleMaxSPRange = " << m_seedFinderCfg.deltaRMiddleMaxSPRange
	    << std::endl << "rMinMiddle = " << m_seedFinderCfg.rMinMiddle
	    << std::endl << "rMaxMiddle = " << m_seedFinderCfg.rMaxMiddle
	    << std::endl << "deltaZMax = " << m_seedFinderCfg.deltaZMax
	    << std::endl << "rMax = " << m_seedFinderCfg.rMax
	    << std::endl << "rMin = " << m_seedFinderCfg.rMin
	    << std::endl << "zMin = " << m_seedFinderCfg.zMin
	    << std::endl << "zMax = " << m_seedFinderCfg.zMax
	    << std::endl << "collisionRegionMin = " << m_seedFinderCfg.collisionRegionMin
	    << std::endl << "collisionRegionMax = " << m_seedFinderCfg.collisionRegionMax
	    << std::endl << "sigmaScattering = " << m_seedFinderCfg.sigmaScattering
	    << std::endl << "maxSeedsPerSpM = " << m_seedFinderCfg.maxSeedsPerSpM
	    << std::endl << "bFieldInZ = " << m_seedFinderOptions.bFieldInZ 
	    << std::endl << "radLengthPerSeed = " << m_seedFinderCfg.radLengthPerSeed
	    << std::endl << "impactMax = " << m_seedFinderCfg.impactMax
	    << std::endl << "zAlign = " << m_seedFinderCfg.zAlign
	    << std::endl << "rAlign = " << m_seedFinderCfg.rAlign
	    << std::endl << "toleranceParam = " << m_seedFinderCfg.toleranceParam
	    << std::endl << "maxPtScattering = " << m_seedFinderCfg.maxPtScattering
	    << std::endl << "sigmaError = " << m_seedFinderCfg.sigmaError
	    << std::endl << "helixcut = " << m_seedFinderCfg.helixcut 
	    << std::endl;

 std::cout << " --------------- SeedFinderOptions ------------------ " << std::endl;
 std::cout << "beamPos = " << m_seedFinderOptions.beamPos.transpose()
	   << std::endl << "bFieldInZ = " << m_seedFinderOptions.bFieldInZ
	   << std::endl << "pTPerHelixRadius = " << m_seedFinderOptions.pTPerHelixRadius
	   << std::endl << "minHelixDiameter2 = " << m_seedFinderOptions.minHelixDiameter2
	   << std::endl << "pT2perRadius = " << m_seedFinderOptions.pT2perRadius
	   << std::endl << "sigmapT2perRadius = " << m_seedFinderOptions.sigmapT2perRadius
	   << std::endl;


   std::cout << " --------------- SpacePointGridConfig ------------------ " << std::endl;

   std::cout << "minpt = " << m_gridCfg.minPt << std::endl
	     << "rMax = " << m_gridCfg.rMax << std::endl
	     << "zMax = " << m_gridCfg.zMax << std::endl
	     << "zMin = " << m_gridCfg.zMin << std::endl
	     << "deltaRMax = " << m_gridCfg.deltaRMax << std::endl
	     << "cotThetaMax = " << m_gridCfg.cotThetaMax << std::endl
	     << "impactMax = " << m_gridCfg.impactMax << std::endl
	     << "phiBinDeflectionCoverage = " << m_gridCfg.phiBinDeflectionCoverage << std::endl
	     << "bFieldInZ = " << m_gridOptions.bFieldInZ << std::endl;

   
}
