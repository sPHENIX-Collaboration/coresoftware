
#include "PHActsKDTreeSeeding.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeed_v1.h>

#include <intt/CylinderGeomIntt.h>

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase/TrackFitUtils.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrCluster.h>            
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>

#include <Acts/Seeding/InternalSeed.hpp>
#include <Acts/Seeding/SeedFilterConfig.hpp>
#include <Acts/Seeding/SeedFinderOrthogonalConfig.hpp>
#include <Acts/Seeding/SpacePointGrid.hpp>
#include <Acts/Utilities/KDTree.hpp>
#include <Acts/Seeding/Seed.hpp>
#include <Acts/Seeding/SeedFilter.hpp>
#include <Acts/Seeding/SeedFinderOrthogonal.hpp>

#include <optional>

namespace
{
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }
}

//____________________________________________________________________________..
PHActsKDTreeSeeding::PHActsKDTreeSeeding(const std::string &name):
 SubsysReco(name)
{
}

//____________________________________________________________________________..
PHActsKDTreeSeeding::~PHActsKDTreeSeeding()
{
}

//____________________________________________________________________________..
int PHActsKDTreeSeeding::Init(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHActsKDTreeSeeding::InitRun(PHCompositeNode* topNode)
{
  int ret = getNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK)
    return ret;
  ret = createNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK)
    return ret;

  m_seedFinderConfig = configureSeedFinder();

  auto beginend = m_geomContainerIntt->get_begin_end();
  int i = 0;
  for(auto iter = beginend.first; iter != beginend.second; ++iter)
    {
      m_nInttLayerRadii[i] = iter->second->get_radius();
      i++;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHActsKDTreeSeeding::process_event(PHCompositeNode* topNode)
{
  if(m_nIteration>0){
    m_iterationMap = findNode::getClass<TrkrClusterIterationMapv1>(topNode, "CLUSTER_ITERATION_MAP");
    if (!m_iterationMap){
      std::cerr << PHWHERE << "Cluster Iteration Map missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  auto seeds = runSeeder();
  
  fillTrackSeedContainer(seeds);

  return Fun4AllReturnCodes::EVENT_OK;
}

SeedContainer PHActsKDTreeSeeding::runSeeder()
{
  Acts::SeedFinderOrthogonal<SpacePoint> finder(m_seedFinderConfig);

  auto spacePoints = getMvtxSpacePoints();

  /// Call acts seeding algo
  SeedContainer seeds = finder.createSeeds(spacePoints);
  if(Verbosity() > 1)
    {
      std::cout << "Acts::OrthogonalSeeder found " << seeds.size() 
		<< " seeds " << std::endl;
    }

  return seeds;
}

void PHActsKDTreeSeeding::fillTrackSeedContainer(SeedContainer& seeds)
{

  for(auto& seed : seeds)
    {
      auto siseed = std::make_unique<TrackSeed_v1>();
      std::map<TrkrDefs::cluskey, Acts::Vector3> positions;

      for (auto& spptr : seed.sp()) 
	{
	  auto ckey = spptr->Id();
	  siseed->insert_cluster_key(ckey);
	   auto globalPosition = m_tGeometry->getGlobalPosition(
                ckey,
		m_clusterMap->findCluster(ckey));
	   positions.insert(std::make_pair(ckey, globalPosition));
	}
      
      siseed->circleFitByTaubin(positions,0,8);
      siseed->lineFit(positions,0,8);

      /// Project to INTT and find matches to add to positions
      findInttMatches(positions, *siseed);

      for(const auto& [key, pos] : positions)
	{
	  if(TrkrDefs::getTrkrId(key) == TrkrDefs::TrkrId::inttId)
	    {
	      siseed->insert_cluster_key(key);
	    }
	}

      /// The acts seed has better z resolution than the circle fit
      siseed->set_Z0(seed.z() / Acts::UnitConstants::cm);
      if(Verbosity() > 2)
	{
	  std::cout << "Found seed" << std::endl;
	  siseed->identify();
	}
      
      m_seedContainer->insert(siseed.get());

    }
}

void PHActsKDTreeSeeding::findInttMatches(
  std::map<TrkrDefs::cluskey, Acts::Vector3>& clusters,
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
      Acts::Vector3 lastclusglob = Acts::Vector3::Zero();
      for(const auto& [key, pos] : clusters)
	{
	  float lastclusglobr = sqrt(square(lastclusglob(0))+square(lastclusglob(1)));
	  float thisr = sqrt(square(pos(0))+square(pos(1)));
	  if(thisr > lastclusglobr) lastclusglob = pos;
	}

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

      if(Verbosity() > 2)
	{
	  std::cout << "Projected point is : " << xProj[layer] << ", "
		    << yProj[layer] << ", " << zProj[layer] << std::endl;
	}
    }
  
  matchInttClusters(clusters, xProj, yProj, zProj);
}


void PHActsKDTreeSeeding::matchInttClusters(
  std::map<TrkrDefs::cluskey, Acts::Vector3>& clusters,
  const double xProj[],
  const double yProj[],
  const double zProj[])
{

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
	    { 
	    
	      continue; 
	    }

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
	     
	      /// Z strip spacing is the entire strip, so because we use fabs
	      /// we divide by two
	      if(fabs(projectionLocal[1] - cluster->getLocalX()) < m_rPhiSearchWin and
		 fabs(projectionLocal[2] - cluster->getLocalY()) < stripZSpacing / 2.)
		{
	
		  /// Cache INTT global positions with seed
		  const auto globalPos = m_tGeometry->getGlobalPosition(
                    cluskey, cluster);
		  clusters.insert(std::make_pair(cluskey, globalPos));

		 		  	      
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
}

std::vector<const SpacePoint*> PHActsKDTreeSeeding::getMvtxSpacePoints()
{
  std::vector<const SpacePoint*> spVec;
  unsigned int numSiliconHits = 0;
 
  for(const auto& hitsetkey : m_clusterMap->getHitSetKeys(TrkrDefs::TrkrId::mvtxId))
    {
      auto range = m_clusterMap->getClusters(hitsetkey);
      for( auto clusIter = range.first; clusIter != range.second; ++clusIter )
	{
	  const auto cluskey = clusIter->first;
	  const auto cluster = clusIter->second;
	  if(m_iterationMap != NULL && m_nIteration > 0 )
	    {
	      if(m_iterationMap->getIteration(cluskey) > 0)
		{
		  continue; 
		}
	    }

	  const auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);
	  const auto surface = m_tGeometry->maps().getSiliconSurface(hitsetkey);
	  if(!surface)
	    { continue; }

	  auto sp = makeSpacePoint(surface, cluskey, cluster).release();
	  spVec.push_back(sp);
	  numSiliconHits++;
	}
    }

  if(Verbosity() > 1)
    {
      std::cout << "Total number of silicon hits to seed find with is "
		<< numSiliconHits << std::endl;
    }

  return spVec;
}


SpacePointPtr PHActsKDTreeSeeding::makeSpacePoint(const Surface& surf,
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
  if(m_clusterVersion==3)
    {
      localCov(0,0) = clus->getActsLocalError(0,0) * Acts::UnitConstants::cm2;
      localCov(1,1) = clus->getActsLocalError(1,1) * Acts::UnitConstants::cm2;
    }
  else if(m_clusterVersion==4)
    {
      auto para_errors = m_clusErrPara.get_si_cluster_error(clus,key);
      localCov(0,0) = para_errors.first * Acts::UnitConstants::cm2;
      localCov(1,1) = para_errors.second * Acts::UnitConstants::cm2;
    }
    
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
    {
      std::cout << "Space point has " 
		<< x << ", " << y << ", " << z << " with local coords "
		<< localPos.transpose() 
		<< " with rphi/z variances " << localCov(0,0) 
		<< ", " << localCov(1,1) << " and rotated variances "
		<< var[0] << ", " << var[1] 
		<< " and cluster key "
		<< key << " and geo id "
		<< surf->geometryId() << std::endl;
    }

  return spPtr;

}


//____________________________________________________________________________..
int PHActsKDTreeSeeding::End(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsKDTreeSeeding::getNodes(PHCompositeNode* topNode)
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

int PHActsKDTreeSeeding::createNodes(PHCompositeNode* topNode)
{

  PHNodeIterator iter(topNode);
  
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHActsKDTreeSeeding::createNodes");
  }

  PHNodeIterator dstIter(dstNode);
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_seedContainer = findNode::getClass<TrackSeedContainer>(topNode, m_trackMapName);
  if(!m_seedContainer)
    {
      m_seedContainer = new TrackSeedContainer_v1;
      PHIODataNode<PHObject> *trackNode = 
	new PHIODataNode<PHObject>(m_seedContainer, m_trackMapName, "PHObject");
      svtxNode->addNode(trackNode);

    }

  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::SeedFinderOrthogonalConfig<SpacePoint> PHActsKDTreeSeeding::configureSeedFinder()
{
  Acts::SeedFinderOrthogonalConfig<SpacePoint> cfg;

  Acts::SeedFilterConfig filterCfg;
  filterCfg.maxSeedsPerSpM = m_maxSeedsPerSpM;

  cfg.seedFilter =
      std::make_unique<Acts::SeedFilter<SpacePoint>>(
          Acts::SeedFilter<SpacePoint>(filterCfg));

  cfg.rMax = m_rMax;
  cfg.deltaRMinTopSP = m_deltaRMinTopSP;
  cfg.deltaRMaxTopSP = m_deltaRMaxTopSP;
  cfg.deltaRMinBottomSP = m_deltaRMinBottomSP;
  cfg.deltaRMaxBottomSP = m_deltaRMaxBottomSP;
  cfg.collisionRegionMin = m_collisionRegionMin;
  cfg.collisionRegionMax = m_collisionRegionMax;
  cfg.zMin = m_zMin;
  cfg.zMax = m_zMax;
  cfg.maxSeedsPerSpM = m_maxSeedsPerSpM;
  cfg.cotThetaMax = m_cotThetaMax;
  cfg.sigmaScattering = m_sigmaScattering;
  cfg.radLengthPerSeed = m_radLengthPerSeed;
  cfg.minPt = m_minPt;
  cfg.bFieldInZ = m_bFieldInZ;
  cfg.beamPos =
      Acts::Vector2(m_beamPosX, m_beamPosY);
  cfg.impactMax = m_impactMax;
  cfg.rMinMiddle = m_rMinMiddle;
  cfg.rMaxMiddle = m_rMaxMiddle;

  // Taken from SeedingOrthogonalAlgorithm.cpp, e.g.
  // https://github.com/acts-project/acts/blob/4b9d9c408158d09bc413e50fc41dbf3277bd751b/Examples/Algorithms/TrackFinding/src/SeedingOrthogonalAlgorithm.cpp
  // calculation of scattering using the highland formula
  // convert pT to p once theta angle is known
  cfg.highland =
      13.6 * std::sqrt(cfg.radLengthPerSeed) *
      (1 + 0.038 * std::log(cfg.radLengthPerSeed));
  float maxScatteringAngle =
      cfg.highland / cfg.minPt;
  cfg.maxScatteringAngle2 =
      maxScatteringAngle * maxScatteringAngle;
  // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and
  // millimeter
  // TODO: change using ACTS units
  cfg.pTPerHelixRadius =
      300. * cfg.bFieldInZ;
  cfg.minHelixDiameter2 =
      std::pow(cfg.minPt * 2 /
                   cfg.pTPerHelixRadius,
               2);

  cfg.pT2perRadius = std::pow(
      cfg.highland / cfg.pTPerHelixRadius,
      2);

  return cfg;
}
