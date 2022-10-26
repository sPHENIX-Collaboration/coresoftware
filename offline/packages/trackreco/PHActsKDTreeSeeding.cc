
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
      /// The acts seed has better z resolution than the circle fit
      siseed->set_Z0(seed.z() / Acts::UnitConstants::cm);

      m_seedContainer->insert(siseed.get());

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
