#include "PHActsSiliconSeeding.h"
#include "SphenixActsDetectorCuts.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#include <Acts/Seeding/BinnedSPGroup.hpp>
#include <Acts/Seeding/InternalSeed.hpp>
#include <Acts/Seeding/InternalSpacePoint.hpp>
#include <Acts/Seeding/Seed.hpp>
#include <Acts/Seeding/SeedFilter.hpp>

PHActsSiliconSeeding::PHActsSiliconSeeding(const std::string& name)
  : SubsysReco(name)
{
}

int PHActsSiliconSeeding::Init(PHCompositeNode *topNode)
{
  m_seedFinderCfg = configureSeeder();
  m_gridCfg = configureSPGrid();
  
  m_bottomBinFinder = 
    std::make_shared<Acts::BinFinder<SpacePoint>>(Acts::BinFinder<SpacePoint>());
  m_topBinFinder = 
    std::make_shared<Acts::BinFinder<SpacePoint>>(Acts::BinFinder<SpacePoint>());

  Acts::SeedFilterConfig sfCfg;
  sfCfg.maxSeedsPerSpM = m_maxSeedsPerSpM;

  SphenixActsDetectorCuts<SpacePoint> detectorCuts =
    SphenixActsDetectorCuts<SpacePoint>();
  m_seedFinderCfg.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
     Acts::SeedFilter<SpacePoint>(sfCfg));

  return Fun4AllReturnCodes::EVENT_OK;
}
int PHActsSiliconSeeding::InitRun(PHCompositeNode *topNode)
{
  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;
  
  if(createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;
  
  return Fun4AllReturnCodes::EVENT_OK;
}
int PHActsSiliconSeeding::process_event(PHCompositeNode *topNode)
{

  if(Verbosity() > -1)
    std::cout << "Processing PHActsSiliconSeeding event "
	      << m_event << std::endl;

  Acts::Seedfinder<SpacePoint> seedFinder(m_seedFinderCfg);
  
  /// Covariance converter tool needed by seed finder
  auto covConverter = [=](const SpacePoint& sp, float, float, float)
    -> Acts::Vector2D { return {sp.m_varianceRphi, sp.m_varianceZ};
  };

  std::vector<const SpacePoint*> spVec;
  unsigned int siliconHits = 0;

  std::cout << "collecting SpacePoints"<<std::endl;

  for(auto &[hitId, sl] : *m_sourceLinks)
    {
      /// collect only source links in silicon
      auto volume = sl.referenceSurface().geometryId().volume();
      
      /// If we run without MMs, volumes are 7, 9, 11 for mvtx, intt, tpc
      /// If we run with MMs, volumes are 10, 12, 14, 16 for mvtx, intt, tpc, mm
      if(volume == 11 or volume > 12)
	continue;

      auto sp = makeSpacePoint(hitId, sl).release();
      spVec.push_back(sp);
      siliconHits++;
    }

  std::unique_ptr<Acts::SpacePointGrid<SpacePoint>> grid = 
    Acts::SpacePointGridCreator::createGrid<SpacePoint>(m_gridCfg);

  auto spGroup = Acts::BinnedSPGroup<SpacePoint>( spVec.begin(),
						  spVec.end(),
						  covConverter,
						  m_bottomBinFinder,
						  m_topBinFinder,
						  std::move(grid),
						  m_seedFinderCfg);

  std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVector;
  auto groupIt = spGroup.begin();
  auto endGroup = spGroup.end();
  std::cout << "Perform seed finding " << std::endl;
  for(; !(groupIt == endGroup); ++groupIt)
    {
      seedVector.push_back(seedFinder.createSeedsForGroup(groupIt.bottom(),
							  groupIt.middle(),
							  groupIt.top()));
    }
  std::cout<<"seed finding done"<<std::endl;
  int numSeeds = 0;
  for(auto& outVec : seedVector)
    numSeeds += outVec.size();

  if(Verbosity() > -1)
    std::cout << "Total number of hits " << siliconHits 
	      << " in " << seedVector.size() << " regions gives " 
	      << numSeeds << " seeds " << std::endl;

  if(Verbosity()> -1)
    std::cout << "Finished PHActsSiliconSeeding process_event"
	      << std::endl;

  m_event++;

  return Fun4AllReturnCodes::EVENT_OK;
}
int PHActsSiliconSeeding::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

SpacePointPtr PHActsSiliconSeeding::makeSpacePoint(const unsigned int& hitId,
						const SourceLink& sl)
{
  Acts::Vector2D localPos(sl.location()(0), sl.location()(1));
  Acts::Vector3D globalPos(0,0,0);
  Acts::Vector3D mom(1,1,1);

  globalPos = sl.referenceSurface().localToGlobal(m_tGeometry->geoContext,
						  localPos, mom);

  auto cov = sl.covariance();
  float x = globalPos.x();
  float y = globalPos.y();
  float z = globalPos.z();
  float r = std::sqrt(x * x + y * y);
  float varianceRphi = cov(0,0);
  float varianceZ = cov(1,1);
  
  SpacePointPtr spPtr(new SpacePoint{sl.hitID(), x, y, z, r, 
	sl.referenceSurface().geometryId(), varianceRphi, varianceZ});

  if(Verbosity() > -1)
    std::cout << "Space point has " 
	      << x << ", " << y << ", " << z
	      << " with variances " << varianceRphi 
	      << ", " << varianceZ 
	      << " and hit id "
	      << sl.hitID() << " and geo id "
	      << sl.referenceSurface().geometryId() << std::endl;
  
  return spPtr;

}

Acts::SpacePointGridConfig PHActsSiliconSeeding::configureSPGrid()
{
  Acts::SpacePointGridConfig config;

  config.bFieldInZ = m_bField;
  config.minPt = m_minSeedPt;
  config.rMax = m_rMax;
  config.zMax = m_zMax;
  config.zMin = m_zMin;
  config.deltaRMax = m_deltaRMax;
  config.cotThetaMax = m_cotThetaMax;

  return config;
}

Acts::SeedfinderConfig<SpacePoint> PHActsSiliconSeeding::configureSeeder()
{
  Acts::SeedfinderConfig<SpacePoint> config;
  
  /// Limiting location of measurements (e.g. detector constraints)
  config.rMax = m_rMax;
  config.rMin = m_rMin;
  config.zMin = m_zMin;
  config.zMax = m_zMax;

  /// Min/max distance between two measurements in one seed
  config.deltaRMin = 1.;
  config.deltaRMax = m_deltaRMax;

  /// Limiting collision region in z
  config.collisionRegionMin = -100.;
  config.collisionRegionMax = 100.;
  config.sigmaScattering = 50.;
  config.maxSeedsPerSpM = m_maxSeedsPerSpM;
  config.cotThetaMax = m_cotThetaMax;
  config.minPt = m_minSeedPt;
  config.bFieldInZ = m_bField;

  /// Average radiation length traversed per seed
  config.radLengthPerSeed = 0.05;

  /// Maximum impact parameter must be smaller than rMin
  config.impactMax = 20;

  return config;
}

int PHActsSiliconSeeding::getNodes(PHCompositeNode *topNode)
{
  m_sourceLinks = findNode::getClass<std::map<unsigned int, SourceLink>>(topNode, "TrkrClusterSourceLinks");
  if(!m_sourceLinks)
    {
      std::cout << PHWHERE << "TrkrClusterSourceLinks node not on node tree. Bailing."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << PHWHERE << "No ActsTrackingGeometry on node tree. Bailing."
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
    throw std::runtime_error("Failed to find DST node in PHActsTracks::createNodes");
  }
  
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }


  return Fun4AllReturnCodes::EVENT_OK;
}
