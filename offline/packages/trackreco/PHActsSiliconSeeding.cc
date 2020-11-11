#include "PHActsSiliconSeeding.h"

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

PHActsSiliconSeeding::PHActsSiliconSeeding(std::string& name)
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
  return Fun4AllReturnCodes::EVENT_OK;
}
int PHActsSiliconSeeding::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
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
  config.rMin = 2. * Acts::UnitConstants::cm;
  config.zMin = m_zMin;
  config.zMax = m_zMax;

  /// Min/max distance between two measurements in one seed
  config.deltaRMin = 0.4 * Acts::UnitConstants::cm;
  config.deltaRMax = m_deltaRMax;

  /// Limiting collision region in z
  config.collisionRegionMin = -15. * Acts::UnitConstants::cm;
  config.collisionRegionMax = 15. * Acts::UnitConstants::cm;
  
  config.maxSeedsPerSpM = m_maxSeedsPerSpM;
  config.cotThetaMax = m_cotThetaMax;
  config.minPt = m_minSeedPt;
  config.bFieldInZ = m_bField;

  /// Average radiation length traversed per seed
  config.radLengthPerSeed = 0.05;

  /// Maximum impact parameter must be smaller than rMin
  config.impactMax = 1.5 * Acts::UnitConstants::cm;

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
