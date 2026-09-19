#include "PHSiliconSeedPruner.h"
#include "PHSiliconSeedPrunerHelper.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>

#include <climits>
#include <iostream>
#include <map>
#include <set>
#include <vector>

PHSiliconSeedPruner::PHSiliconSeedPruner(const std::string &name)
  : SubsysReco(name)
{
  m_rng = gsl_rng_alloc(gsl_rng_mt19937);
  set_random_seed(PHRandomSeed());  // default random seed
}

PHSiliconSeedPruner::~PHSiliconSeedPruner()
{
  gsl_rng_free(m_rng);
}

void PHSiliconSeedPruner::set_random_seed(unsigned int seed)
{
  m_randomSeed = seed;
  gsl_rng_set(m_rng, m_randomSeed);
}

int PHSiliconSeedPruner::InitRun(PHCompositeNode *topNode)
{
  if (m_trackMapName == m_outputTrackMapName)
  {
    std::cout << PHWHERE << " ERROR: Input and output TrackSeedContainer names are both " << m_trackMapName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_siliconSeeds = findNode::getClass<TrackSeedContainer>(topNode, m_trackMapName);
  if (!m_siliconSeeds)
  {
    std::cout << PHWHERE << " ERROR: Can't find " << m_trackMapName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_prunedSiliconSeeds = findNode::getClass<TrackSeedContainer>(topNode, m_outputTrackMapName);
  if (!m_prunedSiliconSeeds)
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      std::cout << PHWHERE << " ERROR: Can't find DST node" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    PHNodeIterator dstIter(dstNode);
    PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "SVTX"));
    if (!svtxNode)
    {
      svtxNode = new PHCompositeNode("SVTX");
      dstNode->addNode(svtxNode);
    }

    m_prunedSiliconSeeds = new TrackSeedContainer_v1;
    auto *outputNode = new PHIODataNode<PHObject>(m_prunedSiliconSeeds, m_outputTrackMapName, "PHObject");
    svtxNode->addNode(outputNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSiliconSeedPruner::process_event(PHCompositeNode * /*topNode*/)
{
  m_prunedSiliconSeeds->Reset();

  std::map<short int, std::vector<size_t>> seedIndicesByCrossing;
  std::set<size_t> selectedSeedIndices;

  size_t inputSeedCount = 0;
  size_t invalidCrossingCount = 0;
  for (size_t seedIndex = 0; seedIndex < m_siliconSeeds->size(); ++seedIndex)
  {
    TrackSeed *seed = m_siliconSeeds->get(seedIndex);
    if (!seed)
      continue;

    ++inputSeedCount;
    const short int crossing = seed->get_crossing();
    if (crossing == SHRT_MAX)
    {
      selectedSeedIndices.insert(seedIndex);
      ++invalidCrossingCount;
      continue;
    }

    seedIndicesByCrossing[crossing].push_back(seedIndex);
  }

  size_t uncertifiedSeedCount = 0;
  for (const auto &[crossing, seedIndices] : seedIndicesByCrossing)
  {
    const PHSiliconSeedPrunerHelper::Result result = PHSiliconSeedPrunerHelper::SelectSeeds(
        *m_siliconSeeds,
        seedIndices,
        m_rng,
        kMvtxLayerCount,
        kMaxRepresentatives,
        kSearchBudget,
        kHeuristicRestarts,
        kDebugMinimumGroupSize,
        Verbosity() > 0,
        Verbosity() > 1);

    selectedSeedIndices.insert(result.selectedSeedIndices.begin(), result.selectedSeedIndices.end());
    uncertifiedSeedCount += result.uncertifiedSeedIndices.size();

    if (Verbosity() > 0)
    {
      std::cout << Name() << ": crossing " << crossing
                << " input seeds " << seedIndices.size()
                << " selected seeds " << result.selectedSeedIndices.size()
                << " uncertified seeds " << result.uncertifiedSeedIndices.size()
                << std::endl;
    }
  }

  for (const size_t seedIndex : selectedSeedIndices)
  {
    m_prunedSiliconSeeds->insert(m_siliconSeeds->get(seedIndex));
  }

  if (Verbosity() > 0)
  {
    std::cout << Name()
              << ": input seeds " << inputSeedCount
              << ", retained seeds " << selectedSeedIndices.size()
              << ", retained invalid-crossing seeds " << invalidCrossingCount
              << ", uncertified selected seeds " << uncertifiedSeedCount
              << ", random seed " << m_randomSeed
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
