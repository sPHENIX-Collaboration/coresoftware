#include "PHSiliconSeedMerger.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

/**
 * @brief Construct a PHSiliconSeedMerger with the given subsystem name.
 *
 * Initializes the PHSiliconSeedMerger and forwards the provided subsystem
 * name to the base SubsysReco constructor.
 *
 * @param name Subsystem name used to register this module in the node tree.
 */
PHSiliconSeedMerger::PHSiliconSeedMerger(const std::string& name)
  : SubsysReco(name)
{
}

/**
 * @brief Default destructor for PHSiliconSeedMerger.
 *
 * Performs default cleanup of the merger object and its owned resources.
 */
PHSiliconSeedMerger::~PHSiliconSeedMerger() = default;

/**
 * @brief Perform module initialization (no operation required).
 *
 * This implementation does not perform any setup and always succeeds.
 *
 * @return int `EVENT_OK` indicating initialization succeeded.
 */
int PHSiliconSeedMerger::Init(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

/**
 * @brief Initializes run-time resources by retrieving required nodes.
 *
 * Calls getNodes(topNode) to locate and cache containers needed for processing this run.
 *
 * @param topNode Root of the node tree from which required nodes are retrieved.
 * @return int `EVENT_OK` on success, `ABORTEVENT` or another non-zero code on failure.
 */
int PHSiliconSeedMerger::InitRun(PHCompositeNode* topNode)
{
  int ret = getNodes(topNode);
  return ret;
}

/**
 * @brief Merge overlapping silicon seed tracks by consolidating MVTX cluster keys.
 *
 * Detects seeds whose MVTX cluster key sets fully overlap (one set equals the
 * intersection) and treats one seed as a duplicate of the other. The merger
 * preserves the seed with the larger MVTX key set; if both seeds share the
 * same MVTX strobe and seed merging is enabled, the smaller seed's MVTX keys
 * are merged into the preserved seed. After consolidation, duplicate seeds are
 * erased from the silicon track container and preserved seeds are updated to
 * include any newly merged MVTX cluster keys.
 *
 * @return Fun4AllReturnCodes::EVENT_OK on successful processing.
 *
 */
int PHSiliconSeedMerger::process_event(PHCompositeNode* /*unused*/)
{
  std::multimap<unsigned int, std::set<TrkrDefs::cluskey>> matches;
  std::set<unsigned int> seedsToDelete;

  if (Verbosity() > 0)
  {
    std::cout << "Silicon seed track container has " << m_siliconTracks->size() << std::endl;
  }

  for (unsigned int track1ID = 0;
       track1ID != m_siliconTracks->size();
       ++track1ID)
  {
    TrackSeed* track1 = m_siliconTracks->get(track1ID);

    if (seedsToDelete.contains(track1ID))
    {
      continue;
    }

    std::set<TrkrDefs::cluskey> mvtx1Keys;
    int track1Strobe = std::numeric_limits<int>::quiet_NaN();
    for (auto iter = track1->begin_cluster_keys();
         iter != track1->end_cluster_keys();
         ++iter)
    {
      TrkrDefs::cluskey ckey = *iter;
      auto trkrid = TrkrDefs::getTrkrId(ckey);
      if (m_mvtxOnly && trkrid == TrkrDefs::TrkrId::inttId)
      {
        continue;
      }
      if (TrkrDefs::getTrkrId(ckey) == TrkrDefs::TrkrId::mvtxId)
      {
        track1Strobe = MvtxDefs::getStrobeId(ckey);
      }
      mvtx1Keys.insert(ckey);
    }

    /// We can speed up the code by only iterating over the track seeds
    /// that are further in the map container from the current track,
    /// since the comparison of e.g. track 1 with track 2 doesn't need
    /// to be repeated with track 2 to track 1.
    for (unsigned int track2ID = track1ID;
         track2ID != m_siliconTracks->size();
         ++track2ID)
    {
      if (track1ID == track2ID)
      {
        continue;
      }

      TrackSeed* track2 = m_siliconTracks->get(track2ID);
      if (track2 == nullptr)
      {
        continue;
      }
      int track2Strobe = std::numeric_limits<int>::quiet_NaN();

      std::set<TrkrDefs::cluskey> mvtx2Keys;
      for (TrackSeed::ConstClusterKeyIter iter = track2->begin_cluster_keys();
           iter != track2->end_cluster_keys();
           ++iter)
      {
        TrkrDefs::cluskey ckey = *iter;
        auto trkrid = TrkrDefs::getTrkrId(ckey);
        if (m_mvtxOnly && trkrid == TrkrDefs::TrkrId::inttId)
        {
          continue;
        }
        mvtx2Keys.insert(ckey);
        if (TrkrDefs::getTrkrId(ckey) == TrkrDefs::TrkrId::mvtxId)
        {
          track2Strobe = MvtxDefs::getStrobeId(ckey);
        }
      }

      std::vector<TrkrDefs::cluskey> intersection;
      std::set_intersection(mvtx1Keys.begin(),
                            mvtx1Keys.end(),
                            mvtx2Keys.begin(),
                            mvtx2Keys.end(),
                            std::back_inserter(intersection));

      /// If the intersection fully encompasses one of the tracks, it is completely duplicated
      if (intersection.size() == mvtx1Keys.size() || intersection.size() == mvtx2Keys.size())
      {
        if (Verbosity() > 2)
        {
          std::cout << "Track " << track1ID << " keys " << std::endl;
          for (const auto& key : mvtx1Keys)
          {
            std::cout << "   ckey: " << key << std::endl;
          }
          std::cout << "Track " << track2ID << " keys " << std::endl;
          for (const auto& key : mvtx2Keys)
          {
            std::cout << "   ckey: " << key << std::endl;
          }
          std::cout << "Intersection keys " << std::endl;
          for (auto& key : intersection)
          {
            std::cout << "   ckey: " << key << std::endl;
          }
        }

        /// one of the tracks is encompassed in the other. Take the larger one
        std::set<TrkrDefs::cluskey> keysToKeep;
        if (mvtx1Keys.size() >= mvtx2Keys.size())
        {
          keysToKeep = mvtx1Keys;
          if (track1Strobe == track2Strobe && m_mergeSeeds)
          {
            keysToKeep.insert(mvtx2Keys.begin(), mvtx2Keys.end());
          }
          matches.insert(std::make_pair(track1ID, mvtx1Keys));
          seedsToDelete.insert(track2ID);
          if (Verbosity() > 2)
          {
            std::cout << "     will delete seed " << track2ID << std::endl;
          }
        }
        else
        {
          keysToKeep = mvtx2Keys;
          if (track1Strobe == track2Strobe && m_mergeSeeds)
          {
            keysToKeep.insert(mvtx1Keys.begin(), mvtx1Keys.end());
          }
          matches.insert(std::make_pair(track2ID, mvtx2Keys));
          seedsToDelete.insert(track1ID);
          if (Verbosity() > 2)
          {
            std::cout << "     will delete seed " << track1ID << std::endl;
          }
        }

        if (Verbosity() > 2)
        {
          std::cout << "Match IDed" << std::endl;
          for (const auto& key : mvtx1Keys)
          {
            std::cout << "  total track keys " << key << std::endl;
          }
        }

        matches.insert(std::make_pair(track1ID, mvtx1Keys));
        seedsToDelete.insert(track2ID);
        break;
      }
    }
  }

  for (const auto& [trackKey, mvtxKeys] : matches)
  {
    auto* track = m_siliconTracks->get(trackKey);
    if (Verbosity() > 2)
    {
      std::cout << "original track: " << std::endl;
      track->identify();
    }

    for (const auto& key : mvtxKeys)
    {
      if (track->find_cluster_key(key) == track->end_cluster_keys())
      {
        track->insert_cluster_key(key);
        if (Verbosity() > 2)
        {
          std::cout << "adding " << key << std::endl;
        }
      }
    }
  }

  for (const auto& key : seedsToDelete)
  {
    if (Verbosity() > 2)
    {
      std::cout << "Erasing track " << key << std::endl;
    }
    m_siliconTracks->erase(key);
  }

  if (Verbosity() > 2)
  {
    for (const auto& seed : *m_siliconTracks)
    {
      if (!seed)
      {
        continue;
      }
      seed->identify();
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

/**
 * @brief Reset per-event state for the merger.
 *
 * This implementation performs no per-event cleanup and always reports success.
 *
 * @return Integer status code: `Fun4AllReturnCodes::EVENT_OK`.
 */
int PHSiliconSeedMerger::ResetEvent(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

/**
 * @brief Perform end-of-run shutdown for the silicon seed merger.
 *
 * @return int EVENT_OK on successful completion.
 */
int PHSiliconSeedMerger::End(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

/**
 * @brief Retrieve required nodes from the top-level node tree and validate availability.
 *
 * Locates the silicon TrackSeedContainer using m_trackMapName and stores it in
 * m_siliconTracks. If the container is not found, the function logs an error
 * message and signals an abort for the current event.
 *
 * @param topNode Root node used to search for the TrackSeedContainer.
 * @return int Fun4AllReturnCodes::EVENT_OK on success, Fun4AllReturnCodes::ABORTEVENT if the silicon TrackSeedContainer is not present.
 */
int PHSiliconSeedMerger::getNodes(PHCompositeNode* topNode)
{
  m_siliconTracks = findNode::getClass<TrackSeedContainer>(topNode, m_trackMapName);
  if (!m_siliconTracks)
  {
    std::cout << PHWHERE << "No silicon track container, can't merge seeds"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


