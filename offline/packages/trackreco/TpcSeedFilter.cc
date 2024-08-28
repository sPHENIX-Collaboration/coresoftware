#include "TpcSeedFilter.h"

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/getClass.h>

#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

int TpcSeedFilter::InitRun(PHCompositeNode* topNode)
{
  _trackseeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!_trackseeds)
  {
    std::cerr << PHWHERE << " ERROR: Can't find TrackSeedContainer " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcSeedFilter::process_event(PHCompositeNode* topNode)
{
  if (topNode == nullptr)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (Verbosity() > 1000)
  {
    topNode->print();
  }

  // iterate through the seeds and remove the ones that fail the cuts
  for (unsigned int iseed = 0; iseed < _trackseeds->size(); ++iseed)
  {
    TrackSeed* seed = _trackseeds->get(iseed);
    if (seed == nullptr)
    {
      continue;
    }

    if (
        (_cut_on_min_pt && (seed->get_pt() < _min_pt)) || (_cut_on_max_pt && (seed->get_pt() > _max_pt)) || (_cut_on_min_eta && (seed->get_eta() < _min_eta)) || (_cut_on_max_eta && (seed->get_eta() > _max_eta)) || (seed->size_cluster_keys() < _nclus_min))
    {
      if (Verbosity() > 4)
      {
        std::cout << " Cutting track seed: id(" << iseed
                  << ")  nclusters: " << static_cast<int>(seed->size_cluster_keys())
                  << "   pt: " << seed->get_pt() << "  eta: " << seed->get_eta() << std::endl;
      }
      if (Verbosity() > 10)
      {
        std::cout << " Cuts are: " << _cut_on_min_pt << " " << _cut_on_max_pt << " " << _cut_on_min_eta << " " << _cut_on_max_eta << " " << _nclus_min << std::endl;
      }
      _trackseeds->erase(iseed);
      continue;
    }
    if (_must_span_sectors)
    {
      // check that there are clusters in at least 2 of the three layers of sectors
      bool in_0 = false;
      bool in_1 = false;
      bool in_2 = false;
      unsigned int sec_cnt = 0;
      for (auto key = seed->begin_cluster_keys();
           key != seed->end_cluster_keys();
           ++key)
      {
        unsigned int layer = TrkrDefs::getLayer(*key);
        if (Verbosity() > 4)
        {
          std::cout << ((int) layer) << " ";
        }

        if (layer < 23 && !in_0)
        {
          in_0 = true;
          sec_cnt += 1;
        }
        else if (layer < 40 && !in_1)
        {
          in_1 = true;
          sec_cnt += 1;
        }
        else if (layer < 49 && !in_2)  // CHECK
        {
          in_2 = true;
          sec_cnt += 1;
        }

        if (sec_cnt >= _min_radial_sectors)
        {
          break;
        }
      }
      if (sec_cnt < _min_radial_sectors)
      {
        _trackseeds->erase(iseed);

        if (Verbosity() > 4)
        {
          std::cout << " Cutting track seed: id(" << static_cast<int>(iseed)
                    << ") :  clusters only in ";
          if (in_0)
          {
            std::cout << " inner ";
          }
          if (in_1)
          {
            std::cout << " middle ";
          }
          if (in_2)
          {
            std::cout << " outer ";
          }
          std::cout << " radial sectors but needs at least "
                    << static_cast<int>(_min_radial_sectors) << std::endl;
        }
        continue;
      }
    }
  }
  if (Verbosity() > 5)
  {
    std::cout << " Done seed on QA for tracks." << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
