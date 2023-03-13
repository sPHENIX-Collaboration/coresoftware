/*!
 * \file Fun4AllDstPileupInputMerger.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "Fun4AllDstPileupMerger.h"

#include "PHG4Hit.h"  // for PHG4Hit
#include "PHG4HitContainer.h"
#include "PHG4Hitv1.h"
#include "PHG4Particle.h"  // for PHG4Particle
#include "PHG4Particlev3.h"
#include "PHG4TruthInfoContainer.h"
#include "PHG4VtxPoint.h"  // for PHG4VtxPoint
#include "PHG4VtxPointv1.h"

#include <phhepmc/PHHepMCGenEvent.h>  // for PHHepMCGenEvent
#include <phhepmc/PHHepMCGenEventMap.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>          // for PHIODataNode
#include <phool/PHNode.h>                // for PHNode
#include <phool/PHNodeIterator.h>        // for PHNodeIterator
#include <phool/PHNodeOperation.h>
#include <phool/PHObject.h>              // for PHObject
#include <phool/getClass.h>

#include <TObject.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#pragma GCC diagnostic pop

#include <climits>
#include <iostream>
#include <iterator>
#include <utility>

// convenient aliases for deep copying nodes
namespace
{
  using PHG4Particle_t = PHG4Particlev3;
  using PHG4VtxPoint_t = PHG4VtxPointv1;
  using PHG4Hit_t = PHG4Hitv1;

  //! utility class to find all PHG4Hit container nodes from the DST node
  class FindG4HitContainer : public PHNodeOperation
  {
   public:
    //! container map alias
    using ContainerMap = std::map<std::string, PHG4HitContainer *>;

    //! get container map
    const ContainerMap &containers() const
    {
      return m_containers;
    }

   protected:
    //! iterator action
    void perform(PHNode *node) override
    {
      // check type name. Only load PHIODataNode
      if (node->getType() != "PHIODataNode") return;

      // cast to IODataNode and check data
      auto ionode = static_cast<PHIODataNode<TObject> *>(node);
      auto data = dynamic_cast<PHG4HitContainer *>(ionode->getData());
      if (data)
      {
        m_containers.insert(std::make_pair(node->getName(), data));
      }
    }

   private:
    //! container map
    ContainerMap m_containers;
  };

}  // namespace

//_____________________________________________________________________________
void Fun4AllDstPileupMerger::load_nodes(PHCompositeNode *dstNode)
{
  // hep mc
  m_geneventmap = findNode::getClass<PHHepMCGenEventMap>(dstNode, "PHHepMCGenEventMap");
  if (!m_geneventmap)
  {
    std::cout << "Fun4AllDstPileupMerger::load_nodes - creating PHHepMCGenEventMap" << std::endl;
    m_geneventmap = new PHHepMCGenEventMap();
    dstNode->addNode(new PHIODataNode<PHObject>(m_geneventmap, "PHHepMCGenEventMap", "PHObject"));
  }

  // find all G4Hit containers under dstNode
  FindG4HitContainer nodeFinder;
  PHNodeIterator(dstNode).forEach(nodeFinder);
  m_g4hitscontainers = nodeFinder.containers();

  // g4 truth info
  m_g4truthinfo = findNode::getClass<PHG4TruthInfoContainer>(dstNode, "G4TruthInfo");
  if (!m_g4truthinfo)
  {
    std::cout << "Fun4AllDstPileupMerger::load_nodes - creating node G4TruthInfo" << std::endl;
    m_g4truthinfo = new PHG4TruthInfoContainer();
    dstNode->addNode(new PHIODataNode<PHObject>(m_g4truthinfo, "G4TruthInfo", "PHObject"));
  }
}

//_____________________________________________________________________________
void Fun4AllDstPileupMerger::copy_background_event(PHCompositeNode *dstNode, double delta_t) const
{
  // copy PHHepMCGenEventMap
  const auto map = findNode::getClass<PHHepMCGenEventMap>(dstNode, "PHHepMCGenEventMap");

  // keep track of new embed id, after insertion as background event
  int new_embed_id = -1;

  if (map && m_geneventmap)
  {
    if (map->size() != 1)
    {
      std::cout << "Fun4AllDstPileupMerger::copy_background_event - cannot merge events that contain more than one PHHepMCGenEventMap" << std::endl;
      return;
    }

    // get event and insert in new map
    auto genevent = map->get_map().begin()->second;
    auto newevent = m_geneventmap->insert_background_event(genevent);

    /*
     * this hack prevents a crash when writting out
     * it boils down to root trying to write deleted items from the HepMC::GenEvent copy if the source has been deleted
     * it does not happen if the source gets written while the copy is deleted
     */
    newevent->getEvent()->swap(*genevent->getEvent());

    // shift vertex time and store new embed id
    newevent->moveVertex(0, 0, 0, delta_t);
    new_embed_id = newevent->get_embedding_id();
  }

  // copy truth container
  // keep track of the correspondance between source index and destination index for vertices, tracks and showers
  using ConversionMap = std::map<int, int>;
  ConversionMap vtxid_map;
  ConversionMap trkid_map;

  const auto container_truth = findNode::getClass<PHG4TruthInfoContainer>(dstNode, "G4TruthInfo");
  if (container_truth && m_g4truthinfo)
  {
    {
      // primary vertices
      auto key = m_g4truthinfo->maxvtxindex();
      const auto range = container_truth->GetPrimaryVtxRange();
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        // clone vertex, insert in map, and add index conversion
        const auto &sourceVertex = iter->second;
        auto newVertex = new PHG4VtxPoint_t(sourceVertex);
        newVertex->set_t(sourceVertex->get_t() + delta_t);
        m_g4truthinfo->AddVertex(++key, newVertex);
        vtxid_map.insert(std::make_pair(sourceVertex->get_id(), key));
      }
    }

    {
      // secondary vertices
      auto key = m_g4truthinfo->minvtxindex();
      const auto range = container_truth->GetSecondaryVtxRange();

      // loop from last to first to preserve order with respect to the original event
      for (
          auto iter = std::reverse_iterator<PHG4TruthInfoContainer::ConstVtxIterator>(range.second);
          iter != std::reverse_iterator<PHG4TruthInfoContainer::ConstVtxIterator>(range.first);
          ++iter)
      {
        // clone vertex, shift time, insert in map, and add index conversion
        const auto &sourceVertex = iter->second;
        auto newVertex = new PHG4VtxPoint_t(sourceVertex);
        newVertex->set_t(sourceVertex->get_t() + delta_t);
        m_g4truthinfo->AddVertex(--key, newVertex);
        vtxid_map.insert(std::make_pair(sourceVertex->get_id(), key));
      }
    }

    {
      // primary particles
      auto key = m_g4truthinfo->maxtrkindex();
      const auto range = container_truth->GetPrimaryParticleRange();
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        const auto &source = iter->second;
        auto dest = new PHG4Particle_t(source);
        m_g4truthinfo->AddParticle(++key, dest);
        dest->set_track_id(key);

        // set parent to zero
        dest->set_parent_id(0);

        // set primary to itself
        dest->set_primary_id(dest->get_track_id());

        // update vertex
        const auto keyiter = vtxid_map.find(source->get_vtx_id());
        if (keyiter != vtxid_map.end())
          dest->set_vtx_id(keyiter->second);
        else
          std::cout << "Fun4AllDstPileupMerger::copy_background_event - vertex id " << source->get_vtx_id() << " not found in map" << std::endl;

        // insert in map
        trkid_map.insert(std::make_pair(source->get_track_id(), dest->get_track_id()));
      }
    }

    {
      // secondary particles
      auto key = m_g4truthinfo->mintrkindex();
      const auto range = container_truth->GetSecondaryParticleRange();

      /*
       * loop from last to first to preserve order with respect to the original event
       * also this ensures that for a given particle its parent has already been converted and thus found in the map
       */
      for (
          auto iter = std::reverse_iterator<PHG4TruthInfoContainer::ConstIterator>(range.second);
          iter != std::reverse_iterator<PHG4TruthInfoContainer::ConstIterator>(range.first);
          ++iter)
      {
        const auto &source = iter->second;
        auto dest = new PHG4Particle_t(source);
        m_g4truthinfo->AddParticle(--key, dest);
        dest->set_track_id(key);

        // update parent id
        auto keyiter = trkid_map.find(source->get_parent_id());
        if (keyiter != trkid_map.end())
          dest->set_parent_id(keyiter->second);
        else
          std::cout << "Fun4AllDstPileupMerger::copy_background_event - track id " << source->get_parent_id() << " not found in map" << std::endl;

        // update primary id
        keyiter = trkid_map.find(source->get_primary_id());
        if (keyiter != trkid_map.end())
          dest->set_primary_id(keyiter->second);
        else
          std::cout << "Fun4AllDstPileupMerger::copy_background_event - track id " << source->get_primary_id() << " not found in map" << std::endl;

        // update vertex
        keyiter = vtxid_map.find(source->get_vtx_id());
        if (keyiter != vtxid_map.end())
          dest->set_vtx_id(keyiter->second);
        else
          std::cout << "Fun4AllDstPileupMerger::copy_background_event - vertex id " << source->get_vtx_id() << " not found in map" << std::endl;

        // insert in map
        trkid_map.insert(std::make_pair(source->get_track_id(), dest->get_track_id()));
      }
    }

    // vertex embed flags
    /* embed flag is stored only for primary vertices, consistently with PHG4TruthEventAction */
    for (const auto &pair : vtxid_map)
    {
      if (pair.first > 0) m_g4truthinfo->AddEmbededVtxId(pair.second, new_embed_id);
    }

    // track embed flags
    /* embed flag is stored only for primary tracks, consistently with PHG4TruthEventAction */
    for (const auto &pair : trkid_map)
    {
      if (pair.first > 0) m_g4truthinfo->AddEmbededTrkId(pair.second, new_embed_id);
    }
  }

  // copy g4hits
  // loop over registered maps
  for (const auto &pair : m_g4hitscontainers)
  {
    // check destination node
    if (!pair.second)
    {
      std::cout << "Fun4AllDstPileupMerger::copy_background_event - invalid destination container " << pair.first << std::endl;
      continue;
    }

    // find source node
    auto container_hit = findNode::getClass<PHG4HitContainer>(dstNode, pair.first);
    if (!container_hit)
    {
      std::cout << "Fun4AllDstPileupMerger::copy_background_event - invalid source container " << pair.first << std::endl;
      continue;
    }
    auto detiter = m_DetectorTiming.find(pair.first);
// apply special  cuts for selected detectors
    if (detiter != m_DetectorTiming.end())
    {
      if (delta_t < detiter->second.first || delta_t > detiter->second.second)
      {
	continue;
      }
    }
    {
      // hits
      const auto range = container_hit->getHits();
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        // clone hit
        const auto &sourceHit = iter->second;
        auto newHit = new PHG4Hit_t(sourceHit);

        // shift time
        newHit->set_t(0, sourceHit->get_t(0) + delta_t);
        newHit->set_t(1, sourceHit->get_t(1) + delta_t);

        // update track id
        const auto keyiter = trkid_map.find(sourceHit->get_trkid());
        if (keyiter != trkid_map.end())
          newHit->set_trkid(keyiter->second);
        else
          std::cout << "Fun4AllDstPileupMerger::copy_background_event - track id " << sourceHit->get_trkid() << " not found in map" << std::endl;

        /*
         * reset shower ids
         * it was decided that showers from the background events will not be copied to the merged event
         * as such we just reset the hits shower id
         */
        newHit->set_shower_id(INT_MIN);

        /*
         * this will generate a new key for the hit and assign it to the hit
         * this ensures that there is no conflict with the hits from the 'main' event
         */
        pair.second->AddHit(newHit->get_detid(), newHit);
      }
    }

    {
      // layers
      const auto range = container_hit->getLayers();
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        pair.second->AddLayer(*iter);
      }
    }
  }
}
