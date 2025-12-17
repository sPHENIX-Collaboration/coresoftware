#include "FillTruthRecoMatchMap.h"

#include <g4tracking/EmbRecoMatch.h>
#include <g4tracking/EmbRecoMatchContainer.h>

#include <trackbase_historic/PHG4ParticleSvtxMap.h>
#include <trackbase_historic/PHG4ParticleSvtxMap_v1.h>
#include <trackbase_historic/SvtxPHG4ParticleMap_v1.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <boost/format.hpp>

#include <iostream>

//____________________________________________________________________________..
FillTruthRecoMatchMap::FillTruthRecoMatchMap(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int FillTruthRecoMatchMap::Init(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    topNode->print();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int FillTruthRecoMatchMap::InitRun(PHCompositeNode *topNode)
{
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int FillTruthRecoMatchMap::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cout << PHWHERE << " DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in FillTruthRecoMatchMap::createNodes");
  }

  m_EmbRecoMatchContainer = findNode::getClass<EmbRecoMatchContainer>(topNode, "TRKR_EMBRECOMATCHCONTAINER");
  if (!m_EmbRecoMatchContainer)
  {
    std::cout << PHWHERE << " Cannot find node TRKR_EMBRECOMATCHCONTAINER on node tree; quitting " << std::endl;
    throw std::runtime_error(" Cannot find node TRKR_EMBRECOMATCHCONTAINER on node tree; quitting");
  }

  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_PHG4ParticleSvtxMap = findNode::getClass<PHG4ParticleSvtxMap>(topNode, "PHG4ParticleToRecoMap");
  if (!m_PHG4ParticleSvtxMap)
  {
    m_PHG4ParticleSvtxMap = new PHG4ParticleSvtxMap_v1;
    PHIODataNode<PHObject> *truthNode =
        new PHIODataNode<PHObject>(m_PHG4ParticleSvtxMap, "PHG4ParticleToRecoMap", "PHObject");
    svtxNode->addNode(truthNode);
  }

  m_SvtxPHG4ParticleMap = findNode::getClass<SvtxPHG4ParticleMap>(topNode, "RecoToPHG4ParticleMap");
  if (!m_SvtxPHG4ParticleMap)
  {
    m_SvtxPHG4ParticleMap = new SvtxPHG4ParticleMap_v1;
    PHIODataNode<PHObject> *recoNode =
        new PHIODataNode<PHObject>(m_SvtxPHG4ParticleMap, "RecoToPHG4ParticleMap", "PHObject");
    svtxNode->addNode(recoNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int FillTruthRecoMatchMap::process_event(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 5)
  {
    std::cout << " FillTruthRecoMatchMap::process_event() " << std::endl;
  }
  // make maps to all the matches for both truth and reco and fill the
  /* auto matches = m_EmbRecoMatchContainer->getMatches(); */

  /* map<unsigned short, set<int>> truth_entries; */
  /* map<unsigned short, set<int>> m_reco_entries; */
  /* int cnt {0}; */
  /* for (auto& match : m_EmbRecoMatchContainer->getMatches()) { */

  /* } */

  SvtxPHG4ParticleMap::Map map_RtoT{};  // RtoT = "Reco to Truth"
  PHG4ParticleSvtxMap::Map map_TtoR{};  // TtoR = "Truth to Reco"

  for (auto &match : m_EmbRecoMatchContainer->getMatches())
  {
    // every match is a unique match from truth to embedded particles
    const int gtrackID = match->idTruthTrack();
    const unsigned int id_reco = match->idRecoTrack();
    const unsigned short n_match = match->nClustersMatched();
    const unsigned short n_truth = match->nClustersTruth();
    const unsigned short n_reco = match->nClustersReco();

    // fill the map_TtoR
    /* PHG4ParticleSvtxMap::WeightedRecoTrackMap*  entry_TtoR; */
    if (!map_TtoR.contains(gtrackID))
    {
      map_TtoR[gtrackID] = PHG4ParticleSvtxMap::WeightedRecoTrackMap{};
    }
    auto &entry_TtoR = map_TtoR[gtrackID];
    float weight_TtoR = (float) n_match + (float) n_truth / 100.;
    if (!entry_TtoR.contains(weight_TtoR))
    {
      entry_TtoR[weight_TtoR] = {id_reco};  // i.e. std::set<unsigned int> { id_reco };
    }
    else
    {
      entry_TtoR[weight_TtoR].insert(id_reco);
    }

    // fill the map_RtoT
    /* SvtxPHG4ParticleMap::WeightedTruthTrackMap*  entry_RtoT; */
    if (!map_RtoT.contains(id_reco))
    {
      map_RtoT[id_reco] = SvtxPHG4ParticleMap::WeightedTruthTrackMap{};
    }
    auto &entry_RtoT = map_RtoT[id_reco];
    float weight_RtoT = (float) n_match + (float) n_reco / 100.;
    if (!entry_RtoT.contains(weight_RtoT))
    {
      entry_RtoT[weight_RtoT] = {gtrackID};  // i.e. std::set<int> { gtrackID }
    }
    else
    {
      entry_RtoT[weight_RtoT].insert(gtrackID);
    }

    if (Verbosity() > 20)
    {
      std::cout << (boost::format("EmbRecoMatch: gtrackID(%2i) id_reco(%2i) nclusters:match(%i),gtrack(%2i),reco(%2i)") % gtrackID % ((int) id_reco) % ((int) n_match) % ((int) n_truth) % ((int) n_reco)).str()
                << std::endl;
      std::cout << (boost::format("   -> in SvtxPHG4ParticleMap {id_reco->{weight->id_true}} = {%2i->{%5.2f->%2i}}") % ((int) id_reco) % weight_RtoT % ((int) gtrackID)).str() << std::endl;
      std::cout << (boost::format("   -> in PHG4ParticleSvtxMap {id_true->{weight->id_reco}} = {%2i->{%5.2f->%2i}}") % gtrackID % weight_TtoR % ((int) id_reco)).str() << std::endl;
    }
  }

  // fill the output maps
  for (auto &entry : map_TtoR)
  {
    m_PHG4ParticleSvtxMap->insert(entry.first, entry.second);
  }
  for (auto &entry : map_RtoT)
  {
    m_SvtxPHG4ParticleMap->insert(entry.first, entry.second);
  }

  if (Verbosity() > 15)
  {
    std::cout << PHWHERE << " Print out: " << std::endl
              << " Contents of SvtxPHG4ParticleMap (node \"RecoToPHG4ParticleMap\")" << std::endl
              << " and PHG4ParticleSvtxMap (node \"PHG4ParticleToRecoMap\")" << std::endl;

    std::cout << " --BEGIN-- Contents of SvtxPHG4ParticleMap: " << std::endl;
    for (auto &iter : *m_SvtxPHG4ParticleMap)
    {
      std::cout << (boost::format("    { %2i ") % ((int) iter.first)).str();  // id_reco
      auto n_matches = iter.second.size();
      long unsigned int cnt_matches = 0;
      for (const auto &matches : iter.second)
      {
        if (cnt_matches == 0)
        {
          std::cout << (boost::format("-> { %5.2f -> ") % matches.first).str();
        }
        else
        {
          std::cout << (boost::format("         -> { %5.2f -> ") % matches.first).str();
        }
        auto size = matches.second.size();
        long unsigned int i = 0;
        for (auto id_true : matches.second)
        {
          if (i < size - 1)
          {
            std::cout << (boost::format("%2i, ") % id_true).str();
          }
          else
          {
            std::cout << (boost::format("%2i }") % id_true).str();
          }
          ++i;
        }  // end matches cnt
        if (cnt_matches < (n_matches - 1))
        {
          std::cout << "," << std::endl;
        }
        else
        {
          std::cout << "}" << std::endl;
        }
        ++cnt_matches;
      }
    }
    std::cout << " --END-- Contents of SvtxPHG4ParticleMap: " << std::endl
              << std::endl;

    std::cout << " --BEGIN-- Contents of PHG4ParticleToRecoMap: " << std::endl;
    for (auto &iter : *m_PHG4ParticleSvtxMap)
    {
      std::cout << (boost::format("    { %2i ") % iter.first).str();  // id_true
      auto n_matches = iter.second.size();
      long unsigned int cnt_matches = 0;
      for (const auto &matches : iter.second)
      {
        if (cnt_matches == 0)
        {
          std::cout << (boost::format("-> { %5.2f -> ") % matches.first).str();
        }
        else
        {
          std::cout << (boost::format("         -> { %5.2f -> ") % matches.first).str();
        }
        auto size = matches.second.size();
        long unsigned int i = 0;
        for (auto id_reco : matches.second)
        {
          if (i < size - 1)
          {
            std::cout << (boost::format("%2i, ") % id_reco).str();
          }
          else
          {
            std::cout << (boost::format("%2i }") % id_reco).str();
          }
          ++i;
        }  // end matches cnt
        if (cnt_matches < (n_matches - 1))
        {
          std::cout << "," << std::endl;
        }
        else
        {
          std::cout << "}" << std::endl;
        }
        ++cnt_matches;
      }
    }
    std::cout << " --END-- Contents of SvtxPHG4ParticleMap: " << std::endl
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int FillTruthRecoMatchMap::End(PHCompositeNode * /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
