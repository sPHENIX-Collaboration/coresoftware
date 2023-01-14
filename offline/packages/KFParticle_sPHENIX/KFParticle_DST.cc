/*****************/
/* Cameron Dean  */
/*   LANL 2020   */
/* cdean@bnl.gov */
/*****************/
/*
 * Class to append reconstructed events to node tree
 */

//Ideas taken from PHRaveVertexing

#include "KFParticle_DST.h"

#include "KFParticle_Container.h"
#include "KFParticle_Tools.h"
#include "KFParticle_truthAndDetTools.h"

#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack
#include <trackbase_historic/SvtxTrackMap.h>  // for SvtxTrackMap, SvtxTr...
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrack_v2.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>

#include <KFParticle.h>  // for KFParticle

#include <cstdlib>   // for exit, size_t, abs
#include <iostream>  // for operator<<, endl
#include <map>       // for map, map<>::mapped_type
#include <utility>   // for pair

KFParticle_Tools kfpTupleTools_DST;
KFParticle_truthAndDetTools kfpTruthTools_DST;

int KFParticle_DST::createParticleNode(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* lowerNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!lowerNode)
  {
    lowerNode = new PHCompositeNode("DST");
    topNode->addNode(lowerNode);
    std::cout << "DST node added" << std::endl;
  }

  std::string baseName;
  std::string trackNodeName;
  std::string particleNodeName;

  if (m_container_name.empty())
    baseName = "reconstructedParticles";
  else
    baseName = m_container_name;

  //Cant have forward slashes in DST or else you make a subdirectory on save!!!
  std::string undrscr = "_";
  std::string nothing = "";
  std::map<std::string, std::string> forbiddenStrings;
  forbiddenStrings["/"] = undrscr;
  forbiddenStrings["("] = undrscr;
  forbiddenStrings[")"] = nothing;
  forbiddenStrings["+"] = "plus";
  forbiddenStrings["-"] = "minus";
  forbiddenStrings["*"] = "star";
  for (auto const& [badString, goodString] : forbiddenStrings)
  {
    size_t pos;
    while ((pos = baseName.find(badString)) != std::string::npos) baseName.replace(pos, 1, goodString);
  }

  trackNodeName = baseName + "_SvtxTrackMap";
  particleNodeName = baseName + "_KFParticle_Container";

  if (m_write_track_container)
  {
    m_recoTrackMap = new SvtxTrackMap_v1();
    PHIODataNode<PHObject>* trackNode = new PHIODataNode<PHObject>(m_recoTrackMap, trackNodeName.c_str(), "PHObject");
    lowerNode->addNode(trackNode);
    std::cout << trackNodeName << " node added" << std::endl;
  }

  if (m_write_particle_container)
  {
    m_recoParticleMap = new KFParticle_Container();
    PHIODataNode<PHObject>* particleNode = new PHIODataNode<PHObject>(m_recoParticleMap, particleNodeName.c_str(), "PHObject");
    lowerNode->addNode(particleNode);
    std::cout << particleNodeName << " node added" << std::endl;
  }

  if (!m_write_track_container && !m_write_particle_container)
  {
    std::cout << "You have asked to put your selection on the node tree but disabled both the SvtxTrackMap and KFParticle_Container\n";
    std::cout << "Check your options" << std::endl;
    exit(0);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void KFParticle_DST::fillParticleNode(PHCompositeNode* topNode, const KFParticle& motherParticle,
                                      const std::vector<KFParticle>& daughters,
                                      const std::vector<KFParticle>& intermediates)
{
  if (m_write_track_container)
  {
    fillParticleNode_Track(topNode, motherParticle, daughters, intermediates);
  }
  if (m_write_particle_container)
  {
    fillParticleNode_Particle(topNode, motherParticle, daughters, intermediates);
  }
}

void KFParticle_DST::fillParticleNode_Track(PHCompositeNode* topNode, const KFParticle& motherParticle,
                                            std::vector<KFParticle> daughters,
                                            std::vector<KFParticle> intermediates)
{
  std::string baseName;
  std::string trackNodeName;

  if (m_container_name.empty())
    baseName = "reconstructedParticles";
  else
    baseName = m_container_name;

  //Cant have forward slashes in DST or else you make a subdirectory on save!!!
  std::string undrscr = "_";
  std::string nothing = "";
  std::map<std::string, std::string> forbiddenStrings;
  forbiddenStrings["/"] = undrscr;
  forbiddenStrings["("] = undrscr;
  forbiddenStrings[")"] = nothing;
  forbiddenStrings["+"] = "plus";
  forbiddenStrings["-"] = "minus";
  forbiddenStrings["*"] = "star";
  for (auto const& [badString, goodString] : forbiddenStrings)
  {
    size_t pos;
    while ((pos = baseName.find(badString)) != std::string::npos) baseName.replace(pos, 1, goodString);
  }

  trackNodeName = baseName + "_SvtxTrackMap";

  m_recoTrackMap = findNode::getClass<SvtxTrackMap>(topNode, trackNodeName.c_str());

  SvtxTrack* m_recoTrack = new SvtxTrack_v2();

  m_recoTrack = buildSvtxTrack(motherParticle);
  m_recoTrackMap->insert(m_recoTrack);
  m_recoTrack->Reset();

  if (m_has_intermediates_DST)
  {
    KFParticle* intermediateArray = &intermediates[0];

    for (unsigned int k = 0; k < intermediates.size(); ++k)
    {
      m_recoTrack = buildSvtxTrack(intermediateArray[k]);
      m_recoTrackMap->insert(m_recoTrack);
      m_recoTrack->Reset();
    }
  }

  SvtxTrackMap* originalTrackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  KFParticle* daughterArray = &daughters[0];

  for (unsigned int k = 0; k < daughters.size(); ++k)
  {
    if (originalTrackMap->size() == 0)
    {
      std::cout << "There was no original track map found, the tracks will have no cluster information!" << std::endl;
      m_recoTrack = buildSvtxTrack(daughterArray[k]);
    }
    else
    {
      m_recoTrack = kfpTruthTools_DST.getTrack(daughterArray[k].Id(), originalTrackMap);
    }

    m_recoTrackMap->insert(m_recoTrack);
  }
}

void KFParticle_DST::fillParticleNode_Particle(PHCompositeNode* topNode, const KFParticle& motherParticle,
                                               std::vector<KFParticle> daughters,
                                               std::vector<KFParticle> intermediates)
{
  std::string baseName;
  std::string particleNodeName;

  if (m_container_name.empty())
    baseName = "reconstructedParticles";
  else
    baseName = m_container_name;

  //Cant have forward slashes in DST or else you make a subdirectory on save!!!
  std::string undrscr = "_";
  std::string nothing = "";
  std::map<std::string, std::string> forbiddenStrings;
  forbiddenStrings["/"] = undrscr;
  forbiddenStrings["("] = undrscr;
  forbiddenStrings[")"] = nothing;
  forbiddenStrings["+"] = "plus";
  forbiddenStrings["-"] = "minus";
  forbiddenStrings["*"] = "star";
  for (auto const& [badString, goodString] : forbiddenStrings)
  {
    size_t pos;
    while ((pos = baseName.find(badString)) != std::string::npos) baseName.replace(pos, 1, goodString);
  }

  particleNodeName = baseName + "_KFParticle_Container";

  m_recoParticleMap = findNode::getClass<KFParticle_Container>(topNode, particleNodeName.c_str());

  m_recoParticleMap->insert(&motherParticle);

  if (m_has_intermediates_DST)
  {
    KFParticle* intermediateArray = &intermediates[0];

    for (unsigned int k = 0; k < intermediates.size(); ++k)
      m_recoParticleMap->insert(&intermediateArray[k]);
  }

  KFParticle* daughterArray = &daughters[0];
  for (unsigned int k = 0; k < daughters.size(); ++k)
    m_recoParticleMap->insert(&daughterArray[k]);
}

SvtxTrack* KFParticle_DST::buildSvtxTrack(KFParticle particle)
{
  SvtxTrack* track = new SvtxTrack_v2();

  track->set_id(std::abs(particle.GetPDG()));
  track->set_charge((int) particle.GetQ());
  track->set_chisq(particle.GetChi2());
  track->set_ndf(particle.GetNDF());

  track->set_x(particle.GetX());
  track->set_y(particle.GetY());
  track->set_z(particle.GetZ());

  track->set_px(particle.GetPx());
  track->set_py(particle.GetPy());
  track->set_pz(particle.GetPz());

  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      track->set_error(i, j, particle.GetCovariance(i, j));

  return track;
}

void KFParticle_DST::printNode(PHCompositeNode* topNode)
{
  std::string baseName;
  std::string trackNodeName;
  std::string particleNodeName;

  if (m_container_name.empty())
    baseName = "reconstructedParticles";
  else
    baseName = m_container_name;

  //Cant have forward slashes in DST or else you make a subdirectory on save!!!
  std::string undrscr = "_";
  std::string nothing = "";
  std::map<std::string, std::string> forbiddenStrings;
  forbiddenStrings["/"] = undrscr;
  forbiddenStrings["("] = undrscr;
  forbiddenStrings[")"] = nothing;
  forbiddenStrings["+"] = "plus";
  forbiddenStrings["-"] = "minus";
  forbiddenStrings["*"] = "star";
  for (auto const& [badString, goodString] : forbiddenStrings)
  {
    size_t pos;
    while ((pos = baseName.find(badString)) != std::string::npos) baseName.replace(pos, 1, goodString);
  }

  if (m_write_track_container)
  {
    trackNodeName = baseName + "_SvtxTrackMap";
    std::cout << "----------------";
    std::cout << " KFParticle_DST: " << trackNodeName << " information ";
    std::cout << "----------------" << std::endl;
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, trackNodeName.c_str());
    for (SvtxTrackMap::Iter iter = trackmap->begin(); iter != trackmap->end(); ++iter)
    {
      SvtxTrack* track = iter->second;
      track->identify();
    }
    std::cout << "--------------------------------------------------------------------------------------------------" << std::endl;
  }

  if (m_write_particle_container)
  {
    particleNodeName = baseName + "_KFParticle_Container";
    std::cout << "----------------";
    std::cout << " KFParticle_DST: " << particleNodeName << " information ";
    std::cout << "----------------" << std::endl;
    KFParticle_Container* particlemap = findNode::getClass<KFParticle_Container>(topNode, particleNodeName.c_str());
    for (KFParticle_Container::Iter iter = particlemap->begin(); iter != particlemap->end(); ++iter)
    {
      KFParticle* particle = iter->second;
      kfpTupleTools_DST.identify(*particle);
    }
    std::cout << "--------------------------------------------------------------------------------------------------" << std::endl;
  }
}
