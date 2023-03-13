
#include "SvtxTruthRecoTableEval.h"
#include "SvtxEvalStack.h"
#include "SvtxTrackEval.h"
#include "SvtxTruthEval.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <trackbase_historic/PHG4ParticleSvtxMap_v1.h>
#include <trackbase_historic/SvtxPHG4ParticleMap_v1.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <CLHEP/Vector/ThreeVector.h>

#include <phool/PHCompositeNode.h>

//____________________________________________________________________________..
SvtxTruthRecoTableEval::SvtxTruthRecoTableEval(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
SvtxTruthRecoTableEval::~SvtxTruthRecoTableEval()
{
}

//____________________________________________________________________________..
int SvtxTruthRecoTableEval::Init(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SvtxTruthRecoTableEval::InitRun(PHCompositeNode *topNode)
{
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SvtxTruthRecoTableEval::process_event(PHCompositeNode *topNode)
{
  if (!m_svtxevalstack)
  {
    m_svtxevalstack = std::make_unique<SvtxEvalStack>(topNode); 
    m_svtxevalstack->set_strict(false);
    m_svtxevalstack->set_verbosity(Verbosity());
    m_svtxevalstack->set_use_initial_vertex(true);
    m_svtxevalstack->set_use_genfit_vertex(false);
    m_svtxevalstack->next_event(topNode);
  }
  else
  {
    m_svtxevalstack->next_event(topNode);
  }
  
  if(Verbosity() > 1)
    {std::cout << "Fill truth map "<< std::endl; }
  fillTruthMap(topNode);

  if(Verbosity() > 1)
    { std::cout << "Fill reco map "<< std::endl; }
  fillRecoMap(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SvtxTruthRecoTableEval::ResetEvent(PHCompositeNode *)
{
  if (Verbosity() > 0)
  {
    std::cout << "Truth track map " << std::endl;
    m_truthMap->identify();
    std::cout << std::endl
              << "Reco track map " << std::endl;
    m_recoMap->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SvtxTruthRecoTableEval::End(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void SvtxTruthRecoTableEval::fillTruthMap(PHCompositeNode *topNode)
{
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  assert(truthinfo);

  SvtxTrackEval *trackeval = m_svtxevalstack->get_track_eval();
  trackeval->set_verbosity(Verbosity());
  assert(trackeval);
  
  PHG4TruthInfoContainer::ConstRange range = truthinfo->GetParticleRange();
  if (m_scanForPrimaries)
  {
    range = truthinfo->GetPrimaryParticleRange();
  }

  for (auto iter = range.first; iter != range.second; ++iter)
  {
    PHG4Particle *g4particle = iter->second;

    const double momentum = CLHEP::
                                Hep3Vector(g4particle->get_px(), g4particle->get_py(), g4particle->get_pz())
                                    .mag();

    // only record particle above minimal momentum requirement.
    if (momentum < m_minMomentumTruthMap) continue;

    int gtrackID = g4particle->get_track_id();
    const std::set<SvtxTrack*>& alltracks = trackeval->all_tracks_from(g4particle);

    // not to record zero associations
    if (alltracks.size() == 0) continue;

    PHG4ParticleSvtxMap::WeightedRecoTrackMap recomap;

    for (std::set<SvtxTrack*>::iterator jter = alltracks.begin();
	 jter != alltracks.end();
	 ++jter)
      {
	SvtxTrack* track = *jter;
	/// We fill the map with a key corresponding to the ncluster contribution.
	/// This weight could in principle be anything we choose
	float clusCont = trackeval->get_nclusters_contribution(track, g4particle);
	
	auto iterator = recomap.find(clusCont);
	if (iterator == recomap.end())
	  {
	    std::set<unsigned int> dumset;
	    dumset.insert(track->get_id());
	    recomap.insert(std::make_pair(clusCont, dumset));
	  }
	else
	  {
	    iterator->second.insert(track->get_id());
	  }
      }
    
    if(Verbosity() > 1)
      { std::cout << " Inserting gtrack id " << gtrackID << " with map size " << recomap.size() << std::endl; }
    
    m_truthMap->insert(gtrackID, recomap);
    
  }
  
  m_truthMap->setProcessed(true);

}

void SvtxTruthRecoTableEval::fillRecoMap(PHCompositeNode *topNode)
{
  SvtxTrackMap *trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  assert(trackMap);

  SvtxTrackEval *trackeval = m_svtxevalstack->get_track_eval();
  assert(trackeval);

  for (const auto &[key, track] : *trackMap)
  {
    const std::set<PHG4Particle *>& allparticles = trackeval->all_truth_particles(track);
    SvtxPHG4ParticleMap::WeightedTruthTrackMap truthmap;
    for (PHG4Particle* g4particle : allparticles)
    {
      float clusCont = trackeval->get_nclusters_contribution(track, g4particle);
      auto iterator = truthmap.find(clusCont);
      if (iterator == truthmap.end())
      {
        std::set<int> dumset;
        dumset.insert(g4particle->get_track_id());
        truthmap.insert(std::make_pair(clusCont, dumset));
      }
      else
      {
        iterator->second.insert(g4particle->get_track_id());
      }
    }
    if(Verbosity() > 1)
      { std::cout << " Inserting track id " << key << " with truth map size " << truthmap.size() << std::endl; }
    m_recoMap->insert(key, truthmap);
  }

  m_recoMap->setProcessed(true);

}

int SvtxTruthRecoTableEval::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in SvtxTruthRecoTableEval::createNodes");
  }

  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_truthMap = findNode::getClass<PHG4ParticleSvtxMap>(topNode, "PHG4ParticleSvtxMap");
  if (!m_truthMap)
  {
    m_truthMap = new PHG4ParticleSvtxMap_v1;
    PHIODataNode<PHObject> *truthNode =
        new PHIODataNode<PHObject>(m_truthMap, "PHG4ParticleSvtxMap", "PHObject");
    svtxNode->addNode(truthNode);
  }

  m_recoMap = findNode::getClass<SvtxPHG4ParticleMap>(topNode, "SvtxPHG4ParticleMap");
  if (!m_recoMap)
  {
    m_recoMap = new SvtxPHG4ParticleMap_v1;
    PHIODataNode<PHObject> *recoNode =
        new PHIODataNode<PHObject>(m_recoMap, "SvtxPHG4ParticleMap", "PHObject");
    svtxNode->addNode(recoNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
