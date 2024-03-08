#include "HepMCCollisionVertex.h"

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <HepMC/SimpleVector.h>

//____________________________________________________________________________..
HepMCCollisionVertex::HepMCCollisionVertex(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int HepMCCollisionVertex::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  PHHepMCGenEventMap *geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!geneventmap)
  {
    geneventmap = new PHHepMCGenEventMap();
    PHIODataNode<PHObject> *newmapnode = new PHIODataNode<PHObject>(geneventmap, "PHHepMCGenEventMap", "PHObject");
    dstNode->addNode(newmapnode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HepMCCollisionVertex::process_event(PHCompositeNode *topNode)
{
  PHHepMCGenEventMap *genevtmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!genevtmap)
  {
    std::cout << "no PHHepMCGenEventMap node" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    std::cout << PHWHERE << " Fatal Error - GlobalVertexMap node is missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (vertexmap->empty())
  {
    std::cout << "no vertex" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  double thisvertex = -9999.;
  GlobalVertex *vtx = vertexmap->begin()->second;
  if (vtx)
  {
    thisvertex = vtx->get_z();
  }
  else
  {
    std::cout << "no vertex in map" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  PHHepMCGenEvent *genevt = new PHHepMCGenEvent();
  HepMC::FourVector collvtx(0, 0, thisvertex, 0);
  genevt->set_collision_vertex(collvtx);
  genevtmap->insert(genevt);
  return Fun4AllReturnCodes::EVENT_OK;
}
