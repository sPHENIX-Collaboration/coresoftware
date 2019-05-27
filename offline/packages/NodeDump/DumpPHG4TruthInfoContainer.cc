#include "DumpPHG4TruthInfoContainer.h"

#include <phool/PHIODataNode.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using namespace std;

typedef PHIODataNode<PHG4TruthInfoContainer> MyNode_t;

DumpPHG4TruthInfoContainer::DumpPHG4TruthInfoContainer(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHG4TruthInfoContainer::process_Node(PHNode *myNode)
{
  PHG4TruthInfoContainer *truthcontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    truthcontainer = thisNode->getData();
  }
  if (truthcontainer)
  {
    *fout << "number of G4 tracks: " << truthcontainer->size() << endl;
    *fout << "number of G4 Vertices: " << truthcontainer->GetNumVertices() << endl;
    PHG4TruthInfoContainer::ConstVtxIterator vtxiter;
    //      std::pair< std::map<int, PHG4VtxPoint *>::const_iterator, std::map<int, PHG4VtxPoint *>::const_iterator > vtxbegin_end = truthcontainer->GetVtxRange();
    PHG4TruthInfoContainer::ConstVtxRange vtxbegin_end = truthcontainer->GetVtxRange();

    for (vtxiter = vtxbegin_end.first; vtxiter != vtxbegin_end.second; vtxiter++)
    {
      *fout << "vtx number: " << vtxiter->first << endl;
      (*vtxiter->second).identify(*fout);
    }
    PHG4TruthInfoContainer::ConstIterator particle_iter;
    PHG4TruthInfoContainer::ConstRange particlebegin_end = truthcontainer->GetParticleRange();
    for (particle_iter = particlebegin_end.first; particle_iter != particlebegin_end.second; particle_iter++)
    {
      *fout << "particle number: " << particle_iter->first << endl;
      (particle_iter->second)->identify(*fout);
    }
  }
  return 0;
}
