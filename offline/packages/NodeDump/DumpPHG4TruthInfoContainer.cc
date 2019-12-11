#include "DumpPHG4TruthInfoContainer.h"

#include <phool/PHIODataNode.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
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
    *fout << "min trk index: " << truthcontainer->mintrkindex() << endl;
    *fout << "max trk index: " << truthcontainer->maxtrkindex() << endl;
    *fout << "number of G4 Vertices: " << truthcontainer->GetNumVertices() << endl;
    *fout << "min vtx index: " << truthcontainer->minvtxindex() << endl;
    *fout << "max vtx index: " << truthcontainer->maxvtxindex() << endl;
    *fout << "number of Showers: " << truthcontainer->shower_size() << endl;
    *fout << "min shower index: " << truthcontainer->minshowerindex() << endl;
    *fout << "max shower index: " << truthcontainer->maxshowerindex() << endl;

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
    PHG4TruthInfoContainer::ConstShowerIterator shower_iter;
    PHG4TruthInfoContainer::ConstShowerRange showerbegin_end = truthcontainer->GetShowerRange();
    for (shower_iter = showerbegin_end.first; shower_iter != showerbegin_end.second; ++shower_iter)
    {
      *fout << "shower " << shower_iter->first << endl;
      *fout << "get_id(): " << shower_iter->second->get_id() << endl;
      *fout << "get_parent_particle_id(): " << shower_iter->second->get_parent_particle_id() << endl;
      *fout << "get_parent_shower_id(): " << shower_iter->second->get_parent_shower_id() << endl;
      *fout << "get_x(): " << shower_iter->second->get_x() << endl;
      *fout << "get_y(): " << shower_iter->second->get_y() << endl;
      *fout << "get_z(): " << shower_iter->second->get_z() << endl;
    }
  }
  return 0;
}
