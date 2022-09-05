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

using MyNode_t = PHIODataNode<PHG4TruthInfoContainer>;

DumpPHG4TruthInfoContainer::DumpPHG4TruthInfoContainer(const std::string &NodeName)
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
    *fout << "number of G4 tracks: " << truthcontainer->size() << std::endl;
    *fout << "min trk index: " << truthcontainer->mintrkindex() << std::endl;
    *fout << "max trk index: " << truthcontainer->maxtrkindex() << std::endl;
    *fout << "number of G4 Vertices: " << truthcontainer->GetNumVertices() << std::endl;
    *fout << "min vtx index: " << truthcontainer->minvtxindex() << std::endl;
    *fout << "max vtx index: " << truthcontainer->maxvtxindex() << std::endl;
    *fout << "number of Showers: " << truthcontainer->shower_size() << std::endl;
    *fout << "min shower index: " << truthcontainer->minshowerindex() << std::endl;
    *fout << "max shower index: " << truthcontainer->maxshowerindex() << std::endl;

    PHG4TruthInfoContainer::ConstVtxIterator vtxiter;
    //      std::pair< std::map<int, PHG4VtxPoint *>::const_iterator, std::map<int, PHG4VtxPoint *>::const_iterator > vtxbegin_end = truthcontainer->GetVtxRange();
    PHG4TruthInfoContainer::ConstVtxRange vtxbegin_end = truthcontainer->GetVtxRange();

    for (vtxiter = vtxbegin_end.first; vtxiter != vtxbegin_end.second; vtxiter++)
    {
      *fout << "vtx number: " << vtxiter->first << std::endl;
      (*vtxiter->second).identify(*fout);
    }
    PHG4TruthInfoContainer::ConstIterator particle_iter;
    PHG4TruthInfoContainer::ConstRange particlebegin_end = truthcontainer->GetParticleRange();
    for (particle_iter = particlebegin_end.first; particle_iter != particlebegin_end.second; particle_iter++)
    {
      *fout << "particle number: " << particle_iter->first << std::endl;
      (particle_iter->second)->identify(*fout);
    }
    PHG4TruthInfoContainer::ConstShowerIterator shower_iter;
    PHG4TruthInfoContainer::ConstShowerRange showerbegin_end = truthcontainer->GetShowerRange();
    for (shower_iter = showerbegin_end.first; shower_iter != showerbegin_end.second; ++shower_iter)
    {
      *fout << "shower " << shower_iter->first << std::endl;
      *fout << "get_id(): " << shower_iter->second->get_id() << std::endl;
      *fout << "get_parent_particle_id(): " << shower_iter->second->get_parent_particle_id() << std::endl;
      *fout << "get_parent_shower_id(): " << shower_iter->second->get_parent_shower_id() << std::endl;
      *fout << "get_x(): " << shower_iter->second->get_x() << std::endl;
      *fout << "get_y(): " << shower_iter->second->get_y() << std::endl;
      *fout << "get_z(): " << shower_iter->second->get_z() << std::endl;
    }
    const std::pair<std::map<int, int>::const_iterator,
                    std::map<int, int>::const_iterator>
        embed_begin_end = truthcontainer->GetEmbeddedVtxIds();
    for (auto embed_iter = embed_begin_end.first; embed_iter != embed_begin_end.second; ++embed_iter)
    {
      *fout << "vtx id " << embed_iter->first << ", embed id: " << embed_iter->second << std::endl;
    }
    const std::pair<std::map<int, int>::const_iterator,
                    std::map<int, int>::const_iterator>
        embed_begin_end1 = truthcontainer->GetEmbeddedTrkIds();
    for (auto embed_iter = embed_begin_end1.first; embed_iter != embed_begin_end1.second; ++embed_iter)
    {
      *fout << "track id " << embed_iter->first << ", embed id: " << embed_iter->second << std::endl;
    }
  }
  return 0;
}
