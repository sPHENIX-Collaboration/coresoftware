#include "DumpPHG4InEvent.h"

#include <phool/PHIODataNode.h>

#include <g4main/PHG4InEvent.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Particle.h>

#include <map>
#include <string>

using namespace std;

typedef PHIODataNode<PHG4InEvent> MyNode_t;

DumpPHG4InEvent::DumpPHG4InEvent(const string &NodeName): DumpObject(NodeName)
{
  return ;
}

int DumpPHG4InEvent::process_Node(PHNode *myNode)
{
  PHG4InEvent *phg4inevent = NULL;
  MyNode_t *thisNode = static_cast <MyNode_t *> (myNode);
  if (thisNode)
    {
      phg4inevent = thisNode->getData();
    }
  if (phg4inevent)
    {
      map<int, PHG4VtxPoint *>::const_iterator vtxiter;
      multimap<int, PHG4Particle *>::const_iterator particle_iter;
      std::pair< std::map<int, PHG4VtxPoint *>::const_iterator, std::map<int, PHG4VtxPoint *>::const_iterator > vtxbegin_end = phg4inevent->GetVertices();

      for (vtxiter = vtxbegin_end.first; vtxiter != vtxbegin_end.second; ++vtxiter)
        {
          *fout << "vtx number: " << vtxiter->first << endl;
          (*vtxiter->second).identify(*fout);
          pair<multimap<int, PHG4Particle *>::const_iterator, multimap<int, PHG4Particle *>::const_iterator > particlebegin_end = phg4inevent->GetParticles(vtxiter->first);
          for (particle_iter = particlebegin_end.first; particle_iter != particlebegin_end.second; ++particle_iter)
            {
              (particle_iter->second)->identify(*fout);
            }
        }
    }
  return 0;
}
