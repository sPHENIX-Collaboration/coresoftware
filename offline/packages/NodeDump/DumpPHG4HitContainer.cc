#include "DumpPHG4HitContainer.h"

#include <phool/PHIODataNode.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <string>

using namespace std;

typedef PHIODataNode<PHG4HitContainer> MyNode_t;

DumpPHG4HitContainer::DumpPHG4HitContainer(const string &NodeName): DumpObject(NodeName)
{
  return ;
}

int DumpPHG4HitContainer::process_Node(PHNode *myNode)
{
  PHG4HitContainer *phg4hitcontainer = NULL;
  MyNode_t *thisNode = static_cast <MyNode_t *> (myNode);
  if (thisNode)
    {
      phg4hitcontainer = thisNode->getData();
    }
  if (phg4hitcontainer)
    {
      PHG4HitContainer::ConstIterator hiter;
      PHG4HitContainer::ConstRange hit_begin_end = phg4hitcontainer->getHits();
      *fout << "size: " << phg4hitcontainer->size() << endl;
      *fout << "num layers: " << phg4hitcontainer->num_layers() << endl;
      for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; hiter++)
        {
          *fout << "id: " << hiter->second->get_hit_id() << endl;
          *fout << "layer: " << hiter->second->get_layer() << endl;
          *fout << "trkid: " << hiter->second->get_trkid() << endl;
          *fout << "edep: " << hiter->second->get_edep() << endl;
          for (int i = 0; i < 2; i++)
            {
              *fout << "x(" << i << "): " << hiter->second->get_x(i) << endl;
              *fout << "y(" << i << "): " << hiter->second->get_y(i) << endl;
              *fout << "z(" << i << "): " << hiter->second->get_z(i) << endl;
              *fout << "t(" << i << "): " << hiter->second->get_t(i) << endl;
            }
        }
    }
  return 0;
}

