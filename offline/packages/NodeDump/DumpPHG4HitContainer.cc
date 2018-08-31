#include "DumpPHG4HitContainer.h"

#include <phool/PHIODataNode.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <string>
#include <stdint.h>

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
	  *fout << "id: 0x" << hex << hiter->second->get_hit_id() << dec << endl;
          *fout << "detid: " << hiter->second->get_detid() << endl;
          *fout << "trkid: " << hiter->second->get_trkid() << endl;
          *fout << "edep: " << hiter->second->get_edep() << endl;
          for (int i = 0; i < 2; i++)
            {
              *fout << "x(" << i << "): " << hiter->second->get_x(i) << endl;
              *fout << "y(" << i << "): " << hiter->second->get_y(i) << endl;
              *fout << "z(" << i << "): " << hiter->second->get_z(i) << endl;
              *fout << "t(" << i << "): " << hiter->second->get_t(i) << endl;
            }
	   for (unsigned char ic = 0; ic < UCHAR_MAX; ic++)
	     {
	       PHG4Hit::PROPERTY prop_id = static_cast<PHG4Hit::PROPERTY> (ic);
	       if (hiter->second->has_property(prop_id))
		 {
		   *fout << "prop id: " << static_cast<unsigned int> (ic);
		   pair<const string, PHG4Hit::PROPERTY_TYPE> property_info = PHG4Hit::get_property_info(prop_id);
		   *fout << ", name " << property_info.first << " value ";
		   switch(property_info.second)
		     {
		     case PHG4Hit::type_int:
		       *fout << hiter->second->get_property_int(prop_id);
		       break;
		     case PHG4Hit::type_uint:
		       *fout << hiter->second->get_property_uint(prop_id);
		       break;
		     case PHG4Hit::type_float:
		       *fout << hiter->second->get_property_float(prop_id);
		       break;
		     default:
		       *fout << " unknown type ";
		     }
		   *fout <<endl;
		 }
	     }
	}
    }
  return 0;
}

