#include "DumpPHG4CylinderCellContainer.h"

#include <phool/PHIODataNode.h>

#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCell.h>

#include <string>

using namespace std;

typedef PHIODataNode<PHG4CylinderCellContainer> MyNode_t;

DumpPHG4CylinderCellContainer::DumpPHG4CylinderCellContainer(const string &NodeName): DumpObject(NodeName)
{
  return ;
}

int DumpPHG4CylinderCellContainer::process_Node(PHNode *myNode)
{
  PHG4CylinderCellContainer *phg4cellcontainer = NULL;
  MyNode_t *thisNode = static_cast <MyNode_t *> (myNode);
  if (thisNode)
    {
      phg4cellcontainer = thisNode->getData();
    }
  if (phg4cellcontainer)
    {
      PHG4CylinderCellContainer::ConstIterator hiter;
      PHG4CylinderCellContainer::ConstRange cell_begin_end = phg4cellcontainer->getCylinderCells();
      *fout << "size: " << phg4cellcontainer->size() << endl;
      for (hiter = cell_begin_end.first; hiter != cell_begin_end.second; hiter++)
        {
          *fout << "id: 0x" << hex << hiter->second->get_cell_id() << dec << endl;
          *fout << "layer: " << hiter->second->get_layer() << endl;
          *fout << "edep: " << hiter->second->get_edep() << endl;
          *fout << "binz: " << hiter->second->get_binz() << endl;
          *fout << "binphi: " << hiter->second->get_binphi() << endl;
          *fout << "bineta: " << hiter->second->get_bineta() << endl;
	   PHG4CylinderCell::EdepConstRange hitedep_begin_end = hiter->second->get_g4hits();
	   for (PHG4CylinderCell::EdepConstIterator iter = hitedep_begin_end.first; iter != hitedep_begin_end.second; ++iter)
	     {
	       *fout << "hit 0x" << hex << iter->first << dec << " edep: " << iter->second << endl;
	     }
	   PHG4CylinderCell::ShowerEdepConstRange shower_begin_end = hiter->second->get_g4showers();
	   for (PHG4CylinderCell::ShowerEdepConstIterator iter = shower_begin_end.first; iter != shower_begin_end.second; ++iter)
	     {
	       *fout << "shower " << iter->first  << " edep: " << iter->second << endl;
	     }
        }
    }
  return 0;
}

