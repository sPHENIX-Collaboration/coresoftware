#include "DumpPHG4CylinderCellContainer.h"

#include <phool/PHIODataNode.h>

#include <g4detectors/PHG4CylinderCell.h>
#include <g4detectors/PHG4CylinderCellContainer.h>

#include <string>

using namespace std;

typedef PHIODataNode<PHG4CylinderCellContainer> MyNode_t;

DumpPHG4CylinderCellContainer::DumpPHG4CylinderCellContainer(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHG4CylinderCellContainer::process_Node(PHNode *myNode)
{
  PHG4CylinderCellContainer *phg4cellcontainer = NULL;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
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
      int tmp = hiter->second->get_binz();
      if (tmp != -1)
      {
        *fout << "binz: " << tmp << endl;
      }
      tmp = hiter->second->get_binphi();
      if (tmp != -1)
      {
        *fout << "binphi: " << tmp << endl;
      }
      tmp = hiter->second->get_bineta();
      if (tmp != -1)
      {
        *fout << "bineta: " << tmp << endl;
      }
      float dtmp = hiter->second->get_light_yield();
      if (isfinite(dtmp))
      {
        *fout << "get_light_yield(): " << dtmp << endl;
      }
      tmp = hiter->second->get_fiber_ID();
      if (tmp != -1)
      {
        *fout << "get_fiber_ID(): " << hiter->second->get_fiber_ID() << endl;
      }
      string tmpstr = hiter->second->get_sensor_index();
      if (!tmpstr.empty())
      {
        *fout << "get_sensor_index(): " << hiter->second->get_sensor_index() << endl;
      }
      tmp = hiter->second->get_ladder_phi_index();
      if (tmp != -9999)
      {
        *fout << "get_ladder_phi_index(): " << hiter->second->get_ladder_phi_index() << endl;
      }
      tmp = hiter->second->get_ladder_z_index();
      if (tmp != -9999)
      {
        *fout << "get_ladder_z_index(): " << hiter->second->get_ladder_z_index() << endl;
      }
      tmp = hiter->second->get_j_index();
      if (tmp != -9999)
      {
        *fout << "get_j_index(): " << hiter->second->get_j_index() << endl;
      }
      tmp = hiter->second->get_k_index();
      if (tmp != -9999)
      {
        *fout << "get_k_index(): " << hiter->second->get_k_index() << endl;
      }
      tmp = hiter->second->get_l_index();
      if (tmp != -9999)
      {
        *fout << "get_l_index(): " << hiter->second->get_l_index() << endl;
      }
      tmp = hiter->second->get_pixel_index();
      if (tmp != -9999)
      {
        *fout << "get_pixel_index(): " << hiter->second->get_pixel_index() << endl;
      }
      tmp = hiter->second->get_chip_index();
      if (tmp != -9999)
      {
        *fout << "get_chip_index(): " << hiter->second->get_chip_index() << endl;
      }
      tmp = hiter->second->get_module_index();
      if (tmp != -9999)
      {
        *fout << "get_module_index(): " << hiter->second->get_module_index() << endl;
      }
      tmp = hiter->second->get_half_stave_index();
      if (tmp != -9999)
      {
        *fout << "get_half_stave_index(): " << hiter->second->get_half_stave_index() << endl;
      }
      tmp = hiter->second->get_stave_index();
      if (tmp != -9999)
      {
        *fout << "get_stave_index(): " << hiter->second->get_stave_index() << endl;
      }

      PHG4CylinderCell::EdepConstRange hitedep_begin_end = hiter->second->get_g4hits();
      for (PHG4CylinderCell::EdepConstIterator iter = hitedep_begin_end.first; iter != hitedep_begin_end.second; ++iter)
      {
        *fout << "hit 0x" << hex << iter->first << dec << " edep: " << iter->second << endl;
      }
      PHG4CylinderCell::ShowerEdepConstRange shower_begin_end = hiter->second->get_g4showers();
      for (PHG4CylinderCell::ShowerEdepConstIterator iter = shower_begin_end.first; iter != shower_begin_end.second; ++iter)
      {
        *fout << "shower " << iter->first << " edep: " << iter->second << endl;
      }
    }
  }
  return 0;
}
