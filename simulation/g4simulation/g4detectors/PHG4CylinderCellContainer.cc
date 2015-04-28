#include "PHG4CylinderCellContainer.h"
#include "PHG4CylinderCell.h"
#include "PHG4CylinderCellv1.h"
#include "PHG4CylinderCellDefs.h"

#include <cstdlib>

using namespace std;

ClassImp(PHG4CylinderCellContainer)

PHG4CylinderCellContainer::PHG4CylinderCellContainer()
{
}

void
PHG4CylinderCellContainer::Reset()
{
   while(cellmap.begin() != cellmap.end())
     {
       delete cellmap.begin()->second;
       cellmap.erase(cellmap.begin());
     }
  return;
}

void
PHG4CylinderCellContainer::identify(ostream& os) const
{
   map<unsigned int,PHG4CylinderCell *>::const_iterator iter;
   os << "Number of cells: " << size() << endl;
   for (iter = cellmap.begin(); iter != cellmap.end(); iter++)
     {
       os << "cell key 0x" << hex << iter->first << dec << endl;
       (iter->second)->identify();
     }
   set<int>::const_iterator siter;
   os << "Number of layers: " << num_layers() << endl;
   for (siter = layers.begin(); siter != layers.end(); siter++)
     {
       os << "layer : " << *siter << endl;
     }
  return;
}

unsigned int
PHG4CylinderCellContainer::genkey(const int detid)
{
  if ((detid >> phg4cylindercelldefs::keybits) > 0)
    {
      cout << " detector id too large: " << detid << endl;
      exit(1);
    }
  int shiftval = detid << phg4cylindercelldefs::cell_idbits;
  int cellid = cellmap.size();
  cellid++;
  int newkey = cellid | shiftval;
  if (cellmap.find(newkey) != cellmap.end())
    {
      cout << " duplicate key: " << newkey << " exiting now" << endl;
      exit(1);
    }
  return newkey;
}

std::map<unsigned int,PHG4CylinderCell *>::const_iterator
PHG4CylinderCellContainer::AddCylinderCell(const int detid, PHG4CylinderCell *newcell)
{
  unsigned int key = genkey(detid);
  layers.insert(newcell->get_layer());
  newcell->set_cell_id(key);
  cellmap[key] = newcell;
  return cellmap.find(key);
}

PHG4CylinderCellContainer::ConstRange PHG4CylinderCellContainer::getCylinderCells(const int detid) const
{
  if ((detid >> phg4cylindercelldefs::keybits) > 0)
    {
      cout << " detector id too large: " << detid << endl;
      exit(1);
    }
  //  unsigned int shiftval = detid << cell_idbits;
  unsigned int keylow = detid << phg4cylindercelldefs::cell_idbits;
  unsigned int keyup = ((detid + 1)<< phg4cylindercelldefs::cell_idbits) -1 ;
//   cout << "keylow: 0x" << hex << keylow << dec << endl;
//   cout << "keyup: 0x" << hex << keyup << dec << endl;
  ConstRange retpair;
  retpair.first = cellmap.lower_bound(keylow);
  retpair.second = cellmap.upper_bound(keyup);
  return retpair;
}

PHG4CylinderCellContainer::ConstRange PHG4CylinderCellContainer::getCylinderCells( void ) const
{ return std::make_pair( cellmap.begin(), cellmap.end() ); }


PHG4CylinderCellContainer::Iterator PHG4CylinderCellContainer::findOrAddCylinderCell(unsigned int key)
{
  PHG4CylinderCellContainer::Iterator it = cellmap.find(key);
  if(it == cellmap.end())
  {
    cellmap[key] = new PHG4CylinderCellv1();
    it = cellmap.find(key);
    PHG4CylinderCell* mcell = it->second;
    mcell->set_cell_id(key);
    layers.insert(mcell->get_layer()); // add layer to our set of layers
  }
  return it;
}

PHG4CylinderCell* PHG4CylinderCellContainer::findCylinderCell(unsigned int key)
{
  PHG4CylinderCellContainer::ConstIterator it = cellmap.find(key);

  if(it != cellmap.end())
    {
      return it->second;
    }

  return 0;
}

double
PHG4CylinderCellContainer::getTotalEdep() const
{
  ConstIterator iter;
  double totalenergy = 0;
  for (iter = cellmap.begin(); iter != cellmap.end(); iter++)
    {
      totalenergy += iter->second->get_edep();
    }
  return totalenergy;
}
