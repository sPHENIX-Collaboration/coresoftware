#include "PHG4CylinderCellContainer.h"

#include "PHG4CylinderCell.h"  // for PHG4CylinderCell
#include "PHG4CylinderCellDefs.h"
#include "PHG4CylinderCellv1.h"

#include <cstdlib>

using namespace std;

void PHG4CylinderCellContainer::Reset()
{
  while (cellmap.begin() != cellmap.end())
  {
    delete cellmap.begin()->second;
    cellmap.erase(cellmap.begin());
  }
  return;
}

void PHG4CylinderCellContainer::identify(ostream& os) const
{
  map<unsigned int, PHG4CylinderCell*>::const_iterator iter;
  os << "Number of cells: " << size() << endl;
  for (iter = cellmap.begin(); iter != cellmap.end(); ++iter)
  {
    os << "cell key 0x" << hex << iter->first << dec << endl;
    (iter->second)->identify();
  }
  set<int>::const_iterator siter;
  os << "Number of layers: " << num_layers() << endl;
  for (siter = layers.begin(); siter != layers.end(); ++siter)
  {
    os << "layer : " << *siter << endl;
  }
  return;
}

PHG4CylinderCellDefs::keytype
PHG4CylinderCellContainer::genkey(const unsigned int detid)
{
  if ((detid >> PHG4CylinderCellDefs::keybits) > 0)
  {
    cout << " detector id too large: " << detid << endl;
    exit(1);
  }
  unsigned int shiftval = detid << PHG4CylinderCellDefs::cell_idbits;
  unsigned int cellid = cellmap.size();
  cellid++;
  PHG4CylinderCellDefs::keytype newkey = cellid | shiftval;
  if (cellmap.find(newkey) != cellmap.end())
  {
    cout << " duplicate key: " << newkey << " exiting now" << endl;
    exit(1);
  }
  return newkey;
}

PHG4CylinderCellContainer::ConstIterator
PHG4CylinderCellContainer::AddCylinderCell(const unsigned int detid, PHG4CylinderCell* newcell)
{
  PHG4CylinderCellDefs::keytype key = genkey(detid);
  layers.insert(newcell->get_layer());
  newcell->set_cell_id(key);
  cellmap[key] = newcell;
  return cellmap.find(key);
}

PHG4CylinderCellContainer::ConstIterator
PHG4CylinderCellContainer::AddCylinderCellSpecifyKey(const PHG4CylinderCellDefs::keytype key, PHG4CylinderCell* newcell)
{
  if (cellmap.find(key) != cellmap.end())
  {
    cout << "PHG4CylinderCellContainer::AddCylinderCellSpecifyKey: duplicate key: " << key << " exiting now" << endl;
    exit(1);
  }
  layers.insert(newcell->get_layer());
  newcell->set_cell_id(key);
  cellmap[key] = newcell;
  return cellmap.find(key);
}

PHG4CylinderCellContainer::ConstRange
PHG4CylinderCellContainer::getCylinderCells(const unsigned int detid) const
{
  if ((detid >> PHG4CylinderCellDefs::keybits) > 0)
  {
    cout << " detector id too large: " << detid << endl;
    exit(1);
  }
  //  unsigned int shiftval = detid << cell_idbits;
  PHG4CylinderCellDefs::keytype keylow = detid << PHG4CylinderCellDefs::cell_idbits;
  PHG4CylinderCellDefs::keytype keyup = ((detid + 1) << PHG4CylinderCellDefs::cell_idbits) - 1;
  //   cout << "keylow: 0x" << hex << keylow << dec << endl;
  //   cout << "keyup: 0x" << hex << keyup << dec << endl;
  ConstRange retpair;
  retpair.first = cellmap.lower_bound(keylow);
  retpair.second = cellmap.upper_bound(keyup);
  return retpair;
}

PHG4CylinderCellContainer::ConstRange
PHG4CylinderCellContainer::getCylinderCells() const
{
  return std::make_pair(cellmap.begin(), cellmap.end());
}

PHG4CylinderCellContainer::Iterator
PHG4CylinderCellContainer::findOrAddCylinderCell(PHG4CylinderCellDefs::keytype key)
{
  PHG4CylinderCellContainer::Iterator it = cellmap.find(key);
  if (it == cellmap.end())
  {
    cellmap[key] = new PHG4CylinderCellv1();
    it = cellmap.find(key);
    PHG4CylinderCell* mcell = it->second;
    mcell->set_cell_id(key);
    layers.insert(mcell->get_layer());  // add layer to our set of layers
  }
  return it;
}

PHG4CylinderCell*
PHG4CylinderCellContainer::findCylinderCell(PHG4CylinderCellDefs::keytype key)
{
  PHG4CylinderCellContainer::ConstIterator it = cellmap.find(key);

  if (it != cellmap.end())
  {
    return it->second;
  }

  return nullptr;
}

double
PHG4CylinderCellContainer::getTotalEdep() const
{
  ConstIterator iter;
  double totalenergy = 0;
  for (iter = cellmap.begin(); iter != cellmap.end(); ++iter)
  {
    totalenergy += iter->second->get_edep();
  }
  return totalenergy;
}
