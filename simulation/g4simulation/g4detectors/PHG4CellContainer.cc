#include "PHG4CellContainer.h"

#include "PHG4Cell.h"  // for PHG4Cell
#include "PHG4CellDefs.h"
#include "PHG4Cellv1.h"

#include <cstdlib>

void PHG4CellContainer::Reset()
{
  while (cellmap.begin() != cellmap.end())
  {
    delete cellmap.begin()->second;
    cellmap.erase(cellmap.begin());
  }
  return;
}

void PHG4CellContainer::identify(std::ostream& os) const
{
  os << "Number of cells: " << size() << std::endl;
  for (auto iter = cellmap.begin(); iter != cellmap.end(); ++iter)
  {
    os << "cell key 0x" << std::hex << iter->first << std::dec << std::endl;
    (iter->second)->identify();
  }
  return;
}

PHG4CellContainer::ConstIterator
PHG4CellContainer::AddCell(PHG4Cell* newcell)
{
  PHG4CellDefs::keytype key = newcell->get_cellid();
  if (cellmap.find(key) != cellmap.end())
  {
    std::cout << "overwriting cell 0x" << std::hex << key << std::dec << std::endl;
    std::cout << "layer: " << PHG4CellDefs::get_detid(key) << std::endl;
  }
  cellmap[key] = newcell;
  return cellmap.find(key);
}

PHG4CellContainer::ConstIterator
PHG4CellContainer::AddCellSpecifyKey(const PHG4CellDefs::keytype key, PHG4Cell* newcell)
{
  if (cellmap.find(key) != cellmap.end())
  {
    std::cout << "PHG4CellContainer::AddCellSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  }
  newcell->set_cellid(key);
  cellmap[key] = newcell;
  return cellmap.find(key);
}

PHG4CellContainer::ConstRange
PHG4CellContainer::getCells(const unsigned short int detid) const
{
  PHG4CellDefs::keytype tmp = detid;
  PHG4CellDefs::keytype keylow = tmp << PHG4CellDefs::bitshift_layer;
  PHG4CellDefs::keytype keyup = ((tmp + 1) << PHG4CellDefs::bitshift_layer) - 1;
  //   std::cout << "keylow: 0x" << std::hex << keylow << std::dec << std::endl;
  //   std::cout << "keyup: 0x" << std::hex << keyup << std::dec << std::endl;
  ConstRange retpair;
  retpair.first = cellmap.lower_bound(keylow);
  retpair.second = cellmap.upper_bound(keyup);
  return retpair;
}

PHG4CellContainer::ConstRange
PHG4CellContainer::getCells() const
{
  return std::make_pair(cellmap.begin(), cellmap.end());
}

PHG4CellContainer::Iterator
PHG4CellContainer::findOrAddCell(PHG4CellDefs::keytype key)
{
  PHG4CellContainer::Iterator it = cellmap.find(key);
  if (it == cellmap.end())
  {
    cellmap[key] = new PHG4Cellv1();
    it = cellmap.find(key);
    PHG4Cell* mcell = it->second;
    mcell->set_cellid(key);
  }
  return it;
}

PHG4Cell*
PHG4CellContainer::findCell(PHG4CellDefs::keytype key)
{
  PHG4CellContainer::ConstIterator it = cellmap.find(key);

  if (it != cellmap.end())
  {
    return it->second;
  }

  return nullptr;
}

double
PHG4CellContainer::getTotalEdep() const
{
  ConstIterator iter;
  double totalenergy = 0;
  for (iter = cellmap.begin(); iter != cellmap.end(); ++iter)
  {
    totalenergy += iter->second->get_edep();
  }
  return totalenergy;
}
