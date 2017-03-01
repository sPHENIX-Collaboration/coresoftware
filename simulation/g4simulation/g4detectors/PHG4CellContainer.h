#ifndef PHG4CELLCONTAINER_H__
#define PHG4CELLCONTAINER_H__

#include "PHG4Cell.h"

#include <phool/PHObject.h>

#include <map>
#include <set>

class PHG4CellContainer: public PHObject
{

  public:
  typedef std::map<PHG4CellDefs::keytype,PHG4Cell *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  PHG4CellContainer();

  virtual ~PHG4CellContainer() {}

  void Reset();

  void identify(std::ostream& os = std::cout) const;

  ConstIterator AddCell(PHG4Cell *newCell);
  ConstIterator AddCellSpecifyKey(const PHG4CellDefs::keytype key, PHG4Cell *newCell);
  
  //! preferred removal method, key is currently the cell id
  void RemoveCell(PHG4CellDefs::keytype key) {
    cellmap.erase(key);
  }

  //! inefficent, use key where possible instead
  void RemoveCell(PHG4Cell *cell)
  {
    Iterator its = cellmap.begin();
    Iterator last = cellmap.end();
    for (; its != last;)
      {
	if (its->second == cell)
	  {
	    cellmap.erase(its++);
	  }
	else
	  {
	    ++its;
	  }
      }
  }


  Iterator findOrAddCell(PHG4CellDefs::keytype key);

  //! return all Cells matching a given detid
  ConstRange getCells(const unsigned short int detid) const;

  //! return all hist
  ConstRange getCells( void ) const;

  PHG4Cell* findCell(PHG4CellDefs::keytype key);

  unsigned int size( void ) const
  { return cellmap.size(); }

  double getTotalEdep() const;

 protected:
  Map cellmap;
  ClassDef(PHG4CellContainer,1)
};

#endif
