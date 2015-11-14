#ifndef NEWGEOMCONTAINER_H__
#define NEWGEOMCONTAINER_H__

#include "RawTowerDefs.h"
#include <phool/PHObject.h>
#include <phool/phool.h>
#include <iostream>
#include <map>

class NewGeom;

class NewGeomContainer : public PHObject
{

 public:

  typedef std::map<RawTowerDefs::keytype ,NewGeom *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  NewGeomContainer( RawTowerDefs::CalorimeterId caloid );
  virtual ~NewGeomContainer();

  void Reset();
  int isValid() const;
  void identify(std::ostream& os=std::cout) const;

  void set_calorimeter_id( RawTowerDefs::CalorimeterId caloid ) { _caloid = caloid; }
  RawTowerDefs::CalorimeterId get_calorimeter_id( ) { return _caloid; }

  ConstIterator add_tower_geometry(NewGeom *geo);
  NewGeom *get_tower_geometry(RawTowerDefs::keytype key);

  //! return all tower geometries
  ConstRange get_tower_geometries( void ) const;
  Range get_tower_geometries( void );

  unsigned int size() const {return _geoms.size();}

 protected:
  RawTowerDefs::CalorimeterId _caloid;
  Map _geoms;

  ClassDef(NewGeomContainer,1)
};

#endif /* NEWGEOMCONTAINER_H__ */
