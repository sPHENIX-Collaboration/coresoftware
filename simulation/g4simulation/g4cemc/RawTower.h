#ifndef RAWTOWER_H_
#define RAWTOWER_H_

#include "RawTowerDefs.h"
#include <g4detectors/PHG4CylinderCellDefs.h>
#include <phool/phool.h>
#include <phool/PHObject.h>
#include <iostream>
#include <map>

class RawTower : public PHObject {

 public:
  typedef std::map<PHG4CylinderCellDefs::keytype, float> CellMap;
  typedef CellMap::iterator CellIterator;
  typedef CellMap::const_iterator CellConstIterator;
  typedef std::pair<CellIterator, CellIterator> CellRange;
  typedef std::pair<CellConstIterator, CellConstIterator> CellConstRange;


  virtual ~RawTower() {}

  virtual void Reset() { PHOOL_VIRTUAL_WARNING; }
  virtual int isValid() const { PHOOL_VIRTUAL_WARN("isValid()"); return 0; }
  virtual void identify(std::ostream& os=std::cout) const { PHOOL_VIRTUAL_WARN("identify()"); }

  virtual RawTowerDefs::keytype get_id() const { PHOOL_VIRTUAL_WARN("get_id()"); return 0; }
  virtual int get_bineta() const { PHOOL_VIRTUAL_WARN("get_ieta()"); return -1; }
  virtual int get_binphi() const { PHOOL_VIRTUAL_WARN("get_iphi()"); return -1; }

  //! energy assigned to the tower. Depending on stage of process and DST node name, it could be energy deposition, light yield or calibrated energies
  virtual double get_energy() const { PHOOL_VIRTUAL_WARN("get_energy()"); return 0.0; }
  //! energy assigned to the tower. Depending on stage of process and DST node name, it could be energy deposition, light yield or calibrated energies
  virtual void set_energy (const double ) { PHOOL_VIRTUAL_WARN("set_energy()"); return  ; }

  virtual CellConstRange get_g4cells() const
  {
    PHOOL_VIRTUAL_WARN("get_g4cells()");
    CellMap dummy;
    return make_pair(dummy.begin(), dummy.end());
  }

  virtual void add_ecell(const PHG4CylinderCellDefs::keytype  g4cellid, const float ecell) {PHOOL_VIRTUAL_WARN("add_ecell(const PHG4CylinderCellDefs::keytype g4cellid, const float ecell)"); return;}

 protected:
  RawTower() {}

  ClassDef(RawTower,1)

};

#endif /* RAWTOWER_H_ */
