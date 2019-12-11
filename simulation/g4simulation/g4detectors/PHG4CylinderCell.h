// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERCELL_H
#define G4DETECTORS_PHG4CYLINDERCELL_H

#include "PHG4Cell.h"

#include "PHG4CylinderCellDefs.h"

#include <g4main/PHG4HitDefs.h>

#include <phool/PHObject.h>

#include <cmath>
#include <map>

class PHG4CylinderCell : public PHG4Cell
{
 public:
/*
  typedef std::map<PHG4HitDefs::keytype, float> EdepMap;
  typedef EdepMap::iterator EdepIterator;
  typedef EdepMap::const_iterator EdepConstIterator;
  typedef std::pair<EdepIterator, EdepIterator> EdepRange;
  typedef std::pair<EdepConstIterator, EdepConstIterator> EdepConstRange;

  typedef std::map<int, float> ShowerEdepMap;
  typedef ShowerEdepMap::iterator ShowerEdepIterator;
  typedef ShowerEdepMap::const_iterator ShowerEdepConstIterator;
  typedef std::pair<ShowerEdepIterator, ShowerEdepIterator> ShowerEdepRange;
  typedef std::pair<ShowerEdepConstIterator, ShowerEdepConstIterator> ShowerEdepConstRange;
*/
  
  virtual ~PHG4CylinderCell(){}

  virtual void identify(std::ostream& os = std::cout) const {
    os << "PHG4CylinderCell base class" << std::endl;
  }
  
  
  
  virtual void set_cell_id(const PHG4CylinderCellDefs::keytype id) {return;}
  virtual void set_layer(const unsigned int i) {return;}

  virtual unsigned int get_layer() const {return 0xFFFFFFFF;}
  virtual PHG4CylinderCellDefs::keytype get_cell_id() const {return 0xFFFFFFFF;}
  virtual int get_binz() const {return -1;}
  virtual int get_binphi() const {return -1;}
  virtual int get_bineta() const {return -1;}

  virtual void set_etabin(const int i) {return;}
  virtual void set_light_yield(float lightYield)  {   return;  }

  virtual void set_fiber_ID(int fiberId) {   return;  }
  virtual int get_fiber_ID() const {return -1;}

  virtual void set_sensor_index(const std::string &si) {return;}
  virtual std::string get_sensor_index() const {return "";}

  virtual void set_ladder_phi_index(const int i) {return;}
  virtual int get_ladder_phi_index() const {return -9999;}

  virtual void set_ladder_z_index(const int i) {return;}
  virtual int get_ladder_z_index() const {return -9999;}

  virtual int get_j_index() const {return -9999;}
  virtual void set_j_index(const int i) {return;}
  virtual int get_k_index() const {return -9999;}
  virtual void set_k_index(const int i) {return;}
  virtual int get_l_index() const {return -9999;}
  virtual void set_l_index(const int i) {return;}

   protected:

  PHG4CylinderCell() {}
  ClassDef(PHG4CylinderCell,2)
};

#endif
