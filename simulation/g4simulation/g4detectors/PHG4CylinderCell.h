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
  
  virtual ~PHG4CylinderCell(){}

  virtual void identify(std::ostream& os = std::cout) const {
    os << "PHG4CylinderCell base class" << std::endl;
  }
  
  virtual EdepConstRange get_g4hits() {
    std::map <PHG4HitDefs::keytype, float> dummy;
    return std::make_pair(dummy.begin(), dummy.end());
  }
  
  virtual void add_edep(const PHG4HitDefs::keytype g4hitid, const float edep) {return;}
  virtual void add_edep(const PHG4HitDefs::keytype g4hitid, const float edep, const float light_yield) {return;}

  virtual ShowerEdepConstRange get_g4showers() {
    std::map <int, float> dummy;
    return std::make_pair(dummy.begin(), dummy.end());
  }
  
  virtual void add_shower_edep(const int g4showerid, const float edep) {return;}
  
  virtual void set_cell_id(const PHG4CylinderCellDefs::keytype id) {return;}
  virtual void set_layer(const unsigned int i) {return;}

  virtual double get_edep() const {return NAN;}
  virtual unsigned int get_layer() const {return 0xFFFFFFFF;}
  virtual PHG4CylinderCellDefs::keytype get_cell_id() const {return 0xFFFFFFFF;}
  virtual int get_binz() const {return -1;}
  virtual int get_binphi() const {return -1;}
  virtual int get_bineta() const {return -1;}
  virtual float  get_light_yield() const  {    return NAN;  }
  virtual int    get_fiber_ID() const {return -1;}

  virtual void set_phibin(const int i) {return;}
  virtual void set_zbin(const int i) {return;}
  virtual void set_etabin(const int i) {return;}
  virtual void  set_light_yield(float lightYield)  {   return;  }
  virtual void  set_fiber_ID(int fiberId) {   return;  }

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

  virtual int get_pixel_index() const {return -9999;}
  virtual int get_chip_index() const {return -9999;}
  virtual int get_module_index() const {return -9999;}
  virtual int get_half_stave_index() const {return -9999;}
  virtual int get_stave_index() const {return -9999;}

  virtual void set_stave_index(const  int si) {return;}
  virtual void set_half_stave_index(const int i) {return;}
  virtual void set_module_index(const int i) {return;}
  virtual void set_chip_index(const int i) {return;}
  virtual void set_pixel_index(const int i) {return;}
  
 protected:

  PHG4CylinderCell() {}
  ClassDef(PHG4CylinderCell,2)
};

#endif
