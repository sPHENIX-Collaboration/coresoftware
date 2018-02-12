#ifndef NEWGEOM_H_
#define NEWGEOM_H_

#include "RawTowerDefs.h"

#include <phool/phool.h>
#include <phool/PHObject.h>
#include <iostream>
#include <map>
#include <cmath>

class RawTowerGeom : public PHObject {

 public:

  virtual ~RawTowerGeom() {}

  virtual void identify(std::ostream& os=std::cout) const;

  virtual void set_id(RawTowerDefs::keytype key)  { PHOOL_VIRTUAL_WARN("set_id()");  }
  virtual RawTowerDefs::keytype get_id() const { PHOOL_VIRTUAL_WARN("get_id()"); return 0; }

  //! x of geometric center of tower
  virtual void set_center_x( double ) { PHOOL_VIRTUAL_WARN("set_center_x()"); return ; }
  //! y of geometric center of tower
  virtual void set_center_y( double ) { PHOOL_VIRTUAL_WARN("set_center_y()"); return ; }
  //! z of geometric center of tower
  virtual void set_center_z( double ) { PHOOL_VIRTUAL_WARN("set_center_z()"); return ; }

  virtual void set_size_x( double ) { PHOOL_VIRTUAL_WARN("set_size_x()"); return ; }
  virtual void set_size_y( double ) { PHOOL_VIRTUAL_WARN("set_size_y()"); return ; }
  virtual void set_size_z( double ) { PHOOL_VIRTUAL_WARN("set_size_z()"); return ; }

  //! x of geometric center of tower
  virtual double get_center_x() const { PHOOL_VIRTUAL_WARN("get_center_x()"); return NAN; }
  //! y of geometric center of tower
  virtual double get_center_y() const { PHOOL_VIRTUAL_WARN("get_center_y()"); return NAN; }
  //! z of geometric center of tower
  virtual double get_center_z() const { PHOOL_VIRTUAL_WARN("get_center_z()"); return NAN; }

  virtual double get_size_x() const { PHOOL_VIRTUAL_WARN("get_size_x()"); return NAN; }
  virtual double get_size_y() const { PHOOL_VIRTUAL_WARN("get_size_y()"); return NAN; }
  virtual double get_size_z() const { PHOOL_VIRTUAL_WARN("get_size_z()"); return NAN; }
  virtual double get_volume() const { PHOOL_VIRTUAL_WARN("get_volume()"); return NAN; }

  virtual double get_center_radius() const { PHOOL_VIRTUAL_WARN("get_center_radius()"); return NAN; }
  virtual double get_eta() const { PHOOL_VIRTUAL_WARN("get_eta()"); return NAN; }
  virtual double get_phi() const { PHOOL_VIRTUAL_WARN("get_phi()"); return NAN; }

  virtual void set_tower_type( int ) { PHOOL_VIRTUAL_WARN("set_tower_type()"); return ; }
  virtual int get_tower_type() const { PHOOL_VIRTUAL_WARN("get_tower_type()"); return NAN; }

  //! x component of the length vector of the tower from front to back, and the amplitude represents the full length of the tower
  virtual void set_tower_lenth_vec_x( double x ) { PHOOL_VIRTUAL_WARN("set_tower_lenth_vec_x()"); return ; }
  //! y component of the length vector of the tower from front to back, and the amplitude represents the full length of the tower
  virtual void set_tower_lenth_vec_y( double y ) { PHOOL_VIRTUAL_WARN("set_tower_lenth_vec_y()"); return ; }
  //! z component of the length vector of the tower from front to back, and the amplitude represents the full length of the tower
  virtual void set_tower_lenth_vec_z( double z ) { PHOOL_VIRTUAL_WARN("set_tower_lenth_vec_z()"); return ; }

  //! x component of the length vector of the tower from front to back, and the amplitude represents the full length of the tower
  virtual double get_tower_lenth_vec_x() const { PHOOL_VIRTUAL_WARN("get_tower_lenth_vec_x()"); return NAN; }
  //! y component of the length vector of the tower from front to back, and the amplitude represents the full length of the tower
  virtual double get_tower_lenth_vec_y() const { PHOOL_VIRTUAL_WARN("get_tower_lenth_vec_y()"); return NAN; }
  //! z component of the length vector of the tower from front to back, and the amplitude represents the full length of the tower
  virtual double get_tower_lenth_vec_z() const { PHOOL_VIRTUAL_WARN("get_tower_lenth_vec_z()"); return NAN; }

  virtual double get_front_x() const { PHOOL_VIRTUAL_WARN("get_front_x()"); return NAN; }
  virtual double get_front_y() const { PHOOL_VIRTUAL_WARN("get_front_y()"); return NAN; }
  virtual double get_front_z() const { PHOOL_VIRTUAL_WARN("get_front_z()"); return NAN; }
 protected:
  RawTowerGeom() {}

  ClassDef(RawTowerGeom,2)

};

#endif /* NEWGEOM_H_ */
