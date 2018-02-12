#ifndef RawTowerGeomv1_h_
#define RawTowerGeomv1_h_

#include "RawTowerGeom.h"

class RawTowerGeomv1 : public RawTowerGeom {

 public:
  RawTowerGeomv1();
  RawTowerGeomv1(RawTowerDefs::keytype id);
  virtual ~RawTowerGeomv1(){}

  virtual void identify(std::ostream& os=std::cout) const;

  void set_id(RawTowerDefs::keytype key) {_towerid = key;}
  RawTowerDefs::keytype get_id() const { return _towerid;}

  //! x of geometric center of tower
  void set_center_x( double x ) { _center_x = x; return ; }
  //! y of geometric center of tower
  void set_center_y( double y ) { _center_y = y; return ; }
  //! z of geometric center of tower
  void set_center_z( double z ) { _center_z = z; return ; }

  //! x of geometric center of tower
  double get_center_x() const { return _center_x; }
  //! y of geometric center of tower
  double get_center_y() const { return _center_y; }
  //! z of geometric center of tower
  double get_center_z() const { return _center_z; }


  double get_center_radius() const;
  double get_eta() const;
  double get_phi() const;


  //! x component of the length vector of the tower from front to back, and the amplitude represents the full length of the tower
  void set_tower_lenth_vec_x( double x ) { _tower_lenth_vec_x = x; return ; }
  //! y component of the length vector of the tower from front to back, and the amplitude represents the full length of the tower
  void set_tower_lenth_vec_y( double y ) { _tower_lenth_vec_y = y; return ; }
  //! z component of the length vector of the tower from front to back, and the amplitude represents the full length of the tower
  void set_tower_lenth_vec_z( double z ) { _tower_lenth_vec_z = z; return ; }

  //! x component of the length vector of the tower from front to back, and the amplitude represents the full length of the tower
  double get_tower_lenth_vec_x() const { return _tower_lenth_vec_x; }
  //! y component of the length vector of the tower from front to back, and the amplitude represents the full length of the tower
  double get_tower_lenth_vec_y() const { return _tower_lenth_vec_y; }
  //! z component of the length vector of the tower from front to back, and the amplitude represents the full length of the tower
  double get_tower_lenth_vec_z() const { return _tower_lenth_vec_z; }

  virtual double get_front_x() const ;
  virtual double get_front_y() const ;
  virtual double get_front_z() const ;


 protected:
  RawTowerDefs::keytype _towerid;

  //! x of geometric center of tower
  double _center_x;
  //! y of geometric center of tower
  double _center_y;
  //! z of geometric center of tower
  double _center_z;

  //! x component of the length vector of the tower from front to back, and the amplitude represents the full length of the tower
  double _tower_lenth_vec_x;
  //! y component of the length vector of the tower from front to back, and the amplitude represents the full length of the tower
  double _tower_lenth_vec_y;
  //! z component of the length vector of the tower from front to back, and the amplitude represents the full length of the tower
  double _tower_lenth_vec_z;

  ClassDef(RawTowerGeomv1,4)
};

#endif /* NEWGEOMV1_H_ */
