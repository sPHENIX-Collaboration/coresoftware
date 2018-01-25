#ifndef RawTowerGeomv1_h_
#define RawTowerGeomv1_h_

#include "RawTowerGeom.h"

class RawTowerGeomv1 : public RawTowerGeom {

 public:
  RawTowerGeomv1();
  RawTowerGeomv1(RawTowerDefs::keytype id);
  virtual ~RawTowerGeomv1(){}

  void identify(std::ostream& os=std::cout) const;

  void set_id(RawTowerDefs::keytype key) {_towerid = key;}
  RawTowerDefs::keytype get_id() const { return _towerid;}

  void set_center_x( double x ) { _center_x = x; return ; }
  void set_center_y( double y ) { _center_y = y; return ; }
  void set_center_z( double z ) { _center_z = z; return ; }

  double get_center_x() const { return _center_x; }
  double get_center_y() const { return _center_y; }
  double get_center_z() const { return _center_z; }


  double get_center_radius() const;
  double get_eta() const;
  double get_phi() const;

 protected:
  RawTowerDefs::keytype _towerid;

  double _center_x;
  double _center_y;
  double _center_z;
  ClassDef(RawTowerGeomv1,4)
};

#endif /* NEWGEOMV1_H_ */
