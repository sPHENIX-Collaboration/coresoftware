#ifndef PHG4CylinderGeomv1_H__
#define PHG4CylinderGeomv1_H__

#include "PHG4CylinderGeom.h"

#include <phool/phool.h>
#include <cmath>

class PHG4CylinderGeomv1: public PHG4CylinderGeom
{
 public:
  PHG4CylinderGeomv1();
  PHG4CylinderGeomv1(const double r, const double zmi, const double zma, const double thickn):
    layer(-1),
    radius(r),
    zmin(zmi),
    zmax(zma),
    thickness(thickn)
      {}

  virtual ~PHG4CylinderGeomv1() {}

  void identify(std::ostream& os = std::cout) const;
  int get_layer() const {return layer;}
  double get_radius() const {return radius;}
  double get_thickness() const {return thickness;}
  double get_zmin() const {return zmin;}
  double get_zmax() const {return zmax;}

  void set_layer(const int i) {layer = i;}
  void set_radius(const double r) {radius = r;}
  void set_thickness(const double t) {thickness = t;}
  void set_zmin(const double z) {zmin = z;}
  void set_zmax(const double z) {zmax = z;}

  //! load parameters from PHParameters, which interface to Database/XML/ROOT files
  virtual void ImportParameters(const PHParameters & param);
  
 protected:
  int layer;
  double radius;
  double zmin;
  double zmax;
  double thickness;

  ClassDef(PHG4CylinderGeomv1,1)
};

#endif
