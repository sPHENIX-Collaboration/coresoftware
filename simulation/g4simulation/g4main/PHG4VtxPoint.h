#ifndef PHG4VTXPOINT_H__
#define PHG4VTXPOINT_H__

#include <phool/PHObject.h>
#include <cmath>

class PHG4VtxPoint: public PHObject
{
 public:
  virtual ~PHG4VtxPoint() {}

  virtual void set_x(const double r) {}
  virtual void set_y(const double r) {}
  virtual void set_z(const double r) {}
  virtual void set_t(const double r) {}

  virtual double get_x() const {return NAN;}
  virtual double get_y() const {return NAN;}
  virtual double get_z() const {return NAN;}
  virtual double get_t() const {return NAN;}

  virtual void identify(std::ostream& os = std::cout) const;

  bool  operator== (const PHG4VtxPoint &) const ;


 protected:
  PHG4VtxPoint(){}
  ClassDef(PHG4VtxPoint,1)

};

#endif
