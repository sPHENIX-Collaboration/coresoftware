// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4VTXPOINT_H
#define G4MAIN_PHG4VTXPOINT_H

#include <phool/PHObject.h>

#include <cmath>
#include <climits>
#include <iostream>

class PHG4VtxPoint: public PHObject
{
 public:
  virtual ~PHG4VtxPoint() {}

  virtual void set_x(const double r) {}
  virtual void set_y(const double r) {}
  virtual void set_z(const double r) {}
  virtual void set_t(const double r) {}
  virtual void set_id(const int i) {}

  virtual double get_x() const {return NAN;}
  virtual double get_y() const {return NAN;}
  virtual double get_z() const {return NAN;}
  virtual double get_t() const {return NAN;}
  virtual int get_id() const {return INT_MIN;}

  virtual void identify(std::ostream& os = std::cout) const;

  //! comparison of vertex value only, not on the id, per algorithm requirement in PHG4TruthInfoContainer::AddPrimaryVertex
  bool  operator== (const PHG4VtxPoint &) const ;


 protected:
  PHG4VtxPoint(){}
  ClassDef(PHG4VtxPoint,1)

};

#endif
