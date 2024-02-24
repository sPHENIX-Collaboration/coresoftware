// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4VTXPOINT_H
#define G4MAIN_PHG4VTXPOINT_H

#include <phool/PHObject.h>

#include <iostream>
#include <limits>

class PHG4VtxPoint : public PHObject
{
 public:
  ~PHG4VtxPoint() override {}

  void identify(std::ostream& os = std::cout) const override;

  virtual void set_x(const double) {}
  virtual void set_y(const double) {}
  virtual void set_z(const double) {}
  virtual void set_t(const double) {}
  virtual void set_id(const int) {}

  virtual double get_x() const { return std::numeric_limits<double>::quiet_NaN(); }
  virtual double get_y() const { return std::numeric_limits<double>::quiet_NaN(); }
  virtual double get_z() const { return std::numeric_limits<double>::quiet_NaN(); }
  virtual double get_t() const { return std::numeric_limits<double>::quiet_NaN(); }
  virtual int get_id() const { return std::numeric_limits<int>::min(); }

  //! comparison of vertex value only, not on the id, per algorithm requirement in PHG4TruthInfoContainer::AddPrimaryVertex
  bool operator==(const PHG4VtxPoint&) const;

 protected:
  PHG4VtxPoint() {}
  ClassDefOverride(PHG4VtxPoint, 1)
};

#endif
