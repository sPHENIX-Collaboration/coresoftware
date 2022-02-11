// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4SCINTILLATORSLAT_H
#define G4DETECTORS_PHG4SCINTILLATORSLAT_H

#include "PHG4ScintillatorSlatDefs.h"

#include <g4main/PHG4HitDefs.h>
#include <phool/PHObject.h>
#include <cmath>
#include <map>
#include <set>

class PHG4ScintillatorSlat : public PHObject
{
 public:
  ~PHG4ScintillatorSlat() override {}

  // from PHObject
  void identify(std::ostream& os = std::cout) const override
  {
    os << "PHG4ScintillatorSlat base class" << std::endl;
  }

  virtual void add_edep(const double /*edep*/, const double /*e*/, const double /*light_yield*/) { return; }

  virtual void set_key(const PHG4ScintillatorSlatDefs::keytype) { return; }
  virtual void add_hit_key(PHG4HitDefs::keytype) { return; }

  virtual short get_column() const { return -1; }
  virtual short get_row() const { return -1; }
  virtual PHG4ScintillatorSlatDefs::keytype get_key() const { return 0xFFFFFFFF; }

  virtual double get_edep() const { return NAN; }
  virtual double get_eion() const { return NAN; }
  virtual double get_light_yield() const { return NAN; }
  virtual std::pair<std::set<PHG4HitDefs::keytype>::const_iterator, std::set<PHG4HitDefs::keytype>::const_iterator> get_hit_ids() const = 0;

 protected:
  PHG4ScintillatorSlat() {}
  ClassDefOverride(PHG4ScintillatorSlat, 1)
};

#endif
