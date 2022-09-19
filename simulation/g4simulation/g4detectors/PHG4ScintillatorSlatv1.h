// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4SCINTILLATORSLATV1_H
#define G4DETECTORS_PHG4SCINTILLATORSLATV1_H

#include "PHG4ScintillatorSlat.h"

#include "PHG4ScintillatorSlatDefs.h"  // for keytype

#include <g4main/PHG4HitDefs.h>

#include <iostream>  // for cout, ostream
#include <set>
#include <utility>  // for make_pair, pair

class PHG4ScintillatorSlatv1 : public PHG4ScintillatorSlat
{
 public:
  PHG4ScintillatorSlatv1() {}
  ~PHG4ScintillatorSlatv1() override {}

  void identify(std::ostream& os = std::cout) const override;

  void add_edep(const double f, const double e, const double ly) override
  {
    edep += f;
    eion += e;
    light_yield += ly;
  }
  void add_hit_key(PHG4HitDefs::keytype i) override { hit_id.insert(i); }

  void set_key(PHG4ScintillatorSlatDefs::keytype i) override { key = i; }
  void set_light_yield(const double lightYield) { light_yield = lightYield; }

  short get_row() const override;
  short get_column() const override;
  PHG4ScintillatorSlatDefs::keytype get_key() const override { return key; }
  double get_edep() const override { return edep; }
  double get_eion() const override { return eion; }
  double get_light_yield() const override { return light_yield; }
  std::pair<std::set<PHG4HitDefs::keytype>::const_iterator, std::set<PHG4HitDefs::keytype>::const_iterator> get_hit_ids() const override { return std::make_pair(hit_id.begin(), hit_id.end()); }

 protected:
  PHG4ScintillatorSlatDefs::keytype key = ~0x0;
  double edep = 0.;
  double eion = 0.;
  double light_yield = 0.;
  std::set<PHG4HitDefs::keytype> hit_id;

  ClassDefOverride(PHG4ScintillatorSlatv1, 1)
};

#endif
