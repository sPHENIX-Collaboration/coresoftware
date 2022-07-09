// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERCELLV1_H
#define G4DETECTORS_PHG4CYLINDERCELLV1_H

#include "PHG4CylinderCell.h"

#include "PHG4Cell.h"              // for PHG4Cell::EdepMap, PHG4Cell::Showe...
#include "PHG4CylinderCellDefs.h"  // for keytype

#include <g4main/PHG4HitDefs.h>

#include <iostream>  // for cout, ostream
#include <map>
#include <utility>  // for make_pair

class PHG4CylinderCellv1 : public PHG4CylinderCell
{
 public:
  PHG4CylinderCellv1();
  ~PHG4CylinderCellv1() override {}

  void identify(std::ostream& os = std::cout) const override;

  EdepConstRange get_g4hits() override
  {
    return std::make_pair(edeps.begin(), edeps.end());
  }

  using PHG4Cell::add_edep;  // avoid warning for not including all overrides of add_edep
  void add_edep(const PHG4HitDefs::keytype g4hitid, const float edep) override;
  void add_edep(const PHG4HitDefs::keytype g4hitid, const float edep, const float light_yield) override;

  ShowerEdepConstRange get_g4showers() override
  {
    return std::make_pair(showeredeps.begin(), showeredeps.end());
  }
  void add_shower_edep(const int g4showerid, const float edep) override;

  void set_cell_id(const PHG4CylinderCellDefs::keytype id) override { cellid = id; }
  void set_layer(const unsigned int i) override { layer = i; }
  double get_edep() const override;
  unsigned int get_layer() const override { return layer; }
  PHG4CylinderCellDefs::keytype get_cell_id() const override { return cellid; }
  int get_binz() const override { return binz; }
  int get_binphi() const override { return binphi; }
  int get_bineta() const override { return get_binz(); }
  float get_light_yield() const override { return light_yield; }

  void set_zbin(const int i) override { binz = i; }
  void set_etabin(const int i) override { set_zbin(i); }
  void set_phibin(const int i) override { binphi = i; }
  void set_light_yield(const float lightYield) override { light_yield = lightYield; }

 protected:
  unsigned int layer;
  PHG4CylinderCellDefs::keytype cellid;
  int binz;
  int binphi;
  EdepMap edeps;
  ShowerEdepMap showeredeps;
  float light_yield;

  ClassDefOverride(PHG4CylinderCellv1, 2)
};

#endif
