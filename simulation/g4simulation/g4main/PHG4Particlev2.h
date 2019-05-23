#ifndef G4MAIN_PHG4PARTICLEV2_H
#define G4MAIN_PHG4PARTICLEV2_H

#include "PHG4Particlev1.h"

#include <iostream>
#include <string>

class PHG4Particle;

class PHG4Particlev2 : public PHG4Particlev1
{
 public:
  PHG4Particlev2();
  PHG4Particlev2(const std::string &name, const int pid, const double px, const double py, const double pz);
  PHG4Particlev2(const PHG4Particle *in);

  virtual ~PHG4Particlev2() {}

  int get_track_id() const { return trkid; }
  int get_vtx_id() const { return vtxid; }
  int get_parent_id() const { return parentid; }
  int get_primary_id() const { return primaryid; }
  double get_e() const { return fe; }

  void set_track_id(const int i) { trkid = i; }
  void set_vtx_id(const int i) { vtxid = i; }
  void set_parent_id(const int i) { parentid = i; }
  void set_primary_id(const int i) { primaryid = i; }
  void set_e(const double e) { fe = e; }

  void identify(std::ostream &os = std::cout) const;

 protected:
  int trkid;
  int vtxid;
  int parentid;
  int primaryid;
  double fe;

  ClassDef(PHG4Particlev2, 2)
};

#endif
