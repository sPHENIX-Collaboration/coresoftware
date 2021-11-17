// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PARTICLE_H
#define G4MAIN_PHG4PARTICLE_H

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>
#include <string>

class PHG4Particle : public PHObject
{
 public:
  PHG4Particle() {}
  ~PHG4Particle() override {}

  void identify(std::ostream &os = std::cout) const override;


  virtual bool isIon() const { return false; }
  virtual int get_pid() const { return 0; }
  virtual std::string get_name() const { return "NONE"; }
  virtual double get_px() const { return NAN; }
  virtual double get_py() const { return NAN; }
  virtual double get_pz() const { return NAN; }
  virtual double get_e() const { return NAN; }

  virtual int get_track_id() const { return -9999; }
  virtual int get_vtx_id() const { return -9999; }
  virtual int get_parent_id() const { return -9999; }
  virtual int get_primary_id() const { return 0xFFFFFFFF; }

  virtual int get_barcode() const { return 0xFFFFFFFF; }

  virtual int get_A() const { return 0xFFFFFFFF; }
  virtual int get_Z() const { return 0xFFFFFFFF; }
  virtual double get_IonCharge() const { return NAN; }
  virtual double get_ExcitEnergy() const { return NAN; }

  virtual void set_track_id(const int) { return; }
  virtual void set_vtx_id(const int) { return; }
  virtual void set_parent_id(const int) { return; }
  virtual void set_primary_id(const int) { return; }
  virtual void set_name(const std::string &) { return; }
  virtual void set_pid(const int) { return; }
  virtual void set_px(const double) { return; }
  virtual void set_py(const double) { return; }
  virtual void set_pz(const double) { return; }
  virtual void set_e(const double) { return; }

  virtual void set_barcode(const int) { return; }

  virtual void set_A(const int) { return; }
  virtual void set_Z(const int) { return; }
  virtual void set_NumCharge(const int) { return; }
  virtual void set_IonCharge(const double) { return; }
  virtual void set_ExcitEnergy(const double) { return; }

  bool operator==(const PHG4Particle &p) const;

 protected:
  ClassDefOverride(PHG4Particle, 1)
};

#endif  // G4MAIN_PHG4PARTICLE_H
