#ifndef __PHG4PARTICLE_H__
#define __PHG4PARTICLE_H__

#include <phool/PHObject.h>
#include <cmath>

class PHG4Particle: public PHObject
{
 public:

  PHG4Particle(){}
  virtual ~PHG4Particle() {}

  virtual int get_pid() const {return 0;}
  virtual std::string get_name() const {return "NONE";}
  virtual double get_px() const {return NAN;}
  virtual double get_py() const {return NAN;}
  virtual double get_pz() const {return NAN;}
  virtual double get_e() const {return NAN;}

  virtual int get_track_id() const {return -9999;}
  virtual int get_vtx_id() const {return -9999;}
  virtual int get_parent_id() const {return -9999;}
  virtual int get_primary_id() const {return 0xFFFFFFFF;}

  virtual int get_barcode() const {return 0xFFFFFFFF;}

  virtual void set_track_id(const int i) {return;}
  virtual void set_vtx_id(const int i) {return;}
  virtual void set_parent_id(const int i) {return;}
  virtual void set_primary_id(const int i) {return;}
  virtual void set_name(const std::string &name) {return;}
  virtual void set_pid(const int i) {return;}
  virtual void set_px(const double x) {return;}
  virtual void set_py(const double x) {return;}
  virtual void set_pz(const double x) {return;}
  virtual void set_e(const double e) {return;}

  virtual void set_barcode(const int barcode) {return;}

  void identify(std::ostream& os = std::cout) const;

  bool  operator== (const PHG4Particle &p) const;

 protected:
  ClassDef(PHG4Particle,1)
};


#endif
