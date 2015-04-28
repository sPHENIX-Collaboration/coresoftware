#ifndef PHG4Hit_H__
#define PHG4Hit_H__

#include <phool/PHObject.h>
#include <cmath>

class PHG4Hit: public PHObject
{
 public:
  PHG4Hit() {}
  virtual ~PHG4Hit() {}

  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Copy(PHG4Hit const &g4hit);
  friend ostream &operator<<(ostream & stream, const PHG4Hit * hit);

  // The indices here represent the entry and exit points of the particle
  virtual float get_x(const int i) const {return NAN;}
  virtual float get_y(const int i) const {return NAN;}
  virtual float get_z(const int i) const {return NAN;}
  virtual float get_px(const int i) const {return NAN;}
  virtual float get_py(const int i) const {return NAN;}
  virtual float get_pz(const int i) const {return NAN;}
  virtual float get_t(const int i) const {return NAN;}
  virtual float get_edep() const {return NAN;}
  virtual float get_eion() const {return NAN;}
  virtual float get_light_yield() const {return NAN;}
  virtual float    get_path_length() const {return NAN;}
  virtual unsigned int get_layer() const {return 0xffffffff;}
  virtual unsigned int get_hit_id() const {return 0xffffffff;}
  virtual int get_scint_id() const {return -9999;}
  virtual int get_trkid() const {return -9999;}
  virtual int get_strip_z_index() const {return -9999;}
  virtual int get_strip_y_index() const {return -9999;}
  virtual int get_ladder_z_index() const {return -9999;}
  virtual int get_ladder_phi_index() const {return -9999;}

  virtual void set_x(const int i, const float f) {return;}
  virtual void set_y(const int i, const float f) {return;}
  virtual void set_z(const int i, const float f) {return;}
  virtual void set_px(const int i, const float f) {return;}
  virtual void set_py(const int i, const float f) {return;}
  virtual void set_pz(const int i, const float f) {return;}
  virtual void set_t(const int i, const float f) {return;}
  virtual void set_edep(const float f) {return;}
  virtual void set_eion(const float f) {return;}
  virtual void    set_light_yield(float lightYield){return;}
  virtual void  set_path_length(float pathLength){return;}
  virtual void set_layer(const unsigned int i) {return;}
  virtual void set_hit_id(const unsigned int i) {return;}
  virtual void set_scint_id(const int i) {return;}
  virtual void set_trkid(const int i) {return;}
  virtual void set_strip_z_index(const int i) {return;}
  virtual void set_strip_y_index(const int i) {return;}
  virtual void set_ladder_z_index(const int i) {return;}
  virtual void set_ladder_phi_index(const int i) {return;}

  virtual float get_avg_x() const;
  virtual float get_avg_y() const;
  virtual float get_avg_z() const;
  virtual float get_avg_t() const;

  virtual void print() const {std::cout<<"PHG4Hit base class - print() not implemented"<<std::endl;}

 protected:
  ClassDef(PHG4Hit,1)
};

inline float PHG4Hit::get_avg_x() const { return 0.5*(get_x(0)+get_x(1)); }
inline float PHG4Hit::get_avg_y() const { return 0.5*(get_y(0)+get_y(1)); }
inline float PHG4Hit::get_avg_z() const { return 0.5*(get_z(0)+get_z(1)); }
inline float PHG4Hit::get_avg_t() const { return 0.5*(get_t(0)+get_t(1)); }

#endif
