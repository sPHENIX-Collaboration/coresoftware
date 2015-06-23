#ifndef __PHG4Hitv1_H__
#define __PHG4Hitv1_H__

#include "PHG4Hit.h"
//#include <Geant4/G4Allocator.hh>
#include <map>
#include <stdint.h>

//#ifndef __CINT__
//class PHG4Hitv1;
//extern G4Allocator<PHG4Hitv1> PHG4Hitv1Allocator;
//#endif

class PHG4Hitv1 : public PHG4Hit
{
 public:
  PHG4Hitv1();
  PHG4Hitv1(PHG4Hit const &g4hit);
  // The indices here represent the entry and exit points of the particle
  float get_x(const int i) const {return x[i];}
  float get_y(const int i) const {return y[i];}
  float get_z(const int i) const {return z[i];}
  float get_t(const int i) const {return t[i];}
  float get_edep() const {return edep;}
  unsigned int get_hit_id() const {return hitid;}
  int get_trkid() const {return trackid;}
  
  void set_x(const int i, const float f) {x[i]=f;}
  void set_y(const int i, const float f) {y[i]=f;}
  void set_z(const int i, const float f) {z[i]=f;}
  void set_t(const int i, const float f) {t[i]=f;}
  void set_edep(const float f) {edep = f;}
  void set_hit_id(const unsigned int i) {hitid=i;}
  void set_trkid(const int i) {trackid=i;}

  virtual void print() const;

  bool  has_property(PROPERTY prop_id) const;
  float get_property_float(PROPERTY prop_id) const;
  int   get_property_int(PROPERTY prop_id) const;
  unsigned int   get_property_uint(PROPERTY prop_id) const;
  void  set_property(PROPERTY prop_id, float value);
  void  set_property(PROPERTY prop_id, int value);
  void  set_property(PROPERTY prop_id, unsigned int value);

  virtual float get_px(const int i) const {return  get_property_float(prop_px);}
  virtual float get_py(const int i) const {return  get_property_float(prop_py);}
  virtual float get_pz(const int i) const {return  get_property_float(prop_pz);}
  virtual float get_eion() const          {return  get_property_float(prop_eion);}
  virtual float get_light_yield() const   {return  get_property_float(prop_light_yield);}
  virtual float    get_path_length() const{return  get_property_float(prop_path_length);}
  virtual unsigned int get_layer() const  {return  get_property_uint(prop_layer);}
  virtual int get_scint_id() const        {return  get_property_int(prop_scint_id);}
  virtual int get_strip_z_index() const   {return  get_property_int(prop_strip_z_index);}
  virtual int get_strip_y_index() const   {return  get_property_int(prop_strip_y_index);}
  virtual int get_ladder_z_index() const  {return  get_property_int(prop_ladder_phi_index);}
  virtual int get_ladder_phi_index() const{return  get_property_int(prop_ladder_phi_index);}
  virtual int get_index_i() const{return  get_property_int(prop_index_i);}
  virtual int get_index_j() const{return  get_property_int(prop_index_j);}
  virtual int get_index_k() const{return  get_property_int(prop_index_k);}
  virtual int get_index_l() const{return  get_property_int(prop_index_l);}

  virtual void set_px(const int i, const float f) {set_property(prop_px,f);}
  virtual void set_py(const int i, const float f) {set_property(prop_py,f);}
  virtual void set_pz(const int i, const float f) {set_property(prop_pz,f);}
  virtual void set_eion(const float f)            {set_property(prop_eion,f);}
  virtual void set_light_yield(float f)           {set_property(prop_light_yield,f);}
  virtual void set_path_length(float f)           {set_property(prop_path_length,f);}
  virtual void set_layer(const unsigned int i)    {set_property(prop_layer,i);}
  virtual void set_scint_id(const int i)          {set_property(prop_scint_id,i);}
  virtual void set_strip_z_index(const int i)     {set_property(prop_strip_z_index,i);}
  virtual void set_strip_y_index(const int i)     {set_property(prop_strip_y_index,i);}
  virtual void set_ladder_z_index(const int i)    {set_property(prop_ladder_phi_index,i);}
  virtual void set_ladder_phi_index(const int i)  {set_property(prop_ladder_phi_index,i);}
  virtual void set_index_i(const int i)  {set_property(prop_index_i,i);}
  virtual void set_index_j(const int i)  {set_property(prop_index_j,i);}
  virtual void set_index_k(const int i)  {set_property(prop_index_k,i);}
  virtual void set_index_l(const int i)  {set_property(prop_index_l,i);}

 protected:
  // Store both the entry and exit points of the particle
  // Remember, particles do not always enter on the inner edge!
  float x[2];
  float y[2];
  float z[2];
  float t[2];
  unsigned int hitid;
  int trackid;
  float edep;

  //! storage types for additional property
  union u_property{
    float fdata;
    int32_t idata;
    uint32_t uidata;
  };
  typedef uint8_t prop_t;
  typedef std::map<prop_t, u_property> prop_map_t;

  //! container for additional property
  prop_map_t prop_map;


  ClassDef(PHG4Hitv1,2)
};

#endif
