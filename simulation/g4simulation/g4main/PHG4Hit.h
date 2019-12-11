// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4HIT_H
#define G4MAIN_PHG4HIT_H

#include "PHG4HitDefs.h"

#include <phool/PHObject.h>

#include <cmath>
#include <climits>
#include <iostream> 
#include <string>
#include <utility>

class PHG4Hit: public PHObject
{
 public:
  PHG4Hit() {}
  virtual ~PHG4Hit() {}

  virtual void identify(std::ostream& os = std::cout) const;
  virtual void CopyFrom(const PHObject *phobj);
  friend std::ostream &operator<<(std::ostream & stream, const PHG4Hit * hit);
  virtual void Reset();

  // The indices here represent the entry and exit points of the particle
  virtual float get_x(const int i) const {return NAN;}
  virtual float get_y(const int i) const {return NAN;}
  virtual float get_z(const int i) const {return NAN;}
  virtual float get_px(const int i) const {return NAN;}
  virtual float get_py(const int i) const {return NAN;}
  virtual float get_pz(const int i) const {return NAN;}
  virtual float get_local_x(const int i) const {return NAN;}
  virtual float get_local_y(const int i) const {return NAN;}
  virtual float get_local_z(const int i) const {return NAN;}
  virtual float get_t(const int i) const {return NAN;}
  virtual float get_edep() const {return NAN;}
  virtual float get_eion() const {return NAN;}
  virtual float get_light_yield() const {return NAN;}
  virtual float get_path_length() const {return NAN;}
  virtual unsigned int get_layer() const {return UINT_MAX;}
  virtual PHG4HitDefs::keytype get_hit_id() const {return ULONG_LONG_MAX;}
  virtual int get_detid() const {return INT_MIN;}
  virtual int get_shower_id() const {return INT_MIN;}
  virtual int get_scint_id() const {return INT_MIN;}
  virtual int get_row() const {return INT_MIN;}
  virtual int get_trkid() const {return INT_MIN;}
  virtual int get_strip_z_index() const {return INT_MIN;}
  virtual int get_strip_y_index() const {return INT_MIN;}
  virtual int get_ladder_z_index() const {return INT_MIN;}
  virtual int get_ladder_phi_index() const {return INT_MIN;}
  virtual int get_index_i() const {return INT_MIN;}
  virtual int get_index_j() const {return INT_MIN;}
  virtual int get_index_k() const {return INT_MIN;}
  virtual int get_index_l() const {return INT_MIN;}
  virtual int get_hit_type() const {return INT_MIN;}

  virtual void set_x(const int i, const float f) {return;}
  virtual void set_y(const int i, const float f) {return;}
  virtual void set_z(const int i, const float f) {return;}
  virtual void set_px(const int i, const float f) {return;}
  virtual void set_py(const int i, const float f) {return;}
  virtual void set_pz(const int i, const float f) {return;}
  virtual void set_local_x(const int i, const float f) {return;}
  virtual void set_local_y(const int i, const float f) {return;}
  virtual void set_local_z(const int i, const float f) {return;}
  virtual void set_t(const int i, const float f) {return;}
  virtual void set_edep(const float f) {return;}
  virtual void set_eion(const float f) {return;}
  virtual void set_light_yield(const float lightYield){return;}
  virtual void set_path_length(const float pathLength){return;}
  virtual void set_layer(const unsigned int i) {return;}
  virtual void set_hit_id(const PHG4HitDefs::keytype i) {return;}
  virtual void set_shower_id(const int i) {return;}
  virtual void set_scint_id(const int i) {return;}
  virtual void set_row(const int i) {return;}
  virtual void set_trkid(const int i) {return;}
  virtual void set_strip_z_index(const int i) {return;}
  virtual void set_strip_y_index(const int i) {return;}
  virtual void set_ladder_z_index(const int i) {return;}
  virtual void set_ladder_phi_index(const int i) {return;}
  virtual void set_index_i(const int i) {return;}
  virtual void set_index_j(const int i) {return;}
  virtual void set_index_k(const int i) {return;}
  virtual void set_index_l(const int i) {return;}
  virtual void set_hit_type(const int i) {return;}

  virtual float get_avg_x() const;
  virtual float get_avg_y() const;
  virtual float get_avg_z() const;
  virtual float get_avg_t() const;

  virtual void print() const {std::cout<<"PHG4Hit base class - print() not implemented"<<std::endl;}


  //! Procedure to add a new PROPERTY tag:
  //! 1.add new tag below with unique value,
  //! 2.add a short name to PHG4Hit::get_property_info
  enum PROPERTY 
  {//

    //-- hit properties: 1 - 10  --
    //! ionizing energy loss
    prop_eion = 1,

    //! for scintillation detectors, the amount of light produced
    prop_light_yield = 2,

    //-- track properties: 10 - 20  --

    //! momentum
    prop_px_0 = 10,
    prop_px_1 = 11,
    prop_py_0 = 12,
    prop_py_1 = 13,
    prop_pz_0 = 14,
    prop_pz_1 = 15,

    //! pathlength
    prop_path_length = 16,

    //! local coordinate
    prop_local_x_0 = 20,
    prop_local_x_1 = 21,
    prop_local_y_0 = 22,
    prop_local_y_1 = 23,
    prop_local_z_0 = 24,
    prop_local_z_1 = 25,

    //-- detector specific IDs: 100+ --

    //! layer ID
    prop_layer = 101,
    //! scintillator ID
    prop_scint_id = 102,
    //! row (mother volume or steel plate id)
    prop_row = 103,

    //! SVX stuff
    prop_strip_z_index = 110,
    prop_strip_y_index = 111,
    prop_ladder_z_index = 112,
    prop_ladder_phi_index = 113,

    // MAPS stuff 
    prop_stave_index = 114,
    prop_half_stave_index = 115,
    prop_module_index = 116,
    prop_chip_index = 117,

    prop_local_pos_x_0 = 118,
    prop_local_pos_y_0 = 119,
    prop_local_pos_z_0 = 120,

    //! generic indexes
    prop_index_i = 121,
    prop_index_j = 122,
    prop_index_k = 123,
    prop_index_l = 124,

    //! hit type
    prop_hit_type = 125,

    prop_local_pos_x_1 = 128,
    prop_local_pos_y_1 = 126,
    prop_local_pos_z_1 = 127,

    //! max limit in order to fit into 8 bit unsigned number
    prop_MAX_NUMBER = UCHAR_MAX
  };

  enum PROPERTY_TYPE 
  {//
    type_int = 1,
    type_uint = 2,
    type_float = 3,
    type_unknown = -1
  };

  virtual bool  has_property(const PROPERTY prop_id) const {return false;}
  virtual float get_property_float(const PROPERTY prop_id) const {return NAN;}
  virtual int   get_property_int(const PROPERTY prop_id) const {return INT_MIN;}
  virtual unsigned int   get_property_uint(const PROPERTY prop_id) const {return UINT_MAX;}
  virtual void  set_property(const PROPERTY prop_id, const float value) {return;}
  virtual void  set_property(const PROPERTY prop_id, const int value) {return;}
  virtual void  set_property(const PROPERTY prop_id, const unsigned int value) {return;}
  static std::pair<const std::string,PROPERTY_TYPE> get_property_info(PROPERTY prop_id);
  static bool check_property(const PROPERTY prop_id, const PROPERTY_TYPE prop_type);
  static std::string get_property_type(const PROPERTY_TYPE prop_type);

 protected:
  virtual unsigned int get_property_nocheck(const PROPERTY prop_id) const {return UINT_MAX;}
  virtual void set_property_nocheck(const PROPERTY prop_id,const unsigned int) {return;}
  ClassDef(PHG4Hit,1)
};


inline float PHG4Hit::get_avg_x() const { return 0.5*(get_x(0)+get_x(1)); }
inline float PHG4Hit::get_avg_y() const { return 0.5*(get_y(0)+get_y(1)); }
inline float PHG4Hit::get_avg_z() const { return 0.5*(get_z(0)+get_z(1)); }
inline float PHG4Hit::get_avg_t() const { return 0.5*(get_t(0)+get_t(1)); }


#endif
