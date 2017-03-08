#ifndef PHG4Cellv1_h__
#define PHG4Cellv1_h__

#include "PHG4Cell.h"
#include "PHG4CellDefs.h"
#ifdef __CINT__
#include <stdint.h>
#else
#include <cstdint>
#endif
#include <iostream>
#include <map>

class PHG4Cellv1: public PHG4Cell
{
 public:
  PHG4Cellv1();
  PHG4Cellv1(const PHG4CellDefs::keytype g4cellid);
  virtual ~PHG4Cellv1();

  void Reset();

  void set_cellid(const PHG4CellDefs::keytype i) {cellid = i;}

  PHG4CellDefs::keytype get_cellid() const {return cellid;}
  bool has_binning(const PHG4CellDefs::CellBinning binning) const;
  short int get_detid() const;

  void add_edep(const PHG4HitDefs::keytype g4hitid, const float edep);
  void add_shower_edep(const int g4showerid, const float edep);


  void add_edep(const float f) {add_property(prop_edep,f);}
  double get_edep() const {return get_property_float(prop_edep);}

  void add_eion(const float f) {add_property(prop_eion,f);}
  double get_eion() const {return get_property_float(prop_eion);}

  void add_light_yield(const float f) {add_property(prop_light_yield,f);}
  float get_light_yield() const {return get_property_float(prop_light_yield);}

  EdepConstRange get_g4hits() {
    return std::make_pair(hitedeps.begin(), hitedeps.end());
  }
  
  ShowerEdepConstRange get_g4showers() {
    return std::make_pair(showeredeps.begin(),showeredeps.end());
  } 

  void print() const;

  bool  has_property(const PROPERTY prop_id) const;
  float get_property_float(const PROPERTY prop_id) const;
  int   get_property_int(const PROPERTY prop_id) const;
  unsigned int   get_property_uint(const PROPERTY prop_id) const;
  void  add_property(const PROPERTY prop_id, const float value);
  void  add_property(const PROPERTY prop_id, const int value);
  void  add_property(const PROPERTY prop_id, const unsigned int value);
  void  set_property(const PROPERTY prop_id, const float value);
  void  set_property(const PROPERTY prop_id, const int value);
  void  set_property(const PROPERTY prop_id, const unsigned int value);


 protected:
  unsigned int get_property_nocheck(const PROPERTY prop_id) const;
  void set_property_nocheck(const PROPERTY prop_id,const unsigned int ui) {prop_map[prop_id]=ui;}

  PHG4CellDefs::keytype cellid;
  EdepMap hitedeps;
  ShowerEdepMap showeredeps;

  //! storage types for additional property
  typedef uint8_t prop_id_t;
  typedef uint32_t prop_storage_t;
  typedef std::map<prop_id_t, prop_storage_t> prop_map_t;

  //! convert between 32bit inputs and storage type prop_storage_t
  union u_property{
    float fdata;
    int32_t idata;
    uint32_t uidata;

    u_property(int32_t in): idata(in) {}
    u_property(uint32_t in): uidata(in) {}
    u_property(float in): fdata(in) {}
    u_property(): uidata(0) {}

    prop_storage_t get_storage() const {return uidata;}
  };

  //! container for additional property
  prop_map_t prop_map;


  ClassDef(PHG4Cellv1,1)
};

#endif
