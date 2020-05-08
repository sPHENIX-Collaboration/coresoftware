// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EICEVENTHEADERV1_H
#define EICEVENTHEADERV1_H

#include "EicEventHeader.h"

#include <climits>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <utility>

class EicEventHeaderv1: public EicEventHeader
{
 public:
  EicEventHeaderv1() {}
  explicit EicEventHeaderv1(const EicEventHeader *eicevt);
  virtual ~EicEventHeaderv1(){}

//  void identify(std::ostream& os  = std::cout) const;
  void Reset();

// property relates methods
  bool  has_property(const PROPERTY prop_id) const;
  float get_property_float(const PROPERTY prop_id) const;
  int   get_property_int(const PROPERTY prop_id) const;
  unsigned int   get_property_uint(const PROPERTY prop_id) const;
  void  set_property(const PROPERTY prop_id, const float value);
  void  set_property(const PROPERTY prop_id, const int value);
  void  set_property(const PROPERTY prop_id, const unsigned int value);


// Generator specific values
  void set_eventgenerator_type(const int i) {set_property( prop_eventgen, i);}
  int get_eventgenerator_type() const {return get_property_int(prop_eventgen);}

// Milou
  void set_milou_weight(const float val) {set_property( prop_milou_weight,val);}
  float get_milou_weight() const {return get_property_float(prop_milou_weight);}
  void set_milou_trueX(const float val) {set_property(prop_milou_truex,val);}
  float get_milou_trueX() const {return get_property_float(prop_milou_truex);}
  void set_milou_trueQ2(const float val) {set_property(prop_milou_trueq2,val);}
  float get_milou_trueQ2() const {return get_property_float(prop_milou_trueq2);}
 

 protected:
  unsigned int get_property_nocheck(const PROPERTY prop_id) const;
  void set_property_nocheck(const PROPERTY prop_id,const unsigned int ui) {prop_map[prop_id]=ui;}

  //! storage types for properties
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

  //! container for  properties
  prop_map_t prop_map;

  ClassDef(EicEventHeaderv1,1)
};

#endif
