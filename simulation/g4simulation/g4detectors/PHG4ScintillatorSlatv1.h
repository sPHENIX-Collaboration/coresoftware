#ifndef PHG4SCINTILLATORSLATV1_H
#define PHG4SCINTILLATORSLATV1_H

#include "PHG4ScintillatorSlat.h"
#include <g4main/PHG4HitDefs.h>

#include <cmath>
#include <set>
#include <map>

class PHG4ScintillatorSlatv1 : public PHG4ScintillatorSlat
{
 public:

  PHG4ScintillatorSlatv1();
  virtual ~PHG4ScintillatorSlatv1(){}

  void identify(std::ostream& os = std::cout) const;

  void add_edep(const double f, const double e, const double ly) {edep+=f; eion+= e; light_yield+=ly;}
  void add_hit_key(PHG4HitDefs::keytype key) {hit_id.insert(key);}
  
  void set_column(const short id) {column = id;}
  void set_key(const unsigned int i) {key = i;}
  void set_row(const short i) {row = i;}
  void set_light_yield(const double lightYield)  {light_yield = lightYield;}

  short get_row() const {return row;}
  short get_column() const {return column;}
  unsigned int get_key() const;
  double get_edep() const {return edep;}
  double get_eion() const {return eion;}
  double get_light_yield() const {return light_yield;}



 protected:
  unsigned int key;
  short row;
  short column;
  double edep;
  double eion;
  double light_yield;
  std::set<PHG4HitDefs::keytype> hit_id;

   
  ClassDef(PHG4ScintillatorSlatv1,1)
};

#endif
