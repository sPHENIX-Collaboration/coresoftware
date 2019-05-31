#ifndef G4DETECTORS_PHG4SCINTILLATORSLATV1_H
#define G4DETECTORS_PHG4SCINTILLATORSLATV1_H

#include "PHG4ScintillatorSlat.h"

#include "PHG4ScintillatorSlatDefs.h"  // for keytype

#include <g4main/PHG4HitDefs.h>

#include <iostream>                    // for cout, ostream
#include <set>
#include <utility>                     // for make_pair, pair

class PHG4ScintillatorSlatv1 : public PHG4ScintillatorSlat
{
 public:

  PHG4ScintillatorSlatv1();
  virtual ~PHG4ScintillatorSlatv1(){}

  void identify(std::ostream& os = std::cout) const;

  void add_edep(const double f, const double e, const double ly) {edep+=f; eion+= e; light_yield+=ly;}
  void add_hit_key(PHG4HitDefs::keytype key) {hit_id.insert(key);}
  
  void set_key(PHG4ScintillatorSlatDefs::keytype i) {key = i;}
  void set_light_yield(const double lightYield)  {light_yield = lightYield;}

  short get_row() const;
  short get_column() const;
  PHG4ScintillatorSlatDefs::keytype get_key() const {return key;}
  double get_edep() const {return edep;}
  double get_eion() const {return eion;}
  double get_light_yield() const {return light_yield;}
  std::pair<std::set<PHG4HitDefs::keytype>::const_iterator, std::set<PHG4HitDefs::keytype>::const_iterator> get_hit_ids() const {return std::make_pair(hit_id.begin(),hit_id.end());}


 protected:
  PHG4ScintillatorSlatDefs::keytype key;
  double edep;
  double eion;
  double light_yield;
  std::set<PHG4HitDefs::keytype> hit_id;

   
  ClassDef(PHG4ScintillatorSlatv1,1)
};

#endif
