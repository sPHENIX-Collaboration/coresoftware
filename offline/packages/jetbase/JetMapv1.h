#ifndef G4JET_JETMAPV1_H
#define G4JET_JETMAPV1_H

#include "JetMap.h"

#include "Jet.h"

#include <cstddef>  // for size_t
#include <iostream>
#include <set>

class PHObject;

class JetMapv1 : public JetMap
{
 public:
  JetMapv1()  = default;
  explicit JetMapv1(const JetMap& jets);
  JetMapv1& operator=(const JetMap& jets);
  ~JetMapv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override { return 1; }
  PHObject* CloneMe() const override  { return new JetMapv1(*this);}

  // map content info ----------------------------------------------------------

  void set_algo(Jet::ALGO algo) override { _algo = algo; }
  Jet::ALGO get_algo() const override { return _algo; }

  void set_par(float par) override { _par = par; }
  float get_par() const override { return _par; }

  // set access to source identifiers ------------------------------------------

  bool empty_src() const override { return _src.empty(); }
  void insert_src(Jet::SRC src) override { _src.insert(src); }

  ConstSrcIter begin_src() const override { return _src.begin(); }
  ConstSrcIter find_src(Jet::SRC src) const override { return _src.find(src); }
  ConstSrcIter end_src() const override { return _src.end(); }

  SrcIter begin_src() override { return _src.begin(); }
  SrcIter find_src(Jet::SRC src) override { return _src.find(src); }
  SrcIter end_src() override { return _src.end(); }

  // map access to jets --------------------------------------------------------

  bool empty() const override { return _map.empty(); }
  size_t size() const override { return _map.size(); }
  size_t count(unsigned int idkey) const override { return _map.count(idkey); }
  void clear() override { Reset(); }

  const Jet* get(unsigned int idkey) const override;
  Jet* get(unsigned int idkey) override;

  /// insert Jet to the map. Once inserted, the JetMap owns the Jet memory
  Jet* insert(Jet* jet) override;
  size_t erase(unsigned int idkey) override { return _map.erase(idkey); }

  ConstIter begin() const override { return _map.begin(); }
  ConstIter find(unsigned int idkey) const override { return _map.find(idkey); }
  ConstIter end() const override { return _map.end(); }

  Iter begin() override { return _map.begin(); }
  Iter find(unsigned int idkey) override { return _map.find(idkey); }
  Iter end() override { return _map.end(); }

  std::vector<Jet*> vec(Jet::SORT sort=Jet::SORT::PT) override; // defaulted to PT in JetMap.h
  std::vector<Jet*> vec(std::function<bool (Jet*, Jet*)> custom_sort) override; // defaulted to PT in JetMap.h

 private:
  Jet::ALGO _algo = Jet::NONE;  //< algorithm used to reconstruct jets
  float _par = NAN;             //< algorithm parameter setting (e.g. radius)
  std::set<Jet::SRC> _src;      //< list of sources (clusters, towers, etc)
  typ_JetMap _map;              //< jet algorithm output storage

  ClassDefOverride(JetMapv1, 1);
};

#endif
