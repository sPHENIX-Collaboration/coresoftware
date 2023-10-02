// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTPLANEINFOV1_H
#define EVENTPLANEINFOV1_H

#include "Eventplaneinfo.h"

#include <cstddef>  // for size_t
#include <iostream>
#include <limits>
#include <map>
#include <utility>  // for pair, make_pair

class PHObject;

class Eventplaneinfov1 : public Eventplaneinfo
{
 public:
    Eventplaneinfov1(const Eventplaneinfo::EPTYPE id = Eventplaneinfo::UNDEFINED);
  ~Eventplaneinfov1() override = default;

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = Eventplaneinfov1(); }
  PHObject* CloneMe() const override { return new Eventplaneinfov1(*this); }
    
    
  unsigned int get_id() const override { return _id; }
 void set_id(unsigned int id) override { _id = id; }

  double get_psi_raw(int order) const override {return _psi[order - 1];}
  void set_psi_raw(unsigned int order, double f) override {_psi[order] = f;}
    
  
  bool empty_ep_ids() const override { return _ep_ids.empty(); }
  size_t size_ep_ids() const override { return _ep_ids.size(); }
  size_t count_ep_ids(Eventplaneinfo::EPTYPE type) const override { return _ep_ids.count(type); }

  void clear_ep_ids() override { _ep_ids.clear(); }
  void insert_ep_ids(Eventplaneinfo::EPTYPE type, unsigned int ep_id) override { _ep_ids.insert(std::make_pair(type, ep_id)); }
  size_t erase_ep_ids(Eventplaneinfo::EPTYPE type) override { return _ep_ids.erase(type); }
  void erase_ep_ids(Eventplaneinfo::EpIter iter) override { _ep_ids.erase(iter); }
  void erase_ep_ids(Eventplaneinfo::EpIter first, Eventplaneinfo::EpIter last) override { _ep_ids.erase(first, last); }

  Eventplaneinfo::ConstEpIter begin_ep_ids() const override { return _ep_ids.begin(); }
  Eventplaneinfo::ConstEpIter find_ep_ids(Eventplaneinfo::EPTYPE type) const override { return _ep_ids.find(type); }
  Eventplaneinfo::ConstEpIter end_ep_ids() const override { return _ep_ids.end(); }

  Eventplaneinfo::EpIter begin_ep_ids() override { return _ep_ids.begin(); }
  Eventplaneinfo::EpIter find_ep_ids(Eventplaneinfo::EPTYPE type) override { return _ep_ids.find(type); }
  Eventplaneinfo::EpIter end_ep_ids() override { return _ep_ids.end(); }

 private:

  unsigned int _id = std::numeric_limits<unsigned int>::max();
  double _psi[24] = {}; 
  std::map<Eventplaneinfo::EPTYPE, unsigned int> _ep_ids;

  ClassDefOverride(Eventplaneinfov1, 1);
};

#endif
