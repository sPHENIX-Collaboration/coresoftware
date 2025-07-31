// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTPLANEINFOMAPV1_H
#define EVENTPLANEINFOMAPV1_H

#include "EventplaneinfoMap.h"

#include "Eventplaneinfo.h"

#include <iostream>
#include <map>

class EventplaneinfoMapv1 : public EventplaneinfoMap
{
 public:
  EventplaneinfoMapv1() = default;
  ~EventplaneinfoMapv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { clear(); }
// cppcheck-suppress [virtualCallInConstructor]
  bool empty() const override {return _map.empty();}
  void clear() override;

  const Eventplaneinfo* get(unsigned int idkey) const override;
  Eventplaneinfo* get(unsigned int idkey) override;

  Eventplaneinfo* insert(Eventplaneinfo* clus, EPTYPE id) override;

  ConstIter begin() const override { return _map.begin(); }
  ConstIter find(unsigned int idkey) const override { return _map.find(idkey); }
  ConstIter end() const override { return _map.end(); }

  Iter begin() override { return _map.begin(); }
  Iter find(unsigned int idkey) override { return _map.find(idkey); }
  Iter end() override { return _map.end(); }

 private:
  std::map<unsigned int, Eventplaneinfo*> _map;

  ClassDefOverride(EventplaneinfoMapv1, 1);
};

#endif  // EVENTPLANEINFOMAPV1_H
