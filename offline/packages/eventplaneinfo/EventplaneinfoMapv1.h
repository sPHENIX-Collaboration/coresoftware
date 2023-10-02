// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTPLANEINFOMAPV1_H
#define EVENTPLANEINFOMAPV1_H

#include "EventplaneinfoMap.h"

#include "Eventplaneinfo.h"

#include <cstddef>  // for size_t
#include <iostream>
#include <map>

class EventplaneinfoMapv1 : public EventplaneinfoMap
{
 public:
  EventplaneinfoMapv1() = default;
  ~EventplaneinfoMapv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { clear(); }

  bool empty() const override { return _map.empty(); }
  size_t size() const override { return _map.size(); }
  size_t count(unsigned int idkey) const override { return _map.count(idkey); }
  void clear() override;

  const Eventplaneinfo* get(unsigned int idkey) const override;
  Eventplaneinfo* get(unsigned int idkey) override;
  Eventplaneinfo* insert(Eventplaneinfo* ep) override;
  size_t erase(unsigned int idkey) override
  {
    delete _map[idkey];
    return _map.erase(idkey);
  }

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
