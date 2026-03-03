// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTPLANEINFOMAP_H
#define EVENTPLANEINFOMAP_H

#include <phool/PHObject.h>
#include <iostream>
#include <map>

class Eventplaneinfo;

class EventplaneinfoMap : public PHObject
{
 public:
  enum EPTYPE
  {
    UNDEFINED = 999,
    sEPDS = 0,
    sEPDN = 1,
    MBDS = 2,
    MBDN = 3,
    sEPDNS = 4,
    MBDNS = 5,
    sEPDRING_SOUTH = 100,
    sEPDRING_NORTH = 200  
  };

  using ConstIter = std::map<unsigned int, Eventplaneinfo*>::const_iterator;
  using Iter = std::map<unsigned int, Eventplaneinfo*>::iterator;

  ~EventplaneinfoMap() override = default;

  void identify(std::ostream& os = std::cout) const override { os << "EventplaneinfoMap base class" << std::endl; }
  virtual bool empty() const {return true;}
  virtual void clear() {}

  virtual const Eventplaneinfo* get(unsigned int /*idkey*/) const { return nullptr; }
  virtual Eventplaneinfo* get(unsigned int /*idkey*/) { return nullptr; }
  virtual Eventplaneinfo* insert(Eventplaneinfo* /*ep*/, EPTYPE /*type*/) { return nullptr; }

  virtual ConstIter begin() const;
  virtual ConstIter find(unsigned int idkey) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(unsigned int idkey);
  virtual Iter end();

 protected:
  EventplaneinfoMap() = default;

  // Rule of Five: Protected to support derived classes
  EventplaneinfoMap(const EventplaneinfoMap&) = default;
  EventplaneinfoMap& operator=(const EventplaneinfoMap&) = default;
  EventplaneinfoMap(EventplaneinfoMap&&) = default;
  EventplaneinfoMap& operator=(EventplaneinfoMap&&) = default;

 private:
  ClassDefOverride(EventplaneinfoMap, 1);
};

#endif  // EVENTPLANEINFOMAP_H
