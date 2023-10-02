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
  typedef std::map<unsigned int, Eventplaneinfo*>::const_iterator ConstIter;
  typedef std::map<unsigned int, Eventplaneinfo*>::iterator Iter;

  ~EventplaneinfoMap() override {}

  void identify(std::ostream& os = std::cout) const override { os << "EventplaneinfoMap base class" << std::endl; }

  virtual bool empty() const { return true; }
  virtual size_t size() const { return 0; }
  virtual size_t count(unsigned int /*idkey*/) const { return 0; }
  virtual void clear() {}

  virtual const Eventplaneinfo* get(unsigned int /*idkey*/) const { return nullptr; }
  virtual Eventplaneinfo* get(unsigned int /*idkey*/) { return nullptr; }
  virtual Eventplaneinfo* insert(Eventplaneinfo* /*ep*/) { return nullptr; }
  virtual size_t erase(unsigned int /*idkey*/) { return 0; }

  virtual ConstIter begin() const;
  virtual ConstIter find(unsigned int idkey) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(unsigned int idkey);
  virtual Iter end();

 protected:
    EventplaneinfoMap() {}

 private:
  ClassDefOverride(EventplaneinfoMap, 1);
};

#endif  // EVENTPLANEINFOMAP_H
