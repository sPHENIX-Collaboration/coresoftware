#ifndef EVENTFLOWINFOMAP_H
#define EVENTFLOWINFOMAP_H

#include <phool/PHObject.h>

#include <iostream>
#include <map>

#include "EventFlowInfo.h"


class EventFlowInfoMap : public PHObject
{

 public:

  typedef std::map<unsigned int, EventFlowInfo*>::const_iterator ConstIter;
  typedef std::map<unsigned int, EventFlowInfo*>::iterator Iter;

  ~EventFlowInfoMap() override {}

  void identify(std::ostream& os = std::cout) const override { os << "EventFlowInfoMap base class" << std::endl; }
  virtual bool empty() const {return true;}
  virtual void clear() {}

  virtual const EventFlowInfo* get(unsigned int /*idkey*/) const { return nullptr; }
  virtual EventFlowInfo* get(unsigned int /*idkey*/) { return nullptr; }
  virtual const EventFlowInfo* get(EventFlowInfo::EventFlowSrc /*type*/) const { return nullptr; }
  virtual EventFlowInfo* get(EventFlowInfo::EventFlowSrc /*type*/) { return nullptr; }

  virtual EventFlowInfo* insert(EventFlowInfo* /*ep*/, EventFlowInfo::EventFlowSrc /*type*/) { return nullptr; }

  virtual ConstIter begin() const { return DummyEventFlowInfoMap.end(); }
  virtual ConstIter find(unsigned int /*idkey*/) const { return DummyEventFlowInfoMap.end(); }
  virtual ConstIter find(EventFlowInfo::EventFlowSrc /*type*/) const { return DummyEventFlowInfoMap.end(); }
  virtual ConstIter end() const { return DummyEventFlowInfoMap.end(); }

  virtual Iter begin() { return DummyEventFlowInfoMap.end(); }
  virtual Iter find(unsigned int /*idkey*/) { return DummyEventFlowInfoMap.end(); }
  virtual Iter find(EventFlowInfo::EventFlowSrc /*type*/) { return DummyEventFlowInfoMap.end(); }
  virtual Iter end() { return DummyEventFlowInfoMap.end(); }

 protected:
  EventFlowInfoMap() {}
  std::map<unsigned int, EventFlowInfo*> DummyEventFlowInfoMap {};

 private:
  ClassDefOverride(EventFlowInfoMap, 1);
};

#endif  // EVENTFLOWINFOMAP_H
