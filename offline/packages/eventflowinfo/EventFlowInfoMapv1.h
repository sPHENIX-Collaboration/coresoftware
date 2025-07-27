// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTFLOWINFOMAPV1_H
#define EVENTFLOWINFOMAPV1_H

#include "EventFlowInfoMap.h"

#include "EventFlowInfo.h"

#include <iostream>
#include <map>

class EventFlowInfoMapv1 : public EventFlowInfoMap
{
 public:
  EventFlowInfoMapv1() = default;
  ~EventFlowInfoMapv1() override {}

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { clear(); }
  bool empty() const override {return _map.empty();}
  void clear() override;

  const EventFlowInfo* get(unsigned int idkey) const override;
  EventFlowInfo* get(unsigned int idkey) override;
  const EventFlowInfo* get(EventFlowInfo::EventFlowSrc type) const override;
  EventFlowInfo* get(EventFlowInfo::EventFlowSrc type) override;

  EventFlowInfo* insert(EventFlowInfo* clus, EventFlowInfo::EventFlowSrc id) override;

  ConstIter begin() const override { return _map.begin(); }
  ConstIter find(unsigned int idkey) const override { return _map.find(idkey); }
  ConstIter find(EventFlowInfo::EventFlowSrc type) const override;
  ConstIter end() const override { return _map.end(); }

  Iter begin() override { return _map.begin(); }
  Iter find(unsigned int idkey) override { return _map.find(idkey); }
  Iter find(EventFlowInfo::EventFlowSrc type) override;
  Iter end() override { return _map.end(); }

 private:
  std::map<unsigned int, EventFlowInfo*> _map;


  ClassDefOverride(EventFlowInfoMapv1, 1);
};

#endif  // EVENTFLOWINFOMAPV1_H
