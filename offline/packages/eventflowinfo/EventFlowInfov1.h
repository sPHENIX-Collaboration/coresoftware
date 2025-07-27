#ifndef EVENTFLOWINFOV1_H
#define EVENTFLOWINFOV1_H

#include "EventFlowInfo.h"

#include <iostream>
#include <limits>
#include <map>

class PHObject;

class EventFlowInfov1 : public EventFlowInfo
{
 public:

  EventFlowInfov1() = default;
  ~EventFlowInfov1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = EventFlowInfov1(); }
  PHObject* CloneMe() const override { return new EventFlowInfov1(*this); }

  void SetPsi(int order, double psi) override { _psi_map[order] = psi; }
  double GetPsi(int order) const override;

  void SetVn(int order, double vn) override { _vn_map[order] = vn; }
  double GetVn(int order) const override;

  void SetSrc(EventFlowSrc src) override { _src = src; }
  EventFlowSrc GetSrc() const override { return _src; }

 private:

  std::map<int, double> _psi_map {};
  std::map<int, double> _vn_map {};
  EventFlowSrc _src = VOID;

  ClassDefOverride(EventFlowInfov1, 1);
};

#endif
     
