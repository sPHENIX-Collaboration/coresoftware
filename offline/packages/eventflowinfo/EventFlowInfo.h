#ifndef EVENTFLOWINFO_H
#define EVENTFLOWINFO_H

#include <phool/PHObject.h>

#include <iostream>
#include <cmath>

class EventFlowInfo : public PHObject
{
 public:

  
  enum EventFlowSrc
  {
    VOID = -1,
    sEPD = 0,
    MBD = 1,
    CALO = 2
  };



  
  ~EventFlowInfo() override {}

  void identify(std::ostream& os = std::cout) const override {
    os << "EventFlowInfo base class" << std::endl;
  }
  PHObject* CloneMe() const override { return nullptr; }



  virtual void SetPsi( int /*order*/, double /*psi*/) { return; }
  virtual double GetPsi( int /*order*/) const { return NAN; }
  virtual void SetVn( int /*order*/, double /*vn*/) { return; }
  virtual double GetVn( int /*order*/) const { return NAN; }
  virtual void SetSrc(EventFlowSrc /*src*/) { return; }
  virtual EventFlowSrc GetSrc() const { return VOID; }


 protected:
  EventFlowInfo() {}

 private:
  ClassDefOverride(EventFlowInfo, 1);
};

#endif
