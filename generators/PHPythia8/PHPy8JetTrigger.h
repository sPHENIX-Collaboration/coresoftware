#ifndef __PHPY8JETTRIGGER_H__
#define __PHPY8JETTRIGGER_H__

#include "PHPy8GenTrigger.h"
#include <string>

namespace Pythia8 {
  class Pythia;
};

class PHPy8JetTrigger : public PHPy8GenTrigger {

 public:

  PHPy8JetTrigger(const std::string &name = "PHPy8JetTrigger");
  virtual ~PHPy8JetTrigger();

  #ifndef __CINT__
  bool Apply(Pythia8::Pythia *pythia);
  #endif

  void SetEtaHighLow(double etaHigh, double etaLow);
  void SetMinJetPt(double minPt);
  void SetJetR(double R);
  void SetMinJetEnergy(double minenergy);
  void SetMinJetConstituents(int numconstituents);
  void PrintConfig();

 private:
  double _minenergy;
  double _theEtaHigh;
  double _theEtaLow;
  double _minPt; 
  double _R; 
  double _nconstituents;

};

#endif	
