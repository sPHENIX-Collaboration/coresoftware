#ifndef __PHPY8FWDJETTRIGGER_H__
#define __PHPY8FWDJETTRIGGER_H__

#include "PHPy8GenTrigger.h"
#include <string>

namespace Pythia8 {
  class Pythia;
};

class PHPy8FwdJetTrigger : public PHPy8GenTrigger {

 public:

  PHPy8FwdJetTrigger(const std::string &name = "PHPy8FwdJetTrigger");
  virtual ~PHPy8FwdJetTrigger();

  #ifndef __CINT__
  bool Apply(Pythia8::Pythia *pythia);
  #endif

  void SetEtaHighLow(double etaHigh, double etaLow);
  void SetMinJetPt(double minPt);
  void SetJetR(double R);

  void PrintConfig();

 private:

  double _theEtaHigh;
  double _theEtaLow;
  double _minPt; 
  double _R; 

};

#endif	
