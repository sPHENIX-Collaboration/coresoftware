#ifndef PHPYTHIA6_PHPY6JETTRIGGER_H
#define PHPYTHIA6_PHPY6JETTRIGGER_H

#include "PHPy6GenTrigger.h"
#include <HepMC/GenEvent.h>
#include <string>

namespace HepMC
{
  class GenEvent;
};


class PHPy6JetTrigger : public PHPy6GenTrigger {

 public:

  PHPy6JetTrigger(const std::string &name = "PHPy6JetTrigger");
  virtual ~PHPy6JetTrigger();

  #ifndef __CINT__
  bool Apply(const HepMC::GenEvent* evt);
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
