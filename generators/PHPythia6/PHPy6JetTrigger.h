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
  void SetMinJetPt(double minPt) { m_minPt = minPt;}
  void SetJetR(double R) {m_R = R;}

  void PrintConfig();

 private:

  double m_theEtaHigh;
  double m_theEtaLow;
  double m_minPt;
  double m_R;

};

#endif	
