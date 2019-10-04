#ifndef PHPYTHIA8_PHPY8JETTRIGGER_H
#define PHPYTHIA8_PHPY8JETTRIGGER_H

#include "PHPy8GenTrigger.h"

#include <string>

namespace Pythia8
{
  class Pythia;
}

class PHPy8JetTrigger : public PHPy8GenTrigger
{
 public:
  PHPy8JetTrigger(const std::string &name = "PHPy8JetTrigger");
  virtual ~PHPy8JetTrigger();

  bool Apply(Pythia8::Pythia *pythia);

  void SetEtaHighLow(double etaHigh, double etaLow);
  void SetMinJetPt(double minPt);
  void SetJetR(double R);
  void SetMinLeadingZ(double minZ);

  void PrintConfig();

 private:
  double _theEtaHigh;
  double _theEtaLow;
  double _minPt;
  double _minZ;
  double _R;
};

#endif
