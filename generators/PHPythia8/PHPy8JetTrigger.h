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
  ~PHPy8JetTrigger() override;

  bool Apply(Pythia8::Pythia *pythia) override;

  void SetEtaHighLow(double etaHigh, double etaLow);
  void SetMinJetPt(double minPt);
  void SetJetR(double R);
  void SetMinLeadingZ(double minZ);
  void SetMinNumConstituents(int nconst);

  void PrintConfig();

 private:
  double _theEtaHigh;
  double _theEtaLow;
  double _minPt;
  double _minZ;
  double _R;
  int _nconst;
};

#endif
