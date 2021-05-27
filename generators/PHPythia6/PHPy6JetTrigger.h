#ifndef PHPYTHIA6_PHPY6JETTRIGGER_H
#define PHPYTHIA6_PHPY6JETTRIGGER_H

#include "PHPy6GenTrigger.h"

#include <string>

namespace HepMC
{
class GenEvent;
}

class PHPy6JetTrigger : public PHPy6GenTrigger
{
 public:
  PHPy6JetTrigger(const std::string& name = "PHPy6JetTrigger");
  ~PHPy6JetTrigger() override;

  bool Apply(const HepMC::GenEvent* evt) override;

  void SetEtaHighLow(double etaHigh, double etaLow);
  void SetMinJetPt(double minPt) { m_minPt = minPt; }
  void SetJetR(double R) { m_R = R; }
  void SetMinNumConstituents(int nconst) { m_nconst = nconst; }

  void PrintConfig();

 private:
  double m_theEtaHigh = 4.;
  double m_theEtaLow = 1.;
  double m_minPt = 10.;
  double m_R = 1.;
  int m_nconst = 0;
};

#endif
