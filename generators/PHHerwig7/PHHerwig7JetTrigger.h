#ifndef PHHERWIG7_PHHERWIG7JETTRIGGER_H
#define PHHERWIG7_PHHERWIG7JETTRIGGER_H

#include "PHHerwig7GenTrigger.h"

#include <string>

namespace HepMC
{
  class GenEvent;
}

class PHHerwig7JetTrigger : public PHHerwig7GenTrigger
{
 public:
  PHHerwig7JetTrigger(const std::string &name = "PHHerwig7JetTrigger");
  ~PHHerwig7JetTrigger() override;

  bool Apply(HepMC::GenEvent *herwig) override;

  void SetEtaHighLow(double etaHigh, double etaLow);
  void SetMinJetPt(double minPt);
  void SetJetR(double R);
  void SetMinLeadingZ(double minZ);
  void SetMinNumConstituents(int nconst);

  void PrintConfig() const;

 private:
  double _theEtaHigh{4.0};
  double _theEtaLow{1.0};
  double _minPt{0.0};
  double _minZ{0.0};
  double _R{1.0};
  int _nconst{0};
};

#endif
