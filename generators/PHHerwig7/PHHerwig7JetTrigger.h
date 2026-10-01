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

  bool Apply(HepMC::GenEvent* herwig) override;

  void SetEtaHighLow(double etaHigh, double etaLow);
  void SetMinJetPt(double minPt);
  void SetJetR(double R);
  void SetMinLeadingZ(double minZ);
  void SetMinNumConstituents(int nconst);

  void PrintConfig() const;

 private:
  double _theEtaHigh;
  double _theEtaLow;
  double _minPt;
  double _minZ;
  double _R;
  int _nconst;
};

#endif
