#ifndef PHPYTHIA8_PHPY8PARTICLETRIGGER_H
#define PHPYTHIA8_PHPY8PARTICLETRIGGER_H

#include "PHPy8GenTrigger.h"

#include <string>
#include <vector>  // for vector

namespace Pythia8
{
  class Pythia;
}

class PHPy8ParticleTrigger : public PHPy8GenTrigger
{
 public:
  PHPy8ParticleTrigger(const std::string &name = "PHPy8ParticleTrigger");
  ~PHPy8ParticleTrigger() override;

  bool Apply(Pythia8::Pythia *pythia) override;

  void AddParticles(const std::string &particles);
  void AddParticles(int particle);
  void AddParticles(std::vector<int> particles);

  void AddParents(const std::string &parents);
  void AddParents(int parent);
  void AddParents(std::vector<int> parents);

  void SetPtHigh(double pt);
  void SetPtLow(double pt);
  void SetPtHighLow(double ptHigh, double ptLow);

  void SetPHigh(double p);
  void SetPLow(double p);
  void SetPHighLow(double pHigh, double pLow);

  //! rapidity cuts
  void SetYHigh(double Y);
  void SetYLow(double Y);
  void SetYHighLow(double YHigh, double YLow);

  void SetEtaHigh(double eta);
  void SetEtaLow(double eta);
  void SetEtaHighLow(double etaHigh, double etaLow);

  void SetAbsEtaHigh(double eta);
  void SetAbsEtaLow(double eta);
  void SetAbsEtaHighLow(double etaHigh, double etaLow);

  void SetPzHigh(double pz);
  void SetPzLow(double pz);
  void SetPzHighLow(double pzHigh, double pzLow);

  //! Whether to apply the criteria to unstable particles in the Pythia records too (default = true)
  void SetStableParticleOnly(bool b) { m_doStableParticleOnly = b; }

  void PrintConfig();

 private:
  std::vector<int> _theParents;
  std::vector<int> _theParticles;

  double _theYHigh{999.9};
  double _theYLow{999.9};
  double _theEtaHigh{999.9};
  double _theEtaLow{999.9};
  double _thePtHigh{999.9};
  double _thePtLow{999.9};
  double _thePHigh{999.9};
  double _thePLow{999.9};
  double _thePzHigh{999.9};
  double _thePzLow{999.9};

  bool _doYHighCut{false};
  bool _doYLowCut{false};
  bool _doBothYCut{false};
  bool _doEtaHighCut{false};
  bool _doEtaLowCut{false};
  bool _doBothEtaCut{false};
  bool _doAbsEtaHighCut{false};
  bool _doAbsEtaLowCut{false};
  bool _doBothAbsEtaCut{false};
  bool _doPtHighCut{false};
  bool _doPtLowCut{false};
  bool _doBothPtCut{false};
  bool _doPHighCut{false};
  bool _doPLowCut{false};
  bool _doBothPCut{false};
  bool _doPzHighCut{false};
  bool _doPzLowCut{false};
  bool _doBothPzCut{false};

  bool m_doStableParticleOnly = true;
};

#endif
