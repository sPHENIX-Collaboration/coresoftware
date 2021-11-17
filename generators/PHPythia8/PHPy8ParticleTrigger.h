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
  void SetYHigh(double eta);
  void SetYLow(double eta);
  void SetYHighLow(double etaHigh, double etaLow);

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

  double _theYHigh = 999, _theYLow = -999;
  double _theEtaHigh, _theEtaLow;
  double _thePtHigh, _thePtLow;
  double _thePHigh, _thePLow;
  double _thePzHigh, _thePzLow;

  bool _doYHighCut = false, _doYLowCut = false, _doBothYCut = false;
  bool _doEtaHighCut, _doEtaLowCut, _doBothEtaCut;
  bool _doAbsEtaHighCut, _doAbsEtaLowCut, _doBothAbsEtaCut;
  bool _doPtHighCut, _doPtLowCut, _doBothPtCut;
  bool _doPHighCut, _doPLowCut, _doBothPCut;
  bool _doPzHighCut, _doPzLowCut, _doBothPzCut;

  bool m_doStableParticleOnly = true;
};

#endif
