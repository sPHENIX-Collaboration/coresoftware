#ifndef PHHERWIG7_PHHERWIG7PARTICLETRIGGER_H
#define PHHERWIG7_PHHERWIG7PARTICLETRIGGER_H

#include "PHHerwig7GenTrigger.h"

#include <limits>
#include <string>
#include <vector>

namespace HepMC
{
  class GenEvent;
}

//! Particle trigger for PHHerwig7, with the same interface as PHPy8ParticleTrigger.
//! An event passes if any particle in the HepMC record satisfies all configured cuts.
class PHHerwig7ParticleTrigger : public PHHerwig7GenTrigger
{
 public:
  PHHerwig7ParticleTrigger(const std::string &name = "PHHerwig7ParticleTrigger");
  ~PHHerwig7ParticleTrigger() override;

  bool Apply(HepMC::GenEvent *herwig) override;

  //! trigger particle PDG ids (OR'd, sign matters)
  void AddParticles(const std::string &particles);
  void AddParticles(int particle);
  void AddParticles(const std::vector<int> &particles);

  //! required immediate parent PDG ids (OR'd, sign ignored); none registered = no parent requirement
  void AddParents(const std::string &parents);
  void AddParents(int parent);
  void AddParents(const std::vector<int> &parents);

  void SetPtHigh(double pt) { m_PtHigh = pt; }
  void SetPtLow(double pt) { m_PtLow = pt; }
  void SetPtHighLow(double ptHigh, double ptLow);

  void SetPHigh(double p) { m_PHigh = p; }
  void SetPLow(double p) { m_PLow = p; }
  void SetPHighLow(double pHigh, double pLow);

  //! rapidity cuts
  void SetYHigh(double Y) { m_YHigh = Y; }
  void SetYLow(double Y) { m_YLow = Y; }
  void SetYHighLow(double YHigh, double YLow);

  void SetEtaHigh(double eta) { m_EtaHigh = eta; }
  void SetEtaLow(double eta) { m_EtaLow = eta; }
  void SetEtaHighLow(double etaHigh, double etaLow);

  void SetAbsEtaHigh(double eta) { m_AbsEtaHigh = eta; }
  void SetAbsEtaLow(double eta) { m_AbsEtaLow = eta; }
  void SetAbsEtaHighLow(double etaHigh, double etaLow);

  void SetPzHigh(double pz) { m_PzHigh = pz; }
  void SetPzLow(double pz) { m_PzLow = pz; }
  void SetPzHighLow(double pzHigh, double pzLow);

  //! true (default): only final-state particles (HepMC status 1); false: all copies in the record
  void SetStableParticleOnly(bool b) { m_StableParticleOnly = b; }

  //! true: reject particles with a hadron as immediate parent (decay photons, for instance),
  //! keeping direct and fragmentation photons. Any hadron parent vetoes, except p/n (beam remnants).
  void SetPromptOnly(bool b) { m_PromptOnly = b; }

  void PrintConfig() const;

 private:
  static constexpr double kNoCut = std::numeric_limits<double>::infinity();

  std::vector<int> m_Particles;
  std::vector<int> m_Parents;

  double m_PtLow{-kNoCut};
  double m_PtHigh{kNoCut};
  double m_PLow{-kNoCut};
  double m_PHigh{kNoCut};
  double m_YLow{-kNoCut};
  double m_YHigh{kNoCut};
  double m_EtaLow{-kNoCut};
  double m_EtaHigh{kNoCut};
  double m_AbsEtaLow{-kNoCut};
  double m_AbsEtaHigh{kNoCut};
  double m_PzLow{-kNoCut};
  double m_PzHigh{kNoCut};

  bool m_StableParticleOnly{true};
  bool m_PromptOnly{false};
};

#endif
