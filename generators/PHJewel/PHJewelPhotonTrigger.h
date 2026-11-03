// PHJewelPhotonTrigger -- truth photon trigger for PHJewel.
//
// Passes if any status-1 photon has pT >= ptmin and |eta| <= etamax.
// Optionally isolated (scalar pT sum of all other status-1 particles except
// neutrinos inside iso_r below iso_frac * pT) and/or prompt (no hadron among
// the photon's immediate parents).
#ifndef PHJEWEL_PHJEWELPHOTONTRIGGER_H
#define PHJEWEL_PHJEWELPHOTONTRIGGER_H

#include "PHJewelGenTrigger.h"

#include <string>

class PHJewelPhotonTrigger : public PHJewelGenTrigger
{
 public:
  explicit PHJewelPhotonTrigger(const std::string &name = "PHJewelPhotonTrigger");
  ~PHJewelPhotonTrigger() override = default;

  bool Apply(const HepMC::GenEvent *evt) override;
  void PrintConfig() const override;

  void SetMinPt(double pt) { m_PtMin = pt; }
  void SetEtaMax(double eta) { m_EtaMax = eta; }
  void SetIsolation(bool on, double r = 0.4, double frac = 0.1)
  {
    m_Iso = on;
    m_IsoR = r;
    m_IsoFrac = frac;
  }
  void SetPrompt(bool b = true) { m_Prompt = b; }

  // pT of the highest photon that passed in the last event (0 if none)
  double GetLeadPhotonPt() const { return m_LeadPhotonPt; }

 private:
  double m_PtMin{0.0};
  double m_EtaMax{1.1};
  bool m_Iso{false};
  double m_IsoR{0.4};
  double m_IsoFrac{0.1};
  bool m_Prompt{false};
  double m_LeadPhotonPt{0.0};
};

#endif
