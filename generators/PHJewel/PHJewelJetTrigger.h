// PHJewelJetTrigger -- leading truth jet trigger for PHJewel.
//
// Leading anti-kT(R) jet built from all status-1 particles except neutrinos,
// |eta_jet| < etamax, pT >= ptmin; optionally the leading jet must also lie in
// [window_lo, window_hi) (pT-hat stitching).  Status-3 entries (JEWEL's
// thermal momenta, WRITESCATCEN T) are not used, but the recoils are, so with
// KEEPRECOILS T the trigger sees UNSUBTRACTED jets.  Those are harder than
// subtracted ones: a threshold at X keeps the events whose subtracted jets are
// above X.
#ifndef PHJEWEL_PHJEWELJETTRIGGER_H
#define PHJEWEL_PHJEWELJETTRIGGER_H

#include "PHJewelGenTrigger.h"

#include <string>

class PHJewelJetTrigger : public PHJewelGenTrigger
{
 public:
  explicit PHJewelJetTrigger(const std::string &name = "PHJewelJetTrigger");
  ~PHJewelJetTrigger() override = default;

  bool Apply(const HepMC::GenEvent *evt) override;
  void PrintConfig() const override;

  void SetJetR(double r) { m_R = r; }
  void SetMinJetPt(double pt) { m_PtMin = pt; }
  void SetEtaMax(double eta) { m_EtaMax = eta; }
  void SetLeadJetWindow(double lo, double hi)
  {
    m_WindowLo = lo;
    m_WindowHi = hi;
  }

  // leading-jet pT of the last event given to Apply (0 if none)
  double GetLeadJetPt() const { return m_LeadJetPt; }

 private:
  double m_R{0.4};
  double m_PtMin{0.0};
  double m_EtaMax{1.1};
  double m_WindowLo{0.0};
  double m_WindowHi{1e9};
  double m_LeadJetPt{0.0};
};

#endif
