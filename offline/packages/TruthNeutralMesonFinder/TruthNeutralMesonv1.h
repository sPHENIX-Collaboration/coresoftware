#ifndef TRUTHNEUTRALMESONV1_H
#define TRUTHNEUTRALMESONV1_H

#include <array>
#include "TruthNeutralMeson.h"

class TruthNeutralMesonv1 : public TruthNeutralMeson
{
 public:
  TruthNeutralMesonv1() = default;
  ~TruthNeutralMesonv1() override = default;

  int get_pid() const override { return m_pid; }
  void set_pid(int pid) override { m_pid = pid; }

  bool is_prompt() const override { return m_is_prompt; }
  void set_prompt(bool isprompt) override { m_is_prompt = isprompt; }

  int get_parent_pid() const override { return m_parent_pid; }
  void set_parent_pid(int pid) override { m_parent_pid = pid; }

  bool mother_is_eta() const override { return m_mother_is_eta; }
  void set_mother_is_eta(bool ismothereta) override { m_mother_is_eta = ismothereta; }

  float get_e() const override { return m_e; }
  float get_pt() const override { return m_pt; }
  float get_eta() const override { return m_eta; }
  float get_phi() const override { return m_phi; }
  float get_p() const override { return m_p; }

  void set_e(float val) override { m_e = val; }
  void set_pt(float val) override { m_pt = val; }
  void set_eta(float val) override { m_eta = val; }
  void set_phi(float val) override { m_phi = val; }
  void set_p(float val) override { m_p = val; }

  int get_n_photons() const override { return m_n_photons; }
  float get_photon_e(int i) const override { return (i >= 0 && i < m_n_photons) ? m_photon_e[i] : 0; }
  float get_photon_pt(int i) const override { return (i >= 0 && i < m_n_photons) ? m_photon_pt[i] : 0; }
  float get_photon_eta(int i) const override { return (i >= 0 && i < m_n_photons) ? m_photon_eta[i] : 0; }
  float get_photon_phi(int i) const override { return (i >= 0 && i < m_n_photons) ? m_photon_phi[i] : 0; }
  float get_photon_p(int i) const override { return (i >= 0 && i < m_n_photons) ? m_photon_p[i] : 0; }
  float get_photon_converted(int i) const override { return (i >= 0 && i < m_n_photons) ? m_photon_converted[i] : false; }

  void add_photon(float e, float pt, float eta, float phi, float p, bool isconv) override;

  bool is_eta_3pi0() const override { return m_is_eta_3pi0; }
  void set_eta_3pi0(bool isdecay) override { m_is_eta_3pi0 = isdecay; }

  bool is_eta_pi0pipm() const override { return m_is_eta_pi0pipm; }
  void set_eta_pi0pipm(bool isdecay) override { m_is_eta_pi0pipm = isdecay; }

  void identify(std::ostream& os = std::cout) const override;
  PHObject* CloneMe() const override { return new TruthNeutralMesonv1(*this); }
  void Reset() override;
  int isValid() const override { return 1; }

 private:
  bool m_is_eta_3pi0{false};
  bool m_is_eta_pi0pipm{false};
  bool m_is_prompt{true};
  bool m_mother_is_eta{false};

  int m_pid{0};
  int m_parent_pid{0};

  float m_e{0};
  float m_eta{0};
  float m_p{0};
  float m_phi{0};
  float m_pt{0};

  int m_n_photons{0};
  std::array<float, 2> m_photon_e{{0, 0}};
  std::array<float, 2> m_photon_eta{{0, 0}};
  std::array<float, 2> m_photon_p{{0, 0}};
  std::array<float, 2> m_photon_phi{{0, 0}};
  std::array<float, 2> m_photon_pt{{0, 0}};
  std::array<bool, 2> m_photon_converted{{false, false}};

  ClassDefOverride(TruthNeutralMesonv1, 1);
};

#endif
