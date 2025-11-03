#ifndef TRUTHNEUTRALMESON_H
#define TRUTHNEUTRALMESON_H

#include <phool/PHObject.h>
#include <array>
#include <iostream>

class TruthNeutralMeson : public PHObject
{
 public:
  TruthNeutralMeson() = default;
  virtual ~TruthNeutralMeson();

  virtual int get_pid() const { return 0; }
  virtual void set_pid(int /*unused*/) {}

  virtual bool is_prompt() const { return true; }
  virtual void set_prompt(bool /*unused*/) {}

  virtual int get_parent_pid() const { return 0; }
  virtual void set_parent_pid(int /*unused*/) {}

  virtual bool mother_is_eta() const { return false; }
  virtual void set_mother_is_eta(bool /*unused*/) {}

  virtual float get_e() const { return 0; }
  virtual float get_pt() const { return 0; }
  virtual float get_eta() const { return 0; }
  virtual float get_phi() const { return 0; }
  virtual float get_p() const { return 0; }

  virtual void set_e(float /*unused*/) {}
  virtual void set_pt(float /*unused*/) {}
  virtual void set_eta(float /*unused*/) {}
  virtual void set_phi(float /*unused*/) {}
  virtual void set_p(float /*unused*/) {}

  virtual int get_n_photons() const { return 0; }

  virtual float get_photon_e(int /*unused*/) const { return 0; }
  virtual float get_photon_pt(int /*unused*/) const { return 0; }
  virtual float get_photon_eta(int /*unused*/) const { return 0; }
  virtual float get_photon_phi(int /*unused*/) const { return 0; }
  virtual float get_photon_p(int /*unused*/) const { return 0; }
  virtual float get_photon_converted(int /*unused*/) const { return false; }

  virtual void add_photon(float e, float pt, float eta, float phi, float p, bool isconv)
  {
    (void) e;
    (void) pt;
    (void) eta;
    (void) phi;
    (void) p;
    (void) isconv;
  }

  virtual bool is_diphoton() const { return get_n_photons() == 2; }

  virtual bool is_eta_3pi0() const { return false; }
  virtual void set_eta_3pi0(bool /*unused*/) {}

  virtual bool is_eta_pi0pipm() const { return false; }
  virtual void set_eta_pi0pipm(bool /*unused*/) {}

  void Reset() override {}
  int isValid() const override { return 1; }
  PHObject* CloneMe() const override { return nullptr; }
  void identify(std::ostream& os = std::cout) const override { os << "TruthNeutralMeson base" << std::endl; }

 protected:
  ClassDefOverride(TruthNeutralMeson, 1);
};

#endif
