#ifndef PHHEPMC_PHHEPMCEVENTv1_H
#define PHHEPMC_PHHEPMCEVENTv1_H

#include "PHHepMCGenEvent.h"

class PHHepMCGenEventv1 : public PHHepMCGenEvent
{
 public:
  PHHepMCGenEventv1();

  PHHepMCGenEventv1(const PHHepMCGenEventv1& event);
  PHHepMCGenEventv1& operator=(const PHHepMCGenEventv1& event);
  ~PHHepMCGenEventv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override
  {
    PHOOL_VIRTUAL_WARNING;
    return 0;
  }
  PHObject* CloneMe() const override { return new PHHepMCGenEventv1(*this); }

  //! boost beta vector for Lorentz Transform, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  const HepMC::ThreeVector& get_boost_beta_vector() const final { return m_boost_beta_vector; }

  //! boost beta vector for Lorentz Transform, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  void set_boost_beta_vector(const HepMC::ThreeVector& v) final { m_boost_beta_vector = v; }

  //! rotation axis vector, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  const HepMC::ThreeVector& get_rotation_vector() const final { return m_rotation_vector; }

  //! rotation axis vector, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  void set_rotation_vector(const HepMC::ThreeVector& v) final { m_rotation_vector = v; }

  //! rotation angle, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  double get_rotation_angle() const final { return m_rotation_angle; }

  //! rotation angle, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  void set_rotation_angle(const double a) final { m_rotation_angle = a; }

  //!LorentzRotation to translate from hepmc event frame to lab frame
  CLHEP::HepLorentzRotation get_LorentzRotation_EvtGen2Lab() const final;

  //!LorentzRotation to translate from lab frame to hepmc event frame
  CLHEP::HepLorentzRotation get_LorentzRotation_Lab2EvtGen() const final;

  //! reaction plane angles thrown by hijing flowAfterburner
  float get_flow_psi(unsigned int n) const final;

  // ! get the reaction plane angle psi_n[n] for all n
  const FlowAfterburner_PsiMap& get_flow_psi_map() const final { return m_psi_n; }

  // ! set the reaction plane angle psi_n[n]
  void set_flow_psi(unsigned int n, float psi) final { m_psi_n[n] = psi; }
  
  // ! set the reaction plane angle psi_n[n] for all n
  void set_flow_psi_map(const FlowAfterburner_PsiMap& psi_map) final { m_psi_n = psi_map; }



 protected:
  //! boost beta vector for Lorentz Transform, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  HepMC::ThreeVector m_boost_beta_vector;

  //! rotation axis vector, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  HepMC::ThreeVector m_rotation_vector;

  //! rotation angle, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  double m_rotation_angle;

  //! reaction plane angles thrown by hijing flowAfterburner
  PHHepMCGenEvent::FlowAfterburner_PsiMap m_psi_n;  // only will be filled if the flowAfterburner enabled and input HepMC is from hijing
  // it is a map of n to psi_n[n]

  ClassDefOverride(PHHepMCGenEventv1, 2)
};

#endif  // PHHEPMC_PHHEPMCEVENTv1_H
