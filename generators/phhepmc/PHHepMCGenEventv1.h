#ifndef PHHEPMC_PHHEPMCEVENTv1_H
#define PHHEPMC_PHHEPMCEVENTv1_H

#include "PHHepMCGenEvent.h"

class PHHepMCGenEventv1 : public PHHepMCGenEvent
{
 public:
  PHHepMCGenEventv1();

  PHHepMCGenEventv1(const PHHepMCGenEventv1& event);
  PHHepMCGenEventv1& operator=(const PHHepMCGenEventv1& event);
  virtual ~PHHepMCGenEventv1();

  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset();
  virtual int isValid() const
  {
    PHOOL_VIRTUAL_WARNING;
    return 0;
  }
  PHObject* CloneMe() const { return new PHHepMCGenEventv1(*this); }

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

 protected:
  //! collision vertex position in the Hall coordinate system, use PHENIX units of cm, ns
  HepMC::FourVector _collisionVertex;

  //! boost beta vector for Lorentz Transform, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  HepMC::ThreeVector m_boost_beta_vector;

  //! rotation axis vector, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  HepMC::ThreeVector m_rotation_vector;

  //! rotation angle, part of composition of a LorentzRotation to translate from hepmc event frame to lab frame
  double m_rotation_angle;

  ClassDef(PHHepMCGenEventv1, 1)
};

#endif  // PHHEPMC_PHHEPMCEVENTv1_H
