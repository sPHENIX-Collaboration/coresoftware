#include "PHHepMCGenEventv1.h"

#include <HepMC/GenEvent.h>
#include <HepMC/SimpleVector.h>  // for FourVector

#include <CLHEP/Vector/Boost.h>
#include <CLHEP/Vector/LorentzRotation.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Vector/Rotation.h>

#include <sstream>
#include <utility>  // for swap

using namespace std;

PHHepMCGenEventv1::PHHepMCGenEventv1()
  : m_boost_beta_vector(0, 0, 0)
  , m_rotation_vector(0, 0, 1)
  , m_rotation_angle(0)
{
}

PHHepMCGenEventv1::PHHepMCGenEventv1(const PHHepMCGenEventv1& event)
  : PHHepMCGenEvent(event)
  , m_boost_beta_vector(event.get_boost_beta_vector())
  , m_rotation_vector(event.get_rotation_vector())
  , m_rotation_angle(event.get_rotation_angle())
{
  return;
}

PHHepMCGenEventv1& PHHepMCGenEventv1::operator=(const PHHepMCGenEventv1& event)
{
  if (&event == this) return *this;

  Reset();

  _embedding_id = event.get_embedding_id();
  _isSimulated = event.is_simulated();
  _theEvt = new HepMC::GenEvent(*event.getEvent());

  return *this;
}

PHHepMCGenEventv1::~PHHepMCGenEventv1()
{
}

void PHHepMCGenEventv1::Reset()
{
  PHHepMCGenEvent::Reset();

  m_boost_beta_vector.set(0, 0, 0);
  m_rotation_vector.set(0, 0, 1);
  m_rotation_angle = 0;
}

//_____________________________________________________________________________
void PHHepMCGenEventv1::identify(std::ostream& os) const
{
  PHHepMCGenEvent::identify(os);

  os << " m_boost_beta_vector = (" << m_boost_beta_vector.x() << "," << m_boost_beta_vector.y() << "," << m_boost_beta_vector.z() << ") " << endl;
  os << " m_rotation_vector = (" << m_rotation_vector.x() << "," << m_rotation_vector.y() << "," << m_rotation_vector.z() << ") by " << m_rotation_angle << " rad" << endl;

  static const CLHEP::HepLorentzVector zp_lightcone(0, 0, 1, 1);
  static const CLHEP::HepLorentzVector zm_lightcone(0, 0, -1, 1);

  os << " HepMC Frame unit light cone vector along +z axis "<<zp_lightcone<<" translate to lab at : "
     << (get_LorentzRotation_EvtGen2Lab() * zp_lightcone)
     << endl;
  os << " HepMC Frame unit light cone vector along -z axis  "<<zm_lightcone<<" translate to lab at : "
     << (get_LorentzRotation_EvtGen2Lab() * zm_lightcone)
     << endl;

  return;
}

CLHEP::HepLorentzRotation PHHepMCGenEventv1::get_LorentzRotation_EvtGen2Lab() const
{
  CLHEP::HepBoost boost(m_boost_beta_vector.x(), m_boost_beta_vector.y(), m_boost_beta_vector.z());
  CLHEP::Hep3Vector axis(m_rotation_vector.x(), m_rotation_vector.y(), m_rotation_vector.z());
  CLHEP::HepRotation rotation(axis, m_rotation_angle);

  return CLHEP::HepLorentzRotation(boost, rotation);
}

CLHEP::HepLorentzRotation PHHepMCGenEventv1::get_LorentzRotation_Lab2EvtGen() const
{
  return get_LorentzRotation_EvtGen2Lab().inverse();
}
