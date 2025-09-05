#include "TruthNeutralMesonv1.h"

#include <phool/phool.h>

#include <iostream>

void TruthNeutralMesonv1::Reset()
{
  m_pid = 0;
  m_parent_pid = 0;
  m_is_prompt = true;
  m_mother_is_eta = false;
  m_is_eta_3pi0 = false;
  m_is_eta_pi0pipm = false;

  m_e = 0;
  m_pt = 0;
  m_eta = 0;
  m_phi = 0;
  m_p = 0;

  m_n_photons = 0;
  for (int i = 0; i < 2; ++i)
  {
    m_photon_e[i] = 0;
    m_photon_pt[i] = 0;
    m_photon_eta[i] = 0;
    m_photon_phi[i] = 0;
    m_photon_p[i] = 0;
    m_photon_converted[i] = false;
  }
}

void TruthNeutralMesonv1::identify(std::ostream& os) const
{
  os << "TruthNeutralMesonv1:"
     << " pid = " << m_pid
     << " origin prompt / feed-down pid = " << m_is_prompt << " / " << m_parent_pid
     << " mother_is_eta = " << m_mother_is_eta
     << " e = " << m_e
     << " pt = " << m_pt
     << " eta = " << m_eta
     << " phi = " << m_phi
     << " p = " << m_p
     << " n_photons = " << m_n_photons
     << std::endl;

  for (int i = 0; i < m_n_photons; ++i)
  {
    os << "  Photon[" << i << "] "
       << " e = " << m_photon_e[i]
       << " pt = " << m_photon_pt[i]
       << " eta = " << m_photon_eta[i]
       << " phi = " << m_photon_phi[i]
       << " p = " << m_photon_p[i]
       << std::endl;
  }
}

void TruthNeutralMesonv1::add_photon(float e, float pt, float eta, float phi, float p, bool isconv)
{
  if (m_n_photons >= 2)
  {
    std::cout << PHWHERE << " Tried to add more than 2 photons!" << std::endl;
    return;
  }

  m_photon_e[m_n_photons] = e;
  m_photon_pt[m_n_photons] = pt;
  m_photon_eta[m_n_photons] = eta;
  m_photon_phi[m_n_photons] = phi;
  m_photon_p[m_n_photons] = p;
  m_photon_converted[m_n_photons] = isconv;

  ++m_n_photons;
}
