#include "PHG4ParticleGenerator.h"

#include "PHG4InEvent.h"
#include "PHG4Particle.h"  // for PHG4Particle
#include "PHG4Particlev2.h"

#include <phool/getClass.h>

#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos

#include <cmath>
#include <iostream>  // for operator<<, endl, basic_ostream, basic_o...
#include <vector>    // for vector, vector<>::iterator

class PHCompositeNode;

PHG4ParticleGenerator::PHG4ParticleGenerator(const std::string &name)
  : PHG4ParticleGeneratorBase(name)
{
  return;
}

void PHG4ParticleGenerator::set_z_range(const double min, const double max)
{
  m_ZMin = min;
  m_ZMax = max;
  return;
}

void PHG4ParticleGenerator::set_eta_range(const double min, const double max)
{
  m_EtaMin = min;
  m_EtaMax = max;
  return;
}

void PHG4ParticleGenerator::set_phi_range(const double min, const double max)
{
  m_PhiMin = min;
  m_PhiMax = max;
  return;
}

void PHG4ParticleGenerator::set_mom_range(const double min, const double max)
{
  m_MomMin = min;
  m_MomMax = max;
  return;
}

int PHG4ParticleGenerator::process_event(PHCompositeNode *topNode)
{
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");

  if (!ReuseExistingVertex(topNode))
  {
    set_vtx_z((m_ZMax - m_ZMin) * gsl_rng_uniform_pos(RandomGenerator()) + m_ZMin);
  }
  int vtxindex = ineve->AddVtx(get_vtx_x(), get_vtx_y(), get_vtx_z(), get_t0());

  std::vector<PHG4Particle *>::iterator iter;
  for (iter = particlelist_begin(); iter != particlelist_end(); ++iter)
  {
    PHG4Particle *particle = new PHG4Particlev2(*iter);
    SetParticleId(particle, ineve);
    double mom = (m_MomMax - m_MomMin) * gsl_rng_uniform_pos(RandomGenerator()) + m_MomMin;
    double eta = (m_EtaMax - m_EtaMin) * gsl_rng_uniform_pos(RandomGenerator()) + m_EtaMin;
    double phi = (m_PhiMax - m_PhiMin) * gsl_rng_uniform_pos(RandomGenerator()) + m_PhiMin;
    double pt = mom / cosh(eta);

    particle->set_e(mom);
    particle->set_px(pt * cos(phi));
    particle->set_py(pt * sin(phi));
    particle->set_pz(pt * sinh(eta));
    // update internal particle list with changed momenta
    // needed for correct printout of particle kinematics
    // in PHG4ParticleGenerator::Print()
    (*iter)->set_e(particle->get_e());
    (*iter)->set_px(particle->get_px());
    (*iter)->set_py(particle->get_py());
    (*iter)->set_pz(particle->get_pz());
    ineve->AddParticle(vtxindex, particle);
  }
  if (Verbosity() > 0)
  {
    ineve->identify();
  }
  return 0;
}

void PHG4ParticleGenerator::Print(const std::string &what) const
{
  std::cout << "PHG4ParticleGenerator settings:" << std::endl;
  std::cout << "ZMin, ZMax: " << m_ZMin << "/" << m_ZMax << std::endl;
  std::cout << "EtaMin, EtaMax: " << m_EtaMin << "/" << m_EtaMax << std::endl;
  std::cout << "PhiMin, PhiMax: " << m_PhiMin << "/" << m_PhiMax << std::endl;
  std::cout << "MomMin, MomMax: " << m_MomMin << "/" << m_MomMax << std::endl;
  PrintParticles(what);
  return;
}
