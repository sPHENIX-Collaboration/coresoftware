#include "PHG4ParticleGenerator.h"
#include "PHG4Particlev2.h"

#include "PHG4InEvent.h"


#include <phool/getClass.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>

#include <gsl/gsl_randist.h>

#include <cmath>

using namespace std;

PHG4ParticleGenerator::PHG4ParticleGenerator(const string &name): 
  PHG4ParticleGeneratorBase(name),
  z_min(-10.0),
  z_max(10.0),
  eta_min(-1.0),
  eta_max(1.0),
  phi_min(-M_PI),
  phi_max(M_PI),
  mom_min(0.0),
  mom_max(10.0)
{
  return;
}

void
PHG4ParticleGenerator::set_z_range(const double min, const double max)
{
  z_min = min;
  z_max = max;
  return;
}

void
PHG4ParticleGenerator::set_eta_range(const double min, const double max)
{
  eta_min = min;
  eta_max = max;
  return;
}

void
PHG4ParticleGenerator::set_phi_range(const double min, const double max)
{
  phi_min = min;
  phi_max = max;
  return;
}

void
PHG4ParticleGenerator::set_mom_range(const double min, const double max)
{
  mom_min = min;
  mom_max = max;
  return;
}

int
PHG4ParticleGenerator::process_event(PHCompositeNode *topNode)
{
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode,"PHG4INEVENT");

  if (! ReuseExistingVertex(topNode))
    {
      vtx_z = (z_max-z_min)*gsl_rng_uniform_pos(RandomGenerator) + z_min;
    }
  int vtxindex = ineve->AddVtx(vtx_x,vtx_y,vtx_z,t0);

  vector<PHG4Particle *>::iterator iter;
  for (iter = particlelist.begin(); iter != particlelist.end(); ++iter)
    {
      PHG4Particle *particle = new PHG4Particlev2(*iter);
      SetParticleId(particle,ineve);
      double mom = (mom_max - mom_min)*gsl_rng_uniform_pos(RandomGenerator) + mom_min;
      double eta = (eta_max - eta_min)*gsl_rng_uniform_pos(RandomGenerator) + eta_min;
      double phi = (phi_max - phi_min)*gsl_rng_uniform_pos(RandomGenerator) + phi_min;
      double pt = mom/cosh(eta);

      particle->set_e(mom);
      particle->set_px(pt*cos(phi));
      particle->set_py(pt*sin(phi));
      particle->set_pz(pt*sinh(eta));
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
  cout << "PHG4ParticleGenerator settings:" << endl;
  cout << "z_min, z_max: " << z_min << "/" << z_max << endl;
  cout << "eta_min, eta_max: " << eta_min << "/" << eta_max << endl;
  cout << "phi_min, phi_max: " << phi_min << "/" << phi_max << endl;
  cout << "mom_min, mom_max: " << mom_min << "/" << mom_max << endl;
  PrintParticles(what);
  return;
}
