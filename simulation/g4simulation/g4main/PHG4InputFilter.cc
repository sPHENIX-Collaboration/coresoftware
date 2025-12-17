#include "PHG4InputFilter.h"
#include "PHG4InEvent.h"
#include "PHG4Particle.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <cmath>
#include <iostream>  // for operator<<, endl, basic_ostream
#include <map>       // for multimap<>::iterator, _Rb_tr...
#include <utility>   // for pair

PHG4InputFilter::PHG4InputFilter(const std::string &name)
  : SubsysReco(name)
{
}

int PHG4InputFilter::process_event(PHCompositeNode *topNode)
{
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!ineve)
  {
    std::cout << PHWHERE << "no PHG4INEVENT node" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  std::pair<std::multimap<int, PHG4Particle *>::iterator, std::multimap<int, PHG4Particle *>::iterator> beginend = ineve->GetParticles_Modify();
  std::multimap<int, PHG4Particle *>::iterator particleiter;
  if (Verbosity() > 0)
  {
    std::cout << "PHG4InputFilter before filter" << std::endl;
    ineve->identify();
  }
  particleiter = beginend.first;
  while (particleiter != beginend.second)
  {
    if (std::isfinite(etamin))
    {
      double eta = get_eta((particleiter->second)->get_px(), (particleiter->second)->get_py(), (particleiter->second)->get_pz());
      if (eta < etamin)
      {
        std::multimap<int, PHG4Particle *>::iterator particleiter_cache = particleiter;
        ++particleiter_cache;
        ineve->DeleteParticle(particleiter);
        particleiter = particleiter_cache;
        continue;
      }
    }
    if (std::isfinite(etamax))
    {
      double eta = get_eta((particleiter->second)->get_px(), (particleiter->second)->get_py(), (particleiter->second)->get_pz());
      if (eta > etamax)
      {
        std::multimap<int, PHG4Particle *>::iterator particleiter_cache = particleiter;
        ++particleiter_cache;
        ineve->DeleteParticle(particleiter);
        particleiter = particleiter_cache;
        continue;
      }
    }
    if (std::isfinite(ptmin))
    {
      double pt = sqrt((particleiter->second)->get_px() * (particleiter->second)->get_px() + (particleiter->second)->get_py() * (particleiter->second)->get_py());
      if (pt < ptmin)
      {
        std::multimap<int, PHG4Particle *>::iterator particleiter_cache = particleiter;
        ++particleiter_cache;
        ineve->DeleteParticle(particleiter);
        particleiter = particleiter_cache;
        continue;
      }
    }
    if (std::isfinite(ptmax))
    {
      double pt = sqrt((particleiter->second)->get_px() * (particleiter->second)->get_px() + (particleiter->second)->get_py() * (particleiter->second)->get_py());
      if (pt > ptmax)
      {
        std::multimap<int, PHG4Particle *>::iterator particleiter_cache = particleiter;
        ++particleiter_cache;
        ineve->DeleteParticle(particleiter);
        particleiter = particleiter_cache;
        continue;
      }
    }
    ++particleiter;
  }
  if (Verbosity() > 0)
  {
    std::cout << "PHG4InputFilter: after filter" << std::endl;
    ineve->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

double
PHG4InputFilter::get_eta(const double x, const double y, const double z)
{
  double eta;
  double radius;
  double theta;
  radius = sqrt(x * x + y * y);
  theta = atan2(radius, z);
  eta = -log(tan(theta / 2.));
  return eta;
}
