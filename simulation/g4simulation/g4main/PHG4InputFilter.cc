#include "PHG4InputFilter.h"
#include "PHG4InEvent.h"
#include "PHG4Particle.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/phool.h>                 // for PHWHERE

#include <cmath>
#include <iostream>                      // for operator<<, endl, basic_ostream
#include <map>                           // for multimap<>::iterator, _Rb_tr...
#include <utility>                       // for pair

using namespace std;

PHG4InputFilter::PHG4InputFilter(const std::string &name):
  SubsysReco(name),
  etamin(NAN),
  etamax(NAN),
  ptmin(NAN),
  ptmax(NAN)
{}

int
PHG4InputFilter::process_event(PHCompositeNode *topNode)
{
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!ineve)
    {
      cout << PHWHERE << "no PHG4INEVENT node" << endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  pair<multimap<int, PHG4Particle *>::iterator, multimap<int, PHG4Particle *>::iterator > beginend =   ineve->GetParticles_Modify();
  multimap<int, PHG4Particle *>::iterator particleiter;
  if (Verbosity() > 0)
    {
      cout << "PHG4InputFilter before filter" << endl;
      ineve->identify();
    }
  particleiter = beginend.first;
  while (particleiter != beginend.second)
    {
      if (isfinite(etamin))
        {
          double eta = get_eta((particleiter->second)->get_px(), (particleiter->second)->get_py(), (particleiter->second)->get_pz());
          if (eta < etamin)
            {
              multimap<int, PHG4Particle *>::iterator particleiter_cache = particleiter  ;
              ++particleiter_cache;
              ineve->DeleteParticle(particleiter);
              particleiter = particleiter_cache;
              continue;
            }
        }
      if (isfinite(etamax))
        {
          double eta = get_eta((particleiter->second)->get_px(), (particleiter->second)->get_py(), (particleiter->second)->get_pz());
          if (eta > etamax)
            {
              multimap<int, PHG4Particle *>::iterator particleiter_cache = particleiter  ;
              ++particleiter_cache;
              ineve->DeleteParticle(particleiter);
              particleiter = particleiter_cache;
              continue;
            }
        }
      if (isfinite(ptmin))
        {
          double pt = sqrt((particleiter->second)->get_px()*(particleiter->second)->get_px()+ (particleiter->second)->get_py()*(particleiter->second)->get_py());
          if (pt < ptmin)
            {
              multimap<int, PHG4Particle *>::iterator particleiter_cache = particleiter  ;
              ++particleiter_cache;
              ineve->DeleteParticle(particleiter);
              particleiter = particleiter_cache;
              continue;
            }
        }
      if (isfinite(ptmax))
        {
          double pt = sqrt((particleiter->second)->get_px()*(particleiter->second)->get_px()+ (particleiter->second)->get_py()*(particleiter->second)->get_py());
          if (pt > ptmax)
            {
              multimap<int, PHG4Particle *>::iterator particleiter_cache = particleiter  ;
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
      cout << "PHG4InputFilter: after filter" << endl;
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
