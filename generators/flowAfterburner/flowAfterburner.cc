// File:  Generators/FlowAfterburnber/AddFlowByShifting.cxx
// Description:
//    This code is used to introduce particle flow
//    to particles from generated events
//
// AuthorList:
// Andrzej Olszewski: Initial Code February 2006
// 11.10.2006: Add predefined flow function by name

// The initialization of this class is trivial.  There's really no
// need for it.  Pass in a pointer to a random generator, algorithm
// selection and an event.

#include "flowAfterburner.h"
#include "AfterburnerAlgo.h"

#include <phool/phool.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>  // for GenParticle
#include <HepMC/GenRanges.h>
#include <HepMC/GenVertex.h>      // for GenVertex, GenVertex::part...
#include <HepMC/HeavyIon.h>       // for HeavyIon
#include <HepMC/IteratorRange.h>  // for children, descendants
#include <HepMC/SimpleVector.h>   // for FourVector

#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Random/MTwistEngine.h>

#include <cmath>
#include <cstdlib>  // for exit
#include <iostream>
#include <map>  // for map
#include <ctime>  // for time
#include <string>
#include <algorithm>  // for max, min

namespace CLHEP
{
  class HepRandomEngine;
}

Afterburner::Afterburner( const std::string &algorithmName,
                          CLHEP::HepRandomEngine *engine,
                          float mineta, float maxeta,
                          float minpt, float maxpt )
  : m_algo(new AfterburnerAlgo(AfterburnerAlgo::getAlgoFromName(algorithmName)))
  , m_engine(engine)
  , m_ownAlgo(true)
  , m_ownEngine(engine == nullptr)
  , m_mineta(mineta)
  , m_maxeta(maxeta)
  , m_minpt(minpt)
  , m_maxpt(maxpt)
  , m_phishift(0.0)
{
    
    if (!m_engine)
    {
      // seed with time if no engine is provided
      m_engine = new CLHEP::MTwistEngine(static_cast<long>(time(nullptr)));
    }

}


Afterburner::~Afterburner() 
{
  if (m_ownAlgo) 
  {
    delete m_algo;
  }
  if (m_ownEngine)
  {
    delete m_engine;
  }
}

Afterburner::Afterburner(Afterburner&& o) noexcept
  : m_algo(o.m_algo)
  , m_engine(o.m_engine)
  , m_ownAlgo(o.m_ownAlgo)
  , m_ownEngine(o.m_ownEngine)
  , m_mineta(o.m_mineta)
  , m_maxeta(o.m_maxeta)
  , m_minpt(o.m_minpt)
  , m_maxpt(o.m_maxpt)
  , m_phishift(o.m_phishift)
{
  std::copy(std::begin(o.m_psi_n), std::end(o.m_psi_n), std::begin(m_psi_n));
  o.m_algo = nullptr;  
  o.m_engine = nullptr;
  o.m_ownAlgo = false; 
  o.m_ownEngine = false;
}


Afterburner& Afterburner::operator=(Afterburner&& o) noexcept 
{
  if (this != &o) {
    if (m_ownAlgo) 
    {
      delete m_algo;
    }
    if (m_ownEngine)
    {
      delete m_engine;
    }
    m_algo = o.m_algo;             
    o.m_algo = nullptr;
    m_engine = o.m_engine;         
    o.m_engine = nullptr;
    m_ownAlgo = o.m_ownAlgo;       
    o.m_ownAlgo = false;
    m_ownEngine = o.m_ownEngine;   
    o.m_ownEngine = false;
    m_mineta = o.m_mineta;
    m_maxeta = o.m_maxeta;
    m_minpt  = o.m_minpt;  
    m_maxpt  = o.m_maxpt;
    m_phishift = o.m_phishift;
    std::copy(std::begin(o.m_psi_n), std::end(o.m_psi_n), std::begin(m_psi_n));
  }
  return *this;
}


void Afterburner::setAlgo(AfterburnerAlgo::flowAfterburnerAlgorithm algo_type) 
{
  if (m_ownAlgo && m_algo)
  {
    delete m_algo;
  } 
  m_algo = new AfterburnerAlgo(algo_type);
  m_ownAlgo = true;
}

void Afterburner::setAlgo(const std::string &name) 
{
  setAlgo(AfterburnerAlgo::getAlgoFromName(name));
}
void Afterburner::setAlgo(AfterburnerAlgo* algo) 
{
  if (m_ownAlgo && m_algo)
  {
    delete m_algo;
  }
  m_algo = algo;
  m_ownAlgo = false; // external ownership
}

void Afterburner::setEngine(CLHEP::HepRandomEngine* engine)
{
  if (m_ownEngine && m_engine)
  {
    delete m_engine;
  }
  m_engine = engine;
  m_ownEngine = false; // external ownership
}

void Afterburner::setEtaRange(float mineta, float maxeta)
{
  m_mineta = std::min(mineta, maxeta);
  m_maxeta = std::max(mineta, maxeta);
}

void Afterburner::setPtRange(float minpt, float maxpt)
{
  m_minpt = std::min(minpt, maxpt);
  m_maxpt = std::max(minpt, maxpt);
}

void Afterburner::setPsiN(unsigned int n, float psi)
{
  if (n < 1 || n > 6)
  {
    std::cout << PHWHERE << ": Flow Afterburner set reaction plane angle psi_n for n=" << n << " which is out of range. Ignoring." << std::endl;
    return;
  }
  m_psi_n[n - 1] = psi;
}

float Afterburner::getPsiN(unsigned int n) const
{
  if (n < 1 || n > 6)
  {
    std::cout << PHWHERE << ": Flow Afterburner requested reaction plane angle psi_n for n=" << n << " which is out of range. Returning 0.0" << std::endl;
    return 0.0;
  }
  return m_psi_n[n - 1];

}

double Afterburner::vn_func(double x, void* params)
{
  const float* par = static_cast<const float*>(params);
  const double phi0 = par[0];
  const float* vn   = par + 1;
  const float* psi  = par + 7;
  double s = 0.0;
  for (int n = 1; n <= 6; ++n) {
    s += vn[n-1] * std::sin(static_cast<double>(n) * (x - psi[n-1])) / static_cast<double>(n);
  }
  return (x + 2.0 * s) - phi0;
}


void Afterburner::throw_psi_n(HepMC::GenEvent *event)
{
  HepMC::HeavyIon *hi = event->heavy_ion();
  if (!hi)
  {
    std::cout << PHWHERE << ": Flow Afterburner needs the Heavy Ion Event Info, GenEvent::heavy_ion() returns NULL" << std::endl;
    exit(1);
  }
  float psi_n[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // reaction plane angles
  for (int i = 0; i < 6; i++)
  {
    // Principal value must be within -PI/n to PI/n
    psi_n[i] = ((CLHEP::RandFlat::shoot(m_engine) - 0.5) * 2 * M_PI ) 
                /( static_cast<float>(i + 1) );
  }

  // The psi2 plane is aligned with the impact parameter
  psi_n[1] = hi->event_plane_angle();

  // Ensure that Psi2 is within [-PI/2,PI/2]
  psi_n[1] = std::atan2(std::sin(2.0 * psi_n[1]), std::cos(2.0 * psi_n[1])) / 2.0;


  for (unsigned int i = 1; i <= 6; ++i)
  {
    setPsiN(i, psi_n[i - 1]);
  }
  return;
}

void Afterburner::AddFlowToParentAndMoveDescendants(HepMC::GenEvent *event, HepMC::GenParticle *parent )
{
  CLHEP::HepLorentzVector momentum(parent->momentum().px(),
                                    parent->momentum().py(),
                                    parent->momentum().pz(),
                                    parent->momentum().e());
  double pt = momentum.perp();
  double eta = momentum.pseudoRapidity();
  double phi_0 = momentum.phi();

  HepMC::HeavyIon* hi = event->heavy_ion();
  if (!hi)
  {
    std::cout << PHWHERE << ": HeavyIon info missing in GenEvent. Cannot apply flow." << std::endl;
    std::exit(1);
  }

  m_algo->set_impact_parameter(hi->impact_parameter());
  m_algo->calc_flow(eta, pt, m_engine); // add engine for fluctuations (if enabled)

  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
  double x_lo = -2 * M_PI;
  double x_hi = 2 * M_PI;
  float params[13] = {};
  params[0]  = static_cast<float>(phi_0);
  for (int i = 0; i < 6; ++i) 
  {
    params[1 + i] = m_algo->get_vn(i+1);
    params[7 + i] = getPsiN(i+1);
  }

  gsl_function F;
  F.function = &Afterburner::vn_func;
  F.params = params;
  gsl_root_fsolver_set(s, &F, x_lo, x_hi);

  int status;
  int iter = 0;
  double phi = 0;
  do
  {
    ++iter;
    gsl_root_fsolver_iterate(s);
    phi  = gsl_root_fsolver_root(s);
    x_lo = gsl_root_fsolver_x_lower(s);
    x_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(x_lo, x_hi, 0, 1e-5);
  } while (status == GSL_CONTINUE && iter < 1000);

  if (iter >= 1000) {
    gsl_root_fsolver_free(s);
    m_phishift = 0.0;
    return; // do not rotate anything on failure
  }

  gsl_root_fsolver_free(s);

  m_phishift = phi - phi_0;
  if (fabs(m_phishift) > 1e-7)
  {
    momentum.rotateZ(m_phishift);  // DPM check units * Gaudi::Units::rad);
    parent->set_momentum(momentum);
  }

  HepMC::GenVertex *endvtx = parent->end_vertex();
  if (endvtx)
  {
    // now rotate descendant vertices
    for (HepMC::GenVertex::vertex_iterator descvtxit = endvtx->vertices_begin(HepMC::descendants);
          descvtxit != endvtx->vertices_end(HepMC::descendants);
          ++descvtxit)
    {
      HepMC::GenVertex *descvtx = (*descvtxit);

      // rotate vertex (magic number?)
      if (fabs(m_phishift) > 1e-7)
      {
        CLHEP::HepLorentzVector position(descvtx->position().x(),
                                          descvtx->position().y(),
                                          descvtx->position().z(),
                                          descvtx->position().t());
        position.rotateZ(m_phishift);  // DPM check units
        descvtx->set_position(position);
      }

      // now rotate their associated particles
      for (HepMC::GenVertex::particle_iterator descpartit = descvtx->particles_begin(HepMC::children);
            descpartit != descvtx->particles_end(HepMC::children);
            ++descpartit)
      {
        HepMC::GenParticle *descpart = (*descpartit);
        CLHEP::HepLorentzVector desmomentum(descpart->momentum().px(),
                                          descpart->momentum().py(),
                                          descpart->momentum().pz(),
                                          descpart->momentum().e());
        // Rotate particle
        if (fabs(m_phishift) > 1e-7)
        {
          desmomentum.rotateZ(m_phishift);  // DPM check units * Gaudi::Units::rad);
          descpart->set_momentum(desmomentum);
        }
      } // end of particles loop
    } // end of vertices loop
  } // end of if (endvtx)


  return ;
}

void Afterburner::readLegacyArguments(
    CLHEP::HepRandomEngine *engine,
    const std::string &algorithmName,
    float mineta, float maxeta,
    float minpt, float maxpt)
{
  if (engine)
  {
    setEngine(engine);
  }
  if (!algorithmName.empty())
  {
    setAlgo(algorithmName);
  }
  if (mineta != -5.0 || maxeta != 5.0)
  {
    setEtaRange(mineta, maxeta);
  }
  if (minpt != 0.0 || maxpt != 100.0)
  {
    setPtRange(minpt, maxpt);
  }

  return;
}

int Afterburner::flowAfterburner(HepMC::GenEvent *event,
                                  CLHEP::HepRandomEngine *engine, 
                                  const std::string &algorithmName,
                                  float mineta, float maxeta,
                                  float minpt, float maxpt )
{
  // check legacy arguments to see if they are set
  readLegacyArguments(engine, algorithmName, mineta, maxeta, minpt, maxpt);

  throw_psi_n(event); // throw reaction plane angles

  HepMC::GenVertex *mainvtx = event->barcode_to_vertex(-1);
  if (!mainvtx)
  {
    std::cout << PHWHERE << ": Flow Afterburner cannot find main vertex in GenEvent. Cannot apply flow." << std::endl;
    return -1;
  }

  // Loop over all children of this vertex
  HepMC::GenVertexParticleRange r(*mainvtx, HepMC::children);

  for (HepMC::GenVertex::particle_iterator it = r.begin(); it != r.end(); it++)
  {
    // Process particles from main vertex
    HepMC::GenParticle *parent = (*it);

    CLHEP::HepLorentzVector momentum(parent->momentum().px(),
                                     parent->momentum().py(),
                                     parent->momentum().pz(),
                                     parent->momentum().e());

    float eta = momentum.pseudoRapidity();
    if (eta < m_mineta || eta > m_maxeta)
    {
      continue;
    }

    float pT = momentum.perp();
    if (pT < m_minpt || pT > m_maxpt)
    {
      continue;
    }

    AddFlowToParentAndMoveDescendants(event, parent);
  } // end of particles loop

  return 0;

}

// Legacy function, use the Afterburner class instead
int flowAfterburner(HepMC::GenEvent *event,
                    CLHEP::HepRandomEngine *engine,
                    const std::string &algorithmName,
                    float mineta, float maxeta,
                    float minpt, float maxpt )
{
  Afterburner afterburner(algorithmName, engine, mineta, maxeta, minpt, maxpt);
  return afterburner.flowAfterburner(event, engine, algorithmName, mineta, maxeta, minpt, maxpt);
}


