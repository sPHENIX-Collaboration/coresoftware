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

#include <cmath>
#include <cstdlib>  // for exit
#include <iostream>
#include <map>  // for map

namespace CLHEP
{
  class HepRandomEngine;
}

flowAfterburnerAlgorithm algorithm;
std::map<std::string, flowAfterburnerAlgorithm> algorithms;

struct loaderObj
{
  loaderObj()
  {
    static bool init = false;
    if (!init)
    {
      algorithms["MINBIAS"] = minbias_algorithm;
      algorithms["MINBIAS_V2_ONLY"] = minbias_v2_algorithm;
      algorithms["CUSTOM"] = custom_algorithm;
    }
  }
};

loaderObj loader;

double
vn_func(double x, void *params)
{
  float *par_float = (float *) params;
  float phi_0 = par_float[0];
  float *vn = par_float + 1;
  float *psi_n = vn + 6;
  double val =
      x + 2 * (vn[0] * sin(1 * (x - psi_n[0])) / 1.0 +
               vn[1] * sin(2 * (x - psi_n[1])) / 2.0 +
               vn[2] * sin(3 * (x - psi_n[2])) / 3.0 +
               vn[3] * sin(4 * (x - psi_n[3])) / 4.0 +
               vn[4] * sin(5 * (x - psi_n[4])) / 5.0 +
               vn[5] * sin(6 * (x - psi_n[5])) / 6.0);
  return val - phi_0;
}

double
vn_func_derivative(double x, void *params)
{
  float *par_float = (float *) params;
  float *vn = par_float + 1;
  float *psi_n = vn + 6;
  double val =
      1 + 2 * (vn[0] * cos(1 * (x - psi_n[0])) / 1.0 +
               vn[1] * cos(2 * (x - psi_n[1])) / 2.0 +
               vn[2] * cos(3 * (x - psi_n[2])) / 3.0 +
               vn[3] * cos(4 * (x - psi_n[3])) / 4.0 +
               vn[4] * cos(5 * (x - psi_n[4])) / 5.0 +
               vn[5] * cos(6 * (x - psi_n[5])) / 6.0);
  return val;
}

float psi_n[6], v1, v2, v3, v4, v5, v6;

void MoveDescendantsToParent(HepMC::GenParticle *parent,
                             double phishift)
{
  // Move the branch of descendant vertices and particles by phishift
  // to parent particle position
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
      if (fabs(phishift) > 1e-7)
      {
        CLHEP::HepLorentzVector position(descvtx->position().x(),
                                         descvtx->position().y(),
                                         descvtx->position().z(),
                                         descvtx->position().t());
        position.rotateZ(phishift);  // DPM check units
        descvtx->set_position(position);
      }

      // now rotate their associated particles
      for (HepMC::GenVertex::particle_iterator descpartit = descvtx->particles_begin(HepMC::children);
           descpartit != descvtx->particles_end(HepMC::children);
           ++descpartit)
      {
        HepMC::GenParticle *descpart = (*descpartit);
        CLHEP::HepLorentzVector momentum(descpart->momentum().px(),
                                         descpart->momentum().py(),
                                         descpart->momentum().pz(),
                                         descpart->momentum().e());
        // Rotate particle
        if (fabs(phishift) > 1e-7)
        {
          momentum.rotateZ(phishift);  // DPM check units * Gaudi::Units::rad);
          descpart->set_momentum(momentum);
        }
      }
    }
  }

  return;
}

float calc_v2(double b, double eta, double pt)
{
  float a1, a2, a3, a4;
  a1 = 0.4397 * exp(-(b - 4.526) * (b - 4.526) / 72.0) + 0.636;
  a2 = 1.916 / (b + 2) + 0.1;
  a3 = 4.79 * 0.0001 * (b - 0.621) * (b - 10.172) * (b - 23) + 1.2;  // this is >0 for b>0
  a4 = 0.135 * exp(-0.5 * (b - 10.855) * (b - 10.855) / 4.607 / 4.607) + 0.0120;

  float temp1 = pow(pt, a1) / (1 + exp((pt - 3.0) / a3));
  float temp2 = pow(pt + 0.1, -a2) / (1 + exp(-(pt - 4.5) / a3));
  float temp3 = 0.01 / (1 + exp(-(pt - 4.5) / a3));

  //v2 = (a4 * (temp1 + temp2) + temp3) * exp (-0.5 * eta * eta / 6.27 / 6.27);

  // Adjust flow rapidity dependence to better match PHOBOS 200 GeV Au+Au data
  // JGL 9/9/2019
  // See JS ToG talk at https://indico.bnl.gov/event/6764/

  v2 = (a4 * (temp1 + temp2) + temp3) * exp(-0.5 * eta * eta / 3.43 / 3.43);

  return v2;
}

// New parameterization for vn
void jjia_minbias_new(double b, double eta, double pt)
{
  v2 = calc_v2(b, eta, pt);

  float fb = 0.97 + 1.06 * exp(-0.5 * b * b / 3.2 / 3.2);
  v3 = pow(fb * sqrt(v2), 3);

  float gb = 1.096 + 1.36 * exp(-0.5 * b * b / 3.0 / 3.0);
  gb = gb * sqrt(v2);
  v4 = pow(gb, 4);
  v5 = pow(gb, 5);
  v6 = pow(gb, 6);
  v1 = 0;
}

// New parameterization for v2
void jjia_minbias_new_v2only(double b, double eta, double pt)
{
  v2 = calc_v2(b, eta, pt);

  v1 = 0;
  v3 = 0;
  v4 = 0;
  v5 = 0;
  v6 = 0;
}

// Custom vn
void custom_vn(double /*b*/, double /*eta*/, double /*pt*/)
{
  v1 = 0.0000;
  v2 = 0.0500;
  v3 = 0.0280;
  v4 = 0.0130;
  v5 = 0.0045;
  v6 = 0.0015;
}

double
AddFlowToParent(HepMC::GenEvent *event, HepMC::GenParticle *parent)
{
  CLHEP::HepLorentzVector momentum(parent->momentum().px(),
                                   parent->momentum().py(),
                                   parent->momentum().pz(),
                                   parent->momentum().e());
  double pt = momentum.perp();
  double eta = momentum.pseudoRapidity();
  double phi_0 = momentum.phi();

  HepMC::HeavyIon *hi = event->heavy_ion();
  double b = hi->impact_parameter();

  v1 = 0, v2 = 0, v3 = 0, v4 = 0, v5 = 0, v6 = 0;

  //Call the appropriate function to set the vn values
  if (algorithm == minbias_algorithm)
  {
    jjia_minbias_new(b, eta, pt);
  }
  else if (algorithm == minbias_v2_algorithm)
  {
    jjia_minbias_new_v2only(b, eta, pt);
  }
  else if (algorithm == custom_algorithm)
  {
    custom_vn(b, eta, pt);
  }

  double phishift = 0;

  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
  double x_lo = -2 * M_PI, x_hi = 2 * M_PI;
  float params[13];
  for (int ipar = 0; ipar < 13; ipar++)
  {
    params[ipar] = 0;
  }
  gsl_function F;
  F.function = &vn_func;
  F.params = &params;
  gsl_root_fsolver_set(s, &F, x_lo, x_hi);
  int iter = 0;
  params[0] = phi_0;
  params[1] = v1;
  params[2] = v2;
  params[3] = v3;
  params[4] = v4;
  params[5] = v5;
  params[6] = v6;
  params[7] = psi_n[0];
  params[8] = psi_n[1];
  params[9] = psi_n[2];
  params[10] = psi_n[3];
  params[11] = psi_n[4];
  params[12] = psi_n[5];
  int status;
  double phi;
  do
  {
    iter++;
    gsl_root_fsolver_iterate(s);
    phi = gsl_root_fsolver_root(s);
    x_lo = gsl_root_fsolver_x_lower(s);
    x_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(x_lo, x_hi, 0, 0.00001);
  } while (status == GSL_CONTINUE && iter < 1000);
  gsl_root_fsolver_free(s);

  if (iter >= 1000)
    return 0;

  phishift = phi - phi_0;

  if (fabs(phishift) > 1e-7)
  {
    momentum.rotateZ(phishift);  // DPM check units * Gaudi::Units::rad);
    parent->set_momentum(momentum);
  }

  return phishift;
}

int flowAfterburner(HepMC::GenEvent *event,
                    CLHEP::HepRandomEngine *engine,
                    std::string algorithmName,
                    float mineta, float maxeta,
                    float minpt, float maxpt)
{
  algorithm = algorithms[algorithmName];
  HepMC::HeavyIon *hi = event->heavy_ion();
  if (!hi)
  {
    std::cout << PHWHERE << ": Flow Afterburner needs the Heavy Ion Event Info, GenEvent::heavy_ion() returns NULL" << std::endl;
    exit(1);
  }
  // Generate the v_n reaction plane angles (some of them may or may
  // not be used later on).
  for (int i = 0; i < 6; i++)
  {
    // Principal value must be within -PI/n to PI/n
    psi_n[i] = (CLHEP::RandFlat::shoot(engine) - 0.5) * 2 * M_PI / (i + 1);
  }

  // The psi2 plane is aligned with the impact parameter
  psi_n[1] = hi->event_plane_angle();

  // Ensure that Psi2 is within [-PI/2,PI/2]
  psi_n[1] = atan2(sin(2 * psi_n[1]), cos(2 * psi_n[1])) / 2.0;

  HepMC::GenVertex *mainvtx = event->barcode_to_vertex(-1);

  // Loop over all children of this vertex
  HepMC::GenVertexParticleRange r(*mainvtx, HepMC::children);

  for (HepMC::GenVertex::particle_iterator it = r.begin(); it != r.end(); it++)
  {
    // Process particles from main vertex
    HepMC::GenParticle *parent = (*it);

    // Skip the "jets" found during the Hijing run itself
    // if (parent->status() == 103)
    // 	{
    // 	  continue;
    // 	}
    CLHEP::HepLorentzVector momentum(parent->momentum().px(),
                                     parent->momentum().py(),
                                     parent->momentum().pz(),
                                     parent->momentum().e());

    float eta = momentum.pseudoRapidity();
    if (eta < mineta || eta > maxeta)
    {
      continue;
    }

    // Skip particle if pT is outside implementation range
    float pT = momentum.perp();
    if (pT < minpt || pT > maxpt)
    {
      continue;
    }

    // Add flow to particles from main vertex
    double phishift = AddFlowToParent(event, parent);
    MoveDescendantsToParent(parent, phishift);
  }

  return 0;
}
