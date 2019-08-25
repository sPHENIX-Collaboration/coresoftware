// $Id: $                                                                                             

/*!
 * \file PHG4HitEval.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4HitEval.h"

#include <cmath>

class PHG4Hit;

PHG4HitEval::PHG4HitEval() :
    eion(NAN), scint_id(-9999), light_yield(NAN), path_length(NAN)

{
  // TODO Auto-generated constructor stub

}

PHG4HitEval::PHG4HitEval(PHG4Hit const &g4hit)
{
  CopyFrom(g4hit);
}

void
PHG4HitEval::CopyFrom(PHObject const &g4hit)
{
  PHG4Hit::CopyFrom(g4hit);

  // fill the branched variables from the property arrays

  const PHG4Hit *g4h = dynamic_cast<const PHG4Hit *> (&g4hit);
  eion = g4h->get_eion();
  scint_id = g4h->get_scint_id();
  light_yield = g4h->get_light_yield();
  path_length = g4h->get_path_length();

}
