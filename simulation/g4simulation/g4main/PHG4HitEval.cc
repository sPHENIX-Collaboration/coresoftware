// $Id: $

/*!
 * \file PHG4HitEval.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4HitEval.h"

#include "PHG4Hit.h"  // for PHG4Hit

#include <phool/PHObject.h>  // for PHObject

#include <cassert>
#include <cmath>

PHG4HitEval::PHG4HitEval()
  : eion(NAN)
  , scint_id(-9999)
  , light_yield(NAN)
  , path_length(NAN)

{
  // TODO Auto-generated constructor stub
}

PHG4HitEval::PHG4HitEval(const PHG4Hit *g4hit)
  : eion(NAN)
  , scint_id(-9999)
  , light_yield(NAN)
  , path_length(NAN)
{
  CopyFrom(g4hit);
}

void PHG4HitEval::CopyFrom(const PHObject *phobj)
{
  const PHG4Hit *g4hit = dynamic_cast<const PHG4Hit *>(phobj);
  assert(g4hit);

  PHG4Hit::CopyFrom(phobj);

  // fill the branched variables from the property arrays

  eion = g4hit->get_eion();
  scint_id = g4hit->get_scint_id();
  light_yield = g4hit->get_light_yield();
  path_length = g4hit->get_path_length();
}
