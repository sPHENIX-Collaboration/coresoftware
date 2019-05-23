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
  Copy(g4hit);
}

void
PHG4HitEval::Copy(PHG4Hit const &g4hit)
{
  PHG4Hitv1::Copy(g4hit);

  // fill the branched variables from the property arrays
  eion = PHG4Hitv1::get_eion();
  scint_id = PHG4Hitv1::get_scint_id();
  light_yield = PHG4Hitv1::get_light_yield();
  path_length = PHG4Hitv1::get_path_length();

}
