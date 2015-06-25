// $Id: $                                                                                             

/*!
 * \file PHG4HitEval.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4HitEval.h"

PHG4HitEval::PHG4HitEval() :
    eion(NAN), scint_id(-9999), light_yield(NAN), path_length(NAN)

{
  // TODO Auto-generated constructor stub

}

PHG4HitEval::~PHG4HitEval()
{
  // TODO Auto-generated destructor stub
}

PHG4HitEval::PHG4HitEval(PHG4Hit const &g4hit)
{
  Copy(g4hit);
}
