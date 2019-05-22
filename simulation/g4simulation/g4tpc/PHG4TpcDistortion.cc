// $Id: $

/*!
 * \file PHG4TpcDistortion.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4TpcDistortion.h"

#include <phool/PHRandomSeed.h>

PHG4TpcDistortion::PHG4TpcDistortion(int verbose)
  : verbosity(verbose)
{
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  gsl_rng_set(RandomGenerator, seed);
}

PHG4TpcDistortion::~PHG4TpcDistortion()
{
  gsl_rng_free(RandomGenerator);
}
