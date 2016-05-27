// $Id: $                                                                                             

/*!
 * \file PHG4TPCDistortion.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4TPCDistortion.h"

#include <phool/PHRandomSeed.h>

PHG4TPCDistortion::PHG4TPCDistortion(int verbose) :
    verbosity(verbose)
{
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned int seed = PHRandomSeed(); // fixed seed is handled in this funtcion
  gsl_rng_set(RandomGenerator, seed);
}

PHG4TPCDistortion::~PHG4TPCDistortion()
{
  gsl_rng_free(RandomGenerator);
}

