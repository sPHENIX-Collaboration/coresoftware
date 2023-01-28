#include "PHEvtGenRandomEngine.hh"

#include <phool/PHRandomSeed.h>

#include <gsl/gsl_rng.h>  // for gsl_rng_alloc, gsl_rng_free, gsl_rng...

PHEvtGenRandomEngine::PHEvtGenRandomEngine()
{
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  m_Seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  gsl_rng_set(m_RandomGenerator, m_Seed);
}

PHEvtGenRandomEngine::~PHEvtGenRandomEngine()
{
  gsl_rng_free(RandomGenerator());
}

double PHEvtGenRandomEngine::random()
{
  return gsl_rng_uniform_pos(RandomGenerator());
}
