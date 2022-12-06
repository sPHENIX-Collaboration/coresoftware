#include "PHEvtGenRandomEngine.hh"
#include <phool/PHRandomSeed.h>
//#include <g4main/PHG4ParticleGeneratorBase.h>

#include <iostream>

PHEvtGenRandomEngine::PHEvtGenRandomEngine()
{
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  m_Seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  gsl_rng_set(m_RandomGenerator, m_Seed);
}

PHEvtGenRandomEngine::~PHEvtGenRandomEngine()
{
  // std::cout << "Detroy the Random Generator" << std::endl;
  gsl_rng_free(RandomGenerator());
  // delete m_RandomGenerator;
}

double PHEvtGenRandomEngine::random()
{
  return gsl_rng_uniform_pos(RandomGenerator());
}
