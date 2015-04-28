#include "boost/foreach.hpp"

#include "HepMC/GenVertex.h"
#include "HepMC/GenRanges.h"

int
main ()
{
  HepMC::GenVertexParticleRange r(HepMC::GenVertex *vertex, HepMC::children);

  for (HepMC::GenVertex::particle_iterator it = r.begin (); it != r.end (); it++);

  BOOST_FOREACH(HepMC::GenVertex::particle_iterator it, r);

  return 0;
}

  
