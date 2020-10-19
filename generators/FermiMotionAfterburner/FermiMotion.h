#ifndef FERMIMOTION_FERMIMOTION_H
#define FERMIMOTION_FERMIMOTION_H

#include <string>
#include <gsl/gsl_rng.h>


namespace CLHEP
{
  class HepRandomEngine;
}
namespace HepMC { 
class GenEvent; 
}




  int FermiMotion (HepMC::GenEvent *event, gsl_rng *RandomGenerator);


#endif
