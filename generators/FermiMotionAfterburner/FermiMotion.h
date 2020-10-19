#ifndef FERMIMOTION_FERMIMOTION_H
#define FERMIMOTION_FERMIMOTION_H

#include <string>



namespace CLHEP
{
  class HepRandomEngine;
}
namespace HepMC { 
class GenEvent; 
}




  int FermiMotion (HepMC::GenEvent *event, CLHEP::HepRandomEngine *engine);


#endif
