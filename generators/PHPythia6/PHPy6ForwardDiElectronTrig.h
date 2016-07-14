#ifndef PHPy6ForwardDiElectronTrig_h
#define PHPy6ForwardDiElectronTrig_h

#include <PHPy6GenTrigger.h>
#include <HepMC/GenEvent.h>

namespace HepMC
{
  class GenEvent;
};

class PHPy6ForwardDiElectronTrig: public PHPy6GenTrigger
{
  public:
  
  //! constructor
  PHPy6ForwardDiElectronTrig(const std::string &name = "PHPy6ForwardDiElectronTrigger");
  
  //! destructor 
  ~PHPy6ForwardDiElectronTrig( void ){}
 
  #ifndef __CINT__ 
  bool Apply(const HepMC::GenEvent* evt);
  #endif
 
  protected:
	
  int ntriggered_forward_dielectron;
  int nconsidered_forward_dielectron;
  
};

#endif	
