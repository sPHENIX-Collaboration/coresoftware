#ifndef PHPy6ForwardElectronTrig_h
#define PHPy6ForwardElectronTrig_h

#include <PHPy6GenTrigger.h>
#include <HepMC/GenEvent.h>

namespace HepMC
{
  class GenEvent;
};

class PHPy6ForwardElectronTrig: public PHPy6GenTrigger
{
  public:
  
  //! constructor
  PHPy6ForwardElectronTrig(const std::string &name = "PHPy6ForwardElectronTrigger");
  
  //! destructor 
  ~PHPy6ForwardElectronTrig( void ){}
 
  #ifndef __CINT__ 
  bool Apply(const HepMC::GenEvent* evt);
  #endif
 
  protected:
	
  int ntriggered_forward_electron;
  int nconsidered_forward_electron;
  
};

#endif	
