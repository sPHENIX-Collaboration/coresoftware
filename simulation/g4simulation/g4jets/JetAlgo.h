#ifndef __JETALGO_H__
#define __JETALGO_H__

#include <phool/PHCompositeNode.h>
#include "Jet.h"
#include <cmath>

class JetAlgo {
  
public:

  virtual ~JetAlgo() {}

  virtual Jet::ALGO get_algo() {return Jet::NONE;}
  virtual float get_par() {return NAN;}
  
  virtual std::vector<Jet*> get_jets(std::vector<Jet*> particles) {
    return std::vector<Jet*>();
  }

protected:
  JetAlgo() {}
  
private:
    
};

#endif
