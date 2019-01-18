#ifndef __JETINPUT_H__
#define __JETINPUT_H__

#include "Jet.h"
#include <vector>
#include <iostream>

class PHCompositeNode;

class JetInput {
  
public:

  virtual ~JetInput() {}

  virtual void identify(std::ostream& os = std::cout) {
    os << "JetInput base class" << std::endl;
  }
  
  virtual Jet::SRC get_src() {return Jet::VOID;}
  
  virtual std::vector<Jet*> get_input(PHCompositeNode *topNode) {
    return std::vector<Jet*>();
  }


protected:
  JetInput() {}
  
private:
    
};

#endif
