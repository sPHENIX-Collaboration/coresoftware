#ifndef G4JET_JETINPUT_H
#define G4JET_JETINPUT_H

#include "Jet.h"

#include <iostream>
#include <vector>

class PHCompositeNode;

class JetInput
{
 public:
  virtual ~JetInput() {}

  virtual void identify(std::ostream& os = std::cout)
  {
    os << "JetInput base class" << std::endl;
  }

  virtual Jet::SRC get_src() { return Jet::VOID; }

  virtual std::vector<Jet*> get_input(PHCompositeNode* topNode)
  {
    return std::vector<Jet*>();
  }

 protected:
  JetInput() {}

 private:
};

#endif
