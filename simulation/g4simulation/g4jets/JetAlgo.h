#ifndef G4JET_JETALGO_H
#define G4JET_JETALGO_H

#include "Jet.h"

#include <cmath>

class JetAlgo
{
 public:
  virtual ~JetAlgo() {}

  virtual void identify(std::ostream& os = std::cout)
  {
    os << "JetAlgo base class" << std::endl;
  }

  virtual Jet::ALGO get_algo() { return Jet::NONE; }
  virtual float get_par() { return NAN; }

  virtual std::vector<Jet*> get_jets(std::vector<Jet*> /* particles*/)
  {
    return std::vector<Jet*>();
  }

 protected:
  JetAlgo() {}

 private:
};

#endif
