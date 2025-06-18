#ifndef JETBASE_JETALGO_H
#define JETBASE_JETALGO_H

#include "Jet.h"

#include <limits>

class JetContainer;
class JetAlgo
{
 public:
  virtual ~JetAlgo() = default;

  virtual void identify(std::ostream& os = std::cout)
  {
    os << "JetAlgo base class" << std::endl;
  }

  virtual Jet::ALGO get_algo() { return Jet::NONE; }
  virtual float get_par() { return std::numeric_limits<float>::quiet_NaN(); }

  // old version -- get jets to fill into JetMap
  virtual std::vector<Jet*> get_jets(std::vector<Jet*> /* particles*/)  // ? Why isn't this passed as a reference?
  {
    return std::vector<Jet*>();
  }

  // new version -- pass JetContainer into clusterFillJets to fill it
  virtual void cluster_and_fill(std::vector<Jet*>& /* particles*/, JetContainer* /*clones*/)
  {
  }

  virtual std::map<Jet::PROPERTY, unsigned int>& property_indices();

 protected:
  JetAlgo() = default;
};

#endif
