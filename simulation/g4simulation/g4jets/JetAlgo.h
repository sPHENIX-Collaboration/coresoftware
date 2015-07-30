#ifndef __JETALGO_H__
#define __JETALGO_H__

#include <phool/PHCompositeNode.h>
#include <fastjet/PseudoJet.hh>

class JetAlgo {
  
public:

  virtual ~JetAlgo() {}

  virtual std::vector<fastjet::PseudoJet> get_jets(const std::vector<fastjet::PseudoJet>& particles) {
    return std::vector<fastjet::PseudoJet>();
  }

protected:
  JetAlgo();
  
private:
    
};

#endif
