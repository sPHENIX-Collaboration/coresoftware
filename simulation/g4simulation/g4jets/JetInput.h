#ifndef __JETINPUT_H__
#define __JETINPUT_H__

#include <phool/PHCompositeNode.h>
#include <fastjet/PseudoJet.hh>
#include <vector>

class JetInput {
  
public:

  virtual ~JetInput() {}

  virtual std::vector<fastjet::PseudoJet> get_input(PHCompositeNode *topNode) {
    return std::vector<fastjet::PseudoJet>();
  }


protected:
  JetInput();
  
private:
    
};

#endif
