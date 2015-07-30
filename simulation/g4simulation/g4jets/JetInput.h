#ifndef __JETINPUT_H__
#define __JETINPUT_H__

#include <phool/PHCompositeNode.h>
#include <TLorentzVector.h>

class JetInput {
  
public:

  virtual ~JetInput() {}

  virtual std::vector<TLorentzVector> get_input(PHCompositeNode *topNode) {
    return std::vector<TLorentzVector>();
  }

protected:
  JetInput();
  
private:
    
  ClassDef(JetInput, 1);
};

#endif
