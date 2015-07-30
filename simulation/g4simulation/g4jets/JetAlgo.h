#ifndef __JETALGO_H__
#define __JETALGO_H__

#include <phool/PHCompositeNode.h>
#include <TLorentzVector.h>

class JetAlgo {
  
public:

  virtual ~JetAlgo() {}

  virtual std::vector<TLorentzVector> get_jets(const std::vector<TLorentzVector>& particles) {
    return std::vector<TLorentzVector>();
  }

protected:
  JetAlgo();
  
private:
    
  ClassDef(JetAlgo, 1);
};

#endif
