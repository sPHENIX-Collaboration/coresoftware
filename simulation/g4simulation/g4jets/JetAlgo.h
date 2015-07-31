#ifndef __JETALGO_H__
#define __JETALGO_H__

#include <phool/PHCompositeNode.h>
#include <Jet.h>

class JetAlgo {
  
public:

  virtual ~JetAlgo() {}

  virtual std::vector<Jet> get_jets(const std::vector<Jet>& particles) {
    return std::vector<Jet>();
  }

protected:
  JetAlgo() {}
  
private:
    
};

#endif
