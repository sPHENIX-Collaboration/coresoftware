#ifndef G4ROOTRAWTOWER_H_
#define G4ROOTRAWTOWER_H_

#include "phool/PHObject.h"

#include <map>

class G4RootRawTower : public PHObject {

 public:
  G4RootRawTower();
  G4RootRawTower(const float ieta, const float iphi, const float e);
  virtual ~G4RootRawTower() {}

  void Reset();
  int isValid() const;
  void identify(std::ostream& os=std::cout) const;

  float get_eta() const { return eta; }
  float get_phi() const { return phi; }
  float get_energy() const {return energy;}

 protected:
  float eta;
  float phi;
  float energy;

  ClassDef(G4RootRawTower,1)
};
 
#endif /* RAWTOWERV1_H_ */
