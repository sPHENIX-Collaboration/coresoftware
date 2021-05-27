#ifndef G4HISTOS_G4ROOTRAWTOWER_H
#define G4HISTOS_G4ROOTRAWTOWER_H

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream

class G4RootRawTower : public PHObject
{
 public:
  G4RootRawTower();
  G4RootRawTower(const float ieta, const float iphi, const float e);
  ~G4RootRawTower() override {}

  void Reset() override;
  int isValid() const override;
  void identify(std::ostream& os = std::cout) const override;

  float get_eta() const { return eta; }
  float get_phi() const { return phi; }
  float get_energy() const { return energy; }

 protected:
  float eta;
  float phi;
  float energy;

  ClassDefOverride(G4RootRawTower, 1)
};

#endif /* G4HISTOS_RAWTOWERV1_H */
