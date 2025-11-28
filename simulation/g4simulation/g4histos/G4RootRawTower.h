#ifndef G4HISTOS_G4ROOTRAWTOWER_H
#define G4HISTOS_G4ROOTRAWTOWER_H

#include <phool/PHObject.h>

#include <iostream>  // for cout, ostream
#include <limits>

class G4RootRawTower : public PHObject
{
 public:
  G4RootRawTower() = default;
  G4RootRawTower(const float ieta, const float iphi, const float e);
  ~G4RootRawTower() override = default;

  void Reset() override;
  int isValid() const override;
  void identify(std::ostream& os = std::cout) const override;

  float get_eta() const { return eta; }
  float get_phi() const { return phi; }
  float get_energy() const { return energy; }

 protected:
  float eta{std::numeric_limits<float>::quiet_NaN()};
  float phi{std::numeric_limits<float>::quiet_NaN()};
  float energy{std::numeric_limits<float>::quiet_NaN()};

  ClassDefOverride(G4RootRawTower, 1)
};

#endif /* G4HISTOS_RAWTOWERV1_H */
