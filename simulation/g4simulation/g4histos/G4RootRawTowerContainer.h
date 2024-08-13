#ifndef G4HISTOS_G4ROOTRAWTOWERCONTAINER_H
#define G4HISTOS_G4ROOTRAWTOWERCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>  // for cout, ostream
#include <limits>

class G4RootRawTower;
class TClonesArray;

class G4RootRawTowerContainer : public PHObject
{
 public:
  G4RootRawTowerContainer();
  ~G4RootRawTowerContainer() override;

  // from PHObject
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;

  G4RootRawTower* AddG4RootRawTower(const G4RootRawTower& g4tower);
  void set_etotal(const float e) { etotal = e; }
  float get_etotal() const { return etotal; }

  void set_event(const int i) { event = i; }
  int get_event() const { return event; }

 protected:
  float etotal{std::numeric_limits<float>::quiet_NaN()};
  int event{0};
  TClonesArray* SnglG4RootRawTowers{nullptr};

  ClassDefOverride(G4RootRawTowerContainer, 1)
};

#endif
