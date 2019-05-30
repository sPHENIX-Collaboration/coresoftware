#ifndef G4HISTOS_G4ROOTRAWTOWERCONTAINER_H
#define G4HISTOS_G4ROOTRAWTOWERCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream

class G4RootRawTower;
class TClonesArray;

class G4RootRawTowerContainer : public PHObject
{
 public:
  G4RootRawTowerContainer();
  virtual ~G4RootRawTowerContainer();

  void Reset();

  G4RootRawTower* AddG4RootRawTower(const G4RootRawTower& g4tower);
  void set_etotal(const float e) { etotal = e; }
  float get_etotal() const { return etotal; }

  void set_event(const int i) { event = i; }
  int get_event() const { return event; }

  void identify(std::ostream& os = std::cout) const;

 protected:
  float etotal;
  int event;
  TClonesArray* SnglG4RootRawTowers;

  ClassDef(G4RootRawTowerContainer, 1)
};

#endif
