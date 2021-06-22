#ifndef G4HISTOS_G4ROOTHITCONTAINER_H
#define G4HISTOS_G4ROOTHITCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream

class PHG4Hit;
class TClonesArray;

class G4RootHitContainer : public PHObject
{
 public:
  G4RootHitContainer();
  ~G4RootHitContainer() override;

// from PHObject
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;

  PHG4Hit* AddHit(const PHG4Hit *g4hit);
  void set_etotal(const float e) { etotal = e; }
  float get_etotal() const { return etotal; }

  void set_eion(const float e) { eion = e; }
  float get_eion() const { return eion; }

  void set_leakage(const float f) { leakage = f; }
  float get_leakage() const { return leakage; }

  void set_event(const int i) { event = i; }
  int get_event() const { return event; }


 protected:
  float etotal;
  float eion;
  float leakage;
  int event;
  TClonesArray* SnglHits;

  ClassDefOverride(G4RootHitContainer, 1)
};

#endif
