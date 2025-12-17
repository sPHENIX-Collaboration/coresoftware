// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4HISTOS_G4ROOTSCINTILLATORTOWER_H
#define G4HISTOS_G4ROOTSCINTILLATORTOWER_H

#include <phool/PHObject.h>

#include <iostream>  // for cout, ostream

/* class TowerInfo; */

class G4RootScintillatorTower : public PHObject
{
 public:
  G4RootScintillatorTower() = default;
  G4RootScintillatorTower(double towerenergy, int ieta, int iphi);
  ~G4RootScintillatorTower() override = default;

  void Reset() override;
  int isValid() const override;
  void identify(std::ostream& os = std::cout) const override;

  int get_row() const { return row; }
  int get_column() const { return column; }
  double get_energy() const { return energy; }

 protected:
  short row{-1};
  short column{-1};
  double energy{-1};

  ClassDefOverride(G4RootScintillatorTower, 1)
};

#endif  // G4HISTOS_G4ROOTSCINTILLATORTOWER_H
