// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4HISTOS_G4ROOTSCINTILLATORTOWER_H
#define G4HISTOS_G4ROOTSCINTILLATORTOWER_H

#include <phool/PHObject.h>

class RawTower;

class G4RootScintillatorTower : public PHObject
{
 public:
  G4RootScintillatorTower();
  G4RootScintillatorTower(const RawTower& tower);
  virtual ~G4RootScintillatorTower() {}
  void Reset();
  int isValid() const;
  void identify(std::ostream& os = std::cout) const;

  int get_row() const { return row; }
  int get_column() const { return column; }
  double get_energy() const { return energy; }
 protected:
  short row;
  short column;
  double energy;

  ClassDef(G4RootScintillatorTower, 1)
};

#endif  // G4HISTOS_G4ROOTSCINTILLATORTOWER_H
