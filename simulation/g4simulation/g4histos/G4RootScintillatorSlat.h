#ifndef G4HISTOS_G4ROOTSCINTILLATORSLAT_H
#define G4HISTOS_G4ROOTSCINTILLATORSLAT_H

#include <phool/PHObject.h>

#include <iostream>          // for cout, ostream

class PHG4ScintillatorSlat;

class G4RootScintillatorSlat : public PHObject
{
 public:
  G4RootScintillatorSlat();
  G4RootScintillatorSlat(const PHG4ScintillatorSlat& slat);
  ~G4RootScintillatorSlat() override {}

  void Reset() override;
  int isValid() const override;
  void identify(std::ostream& os = std::cout) const override;

  int get_row() const { return row; }
  int get_column() const { return column; }

  double get_edep() const { return edep; }
  double get_eion() const { return eion; }
  double get_light_yield() const { return light_yield; }

 protected:
  short row;
  short column;
  double edep;
  double eion;
  double light_yield;

  ClassDefOverride(G4RootScintillatorSlat, 1)
};

#endif /* G4HISTOS_G4ROOTSCINTILLATORSLAT_H */
