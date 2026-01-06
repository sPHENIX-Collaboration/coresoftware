#ifndef G4HISTOS_G4ROOTSCINTILLATORSLAT_H
#define G4HISTOS_G4ROOTSCINTILLATORSLAT_H

#include <phool/PHObject.h>

#include <iostream>  // for cout, ostream

class PHG4ScintillatorSlat;

class G4RootScintillatorSlat : public PHObject
{
 public:
  G4RootScintillatorSlat() = default;
  G4RootScintillatorSlat(const PHG4ScintillatorSlat& slat);
  ~G4RootScintillatorSlat() override = default;

  void Reset() override;
  int isValid() const override;
  void identify(std::ostream& os = std::cout) const override;

  int get_row() const { return row; }
  int get_column() const { return column; }

  double get_edep() const { return edep; }
  double get_eion() const { return eion; }
  double get_light_yield() const { return light_yield; }

 protected:
  short row{-1};
  short column{-1};
  double edep{0};
  double eion{0};
  double light_yield{0};

  ClassDefOverride(G4RootScintillatorSlat, 1)
};

#endif /* G4HISTOS_G4ROOTSCINTILLATORSLAT_H */
