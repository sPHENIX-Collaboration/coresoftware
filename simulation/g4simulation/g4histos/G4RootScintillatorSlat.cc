#include "G4RootScintillatorSlat.h"

#include <g4detectors/PHG4ScintillatorSlat.h>

#include <iostream>

using namespace std;

G4RootScintillatorSlat::G4RootScintillatorSlat()
  : row(-1)
  , column(-1)
  , edep(0)
  , eion(0)
  , light_yield(0)
{
}

G4RootScintillatorSlat::G4RootScintillatorSlat(const PHG4ScintillatorSlat& slat)
  : row(slat.get_row())
  , column(slat.get_column())
  , edep(slat.get_edep())
  , eion(slat.get_eion())
  , light_yield(slat.get_light_yield())
{
}

void G4RootScintillatorSlat::Reset()
{
  row = -1;
  column = -1;
  eion = 0;
  edep = 0;
  light_yield = 0;
}

int G4RootScintillatorSlat::isValid() const
{
  return (row >= 0);
}

void G4RootScintillatorSlat::identify(std::ostream& os) const
{
  os << "G4RootScintillatorSlat: row: " << row << ", column: " << column
     << " energy=" << edep << std::endl;
}
