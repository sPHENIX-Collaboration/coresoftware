#include "PHG4ScintillatorSlatv1.h"

using namespace std;

PHG4ScintillatorSlatv1::PHG4ScintillatorSlatv1():
  row(-1),
  column(-1),
  edep(0),
  eion(0),
  light_yield(0)
{}

unsigned int
PHG4ScintillatorSlatv1::get_key() const
{
  unsigned int i = 0;
  i|=column; // lower 16 bits
  i|=(row << 16);
  return i;
}

void
PHG4ScintillatorSlatv1::identify(std::ostream& os) const
{
  os << "row " << row << " ";
  os << " column " << column << " ";
  os << " energy deposition " << get_edep();
  os << " ionization energy " << get_eion();
  os << " light yield " << get_light_yield();
  os << endl;
}
