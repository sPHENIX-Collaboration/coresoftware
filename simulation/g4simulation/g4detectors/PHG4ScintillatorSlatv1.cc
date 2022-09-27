#include "PHG4ScintillatorSlatv1.h"

#include <g4main/PHG4HitDefs.h>  // for keytype

short PHG4ScintillatorSlatv1::get_row() const
{
  return (key & 0xFFFF);
}

short PHG4ScintillatorSlatv1::get_column() const
{
  return (key >> 16);
}

void PHG4ScintillatorSlatv1::identify(std::ostream& os) const
{
  os << "row " << get_row() << " ";
  os << " column " << get_column() << " ";
  os << " energy deposition " << get_edep();
  os << " ionization energy " << get_eion();
  os << " light yield " << get_light_yield();
  os << std::endl;
}
