#include "PHG4ScintillatorSlatDefs.h"

PHG4ScintillatorSlatDefs::keytype
PHG4ScintillatorSlatDefs::genkey(const short irow, const short icolumn)
{
  keytype key = icolumn;
  key |= (irow << columnbits);
  return key;
}

std::pair<short,short>
PHG4ScintillatorSlatDefs::getrowcol(const keytype key)
{
  short irow = key&0xFFFF;
  short icolumn = key >> columnbits;
  return std::make_pair(irow,icolumn);
}
