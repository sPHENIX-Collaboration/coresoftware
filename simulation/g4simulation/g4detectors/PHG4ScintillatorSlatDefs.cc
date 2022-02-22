#include "PHG4ScintillatorSlatDefs.h"

PHG4ScintillatorSlatDefs::keytype
PHG4ScintillatorSlatDefs::genkey(const short irow, const short icolumn)
{
  keytype key = irow;              // lower bits used by row
  key |= (icolumn << columnbits);  // upper bits used by column, so we can easily extract
                                   // slats by column which are combined to towers
  return key;
}

std::pair<short, short>
PHG4ScintillatorSlatDefs::getrowcol(const keytype key)
{
  short irow = key & 0xFFFF;
  short icolumn = key >> columnbits;
  return std::make_pair(irow, icolumn);
}
