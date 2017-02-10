#include "PHG4CellDefs.h"

#include <cmath>
#include <iostream>

using namespace std;

PHG4CellDefs::keytype
PHG4CellDefs::genkey_scintillator_slat(const short irow, const short icolumn)
{
  keytype tmp = scintillatorslatbinning;
  keytype key = tmp << cell_idbits; // binning method used to decode the key
  key |= (icolumn << scintillatorslat_bits); // upper bits used by column, so we can easily extract 
                                  // slats by column which are combined to towers
  key += irow;
  return key;
}

int
PHG4CellDefs::get_row(PHG4CellDefs::keytype key)
{
  // check correct binning first
 keytype tmp = scintillatorslatbinning;
 tmp = (tmp << cell_idbits);
 if ((key & tmp) == tmp)
   {
     int row = (key & 0x1FFF);
     return row;
   }
 cout << "row could not decode 0x" << hex << key << dec << endl; 
 exit(1);
}

int
PHG4CellDefs::get_column(PHG4CellDefs::keytype key)
{
  // check correct binning first
 keytype tmp = scintillatorslatbinning;
 keytype keysave = key; 
 tmp = (tmp << cell_idbits);
 if ((key & tmp) == tmp)
   {
     key = key >> scintillatorslat_bits;
     int column = (key & 0x1FFF);
     return column;
   }
 cout << "col could not decode 0x" << hex << keysave << dec << endl; 
 exit(1);
}

bool
PHG4CellDefs::has_binning(PHG4CellDefs::keytype key, PHG4CellDefs::CellBinning binning)
{
 keytype tmp = binning;
 tmp = (tmp << cell_idbits);
 if ( (key & tmp) == tmp)
   {
     return true;
   }
 return false;
}
