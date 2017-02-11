#include "PHG4CellDefs.h"

#include <cmath>
#include <iostream>

using namespace std;

PHG4CellDefs::keytype
PHG4CellDefs::genkey_scintillator_slat(const unsigned short layer, const unsigned short irow, const unsigned short icolumn)
{
  keytype tmp = layer;
  keytype key = tmp << (cell_idbits+keybits); // binning method used to decode the key
  tmp = scintillatorslatbinning;
  key |= (tmp << cell_idbits); // binning method used to decode the key
  key |= (icolumn << scintillatorslat_bits); // upper bits used by column, so we can easily extract 
                                  // slats by column which are combined to towers
  key += irow;
  return key;
}

unsigned short int
PHG4CellDefs::get_row(PHG4CellDefs::keytype key)
{
  // check correct binning first
 keytype tmp = scintillatorslatbinning;
 tmp = (tmp << cell_idbits);
 if ((key & tmp) == tmp)
   {
     unsigned short int row = (key & 0xFFFF);
     return row;
   }
 cout << "row could not decode 0x" << hex << key << dec << endl; 
 exit(1);
}

unsigned short int
PHG4CellDefs::get_column(PHG4CellDefs::keytype key)
{
  // check correct binning first
 keytype tmp = scintillatorslatbinning;
 keytype keysave = key; 
 tmp = (tmp << cell_idbits);
 if ((key & tmp) == tmp)
   {
     key = key >> scintillatorslat_bits;
     unsigned short int column = (key & 0x1FFF);
     return column;
   }
 cout << "col could not decode 0x" << hex << keysave << dec << endl; 
 exit(1);
}

PHG4CellDefs::keytype
PHG4CellDefs::genkey_eta_phi(const unsigned short layer, const unsigned short eta, const unsigned short phi)
{
  if (eta > 0x1FFF)
    {
      cout << "etabin " << eta << " exceeds range of " << 0x1FFF << endl;
      exit(1);
    }
  if (phi > 0x1FFF)
    {
      cout << "phibin " << phi << " exceeds range of " << 0x1FFF << endl;
      exit(1);
    }
  keytype tmp = etaphibinning;
  keytype key = tmp << cell_idbits; // binning method used to decode the key
  key |= (phi << scintillatorslat_bits); // upper bits used by column, so we can easily extract 
                                  // slats by column which are combined to towers
  key += eta;
  return key;
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
