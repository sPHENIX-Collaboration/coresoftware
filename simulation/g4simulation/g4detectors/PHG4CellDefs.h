#ifndef PHG4CELLDEFS_H
#define PHG4CELLDEFS_H

namespace PHG4CellDefs
{
  typedef unsigned int keytype;

  static int keybits = 6;
  static int cell_idbits = 32 - keybits;
  static int scintillatorslat_bits
#ifndef __CINT__
 __attribute__((unused)) 
#endif
= 13;
  enum CellBinning {undefined = 0, sizebinning = 1, etaphibinning = 2, etaslatbinning = 3, spacalbinning = 4,
        scintillatorslatbinning = 5};
  keytype genkey_scintillator_slat(const short irow, const short icolumn);
  int get_row(PHG4CellDefs::keytype key);
  int get_column(PHG4CellDefs::keytype key);
  bool has_binning(PHG4CellDefs::keytype key, PHG4CellDefs::CellBinning binning);
}

#endif
