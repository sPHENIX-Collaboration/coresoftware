#ifndef PHG4CELLDEFS_H
#define PHG4CELLDEFS_H

namespace PHG4CellDefs
{
  typedef unsigned long long keytype;

  static unsigned int layer_bits = 16;
  static unsigned int keybits = 16;
  static unsigned int cell_idbits = 64 - (keybits+layer_bits);
  // __attribute__((unused)) prevents warnings about unused variables
#ifndef __CINT__
  static unsigned int scintillatorslat_bits __attribute__((unused)) = 16;
  static unsigned int scintillatorslat_mask __attribute__((unused)) = 0xFFFF;
  static unsigned int etaphi_bits __attribute__((unused)) = 16;
#endif

  enum CellBinning {undefined = 0, sizebinning = 1, etaphibinning = 2, etaslatbinning = 3, spacalbinning = 4,
        scintillatorslatbinning = 5};
  bool has_binning(PHG4CellDefs::keytype key, PHG4CellDefs::CellBinning binning);
  keytype genkey_scintillator_slat(const unsigned short layer, const unsigned short irow, const unsigned short icolumn);
  unsigned short int get_row(PHG4CellDefs::keytype key);
  unsigned short int get_column(PHG4CellDefs::keytype key);
  keytype genkey_eta_phi(const unsigned short layer, const unsigned short etabin, const unsigned short phibin);
  unsigned short int get_etabin(PHG4CellDefs::keytype key);
  unsigned short int get_phibin(PHG4CellDefs::keytype key);
}

#endif
