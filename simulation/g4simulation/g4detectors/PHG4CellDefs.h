#ifndef PHG4CELLDEFS_H
#define PHG4CELLDEFS_H

namespace PHG4CellDefs
{

  typedef unsigned long long keytype;

  // CINT does not know the __attribute__((unused))
#ifndef __CINT__
  // key layout
  // bit
  // 48-64 detector id (scintillator slat id, layer,...)
  // 32-48 binning method
  // 
  // 16-32 binning dependant 1st key
  // 0-16 binning dependant 2nd key

  /* static unsigned int layer_bits = 16; */
  /* static unsigned int keybits = 16; */
  // __attribute__((unused)) prevents warnings about unused variables
  // common upper 32 bits 
  static unsigned int bitshift_layer __attribute__((unused)) = 64 - 16;
  static unsigned int bitshift_binning __attribute__((unused)) = bitshift_layer - 16;
  // binning dependeant bit shifts for lower 32 bits
  static unsigned int bitshift_upperkey __attribute__((unused)) = 16;
  static unsigned int bitshift_row __attribute__((unused)) = 16;
  static unsigned int bitshift_phi __attribute__((unused)) = 16;
#endif

  enum CellBinning {undefined = 0, sizebinning = 1, etaphibinning = 2, etaslatbinning = 3, spacalbinning = 4, scintillatorslatbinning = 5, etaxsizebinning = 6, mapsbinning = 7, tpcbinning = 8};
  bool has_binning(PHG4CellDefs::keytype key, PHG4CellDefs::CellBinning binning);
  short get_binning(PHG4CellDefs::keytype key);
  short int get_detid(const PHG4CellDefs::keytype key);

  namespace SizeBinning
  {
    keytype genkey(const unsigned short layer, const unsigned short zbin, const unsigned short iphibin);
    unsigned short int get_zbin(const PHG4CellDefs::keytype key);
    unsigned short int get_phibin(const PHG4CellDefs::keytype key);
  };

  namespace EtaPhiBinning
  {
    keytype genkey(const unsigned short layer, const unsigned short etabin, const unsigned short phibin);
    unsigned short int get_etabin(const PHG4CellDefs::keytype key);
    unsigned short int get_phibin(const PHG4CellDefs::keytype key);
  };

  namespace SpacalBinning
  {
    keytype genkey(const unsigned short etabin, const unsigned short phibin, const unsigned short fiberid);
    unsigned short get_etabin(PHG4CellDefs::keytype key);
    unsigned short get_phibin(PHG4CellDefs::keytype key);
    unsigned short get_fiberid(PHG4CellDefs::keytype key);
  };

  namespace ScintillatorSlatBinning 
  {
    keytype genkey(const unsigned short layer, const unsigned short irow, const unsigned short icolumn);
    unsigned short int get_row(PHG4CellDefs::keytype key);
    unsigned short int get_column(PHG4CellDefs::keytype key);
  };

  namespace EtaXsizeBinning
  {
    keytype genkey(const unsigned short layer, const unsigned short etabin, const unsigned short xbin);
    unsigned short int get_etabin(const PHG4CellDefs::keytype key);
    unsigned short int get_xsizebin(const PHG4CellDefs::keytype key);
  };

  namespace MapsBinning
  {
    keytype genkey(const unsigned short layer, const unsigned int bit32_index);
  };

  namespace TPCBinning
  {
    keytype genkey(const unsigned short lyr, const unsigned short mod, const unsigned short pad);
  };


}

#endif
