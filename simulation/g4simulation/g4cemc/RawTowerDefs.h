#ifndef RAWTOWERDEFS_H
#define RAWTOWERDEFS_H

namespace RawTowerDefs
{
  typedef unsigned int keytype;
  static unsigned int calo_idbits = 8;
  static unsigned int tower_idbits = sizeof(keytype)*8 - calo_idbits;
  static unsigned int eta_idbits = tower_idbits/2;
}

#endif
