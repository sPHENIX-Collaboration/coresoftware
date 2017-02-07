#ifndef PHG4CELLDEFS_H
#define PHG4CELLDEFS_H

namespace PHG4CellDefs
{
  typedef unsigned int keytype;
  static int keybits = 8;
  static int cell_idbits = 32 - keybits;
  enum {undefined = 0, sizebinning = 1, etaphibinning = 2, etaslatbinning = 3, spacalbinning = 4};
}

#endif
