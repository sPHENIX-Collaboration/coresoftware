#ifndef PHG4HITDEFS_H
#define PHG4HITDEFS_H

namespace PHG4HitDefs
{
  typedef unsigned long long keytype;
  static int keybits = 8;
  static int hit_idbits = sizeof(keytype)*8-keybits;
}

#endif


