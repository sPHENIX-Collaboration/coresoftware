#ifndef PHG4HITDEFS_H
#define PHG4HITDEFS_H

#include <string>

namespace PHG4HitDefs
{
  typedef unsigned long long keytype;
  static const unsigned int keybits = 32;
  static const unsigned int hit_idbits = sizeof(keytype)*8-keybits;

  //! convert PHG4HitContainer node names in to ID number for the container.
  //! used in indexing volume ID in PHG4Shower
  int get_volume_id(const std::string & nodename);

}

#endif


