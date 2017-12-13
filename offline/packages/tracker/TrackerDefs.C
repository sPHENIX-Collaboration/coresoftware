#include "TrackerDefs.h"

char TrackerDefs::get_trackerid(const TrackerDefs::keytype key)
{
  keytype tmp = (key >> bitshift_trackerid);
  return tmp;
}

char TrackerDefs::get_layer(const TrackerDefs::keytype key)
{
  keytype tmp = (key >> bitshift_layer);
  return tmp;
}

char TrackerDefs::get_ladder(const TrackerDefs::keytype key)
{
  keytype tmp = (key >> bitshift_ladder);
  return tmp;
}

char TrackerDefs::get_chip(const TrackerDefs::keytype key)
{
  keytype tmp = (key >> bitshift_chip);
  return tmp;
}

unsigned int
TrackerDefs::get_index(const TrackerDefs::keytype key)
{
  return key;
}

TrackerDefs::keytype
TrackerDefs::genkey(const char trackerid, const char layer,
                    const char ladder, const char chip,
                    const unsigned int bit32_index)
{
  TrackerDefs::keytype tmp = trackerid;
  TrackerDefs::keytype key = tmp << TrackerDefs::bitshift_trackerid;  // detector id
  tmp = layer;
  key |= (tmp << TrackerDefs::bitshift_layer);  // layer
  tmp = ladder;
  key |= (tmp << TrackerDefs::bitshift_ladder);  // ladder
  tmp = chip;
  key |= (tmp << TrackerDefs::bitshift_chip);  // chip
  key |= bit32_index;                          // lower 32 bit index
  return key;
}