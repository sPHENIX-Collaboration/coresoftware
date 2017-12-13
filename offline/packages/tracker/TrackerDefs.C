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

char TrackerDefs::MVTXBinning::get_ladder(const TrackerDefs::keytype key)
{
  keytype tmp = (key >> bitshift_ladder);
  return tmp;
}

char TrackerDefs::MVTXBinning::get_chip(const TrackerDefs::keytype key)
{
  keytype tmp = (key >> bitshift_chip);
  return tmp;
}

long
TrackerDefs::get_index(const TrackerDefs::keytype key)
{
  return key;
}

TrackerDefs::keytype
TrackerDefs::MVTXBinning::genhitkey(const char trackerid, const char layer,
                    const char ladder, const char chip,
                    const unsigned int row, const unsigned int col)
{
  TrackerDefs::keytype tmp = trackerid;
  TrackerDefs::keytype key = tmp << TrackerDefs::bitshift_trackerid;  // detector id
  tmp = layer;
  key |= (tmp << TrackerDefs::bitshift_layer);  // layer
  tmp = ladder;
  key |= (tmp << TrackerDefs::bitshift_ladder);  // ladder
  tmp = chip;
  key |= (tmp << TrackerDefs::bitshift_chip);  // chip
  tmp = row;
  key |= (tmp << TrackerDefs::bitshift_row);  // chip
  tmp = col;
  key |= (tmp << TrackerDefs::bitshift_col);  // chip
  return key;
}

TrackerDefs::keytype
TrackerDefs::MVTXBinning::gencluskey(const char trackerid, const char layer,
                    const char ladder, const char chip,
                    const long bit32_index)
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