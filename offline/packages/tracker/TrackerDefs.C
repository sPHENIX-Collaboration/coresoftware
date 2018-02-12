#include "TrackerDefs.h"
#include <bitset>

char TrackerDefs::get_trackerid(const TrackerDefs::hitkeytype key)
{
  hitkeytype tmp = (key >> bitshift_trackerid);
  return tmp;
}

char TrackerDefs::get_trackerid(const TrackerDefs::keytype key)
{
  hitkeytype tmp = (key >> 32);
  return get_trackerid(tmp);
}


char TrackerDefs::get_layer(const TrackerDefs::hitkeytype key)
{
  hitkeytype tmp = (key >> bitshift_layer);
  return tmp;
}

char TrackerDefs::get_layer(const TrackerDefs::keytype key)
{
  hitkeytype tmp = (key >> 32);
  return get_layer(tmp);
}

long
TrackerDefs::get_index(const TrackerDefs::keytype key)
{
  return key;
}

void 
TrackerDefs::print_bits(const TrackerDefs::hitkeytype key, std::ostream& os)
{
  os << "key: " << std::bitset<32>(key) << std::endl;
}

void 
TrackerDefs::print_bits(const TrackerDefs::keytype key, std::ostream& os)
{
  os << "key: " << std::bitset<64>(key) << std::endl;
}

//----------------------
// MVTXBinning Functions
//----------------------
char TrackerDefs::MVTXBinning::get_ladder(const TrackerDefs::hitkeytype key)
{
  hitkeytype tmp = (key >> bitshift_ladder);
  return tmp;
}

char TrackerDefs::MVTXBinning::get_ladder(const TrackerDefs::keytype key)
{
  hitkeytype tmp = (key >> 32);
  return get_ladder(tmp);
}

char TrackerDefs::MVTXBinning::get_chip(const TrackerDefs::hitkeytype key)
{
  hitkeytype tmp = (key >> bitshift_chip);
  return tmp;
}

char TrackerDefs::MVTXBinning::get_chip(const TrackerDefs::keytype key)
{
  hitkeytype tmp = (key >> 32);
  return get_chip(tmp);
}

unsigned short TrackerDefs::MVTXBinning::get_row(const TrackerDefs::keytype key)
{
  keytype tmp = (key >> bitshift_row);
  return tmp;
}

unsigned short TrackerDefs::MVTXBinning::get_col(const TrackerDefs::keytype key)
{
  keytype tmp = (key >> bitshift_col);
  return tmp;
}

TrackerDefs::hitkeytype
TrackerDefs::MVTXBinning::genhitkey(const char trackerid, const char layer,
                    const char ladder, const char chip)
{
  TrackerDefs::hitkeytype tmp = trackerid;
  TrackerDefs::hitkeytype key = tmp << TrackerDefs::bitshift_trackerid;  // detector id
  tmp = layer;
  key |= (tmp << TrackerDefs::bitshift_layer);  // layer
  tmp = ladder;
  key |= (tmp << TrackerDefs::bitshift_ladder);  // ladder
  tmp = chip;
  key |= (tmp << TrackerDefs::bitshift_chip);  // chip
  return key;
}

TrackerDefs::keytype
TrackerDefs::MVTXBinning::gencluskey(const char trackerid, const char layer,
                    const char ladder, const char chip,
                    const unsigned int bit32_index)
{
  TrackerDefs::keytype hitkey = genhitkey(trackerid, layer, ladder, chip);
  TrackerDefs::keytype key = hitkey << 32;  // detector id
  key |= bit32_index;                          // lower 32 bit index
  return key;
}