#include "MvtxRawDefs.h"

uint8_t MvtxRawDefs::getStaveIndex( const uint8_t& lyrId, const uint8_t& stvId )
{
  return firstStaveIndex[lyrId] + stvId;
};

std::pair<uint8_t, uint8_t> const& MvtxRawDefs::get_flx_endpoint( const uint8_t& lyrId, const uint8_t& stvId )
{
  return stave_felix_map.at( getStaveIndex(lyrId, stvId) );
};

MvtxRawDefs::linkId_t MvtxRawDefs::decode_feeid( const uint16_t feeid )
{
  linkId_t ret = {};
  ret.layer = (feeid >> 12) & 0x7;
  ret.stave = feeid & 0x1F;
  ret.gbtid = (feeid >> 8) & 0x3;
  return ret;
};
