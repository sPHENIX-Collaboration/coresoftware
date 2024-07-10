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
// the static_cast< uint16_t> is needed to because the result of (feeid >> 12U)
// is promoted to int which then triggers a (correct) clang-tidy warning that
// a bitwise operation is performed on a signed integer
  ret.layer = static_cast<uint16_t>(feeid >> 12U) & 0x7U;
  ret.stave = feeid & 0x1FU;
  ret.gbtid = static_cast<uint16_t>(feeid >> 8U) & 0x3U;
  return ret;
};
