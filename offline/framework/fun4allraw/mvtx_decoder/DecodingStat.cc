// @file ChipStat.cxx
// @brief Alpide Chip decoding statistics
// @sa <O2/Detectors/ITSMFT/common/reconstruction/src/DecodingStat.cxx>
//     <3cbcf82df>

#include "mvtx_decoder/DecodingStat.h"
#include "mvtx_decoder/PixelData.h"
#include "mvtx_decoder/mvtx_utils.h"

#include <bitset>
#include <sstream>

using namespace mvtx;

#if (__cplusplus >= CXX_17)
constexpr std::array<std::string_view, ChipStat::NErrorsDefined> ChipStat::ErrNames;
#endif

constexpr std::array<uint32_t, ChipStat::NErrorsDefined> ChipStat::ErrActions;

//________________________________________________________________________________
uint32_t ChipStat::getNErrors() const
{
  uint32_t nerr = 0;
  for (int i = NErrorsDefined; i--;)
  {
    nerr += errorCounts[i];
  }
  return nerr;
}

//________________________________________________________________________________
/// print link decoding statistics
uint32_t ChipStat::addErrors(uint32_t mask, uint16_t chID, int verbosity)
{
  uint32_t res = 0;
  if (mask)
  {
    for (uint32_t i = NErrorsDefined; i--;)
    {
      if (mask & (0x1U << i))
      {
        res |= ErrActions[i] & ErrActPropagate;
        if (verbosity > -1 && (!errorCounts[i] || verbosity > 1))
        {
          log_error << "New error registered on the FEEID: " << feeID << ": chip: " << chID << " " << ErrNames[i] << std::endl;
          res |= ErrActions[i] & ErrActDump;
        }
        errorCounts[i]++;
      }
    }
  }
  return res;
}

//________________________________________________________________________________
/// print link decoding statistics
uint32_t ChipStat::addErrors(const ChipPixelData& d, int verbosity)
{
  uint32_t res = 0;
  if (d.getErrorFlags())
  {
    for (uint32_t i = NErrorsDefined; i--;)
    {
      if (d.getErrorFlags() & (0x1U << i))
      {
        res |= ErrActions[i] & ErrActPropagate;
        if (verbosity > -1 && (!errorCounts[i] || verbosity > 1))
        {
          //          LOGP(info, "New error registered at bc/orbit {}/{} on the FEEID:{:#04x} chip#{}: {}{}",
          //               d.getInteractionRecord().bc, d.getInteractionRecord().orbit,
          //               feeID, int16_t(d.getChipID()), ErrNames[i], d.getErrorDetails(i));
          res |= ErrActions[i] & ErrActDump;
        }
        errorCounts[i]++;
      }
    }
  }
  return res;
}

//________________________________________________________________________________
/// print chip decoding statistics
void ChipStat::print(bool skipNoErr, const std::string& pref) const
{
  uint32_t nErr = 0;
  for (int i = NErrorsDefined; i--;)
  {
    nErr += errorCounts[i];
  }
  if (!skipNoErr || nErr)
  {
    std::ostringstream rep(pref.c_str());
    rep << feeID << " NHits: " << nHits << " errors: " << nErr;
    for (uint32_t i = 0; i < NErrorsDefined; i++)
    {
      if (!skipNoErr || errorCounts[i])
      {
        rep << " | Err.: " << ErrNames[i].data() << " : " << errorCounts[i];
      }
    }
    log_error << rep.str().c_str() << std::endl;
  }
}

//________________________________________________________________________________
/// print link decoding statistics
void GBTLinkDecodingStat::print(bool skipNoErr) const
{
  int nErr = 0;
  for (int i = NErrorsDefined; i--;)
  {
    nErr += errorCounts[i];
  }
  if (!skipNoErr || nErr)
  {
    //    std::string rep = fmt::format("FEEID#{:#04x} Packet States Statistics (total packets: {}, triggers: {})", feeID, nPackets, nTriggers);
    //    bool countsSeen = false;
    //    for (int i = 0; i < GBTDataTrailer::MaxStateCombinations; i++) {
    //      if (packetStates[i]) {
    //        if (!countsSeen) {
    //          countsSeen = true;
    //          rep += " | counts for triggers: ";
    //        } else {
    //          rep += ", ";
    //        }
    //        std::bitset<GBTDataTrailer::NStatesDefined> patt(i);
    //        rep += fmt::format("b{:s}: {}", patt.to_string().c_str(), packetStates[i]);
    //      }
    //    }
    //    rep += fmt::format(" | Decoding errors: {}", nErr);
    //    for (int i = 0; i < NErrorsDefined; i++) {
    //      if (!skipNoErr || errorCounts[i]) {
    //        rep += fmt::format(" [{}: {}]", ErrNames[i].data(), errorCounts[i]);
    //      }
    //    }
    //    LOG(info) << rep;
  }
}
