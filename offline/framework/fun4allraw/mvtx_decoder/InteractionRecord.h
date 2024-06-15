// @brief  Interaction record encoding BC, orbit, time
// @sa <O2/DataFormats/common/include/CommonDataFormat/InteractionRecord.h>
//     <29947a45e>

#ifndef MVTXDECODER_INTERACTIONRECORD_H
#define MVTXDECODER_INTERACTIONRECORD_H

#include <cstdint>
#include <string>

namespace mvtx
{
namespace lhcConstants
{
constexpr int LHCMaxBunches = 3564;                              // max N bunches
/*
constexpr double LHCRFFreq = 400.789e6;                          // LHC RF frequency in MHz
constexpr double LHCBunchSpacingNS = 10 * 1.e9 / LHCRFFreq;      // bunch spacing in ns (10 RFbuckets)
constexpr double LHCOrbitNS = LHCMaxBunches * LHCBunchSpacingNS; // orbit duration in ns
constexpr double LHCRevFreq = 1.e9 / LHCOrbitNS;                 // revolution frequency
constexpr double LHCBunchSpacingMUS = LHCBunchSpacingNS * 1e-3;  // bunch spacing in \mus (10 RFbuckets)
constexpr double LHCOrbitMUS = LHCOrbitNS * 1e-3;                // orbit duration in \mus

constexpr int RHICMaxBunches = 120;                              // max N bunches
constexpr double RHICFFreq = 9.839e6;                            // RHIC RF frequency in MHz
constexpr double RHICOrbitNS = 1.e9 / RHICFFreq;                 // RHIC BCO in nS
*/
} // namespace lhcConstants

//!< TODO: Add RHIC constants

struct InteractionRecord
{
  // information about bunch crossing and orbit
  static constexpr uint64_t DummyOrbit = 0xffffffffff;
  static constexpr uint16_t DummyBC = 0xffff;

  uint64_t orbit = DummyOrbit; ///< RHIC orbit (BCO)
  uint16_t bc = DummyBC;       ///< bunch crossing ID of interaction

  InteractionRecord() = default;

  InteractionRecord( uint64_t orb, uint16_t b ) : orbit(orb), bc(b){};

  InteractionRecord(const InteractionRecord& src) = default;

  InteractionRecord& operator=(const InteractionRecord& src) = default;

  void clear()
  {
    orbit = DummyOrbit;
    bc = DummyBC;
  }

  bool isDummy() const
  {
    return bc > mvtx::lhcConstants::LHCMaxBunches;
  }

  bool operator==(const InteractionRecord& other) const
  {
    return (orbit == other.orbit) && (bc == other.bc);
  }

  bool operator!=(const InteractionRecord& other) const
  {
    return (orbit != other.orbit) || (bc != other.bc);
  }

  void print() const;
  std::string asString() const;
  friend std::ostream& operator<<(std::ostream& stream, InteractionRecord const& ir);
// ClassDefNV(InteractionRecord, 3);
};

}
#endif
