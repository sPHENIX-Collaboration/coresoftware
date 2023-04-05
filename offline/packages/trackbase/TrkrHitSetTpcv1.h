#ifndef TRACKBASE_TrkrHitSetTpcv1_H
#define TRACKBASE_TrkrHitSetTpcv1_H

#include "TrkrHitSetTpc.h"

#include <utility>  // for pair
#include <vector>

// forward declaration
class TrkrHit;

//! Vectorized TPC data time frame storage
class TrkrHitSetTpcv1 final : public TrkrHitSetTpc
{
 public:
  TrkrHitSetTpcv1() = default;

  TrkrHitSetTpcv1(const uint16_t n_pad, const uint16_t n_tbin)
    : TrkrHitSetTpc(n_pad, n_tbin)
  {
    Resize(n_pad, n_tbin);
  }

  ~TrkrHitSetTpcv1() override
  {
  }

  void identify(std::ostream& os = std::cout) const override;

  void Resize(const uint16_t n_pad, const uint16_t n_tbin) override;

  //! For ROOT TClonesArray end of event Operation
  void Clear(Option_t* /*option*/ = "") override { Reset(); }

  void Reset() override;

  void setHitSetKey(const TrkrDefs::hitsetkey key) override
  {
    m_hitSetKey = key;
  }

  TrkrDefs::hitsetkey getHitSetKey() const override
  {
    return m_hitSetKey;
  }

  // legacy TrkrDefs::hitkey accesses
  using TrkrHitSetTpc::getTpcADC;

  TpcDefs::ADCDataType& getTpcADC(const uint16_t local_pad, const uint16_t local_tbin) override;

  const TpcDefs::ADCDataType& getTpcADC(const uint16_t local_pad, const uint16_t local_tbin) const override;

  ConstIterator addHitSpecificKey(const TrkrDefs::hitkey, TrkrHit*) override;

  void removeHit(TrkrDefs::hitkey) override;

  TrkrHit* getHit(const TrkrDefs::hitkey) const override;

  ConstRange getHits() const override;

  unsigned int size() const override
  {
    return m_nPads * n_tBins;
  }

  uint16_t getNPads() const override
  {
    return m_nPads;
  }

  void setNPads(uint16_t nPads = 0) override
  {
    m_nPads = nPads;
  }

  uint16_t getTBins() const override
  {
    return n_tBins;
  }

  void setTBins(uint16_t tBins = 0) override
  {
    n_tBins = tBins;
  }

  uint16_t getPadIndexStart() const override
  {
    return m_padIndexStart;
  }

  void setPadIndexStart(uint16_t padIndexStart) override
  {
    m_padIndexStart = padIndexStart;
  }

  TpcDefs::BCODataType getStartingBco() const override
  {
    return m_StartingBCO;
  }

  void setStartingBco(TpcDefs::BCODataType startingBco) override
  {
    m_StartingBCO = startingBco;
  }

  uint16_t getTBinIndexStart() const override
  {
    return m_tBinIndexStart;
  }

  void setTBinIndexStart(uint16_t tBinIndexStart) override
  {
    m_tBinIndexStart = tBinIndexStart;
  }

  const TimeFrameADCDataType& getTimeFrameAdcData() const override
  {
    return m_timeFrameADCData;
  }

  TimeFrameADCDataType& getTimeFrameAdcData() override
  {
    return m_timeFrameADCData;
  }

  void setTimeFrameAdcData(const TimeFrameADCDataType& timeFrameAdcData) override
  {
    m_timeFrameADCData = timeFrameAdcData;
  }

 private:
  /// unique key for this object
  TrkrDefs::hitsetkey m_hitSetKey = TrkrDefs::HITSETKEYMAX;

  /// vector storage of TPC timeframe without zero suppression
  // Top level indexes are vectors of pads
  // Lower level indexes are vectors of time bin
  TimeFrameADCDataType m_timeFrameADCData;

  //! beam collision (BCO) clock counter at the start of the time frame
  TpcDefs::BCODataType m_StartingBCO = 0;

  //! local to global index conversion
  uint16_t m_padIndexStart = 0;
  uint16_t m_tBinIndexStart = 0;

  //! size of the local index
  uint16_t m_nPads = 0;
  uint16_t n_tBins = 0;

  ClassDefOverride(TrkrHitSetTpcv1, 1);
};

#endif  // TRACKBASE_TrkrHitSetTpcv1_H
