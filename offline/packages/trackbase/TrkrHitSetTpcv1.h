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

  TrkrHitSetTpcv1(const unsigned int n_pad, const unsigned int n_tbin)
    : TrkrHitSetTpc(n_pad, n_tbin)
  {
    Resize(n_pad, n_tbin);
  }

  ~TrkrHitSetTpcv1() override
  {
  }

  void identify(std::ostream& os = std::cout) const override;

  void Resize(const unsigned int n_pad, const unsigned int n_tbin) override;

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

  TpcDefs::ADCDataType& getTpcADC(const uint16_t pad, const uint16_t tbin) override;

  const TpcDefs::ADCDataType& getTpcADC(const uint16_t pad, const uint16_t tbin) const override;

  ConstIterator addHitSpecificKey(const TrkrDefs::hitkey, TrkrHit*) override;

  void removeHit(TrkrDefs::hitkey) override;

  TrkrHit* getHit(const TrkrDefs::hitkey) const override;

  ConstRange getHits() const override;

  unsigned int size() const override
  {
    return m_nPads * n_tBins;
  }

  unsigned int getNPads() const override
  {
    return m_nPads;
  }

  void setNPads(unsigned int nPads = 0) override
  {
    m_nPads = nPads;
  }

  const TimeFrameADCDataType& getTimeFrameAdcData() const override
  {
    return m_timeFrameADCData;
  }

  void setTimeFrameAdcData(const TimeFrameADCDataType& timeFrameAdcData) override
  {
    m_timeFrameADCData = timeFrameAdcData;
  }

  unsigned int getTBins() const override
  {
    return n_tBins;
  }

  void setTBins(unsigned int tBins = 0) override
  {
    n_tBins = tBins;
  }

 private:
  /// unique key for this object
  TrkrDefs::hitsetkey m_hitSetKey = TrkrDefs::HITSETKEYMAX;

  /// vector storage of TPC timeframe without zero suppression
  // Top level indexes are vectors of pads
  // Lower level indexes are vectors of time bin
  TimeFrameADCDataType m_timeFrameADCData;

  unsigned int m_nPads = 0;
  unsigned int n_tBins = 0;

  ClassDefOverride(TrkrHitSetTpcv1, 1);
};

#endif  // TRACKBASE_TrkrHitSetTpcv1_H
