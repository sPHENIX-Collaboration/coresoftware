#ifndef TRACKBASE_TrkrHitSetTpcv1_H
#define TRACKBASE_TrkrHitSetTpcv1_H

#include "TpcDefs.h"
#include "TrkrDefs.h"
#include "TrkrHitSet.h"

#include <utility>  // for pair
#include <vector>

// forward declaration
class TrkrHit;

//! Vectorized TPC data time frame storage
class TrkrHitSetTpcv1 : public TrkrHitSet
{
 public:
  TrkrHitSetTpcv1() = default;

  TrkrHitSetTpcv1(const unsigned int n_pad, const unsigned int n_tbin) { Resize(n_pad, n_tbin); }

  ~TrkrHitSetTpcv1() override
  {
  }

  void identify(std::ostream& os = std::cout) const override;

  void Resize(const unsigned int n_pad, const unsigned int n_tbin);

  void Reset() override;

  void setHitSetKey(const TrkrDefs::hitsetkey key) override
  {
    m_hitSetKey = key;
  }

  TrkrDefs::hitsetkey getHitSetKey() const override
  {
    return m_hitSetKey;
  }

  TpcDefs::ADCDataType& getTpcADC(const TrkrDefs::hitkey);

  const TpcDefs::ADCDataType& getTpcADC(const TrkrDefs::hitkey) const;

  TpcDefs::ADCDataType& getTpcADC(const uint16_t pad, const uint16_t tbin);

  const TpcDefs::ADCDataType& getTpcADC(const uint16_t pad, const uint16_t tbin) const;

  ConstIterator addHitSpecificKey(const TrkrDefs::hitkey, TrkrHit*) override;

  void removeHit(TrkrDefs::hitkey) override;

  TrkrHit* getHit(const TrkrDefs::hitkey) const override;

  ConstRange getHits() const override;

  unsigned int size() const override
  {
    return m_hits.size();
  }

  unsigned int getNPads() const
  {
    return m_nPads;
  }

  void setNPads(unsigned int nPads = 0)
  {
    m_nPads = nPads;
  }

  const std::vector<std::vector<TpcDefs::ADCDataType> >& getTimeFrameAdcData() const
  {
    return m_timeFrameADCData;
  }

  void setTimeFrameAdcData(const std::vector<std::vector<TpcDefs::ADCDataType> >& timeFrameAdcData)
  {
    m_timeFrameADCData = timeFrameAdcData;
  }

  unsigned int getTBins() const
  {
    return n_tBins;
  }

  void setTBins(unsigned int tBins = 0)
  {
    n_tBins = tBins;
  }

 private:
  /// unique key for this object
  TrkrDefs::hitsetkey m_hitSetKey = TrkrDefs::HITSETKEYMAX;

  /// vector storage of TPC timeframe without zero suppression
  // Top level indexes are vectors of pads
  // Lower level indexes are vectors of time bin
  std::vector<std::vector<TpcDefs::ADCDataType> > m_timeFrameADCData;

  unsigned int m_nPads = 0;
  unsigned int n_tBins = 0;

  ClassDefOverride(TrkrHitSetTpcv1, 1);
};

#endif  // TRACKBASE_TrkrHitSetTpcv1_H
