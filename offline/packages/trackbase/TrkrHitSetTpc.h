#ifndef TRACKBASE_TrkrHitSetTpc_H
#define TRACKBASE_TrkrHitSetTpc_H

#include "TpcDefs.h"
#include "TrkrDefs.h"
#include "TrkrHitSet.h"

#include <utility>  // for pair
#include <vector>

// forward declaration
class TrkrHit;

//! Generalized interface class for vectorized TPC data time frame storage
class TrkrHitSetTpc : public TrkrHitSet
{
 public:
  typedef std::vector<std::vector<TpcDefs::ADCDataType> > TimeFrameADCDataType;

  TrkrHitSetTpc() = default;

  TrkrHitSetTpc(const unsigned int /*n_pad*/, const unsigned int /*n_tbin*/) {}

  ~TrkrHitSetTpc() override
  {
  }

  virtual void identify(std::ostream& os = std::cout) const override;

  virtual void Resize(const unsigned int /*n_pad*/, const unsigned int /*n_tbin*/) {}

  TpcDefs::ADCDataType& getTpcADC(const TrkrDefs::hitkey);

  const TpcDefs::ADCDataType& getTpcADC(const TrkrDefs::hitkey) const;

  virtual TpcDefs::ADCDataType& getTpcADC(const uint16_t /*pad*/, const uint16_t /*tbin*/)
  {
    static TpcDefs::ADCDataType v = 0;
    return v;
  };

  virtual const TpcDefs::ADCDataType& getTpcADC(const uint16_t /*pad*/, const uint16_t /*tbin*/) const
  {
    static TpcDefs::ADCDataType v = 0;
    return v;
  };

  virtual unsigned int getNPads() const
  {
    return 0;
  }

  virtual void setNPads(unsigned int /*nPads*/)
  {
  }

  virtual const TimeFrameADCDataType& getTimeFrameAdcData() const
  {
    static TimeFrameADCDataType tmp;

    return tmp;
  }

  virtual void setTimeFrameAdcData(const TimeFrameADCDataType&)
  {
  }

  virtual unsigned int getTBins() const
  {
    return 0;
  }

  virtual void setTBins(unsigned int /*tBins*/)
  {
  }

 private:
  ClassDefOverride(TrkrHitSetTpc, 1);
};

#endif  // TRACKBASE_TrkrHitSetTpc_H
