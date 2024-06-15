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

  TrkrHitSetTpc(const uint16_t /*n_pad*/, const uint16_t /*n_tbin*/) {}

  ~TrkrHitSetTpc() override
  {
  }

  virtual void identify(std::ostream& os = std::cout) const override;

  //! resize to fit getNPads and getNTBins
  virtual void Resize() {}

  //! global -> local conversion
  std::pair<uint16_t, uint16_t> getLocalPhiTBin(const TrkrDefs::hitkey) const;

  //! local -> global conversion
  TrkrDefs::hitkey getHitKeyfromLocalBin(const uint16_t /*local_pad*/, const uint16_t /*local_tbin*/) const;

  TpcDefs::ADCDataType& getTpcADC(const TrkrDefs::hitkey);

  const TpcDefs::ADCDataType& getTpcADC(const TrkrDefs::hitkey) const;

  virtual TpcDefs::ADCDataType& getTpcADC(const uint16_t /*local_pad*/, const uint16_t /*local_tbin*/)
  {
    static TpcDefs::ADCDataType v = 0;
    return v;
  };

  virtual const TpcDefs::ADCDataType& getTpcADC(const uint16_t /*local_pad*/, const uint16_t /*local_tbin*/) const
  {
    static TpcDefs::ADCDataType v = 0;
    return v;
  };

  virtual uint16_t getNPads() const
  {
    return 0;
  }

  virtual void setNPads(uint16_t /*nPads*/)
  {
  }

  virtual const TimeFrameADCDataType& getTimeFrameAdcData() const
  {
    static TimeFrameADCDataType tmp;

    return tmp;
  }
  virtual  TimeFrameADCDataType& getTimeFrameAdcData()
  {
    static TimeFrameADCDataType tmp;

    return tmp;
  }

  virtual void setTimeFrameAdcData(const TimeFrameADCDataType&)
  {
  }

  virtual uint16_t getNTBins() const
  {
    return 0;
  }

  virtual void setNTBins(uint16_t /*tBins*/)
  {
  }

  virtual uint16_t getPadIndexStart() const
  {
    return 0;
  }

  virtual void setPadIndexStart(uint16_t )
  {
  }

  virtual TpcDefs::BCODataType getStartingBco() const
  {
    return 0;
  }

  virtual void setStartingBco(TpcDefs::BCODataType )
  {
  }

  virtual uint16_t getTBinIndexStart() const
  {
    return 0;
  }

  virtual void setTBinIndexStart(uint16_t)
  {
  }

 private:
  ClassDefOverride(TrkrHitSetTpc, 1);
};

#endif  // TRACKBASE_TrkrHitSetTpc_H
