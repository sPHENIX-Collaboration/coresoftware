/**
 * @file trackbase/TpcTpotEventInfo.h
 * @author Thomas Marshall
 * @date September 2023
 * @brief Base class for TPC and TPOT event level information
 */
#ifndef TRACKBASE_TPCTPOTEVENTINFO_H
#define TRACKBASE_TPCTPOTEVENTINFO_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <climits>
#include <cmath>
#include <iostream>
#include <memory>

/**
 * @brief Base class for TPC and TPOT event level information
 *
 * Virtual base class for TPC and TPOT event level information for
 * use in error checking and debugging
 */
class TpcTpotEventInfo : public PHObject
{
 public:
  //! dtor
  ~TpcTpotEventInfo() override = default;

  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TpcTpotEventInfo base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;

  //! copy content from base class
  virtual void CopyFrom(const TpcTpotEventInfo&)
  {
  }

  //! copy content from base class
  virtual void CopyFrom(TpcTpotEventInfo*)
  {
  }

  //
  // event tagger info
  //

  enum SectorID
  {
    kTPCSector0 = 0,
    kTPCSector1,
    kTPCSector2,
    kTPCSector3,
    kTPCSector4,
    kTPCSector5,
    kTPCSector6,
    kTPCSector7,
    kTPCSector8,
    kTPCSector9,
    kTPCSector10,
    kTPCSector11,
    kTPCSector12,
    kTPCSector13,
    kTPCSector14,
    kTPCSector15,
    kTPCSector16,
    kTPCSector17,
    kTPCSector18,
    kTPCSector19,
    kTPCSector20,
    kTPCSector21,
    kTPCSector22,
    kTPCSector23,
    kTPOT = 24
  };
  // TPC has sectors 0..23, TPOT uses index 24 in this storage
  enum PCIeEndPointID
  {
    kEndPoint0 = 0,
    kEndPoint1 = 1
  };
  enum TaggerID
  {
    kLVL1Tagger = 0,
    kEnDatTagger = 1
  };

  virtual void checkIndexes(SectorID, PCIeEndPointID, TaggerID)
  {
  }

  virtual uint64_t getBCO(SectorID, PCIeEndPointID, TaggerID) const { return UINT64_MAX; }
  virtual void setBCO(uint64_t, SectorID, PCIeEndPointID, TaggerID) {}
  virtual uint32_t getLevel1Count(SectorID, PCIeEndPointID, TaggerID) const { return UINT32_MAX; }
  virtual void setLevel1Count(uint32_t, SectorID, PCIeEndPointID, TaggerID) {}
  virtual uint32_t getEnDatCount(SectorID, PCIeEndPointID, TaggerID) const { return UINT32_MAX; }
  virtual void setEnDatCount(uint32_t, SectorID, PCIeEndPointID, TaggerID) {}
  virtual uint64_t getLastBCO(SectorID, PCIeEndPointID, TaggerID) const { return UINT64_MAX; }
  virtual void setLastBCO(uint64_t, SectorID, PCIeEndPointID, TaggerID) {}
  virtual uint8_t getModebits(SectorID, PCIeEndPointID, TaggerID) const { return UINT8_MAX; }
  virtual void setModebits(uint8_t, SectorID, PCIeEndPointID, TaggerID) {}

 protected:
  TpcTpotEventInfo() = default;
  ClassDefOverride(TpcTpotEventInfo, 1)
};

#endif  // TRACKBASE_TRKRCLUSTER_H
