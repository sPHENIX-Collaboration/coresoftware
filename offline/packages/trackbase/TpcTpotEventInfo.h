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
  virtual void CopyFrom( const TpcTpotEventInfo& ) 
  {}

  //! copy content from base class
  virtual void CopyFrom( TpcTpotEventInfo* ) 
  {}

  //
  // event tagger info
  //

  virtual uint16_t getTaggerType(int, int, int) const { return UINT16_MAX; } 
  virtual void setTaggerType(uint16_t, int, int, int) {} 
  virtual uint8_t getEnDat(int, int, int) const { return UINT8_MAX; } 
  virtual void setEnDat(uint8_t, int, int, int) {} 
  virtual uint8_t getIsLevel1(int, int, int) const { return UINT8_MAX; } 
  virtual void setIsLevel1(uint8_t, int, int, int) {} 
  virtual uint64_t getBCO(int, int, int) const { return UINT64_MAX; } 
  virtual void setBCO(uint64_t, int, int, int) {} 
  virtual uint32_t getLevel1Count(int, int, int) const { return UINT32_MAX; } 
  virtual void setLevel1Count(uint32_t, int, int, int) {} 
  virtual uint32_t getEnDatCount(int, int, int) const { return UINT32_MAX; } 
  virtual void setEnDatCount(uint32_t, int, int, int) {} 
  virtual uint64_t getLastBCO(int, int, int) const { return UINT64_MAX; } 
  virtual void setLastBCO(uint64_t, int, int, int) {} 
  virtual uint8_t getModebits(int, int, int) const { return UINT8_MAX; } 
  virtual void setModebits(uint8_t, int, int, int) {} 

 protected:
  TpcTpotEventInfo() = default;
  ClassDefOverride(TpcTpotEventInfo, 1)
};

#endif //TRACKBASE_TRKRCLUSTER_H
