/**
 * @file trackbase/TpcTpotEventInfov1.h
 * @author Thomas Marshall
 * @date September 2023
 * @brief Version 1 of TpcTpotEventInfo
 */
#ifndef TRACKBASE_TPCTPOTEVENTINFOV1_H
#define TRACKBASE_TPCTPOTEVENTINFOV1_H

#include "TpcTpotEventInfo.h"
#include "TrkrDefs.h"

#include <iostream>

class PHObject;

/**
 * @brief Version 1 of TpcTpotEventInfo
 *
 */
class TpcTpotEventInfov1 : public TpcTpotEventInfo
{
 public:
  //! ctor
  TpcTpotEventInfov1();

  //!dtor
  ~TpcTpotEventInfov1() override = default;
  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override 
  {
    std::fill_n(&m_bco[0][0][0],100,UINT64_MAX);
    std::fill_n(&m_lvl1_count[0][0][0],100,UINT32_MAX);
    std::fill_n(&m_endat_count[0][0][0],100,UINT32_MAX);
    std::fill_n(&m_last_bco[0][0][0],100,UINT64_MAX);
    std::fill_n(&m_modebits[0][0][0],100,UINT8_MAX);
  }
  int isValid() const override;
  PHObject* CloneMe() const override { return new TpcTpotEventInfov1(*this); }
 
  //! copy content from base class
  void CopyFrom( const TpcTpotEventInfo& ) override;

  //! copy content from base class
  void CopyFrom( TpcTpotEventInfo* source ) override
  { CopyFrom( *source ); }

  //
  // event tagger info
  //

  void checkIndexes(SectorID, PCIeEndPointID, TaggerID) override;

  uint64_t getBCO(SectorID sector, PCIeEndPointID PCIe, TaggerID tagger) const override { return m_bco[sector][PCIe][tagger]; }
  void setBCO(uint64_t bco, SectorID sector, PCIeEndPointID PCIe, TaggerID tagger) override { m_bco[sector][PCIe][tagger] = bco; }
  uint32_t getLevel1Count(SectorID sector, PCIeEndPointID PCIe, TaggerID tagger) const override { return m_lvl1_count[sector][PCIe][tagger]; }
  void setLevel1Count(uint32_t lvl1count, SectorID sector, PCIeEndPointID PCIe, TaggerID tagger) override { m_lvl1_count[sector][PCIe][tagger] = lvl1count; }
  uint32_t getEnDatCount(SectorID sector, PCIeEndPointID PCIe, TaggerID tagger) const override { return m_endat_count[sector][PCIe][tagger]; }
  void setEnDatCount(uint32_t endatcount, SectorID sector, PCIeEndPointID PCIe, TaggerID tagger) override { m_endat_count[sector][PCIe][tagger] = endatcount; }
  uint64_t getLastBCO(SectorID sector, PCIeEndPointID PCIe, TaggerID tagger) const override { return m_last_bco[sector][PCIe][tagger]; }
  void setLastBCO(uint64_t lastbco, SectorID sector, PCIeEndPointID PCIe, TaggerID tagger) override { m_last_bco[sector][PCIe][tagger] = lastbco; }
  uint8_t getModebits(SectorID sector, PCIeEndPointID PCIe, TaggerID tagger) const override { return m_modebits[sector][PCIe][tagger]; }
  void setModebits(uint8_t modebits, SectorID sector, PCIeEndPointID PCIe, TaggerID tagger) override { m_modebits[sector][PCIe][tagger] = modebits; }

 protected:

  uint64_t m_bco[25][2][2]{};
  uint32_t m_lvl1_count[25][2][2]{};
  uint32_t m_endat_count[25][2][2]{};
  uint64_t m_last_bco[25][2][2]{};
  uint8_t m_modebits[25][2][2]{};

  ClassDefOverride(TpcTpotEventInfov1, 1)
};

#endif //TRACKBASE_TPCTPOTEVENTINFOV1_H
