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
    memset(m_tagger_type,0,sizeof(m_tagger_type));
    memset(m_is_endat,0,sizeof(m_is_endat));
    memset(m_is_lvl1,0,sizeof(m_is_lvl1));
    memset(m_bco,0,sizeof(m_bco));
    memset(m_lvl1_count,0,sizeof(m_lvl1_count));
    memset(m_endat_count,0,sizeof(m_endat_count));
    memset(m_last_bco,0,sizeof(m_last_bco));
    memset(m_modebits,0,sizeof(m_modebits));
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

  uint16_t getTaggerType(int sector, int PCIe, int tagger) const override { return m_tagger_type[sector][PCIe][tagger]; }
  void setTaggerType(uint16_t type, int sector, int PCIe, int tagger) override { m_tagger_type[sector][PCIe][tagger] = type; }
  uint8_t getEnDat(int sector, int PCIe, int tagger) const override { return m_is_endat[sector][PCIe][tagger]; }
  void setEnDat(uint8_t endat, int sector, int PCIe, int tagger) override { m_is_endat[sector][PCIe][tagger] = endat; }
  uint8_t getIsLevel1(int sector, int PCIe, int tagger) const override { return m_is_lvl1[sector][PCIe][tagger]; }
  void setIsLevel1(uint8_t islvl1, int sector, int PCIe, int tagger) override { m_is_lvl1[sector][PCIe][tagger] = islvl1; }
  uint64_t getBCO(int sector, int PCIe, int tagger) const override { return m_bco[sector][PCIe][tagger]; }
  void setBCO(uint64_t bco, int sector, int PCIe, int tagger) override { m_bco[sector][PCIe][tagger] = bco; }
  uint32_t getLevel1Count(int sector, int PCIe, int tagger) const override { return m_lvl1_count[sector][PCIe][tagger]; }
  void setLevel1Count(uint32_t lvl1count, int sector, int PCIe, int tagger) override { m_lvl1_count[sector][PCIe][tagger] = lvl1count; }
  uint32_t getEnDatCount(int sector, int PCIe, int tagger) const override { return m_endat_count[sector][PCIe][tagger]; }
  void setEnDatCount(uint32_t endatcount, int sector, int PCIe, int tagger) override { m_endat_count[sector][PCIe][tagger] = endatcount; }
  uint64_t getLastBCO(int sector, int PCIe, int tagger) const override { return m_last_bco[sector][PCIe][tagger]; }
  void setLastBCO(uint64_t lastbco, int sector, int PCIe, int tagger) override { m_last_bco[sector][PCIe][tagger] = lastbco; }
  uint8_t getModebits(int sector, int PCIe, int tagger) const override { return m_modebits[sector][PCIe][tagger]; }
  void setModebits(uint8_t modebits, int sector, int PCIe, int tagger) override { m_modebits[sector][PCIe][tagger] = modebits; }

 protected:

  uint16_t m_tagger_type[25][2][2]{};
  uint8_t m_is_endat[25][2][2]{};
  uint8_t m_is_lvl1[25][2][2]{};
  uint64_t m_bco[25][2][2]{};
  uint32_t m_lvl1_count[25][2][2]{};
  uint32_t m_endat_count[25][2][2]{};
  uint64_t m_last_bco[25][2][2]{};
  uint8_t m_modebits[25][2][2]{};

  ClassDefOverride(TpcTpotEventInfov1, 1)
};

#endif //TRACKBASE_TPCTPOTEVENTINFOV1_H
