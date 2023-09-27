/**
 * @file trackbase/TpcTpotEventInfov1.cc
 * @author Thomas Marshall
 * @date September 2023
 * @brief Implementation of TpcTpotEventInfov1
 */
#include "TpcTpotEventInfov1.h"

#include <cmath>
#include <utility>          // for swap
#include <cassert>

TpcTpotEventInfov1::TpcTpotEventInfov1()
{
  std::fill_n(&m_bco[0][0][0],100,UINT64_MAX);
  std::fill_n(&m_lvl1_count[0][0][0],100,UINT32_MAX);
  std::fill_n(&m_endat_count[0][0][0],100,UINT32_MAX);
  std::fill_n(&m_last_bco[0][0][0],100,UINT64_MAX);
  std::fill_n(&m_modebits[0][0][0],100,UINT8_MAX);
}

void TpcTpotEventInfov1::identify(std::ostream& os) const
{
  os << "---TpcTpotEventInfov1--------------------" << std::endl;

  for (int i = 0; i < 25; i++)
  {
    os << "Sector " << i << ": (PCIe 0 lvl1, PCIe 0 EnDat, PCIe 1 lvl1, PCIe 1 EnDat)" << std::endl;
    os << "BCO - " << m_bco[i][0][0] << ", " << m_bco[i][0][1] << ", " << m_bco[i][1][0] << ", " << m_bco[i][1][1] << std::endl;
    os << "Level 1 Count - " << m_lvl1_count[i][0][0] << ", " << m_lvl1_count[i][0][1] << ", " << m_lvl1_count[i][1][0] << ", " << m_lvl1_count[i][1][1] << std::endl;
    os << "EnDat Count - " << m_endat_count[i][0][0] << ", " << m_endat_count[i][0][1] << ", " << m_endat_count[i][1][0] << ", " << m_endat_count[i][1][1] << std::endl;
    os << "Last BCO - " << m_last_bco[i][0][0] << ", " << m_last_bco[i][0][1] << ", " << m_last_bco[i][1][0] << ", " << m_last_bco[i][1][1] << std::endl;
    os << "Modebits - " << m_modebits[i][0][0] << ", " << m_modebits[i][0][1] << ", " << m_modebits[i][1][0] << ", " << m_modebits[i][1][1] << std::endl;
  }

  os << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int TpcTpotEventInfov1::isValid() const
{
  return 1;
}

void TpcTpotEventInfov1::CopyFrom( const TpcTpotEventInfo& source )
{
  // do nothing if copying onto oneself
  if( this == &source ) return;
 
  // parent class method
  TpcTpotEventInfo::CopyFrom( source );
  
  for (int i = 0; i < 25; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      for (int k = 0; k < 2; k++)
      {
        setBCO(source.getBCO(static_cast<SectorID>(i), static_cast<PCIeEndPointID>(j), static_cast<TaggerID>(k)),static_cast<SectorID>(i), static_cast<PCIeEndPointID>(j), static_cast<TaggerID>(k));
        setLevel1Count(source.getLevel1Count(static_cast<SectorID>(i), static_cast<PCIeEndPointID>(j), static_cast<TaggerID>(k)),static_cast<SectorID>(i), static_cast<PCIeEndPointID>(j), static_cast<TaggerID>(k));
        setEnDatCount(source.getEnDatCount(static_cast<SectorID>(i), static_cast<PCIeEndPointID>(j), static_cast<TaggerID>(k)), static_cast<SectorID>(i), static_cast<PCIeEndPointID>(j), static_cast<TaggerID>(k));
        setLastBCO(source.getLastBCO(static_cast<SectorID>(i), static_cast<PCIeEndPointID>(j), static_cast<TaggerID>(k)),static_cast<SectorID>(i), static_cast<PCIeEndPointID>(j), static_cast<TaggerID>(k));
        setModebits(source.getModebits(static_cast<SectorID>(i), static_cast<PCIeEndPointID>(j), static_cast<TaggerID>(k)),static_cast<SectorID>(i), static_cast<PCIeEndPointID>(j), static_cast<TaggerID>(k));
      } 
    }
  }
}

void TpcTpotEventInfov1::checkIndexes(SectorID sector, PCIeEndPointID PCIe, TaggerID tagger)
{
  assert( sector>=kTPCSector0 );
  assert( sector<=kTPOT );
  assert( PCIe>=kEndPoint0 );
  assert( PCIe<=kEndPoint1 );
  assert( tagger>=kLVL1Tagger );
  assert( tagger<=kEnDatTagger );
}
