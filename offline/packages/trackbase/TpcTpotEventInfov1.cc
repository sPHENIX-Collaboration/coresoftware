/**
 * @file trackbase/TpcTpotEventInfov1.cc
 * @author Thomas Marshall
 * @date September 2023
 * @brief Implementation of TpcTpotEventInfov1
 */
#include "TpcTpotEventInfov1.h"

#include <cmath>
#include <utility>          // for swap

TpcTpotEventInfov1::TpcTpotEventInfov1()
{}

void TpcTpotEventInfov1::identify(std::ostream& os) const
{
  os << "---TpcTpotEventInfov1--------------------" << std::endl;

  for (int i = 0; i < 25; i++)
  {
    os << "Sector " << i << ": (PCIe 0 lvl1, PCIe 0 EnDat, PCIe 1 lvl1, PCIe 1 EnDat)" << std::endl;
    os << "Tagger Type - " << m_tagger_type[i][0][0] << ", " << m_tagger_type[i][0][1] << ", " << m_tagger_type[i][1][0] << ", " << m_tagger_type[i][1][1] << std::endl;
    os << "Is EnDat - " << m_is_endat[i][0][0] << ", " << m_is_endat[i][0][1] << ", " << m_is_endat[i][1][0] << ", " << m_is_endat[i][1][1] << std::endl;
    os << "Is Level 1 - " << m_is_lvl1[i][0][0] << ", " << m_is_lvl1[i][0][1] << ", " << m_is_lvl1[i][1][0] << ", " << m_is_lvl1[i][1][1] << std::endl;
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
        setTaggerType(source.getTaggerType(i, j, k), i, j, k);
        setEnDat(source.getEnDat(i, j, k), i, j, k);
        setIsLevel1(source.getIsLevel1(i, j, k), i, j, k);
        setBCO(source.getBCO(i, j, k), i, j, k);
        setLevel1Count(source.getLevel1Count(i, j, k), i, j, k);
        setEnDatCount(source.getEnDatCount(i, j, k), i, j, k);
        setLastBCO(source.getLastBCO(i, j, k), i, j, k);
        setModebits(source.getModebits(i, j, k), i, j, k);
      } 
    }
  }
}

