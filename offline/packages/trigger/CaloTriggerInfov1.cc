#include "CaloTriggerInfov1.h"

#include <ostream>

void CaloTriggerInfov1::identify(std::ostream& os) const
{
  os << "CaloTriggerInfo: highest EMCal 2x2 eta/phi = " << m_EMCAL_2x2_BEST_ETA << " / " << m_EMCAL_2x2_BEST_PHI << ", E = " << m_EMCAL_2x2_BEST_E << std::endl;
  os << "CaloTriggerInfo: highest EMCal 4x4 eta/phi = " << m_EMCAL_4x4_BEST_ETA << " / " << m_EMCAL_4x4_BEST_PHI << ", E = " << m_EMCAL_4x4_BEST_E << std::endl;
  os << "CaloTriggerInfo: 2nd highest EMCal 4x4 eta/phi = " << m_EMCAL_4x4_BEST2_ETA << " / " << m_EMCAL_4x4_BEST2_PHI << ", E = " << m_EMCAL_4x4_BEST2_E << std::endl;
  os << "CaloTriggerInfo: highest FullCalo 0.2x0.2 eta/phi = " << m_FULLCALO_0p2x0p2_BEST_ETA << " / " << m_FULLCALO_0p2x0p2_BEST_PHI << ", E = " << m_FULLCALO_0p2x0p2_BEST_E << std::endl;
  os << "CaloTriggerInfo: highest FullCalo 0.4x0.4 eta/phi = " << m_FULLCALO_0p4x0p4_BEST_ETA << " / " << m_FULLCALO_0p4x0p4_BEST_PHI << ", E = " << m_FULLCALO_0p4x0p4_BEST_E << std::endl;
  os << "CaloTriggerInfo: highest FullCalo 0.6x0.6 eta/phi = " << m_FULLCALO_0p6x0p6_BEST_ETA << " / " << m_FULLCALO_0p6x0p6_BEST_PHI << ", E = " << m_FULLCALO_0p6x0p6_BEST_E << std::endl;
  os << "CaloTriggerInfo: highest FullCalo 0.8x0.8 eta/phi = " << m_FULLCALO_0p8x0p8_BEST_ETA << " / " << m_FULLCALO_0p8x0p8_BEST_PHI << ", E = " << m_FULLCALO_0p8x0p8_BEST_E << std::endl;
  os << "CaloTriggerInfo: highest FullCalo 1.0x1.0 eta/phi = " << m_FULLCALO_1p0x1p0_BEST_ETA << " / " << m_FULLCALO_1p0x1p0_BEST_PHI << ", E = " << m_FULLCALO_1p0x1p0_BEST_E << std::endl;

  return;
}
