#include "CaloTriggerInfo_v1.h"


using namespace std;

CaloTriggerInfo_v1::CaloTriggerInfo_v1()
{
  _EMCAL_2x2_BEST_E = 0;
  _EMCAL_2x2_BEST_ETA = 0;
  _EMCAL_2x2_BEST_PHI = 0;

  _EMCAL_4x4_BEST_E = 0;
  _EMCAL_4x4_BEST_ETA = 0;
  _EMCAL_4x4_BEST_PHI = 0;

  _EMCAL_4x4_BEST2_E = 0;
  _EMCAL_4x4_BEST2_ETA = 0;
  _EMCAL_4x4_BEST2_PHI = 0;

  _FULLCALO_0p2x0p2_BEST_E = 0;
  _FULLCALO_0p2x0p2_BEST_ETA = 0;
  _FULLCALO_0p2x0p2_BEST_PHI = 0;

  _FULLCALO_0p4x0p4_BEST_E = 0;
  _FULLCALO_0p4x0p4_BEST_ETA = 0;
  _FULLCALO_0p4x0p4_BEST_PHI = 0;

  _FULLCALO_0p6x0p6_BEST_E = 0;
  _FULLCALO_0p6x0p6_BEST_ETA = 0;
  _FULLCALO_0p6x0p6_BEST_PHI = 0;

  _FULLCALO_0p8x0p8_BEST_E = 0;
  _FULLCALO_0p8x0p8_BEST_ETA = 0;
  _FULLCALO_0p8x0p8_BEST_PHI = 0;

  _FULLCALO_1p0x1p0_BEST_E = 0;
  _FULLCALO_1p0x1p0_BEST_ETA = 0;
  _FULLCALO_1p0x1p0_BEST_PHI = 0;

}

CaloTriggerInfo_v1::~CaloTriggerInfo_v1()
{
}

void CaloTriggerInfo_v1::identify(ostream& os) const
{
  os << "CaloTriggerInfo: highest EMCal 2x2 eta/phi = " << _EMCAL_2x2_BEST_ETA << " / " << _EMCAL_2x2_BEST_PHI << ", E = " << _EMCAL_2x2_BEST_E << std::endl;
  os << "CaloTriggerInfo: highest EMCal 4x4 eta/phi = " << _EMCAL_4x4_BEST_ETA << " / " << _EMCAL_4x4_BEST_PHI << ", E = " << _EMCAL_4x4_BEST_E << std::endl;
  os << "CaloTriggerInfo: 2nd highest EMCal 4x4 eta/phi = " << _EMCAL_4x4_BEST2_ETA << " / " << _EMCAL_4x4_BEST2_PHI << ", E = " << _EMCAL_4x4_BEST2_E << std::endl;
  os << "CaloTriggerInfo: highest FullCalo 0.2x0.2 eta/phi = " << _FULLCALO_0p2x0p2_BEST_ETA << " / " << _FULLCALO_0p2x0p2_BEST_PHI << ", E = " << _FULLCALO_0p2x0p2_BEST_E << std::endl;
  os << "CaloTriggerInfo: highest FullCalo 0.4x0.4 eta/phi = " << _FULLCALO_0p4x0p4_BEST_ETA << " / " << _FULLCALO_0p4x0p4_BEST_PHI << ", E = " << _FULLCALO_0p4x0p4_BEST_E << std::endl;
  os << "CaloTriggerInfo: highest FullCalo 0.6x0.6 eta/phi = " << _FULLCALO_0p6x0p6_BEST_ETA << " / " << _FULLCALO_0p6x0p6_BEST_PHI << ", E = " << _FULLCALO_0p6x0p6_BEST_E << std::endl;
  os << "CaloTriggerInfo: highest FullCalo 0.8x0.8 eta/phi = " << _FULLCALO_0p8x0p8_BEST_ETA << " / " << _FULLCALO_0p8x0p8_BEST_PHI << ", E = " << _FULLCALO_0p8x0p8_BEST_E << std::endl;
  os << "CaloTriggerInfo: highest FullCalo 1.0x1.0 eta/phi = " << _FULLCALO_1p0x1p0_BEST_ETA << " / " << _FULLCALO_1p0x1p0_BEST_PHI << ", E = " << _FULLCALO_1p0x1p0_BEST_E << std::endl;

  return;
}
