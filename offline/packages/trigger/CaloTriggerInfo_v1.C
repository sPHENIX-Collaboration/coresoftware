#include "CaloTriggerInfo_v1.h"

ClassImp(CaloTriggerInfo_v1)

CaloTriggerInfo_v1::CaloTriggerInfo_v1()
{
  _EMCAL_4x4_BEST_E = 0;
  _EMCAL_4x4_BEST_ETA = 0;
  _EMCAL_4x4_BEST_PHI = 0;
}

CaloTriggerInfo_v1::~CaloTriggerInfo_v1()
{
}

void CaloTriggerInfo_v1::identify(ostream& os) const
{
  os << "CaloTriggerInfo: highest 2x2 eta/phi = " << _EMCAL_4x4_BEST_ETA << " / " << _EMCAL_4x4_BEST_PHI << ", E = " << _EMCAL_4x4_BEST_E << std::endl;

  return;
}
