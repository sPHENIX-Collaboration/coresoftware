#include "CaloTriggerInfo.h"

ClassImp(CaloTriggerInfo)

CaloTriggerInfo::CaloTriggerInfo()
{
  _EMCAL_4x4_BEST_E = 0;
  _EMCAL_4x4_BEST_ETA = 0;
  _EMCAL_4x4_BEST_PHI = 0;
}

CaloTriggerInfo::~CaloTriggerInfo()
{
}

void CaloTriggerInfo::identify(ostream& os) const
{
  os << "CaloTriggerInfo: highest 2x2 eta/phi = " << _EMCAL_4x4_BEST_ETA << " / " << _EMCAL_4x4_BEST_PHI << ", E = " << _EMCAL_4x4_BEST_E << std::endl;

  return;
}
