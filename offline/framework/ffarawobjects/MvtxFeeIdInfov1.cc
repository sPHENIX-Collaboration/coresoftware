#include "MvtxFeeIdInfov1.h"

MvtxFeeIdInfov1::MvtxFeeIdInfov1(MvtxFeeIdInfo *feeIdInfo)
{
  set_feeId(feeIdInfo->get_feeId());
  set_detField(feeIdInfo->get_detField());
  set_bco(feeIdInfo->get_bco());
}

void MvtxFeeIdInfov1::identify(std::ostream &os) const
{
  os << "FeeId: " << m_feeId << std::endl;
  os << "BCO: 0x" << std::hex << m_bco << std::dec << std::endl;
  os << "detField: " << m_detField << std::endl;
}
