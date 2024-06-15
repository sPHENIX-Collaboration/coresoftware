#include "MvtxRawEvtHeaderv1.h"

void MvtxRawEvtHeaderv1::Reset()
{
  m_MvtxFeeIdSet.clear();
  m_MvtxL1TrgSet.clear();
}

void MvtxRawEvtHeaderv1::identify(std::ostream& os) const
{
  os << "MvtxRawEvtHeaderv1" << std::endl;
  os << "Feeids in the event:  " << m_MvtxFeeIdSet.size() << std::endl;
  os << "Mvtx L1 triggers in the event:  " << m_MvtxL1TrgSet.size() << std::endl;
}

int MvtxRawEvtHeaderv1::isValid() const
{
  return m_MvtxFeeIdSet.size();
}

void MvtxRawEvtHeaderv1::AddFeeId(const std::set<uint16_t>& mvtxFeeIds)
{
  m_MvtxFeeIdSet.insert(mvtxFeeIds.cbegin(), mvtxFeeIds.cend());
}

void MvtxRawEvtHeaderv1::AddL1Trg(const std::set<uint64_t>& mvtxL1TrgSet)
{
  m_MvtxL1TrgSet.insert(mvtxL1TrgSet.cbegin(), mvtxL1TrgSet.cend());
}
